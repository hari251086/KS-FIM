c     ============================================================
c     test_fim_cartesian.F -- sanity check for fim_build_cartesian
c     against a hand-computable orthogonal-LOS geometry.
c     ============================================================
      program test_fim_cartesian
      implicit double precision (a-h, o-z)

      double precision r(3), stat_pos(3,3), sigma(3), F(3,3)
      double precision expect_det, got_det, tol
      integer ierr, npass, nfail

      npass = 0
      nfail = 0

c     Satellite at origin. Three stations placed so the unit
c     line-of-sight vectors are exactly the +x, +y, +z axes:
c     station k sits at distance dk along axis k, sigma_k arbitrary
c     but distinct, so F = diag(1/sigma1^2, 1/sigma2^2, 1/sigma3^2)
c     exactly -- det(F) = 1/(s1^2 s2^2 s3^2) is hand-checkable.
      r(1) = 0.d0
      r(2) = 0.d0
      r(3) = 0.d0

      stat_pos(1,1) = -500.d0
      stat_pos(2,1) = 0.d0
      stat_pos(3,1) = 0.d0
      sigma(1) = 2.d0

      stat_pos(1,2) = 0.d0
      stat_pos(2,2) = -800.d0
      stat_pos(3,2) = 0.d0
      sigma(2) = 3.d0

      stat_pos(1,3) = 0.d0
      stat_pos(2,3) = 0.d0
      stat_pos(3,3) = -1200.d0
      sigma(3) = 5.d0

      call fim_build_cartesian(r, 3, stat_pos, sigma, F, ierr)

      write(*,*) '=== test_fim_cartesian ==='

c     Z1: off-diagonal terms must vanish for this orthogonal geometry
      tol = 1.d-12
      if (dabs(F(1,2)) .lt. tol .and. dabs(F(1,3)) .lt. tol .and.
     &    dabs(F(2,3)) .lt. tol) then
         write(*,*) '[PASS] Z1 off-diagonals ~0 for orthogonal LOS'
         npass = npass + 1
      else
         write(*,*) '[FAIL] Z1 off-diagonals ~0 for orthogonal LOS'
         nfail = nfail + 1
      end if

c     Z2: diagonal terms match 1/sigma_k^2 exactly
      if (dabs(F(1,1) - 1.d0/(sigma(1)*sigma(1))) .lt. tol .and.
     &    dabs(F(2,2) - 1.d0/(sigma(2)*sigma(2))) .lt. tol .and.
     &    dabs(F(3,3) - 1.d0/(sigma(3)*sigma(3))) .lt. tol) then
         write(*,*) '[PASS] Z2 diagonal = 1/sigma^2'
         npass = npass + 1
      else
         write(*,*) '[FAIL] Z2 diagonal = 1/sigma^2'
         nfail = nfail + 1
      end if

c     Z3: determinant matches the hand-computed product
      expect_det = 1.d0/(sigma(1)**2 * sigma(2)**2 * sigma(3)**2)
      got_det = F(1,1)*F(2,2)*F(3,3)
      if (dabs(got_det - expect_det)/expect_det .lt. 1.d-10) then
         write(*,*) '[PASS] Z3 det(F) matches hand-computed value'
         npass = npass + 1
      else
         write(*,*) '[FAIL] Z3 det(F) matches hand-computed value'
         nfail = nfail + 1
      end if

      write(*,'(A,I2,A,I2)') ' Total pass=', npass, ' fail=', nfail
      if (nfail .gt. 0) then
         write(*,*) 'TESTS FAILED'
         stop 1
      end if
      write(*,*) 'ALL TESTS PASSED'

      end
