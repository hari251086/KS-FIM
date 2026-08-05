c     ============================================================
c     test_kalman.F -- Phase 6 (issue #5 EKF study): validates
c     kalman_scalar_update (src/kalman.F) against a hand-computable
c     1D case and basic sanity properties (variance never increases,
c     symmetry preserved) before it's trusted for the EKF study.
c     ============================================================
      program test_kalman
      implicit double precision (a-h, o-z)

      double precision P1(1,1), H1(1), sigma_prior, sigma_meas
      double precision expected, reldiff
      double precision P3(3,3), H3(3), tr_before, tr_after
      double precision asym
      integer i, j, k, npass, nfail

      npass = 0
      nfail = 0

      write(*,*) '=== test_kalman ==='

c     Z_a: 1D case -- standard information-additivity formula
c     1/P_post = 1/P_prior + 1/sigma_meas^2
      sigma_prior = 2.d0
      sigma_meas = 1.d0
      P1(1,1) = sigma_prior*sigma_prior
      H1(1) = 1.d0
      call kalman_scalar_update(P1, 1, H1, sigma_meas)
      expected = 1.d0/(1.d0/(sigma_prior*sigma_prior)
     &   + 1.d0/(sigma_meas*sigma_meas))
      reldiff = dabs(P1(1,1)-expected)/expected
      write(*,'(A,E16.8,A,E16.8,A,E10.2)') ' 1D: P_post=', P1(1,1),
     &   '  expected=', expected, '  reldiff=', reldiff
      if (reldiff .lt. 1.d-12) then
         npass = npass + 1
      else
         write(*,*) ' [FAIL] 1D scalar update'
         nfail = nfail + 1
      end if

c     Z_b: 3D case -- variance along the observed direction must not
c     increase, matrix stays symmetric
      do i = 1, 3
         do j = 1, 3
            P3(i,j) = 0.d0
         end do
         P3(i,i) = 4.d0
      end do
      H3(1) = 1.d0
      H3(2) = 0.5d0
      H3(3) = 0.d0

      tr_before = P3(1,1)+P3(2,2)+P3(3,3)
      call kalman_scalar_update(P3, 3, H3, 1.d0)
      tr_after = P3(1,1)+P3(2,2)+P3(3,3)

      write(*,'(A,E16.8,A,E16.8)') ' 3D: trace before=', tr_before,
     &   '  after=', tr_after
      if (tr_after .lt. tr_before) then
         npass = npass + 1
      else
         write(*,*) ' [FAIL] 3D update did not reduce trace'
         nfail = nfail + 1
      end if

      asym = 0.d0
      do i = 1, 3
         do j = 1, 3
            asym = dmax1(asym, dabs(P3(i,j)-P3(j,i)))
         end do
      end do
      write(*,'(A,E12.4)') ' 3D: max asymmetry=', asym
      if (asym .lt. 1.d-13) then
         npass = npass + 1
      else
         write(*,*) ' [FAIL] 3D update broke symmetry'
         nfail = nfail + 1
      end if

      write(*,'(A,I3,A,I3)') ' Total pass=', npass, ' fail=', nfail
      if (nfail .gt. 0) then
         write(*,*) 'TESTS FAILED'
         stop 1
      end if
      write(*,*) 'ALL TESTS PASSED'

      end
