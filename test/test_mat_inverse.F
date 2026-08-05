c     ============================================================
c     test_mat_inverse.F -- Phase 6 (issue #5 EKF study): validates
c     mat_inverse (src/linalg.F, Gauss-Jordan with partial pivoting)
c     against A*A^-1=I for 3x3 and 4x4 cases before it's trusted for
c     the EKF study's Schur-complement position marginalization.
c     ============================================================
      program test_mat_inverse
      implicit double precision (a-h, o-z)

      double precision A3(3,3), Ainv3(3,3), I3(3,3)
      double precision A4(4,4), Ainv4(4,4), I4(4,4)
      double precision maxdiff
      integer ierr, i, j, k, npass, nfail

      npass = 0
      nfail = 0

      write(*,*) '=== test_mat_inverse ==='

c     3x3 SPD test case
      A3(1,1) = 4.d0
      A3(1,2) = 1.d0
      A3(1,3) = 0.5d0
      A3(2,1) = 1.d0
      A3(2,2) = 3.d0
      A3(2,3) = 0.7d0
      A3(3,1) = 0.5d0
      A3(3,2) = 0.7d0
      A3(3,3) = 2.d0

      call mat_inverse(A3, 3, 3, Ainv3, ierr)

      do i = 1, 3
         do j = 1, 3
            I3(i,j) = 0.d0
            do k = 1, 3
               I3(i,j) = I3(i,j) + A3(i,k)*Ainv3(k,j)
            end do
         end do
      end do

      maxdiff = 0.d0
      do i = 1, 3
         do j = 1, 3
            if (i .eq. j) then
               maxdiff = dmax1(maxdiff, dabs(I3(i,j)-1.d0))
            else
               maxdiff = dmax1(maxdiff, dabs(I3(i,j)))
            end if
         end do
      end do

      write(*,'(A,I2,A,E12.4)') ' ierr=', ierr,
     &   '  3x3 A*Ainv vs I maxdiff=', maxdiff
      if (ierr .eq. 0 .and. maxdiff .lt. 1.d-10) then
         npass = npass + 1
      else
         write(*,*) ' [FAIL] 3x3 inverse'
         nfail = nfail + 1
      end if

c     4x4 SPD test case
      A4(1,1) = 5.d0
      A4(1,2) = 1.d0
      A4(1,3) = 0.3d0
      A4(1,4) = 0.2d0
      A4(2,1) = 1.d0
      A4(2,2) = 4.d0
      A4(2,3) = 0.6d0
      A4(2,4) = 0.1d0
      A4(3,1) = 0.3d0
      A4(3,2) = 0.6d0
      A4(3,3) = 3.d0
      A4(3,4) = 0.8d0
      A4(4,1) = 0.2d0
      A4(4,2) = 0.1d0
      A4(4,3) = 0.8d0
      A4(4,4) = 2.d0

      call mat_inverse(A4, 4, 4, Ainv4, ierr)

      do i = 1, 4
         do j = 1, 4
            I4(i,j) = 0.d0
            do k = 1, 4
               I4(i,j) = I4(i,j) + A4(i,k)*Ainv4(k,j)
            end do
         end do
      end do

      maxdiff = 0.d0
      do i = 1, 4
         do j = 1, 4
            if (i .eq. j) then
               maxdiff = dmax1(maxdiff, dabs(I4(i,j)-1.d0))
            else
               maxdiff = dmax1(maxdiff, dabs(I4(i,j)))
            end if
         end do
      end do

      write(*,'(A,I2,A,E12.4)') ' ierr=', ierr,
     &   '  4x4 A*Ainv vs I maxdiff=', maxdiff
      if (ierr .eq. 0 .and. maxdiff .lt. 1.d-10) then
         npass = npass + 1
      else
         write(*,*) ' [FAIL] 4x4 inverse'
         nfail = nfail + 1
      end if

      write(*,'(A,I3,A,I3)') ' Total pass=', npass, ' fail=', nfail
      if (nfail .gt. 0) then
         write(*,*) 'TESTS FAILED'
         stop 1
      end if
      write(*,*) 'ALL TESTS PASSED'

      end
