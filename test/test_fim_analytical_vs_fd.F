c     ============================================================
c     test_fim_analytical_vs_fd.F -- issue #4: does the analytical
c     KS-space gradient (range_grad_ks_analytical, textbook dr=2L(u)du
c     identity) agree with the finite-difference one (range_grad_ks)
c     already used for Phase 1's findings?
c
c     If yes: Phase 1's conclusions (rank-deficiency, constant
c     rank-corrected ratio) are independently confirmed via a second,
c     analytically-derived method, not just numerical differencing.
c     ============================================================
      program test_fim_analytical_vs_fd
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      double precision a_kep(7), x0(3), xd0(3)
      double precision u(4), us(4), w
      double precision stat_pos(3,4), sigma(4)
      double precision grad_fd(4), grad_an(4), stat_r(3)
      double precision Ffd(4,4), Fan(4,4)
      double precision eigfd(4), eigan(4)
      double precision reddet_fd, reddet_an
      integer rank_fd, rank_an, ierr, ierr2, npass, nfail
      integer icase, k

      call init_constants()
      npass = 0
      nfail = 0

      write(*,*) '=== test_fim_analytical_vs_fd ==='

      do icase = 1, 3

         if (icase .eq. 1) then
            a_kep(1) = 7000.d0
            a_kep(2) = 0.01d0
         else if (icase .eq. 2) then
            a_kep(1) = 24500.d0
            a_kep(2) = 0.72d0
         else
            a_kep(1) = 26000.d0
            a_kep(2) = 0.85d0
         end if
         a_kep(3) = 15.d0
         a_kep(4) = 40.d0
         a_kep(5) = 60.d0
         a_kep(6) = 30.d0
         a_kep(7) = 0.d0

         call oe2car(a_kep, x0, xd0, 1.d-12)
         w = dsqrt(-0.5d0*(amue/a_kep(1)))
         call car2ks(x0, xd0, u, us, w)

c        Z1: per-station gradient agreement (single station, direct
c        component-by-component check)
         stat_r(1) = R_earth*0.6d0
         stat_r(2) = -R_earth*0.3d0
         stat_r(3) = R_earth*0.4d0

         call range_grad_ks(u, w, stat_r, grad_fd, ierr)
         call range_grad_ks_analytical(u, w, stat_r, grad_an, ierr2)

         gmag = dsqrt(grad_fd(1)**2+grad_fd(2)**2+grad_fd(3)**2+
     &                grad_fd(4)**2)
         diffmag = 0.d0
         do k = 1, 4
            diffmag = diffmag + (grad_fd(k)-grad_an(k))**2
         end do
         diffmag = dsqrt(diffmag)

         write(*,'(A,I1,A,E12.4)') ' case ', icase,
     &      ': |grad_fd - grad_analytical| / |grad_fd| = ',
     &      diffmag/gmag

         if (diffmag/gmag .lt. 1.d-5) then
            npass = npass + 1
         else
            write(*,'(A,I1)') ' [FAIL] Z1 gradient agreement, case ',
     &         icase
            nfail = nfail + 1
         end if

c        Z2: full 4-station FIM agreement -- rank and reduced
c        determinant should match between the two methods
         do k = 1, 4
            ang = dble(k)*1.57d0 + 0.4d0*dble(icase)
            stat_pos(1,k) = R_earth*dcos(ang)
            stat_pos(2,k) = R_earth*dsin(ang)*dcos(0.3d0*dble(k))
            stat_pos(3,k) = R_earth*dsin(0.3d0*dble(k))
            sigma(k) = 0.01d0*dble(k)
         end do

         call fim_build_ks(u, w, 4, stat_pos, sigma, Ffd, ierr)
         call fim_reduced_det_ks(Ffd, reddet_fd, rank_fd, eigfd,
     &                           ierr)

         call fim_build_ks_analytical(u, w, 4, stat_pos, sigma, Fan,
     &                                  ierr)
         call fim_reduced_det_ks(Fan, reddet_an, rank_an, eigan,
     &                           ierr2)

         write(*,'(A,I1,A,I1,A,I1)') ' case ', icase,
     &      ': rank_fd=', rank_fd, ' rank_analytical=', rank_an
         write(*,'(A,E16.8,A,E16.8)') '   reduced_det_fd=', reddet_fd,
     &      '  reduced_det_analytical=', reddet_an

         if (rank_fd .eq. rank_an) then
            npass = npass + 1
         else
            write(*,'(A,I1)') ' [FAIL] Z2 rank agreement, case ',
     &         icase
            nfail = nfail + 1
         end if

         reldiff = dabs(reddet_fd-reddet_an)/reddet_fd
         if (reldiff .lt. 1.d-4) then
            npass = npass + 1
         else
            write(*,'(A,I1,A,E12.4)') ' [FAIL] Z3 reduced-det '//
     &         'agreement, case ', icase, ' reldiff=', reldiff
            nfail = nfail + 1
         end if
      end do

      write(*,'(A,I2,A,I2)') ' Total pass=', npass, ' fail=', nfail
      if (nfail .gt. 0) then
         write(*,*) 'TESTS FAILED'
         stop 1
      end if
      write(*,*) 'ALL TESTS PASSED'

      end
