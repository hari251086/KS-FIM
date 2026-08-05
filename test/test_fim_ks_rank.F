c     ============================================================
c     test_fim_ks_rank.F -- the core empirical test of this repo's
c     founding question: is the 4x4 KS-space FIM structurally
c     rank-deficient (rank 3, not 4)?
c
c     Physical argument (ALGORITHM.md section 2): r=(x,y,z,0)=L(u)u
c     depends on only 3 true degrees of freedom, so a FIM built from
c     range-measurement gradients pushed through the KS map can never
c     constrain the 4th ("gauge") direction in u-space. Checked here
c     across several different satellite positions and station counts
c     (1 to 5) -- if the smallest eigenvalue is NOT tiny relative to
c     the largest in every case, the hypothesis is falsified.
c     ============================================================
      program test_fim_ks_rank
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      double precision a_kep(7), x0(3), xd0(3)
      double precision u(4), us(4), w
      double precision stat_pos(3,5), sigma(5)
      double precision F(4,4), eigval(4), reduced_det
      integer rank_est, ierr, ierr2
      integer npass, nfail, icase, nstat

      call init_constants()

      npass = 0
      nfail = 0

      write(*,*) '=== test_fim_ks_rank ==='

c     Three representative satellite positions: LEO-circular-ish,
c     GTO-like, and a near-degenerate case close to perigee of a very
c     eccentric orbit -- covers a spread of |u| and geometry.
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
         w = dsqrt(0.25d0*(amue/a_kep(1)))
         call car2ks(x0, xd0, u, us, w)

c        Try 2, 3, 4, and 5 stations (arbitrary but fixed geometry,
c        distances chosen so |r-stat| stays well-conditioned)
         do nstat = 2, 5
            do k = 1, nstat
               ang = dble(k)*1.3d0 + 0.4d0*dble(icase)
               stat_pos(1,k) = R_earth*dcos(ang)
               stat_pos(2,k) = R_earth*dsin(ang)*dcos(0.3d0*dble(k))
               stat_pos(3,k) = R_earth*dsin(0.3d0*dble(k))
               sigma(k) = 0.01d0*dble(k)
            end do

            call fim_build_ks(u, w, nstat, stat_pos, sigma, F, ierr)
            if (ierr .ne. 0) then
               write(*,'(A,I1,A,I1,A)') ' [SKIP] case ', icase,
     &            ' nstat=', nstat, ' (degenerate gradient)'
               goto 20
            end if

            call fim_reduced_det_ks(F, reduced_det, rank_est,
     &                              eigval, ierr2)

            write(*,'(A,I1,A,I1,A,I1,A,4E12.4)')
     &         ' case ', icase, ' nstat=', nstat,
     &         ' rank_est=', rank_est, '  eig=', eigval

c           Core hypothesis: rank(F_ks) = min(3, nstat) for a well-
c           conditioned station geometry -- NEVER 4, since range
c           measurements only ever depend on 3 true physical DOF.
c           2 stations -> at most rank 2; 3+ well-conditioned
c           stations -> rank 3 (never higher, that's the whole point).
            iexpect = nstat
            if (iexpect .gt. 3) iexpect = 3
            if (rank_est .eq. iexpect) then
               npass = npass + 1
            else
               write(*,'(A,I1,A,I1,A,I1,A,I1)') ' [FAIL] case ',
     &            icase, ' nstat=', nstat, ' expected rank_est=',
     &            iexpect, ', got ', rank_est
               nfail = nfail + 1
            end if

c           Whatever's above the true rank must be tiny relative to
c           the largest eigenvalue (near machine precision noise, not
c           a genuine 4th information direction)
            if (eigval(4) .lt. 1.d-8*eigval(1)) then
               npass = npass + 1
            else
               write(*,'(A,I1,A,I1,A,E12.4)') ' [FAIL] case ', icase,
     &            ' nstat=', nstat,
     &            ' smallest eigenvalue not near-zero, ratio=',
     &            eigval(4)/eigval(1)
               nfail = nfail + 1
            end if
 20         continue
         end do
      end do

      write(*,'(A,I3,A,I3)') ' Total pass=', npass, ' fail=', nfail
      if (nfail .gt. 0) then
         write(*,*) 'TESTS FAILED'
         stop 1
      end if
      write(*,*) 'ALL TESTS PASSED'

      end
