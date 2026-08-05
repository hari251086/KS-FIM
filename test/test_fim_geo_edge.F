c     ============================================================
c     test_fim_geo_edge.F -- issue #1: extreme low-station-count edge
c     case (1 and 2 stations, near-GEO orbit) -- does rank(F_ks) still
c     track rank(F_cartesian) exactly down to rank 1, and does the raw
c     4x4 determinant correctly go to (numerical) zero whenever the
c     true rank is below 4, as ALGORITHM.md section 5 predicts?
c
c     Also directly relevant to the source presentation's own noted
c     caveat ("if two stations are very near, the observability index
c     can go to zero, misconstrued as a single observing station") --
c     this test shows the raw determinant is zero at low station
c     count generally, not just for near-degenerate close stations.
c     ============================================================
      program test_fim_geo_edge
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      double precision a_kep(7), x0(3), xd0(3)
      double precision u(4), us(4), w
      double precision stat_pos(3,2), sigma(2)
      double precision Fc(3,3), Fk(4,4)
      double precision eigval(4), eigval_c(3)
      double precision det_c, det_k_raw, reduced_det
      integer rank_est, rank_c, ierr, ierr2, npass, nfail, nstat, k

      call init_constants()
      npass = 0
      nfail = 0

      write(*,*) '=== test_fim_geo_edge ==='

c     Near-GEO orbit (same regime as the case-study GEO addition)
      a_kep(1) = 42164.14d0
      a_kep(2) = 0.0002d0
      a_kep(3) = 0.05d0
      a_kep(4) = 90.d0
      a_kep(5) = 0.d0
      a_kep(6) = 45.d0
      a_kep(7) = 0.d0

      call oe2car(a_kep, x0, xd0, 1.d-12)
      w = dsqrt(-0.5d0*(amue/a_kep(1)))
      call car2ks(x0, xd0, u, us, w)

      stat_pos(1,1) = R_earth*0.6d0
      stat_pos(2,1) = -R_earth*0.7d0
      stat_pos(3,1) = R_earth*0.3d0
      sigma(1) = 10.d0

      stat_pos(1,2) = -R_earth*0.4d0
      stat_pos(2,2) = R_earth*0.8d0
      stat_pos(3,2) = R_earth*0.2d0
      sigma(2) = 10.d0

      do nstat = 1, 2
         call fim_build_cartesian(x0, nstat, stat_pos, sigma, Fc,
     &                             ierr)
         call fim_reduced_det_cartesian(Fc, det_c, rank_c, eigval_c,
     &                                    ierr)

         call fim_build_ks(u, w, nstat, stat_pos, sigma, Fk, ierr)
         call fim_reduced_det_ks(Fk, reduced_det, rank_est, eigval,
     &                            ierr2)
         det_k_raw = eigval(1)*eigval(2)*eigval(3)*eigval(4)

         write(*,'(A,I1,A,I1,A,I1,A,E12.4)') ' nstat=', nstat,
     &      ': rank_c=', rank_c, ' rank_ks=', rank_est,
     &      ' raw_det_ks=', det_k_raw

c        Z1: rank(F_ks) must equal rank(F_cartesian), all the way
c        down to nstat=1 (rank 1)
         if (rank_est .eq. rank_c .and. rank_est .eq. nstat) then
            npass = npass + 1
         else
            write(*,'(A,I1)') ' [FAIL] Z1 rank tracking, nstat=',
     &         nstat
            nfail = nfail + 1
         end if

c        Z2: the raw 4x4 determinant must be numerically zero
c        whenever true rank < 4 (always true here, nstat<3)
         if (dabs(det_k_raw) .lt. 1.d-6*eigval(1)) then
            npass = npass + 1
         else
            write(*,'(A,I1)') ' [FAIL] Z2 raw det ~0, nstat=', nstat
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
