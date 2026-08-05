c     ============================================================
c     test_fim_reduced_consistency.F -- is the "KS observability
c     advantage" a real information gain, or a predictable Jacobian-
c     scaling artifact?
c
c     Analytical prediction (ALGORITHM.md section 2): since
c     F_ks = J^T F_cartesian J with J = dr/du evaluated at the
c     satellite's current u, and the KS matrix L(u) satisfies
c     L(u)^T L(u) = |u|^2 * I (a standard KS identity -- the map is
c     conformal, a uniformly-scaled orthogonal transformation), the
c     ratio (product of F_ks's top-3 eigenvalues) / det(F_cartesian)
c     should be a CONSTANT depending only on the satellite's current
c     position (via |u|^2=r), not on which stations are used or how
c     many. If that ratio holds constant across very different
c     station geometries for the same satellite position, the
c     "advantage" is a coordinate-scaling artifact, not new physical
c     information about the satellite's observability.
c     ============================================================
      program test_fim_reduced_consistency
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      double precision a_kep(7), x0(3), xd0(3)
      double precision u(4), us(4), w
      double precision stat_pos(3,4), sigma(4)
      double precision Fc(3,3), Fk(4,4), eigval(4), eigval3(3)
      double precision det_c, reduced_det, ratio(5), ratio_mean
      double precision ratio_max_dev
      integer rank_est, rank_c, ierr, ierr2, npass, nfail, ncfg, k

      call init_constants()
      npass = 0
      nfail = 0

      write(*,*) '=== test_fim_reduced_consistency ==='

c     One fixed satellite position (GTO-like, matches the regime the
c     COSPAR presentation's HEO case studies used).
      a_kep(1) = 24500.d0
      a_kep(2) = 0.72d0
      a_kep(3) = 15.d0
      a_kep(4) = 40.d0
      a_kep(5) = 60.d0
      a_kep(6) = 30.d0
      a_kep(7) = 0.d0

      call oe2car(a_kep, x0, xd0, 1.d-12)
      w = dsqrt(-0.5d0*(amue/a_kep(1)))
      call car2ks(x0, xd0, u, us, w)

c     Five very different station configurations at the SAME satellite
c     position/epoch -- different count, different geometry, different
c     sigma values, deliberately not matched to each other at all.
      ncfg = 5
      do icfg = 1, ncfg

         if (icfg .eq. 1) then
            nstat = 3
            stat_pos(1,1)=R_earth*0.7d0
            stat_pos(2,1)=R_earth*0.3d0
            stat_pos(3,1)=R_earth*0.2d0
            sigma(1)=0.01d0
            stat_pos(1,2)=-R_earth*0.5d0
            stat_pos(2,2)=R_earth*0.6d0
            stat_pos(3,2)=-R_earth*0.4d0
            sigma(2)=0.05d0
            stat_pos(1,3)=R_earth*0.1d0
            stat_pos(2,3)=-R_earth*0.2d0
            stat_pos(3,3)=R_earth*0.9d0
            sigma(3)=0.08d0
         else if (icfg .eq. 2) then
            nstat = 3
            stat_pos(1,1)=R_earth
            stat_pos(2,1)=0.d0
            stat_pos(3,1)=0.d0
            sigma(1)=0.02d0
            stat_pos(1,2)=0.d0
            stat_pos(2,2)=R_earth
            stat_pos(3,2)=0.d0
            sigma(2)=0.02d0
            stat_pos(1,3)=0.d0
            stat_pos(2,3)=0.d0
            stat_pos(3,3)=R_earth
            sigma(3)=0.02d0
         else if (icfg .eq. 3) then
            nstat = 4
            do k = 1, 4
               ang = dble(k)*1.57d0
               stat_pos(1,k) = R_earth*dcos(ang)
               stat_pos(2,k) = R_earth*dsin(ang)
               stat_pos(3,k) = R_earth*0.1d0*dble(k)
               sigma(k) = 0.001d0*dble(k)
            end do
         else if (icfg .eq. 4) then
            nstat = 3
            stat_pos(1,1)=R_earth*0.9d0
            stat_pos(2,1)=-R_earth*0.1d0
            stat_pos(3,1)=R_earth*0.05d0
            sigma(1)=0.1d0
            stat_pos(1,2)=R_earth*0.85d0
            stat_pos(2,2)=-R_earth*0.15d0
            stat_pos(3,2)=R_earth*0.06d0
            sigma(2)=0.1d0
            stat_pos(1,3)=-R_earth*0.3d0
            stat_pos(2,3)=R_earth*0.7d0
            stat_pos(3,3)=-R_earth*0.2d0
            sigma(3)=0.04d0
         else
            nstat = 3
            stat_pos(1,1)=-R_earth*0.6d0
            stat_pos(2,1)=-R_earth*0.6d0
            stat_pos(3,1)=R_earth*0.3d0
            sigma(1)=0.03d0
            stat_pos(1,2)=R_earth*0.2d0
            stat_pos(2,2)=-R_earth*0.9d0
            stat_pos(3,2)=-R_earth*0.1d0
            sigma(2)=0.07d0
            stat_pos(1,3)=R_earth*0.4d0
            stat_pos(2,3)=R_earth*0.4d0
            stat_pos(3,3)=R_earth*0.8d0
            sigma(3)=0.015d0
         end if

         call fim_build_cartesian(x0, nstat, stat_pos, sigma, Fc,
     &                             ierr)
         call fim_reduced_det_cartesian(Fc, det_c, rank_c, eigval3,
     &                                   ierr)

         call fim_build_ks(u, w, nstat, stat_pos, sigma, Fk, ierr)
         call fim_reduced_det_ks(Fk, reduced_det, rank_est, eigval,
     &                           ierr2)

         ratio(icfg) = reduced_det/det_c
         write(*,'(A,I1,A,I1,A,E16.8)') ' cfg=', icfg, ' nstat=',
     &      nstat, '  reduced_det_ks/det_cartesian = ', ratio(icfg)
      end do

      ratio_mean = 0.d0
      do icfg = 1, ncfg
         ratio_mean = ratio_mean + ratio(icfg)
      end do
      ratio_mean = ratio_mean/dble(ncfg)

      ratio_max_dev = 0.d0
      do icfg = 1, ncfg
         dev = dabs(ratio(icfg)-ratio_mean)/ratio_mean
         if (dev .gt. ratio_max_dev) ratio_max_dev = dev
      end do

      write(*,'(A,E16.8)') ' mean ratio    = ', ratio_mean
      write(*,'(A,E16.8)') ' max rel. dev  = ', ratio_max_dev

c     Z1: ratio should be constant (same satellite position) across
c     wildly different station geometries/counts/sigmas -- a few
c     percent tolerance for finite-difference/eigensolver noise.
      if (ratio_max_dev .lt. 1.d-3) then
         write(*,*) '[PASS] Z1 ratio constant across station configs'
         npass = npass + 1
      else
         write(*,*) '[FAIL] Z1 ratio constant across station configs'
         nfail = nfail + 1
      end if

      write(*,'(A,I2,A,I2)') ' Total pass=', npass, ' fail=', nfail
      if (nfail .gt. 0) then
         write(*,*) 'TESTS FAILED'
         stop 1
      end if
      write(*,*) 'ALL TESTS PASSED'

      end
