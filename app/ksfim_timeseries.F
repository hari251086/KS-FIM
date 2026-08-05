c     ============================================================
c     ksfim_timeseries.F -- issue #2: multi-revolution observability
c     tracking, matching the source presentation's own "log(det) vs.
c     time/revolution" plots, but reporting the rank-corrected metric
c     alongside the raw one, and checking the exact closed-form
c     prediction (ALGORITHM.md section 5 / test_fim_ratio_formula)
c     `reduced_det_ks/det_cartesian = 64*r^3` at every sampled point
c     across 2 full revolutions of the 35497 case-study object.
c
c     Simple two-body Keplerian propagation (mean anomaly advanced
c     linearly at the orbit's own mean motion, fed to KSROP's oe2car
c     each step) -- sufficient for this repo's purpose, since the
c     structural rank-deficiency/closed-form-ratio findings are
c     properties of the KS<->Cartesian coordinate map evaluated at
c     whatever position the satellite currently occupies, independent
c     of which force model produced that position (see ALGORITHM.md
c     section 9). Station geometry held fixed at the object's real
c     epoch GMST (not re-rotated per sample -- 2 revolutions of a
c     ~12h-period HEO object is ~day-scale, station rotation over that
c     span is a secondary effect already shown (issue #3) not to
c     affect the corrected ratio at all).
c
c     Output: printed summary + output/ksfim_timeseries_35497.csv.
c     ============================================================
      program ksfim_timeseries
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      double precision a_kep(7), x0(3), xd0(3)
      double precision u(4), us(4), w
      parameter(nstat = 4)
      double precision stat_lat(nstat), stat_lon(nstat)
      double precision stat_alt(nstat), sigma(nstat)
      double precision stat_pos(3,nstat)
      double precision Fc(3,3), Fk(4,4)
      double precision eigval(4), eigval_c(3)
      double precision det_c, reduced_det, rmag
      double precision actual_ratio, predicted_ratio, reldiff
      double precision rmin, rmax, reldiff_max
      integer rank_est, rank_c, ierr, ierr2, k

      parameter(nrev = 2, nstep_per_rev = 72)
      parameter(nstep = nrev*nstep_per_rev)
      double precision m0, mstep, mnow
      integer istep

      double precision epoch_cal(6), jd_epoch, gmst

      call init_constants()

c     35497 (Ariane 5 ESC-A), same elements as ksfim_case_study.F's
c     case 1 -- picked as the representative HEO object for this
c     time-series check.
      a_kep(1) = (R_earth+243.13d0)/(1.d0-0.6049d0)
      a_kep(2) = 0.6049d0
      a_kep(3) = 5.7825d0
      a_kep(4) = 346.1476d0
      a_kep(5) = 154.0352d0
      m0 = 331.2056d0
      a_kep(7) = 0.d0

      mstep = 360.d0*dble(nrev)/dble(nstep)

      epoch_cal(1) = 2015.d0
      epoch_cal(2) = 6.d0
      epoch_cal(3) = 6.d0
      epoch_cal(4) = 14.d0
      epoch_cal(5) = 7.d0
      epoch_cal(6) = 22.d0
      call cal2jd(epoch_cal, jd_epoch)
      call gmst_deg(jd_epoch, gmst)

      stat_lat(1) = 8.717d0
      stat_lon(1) = 167.733d0
      stat_alt(1) = 0.01d0
      sigma(1) = 10.d0
      stat_lat(2) = 7.972d0
      stat_lon(2) = -14.401d0
      stat_alt(2) = 0.01d0
      sigma(2) = 10.d0
      stat_lat(3) = 13.719d0
      stat_lon(3) = 80.23d0
      stat_alt(3) = 0.01d0
      sigma(3) = 10.d0
      stat_lat(4) = 28.3d0
      stat_lon(4) = -16.509d0
      stat_alt(4) = 2.39d0
      sigma(4) = 10.d0

      do k = 1, nstat
         call geodetic_to_eci(stat_lat(k), stat_lon(k), stat_alt(k),
     &      R_earth, gmst, stat_pos(1,k))
      end do

      write(*,*) '======================================'
      write(*,*) ' KS-FIM multi-revolution time series'
      write(*,*) ' 35497 Ariane 5 ESC-A, 2 revolutions,'
      write(*,'(A,I3,A)') '  ', nstep, ' sample points'
      write(*,*) '======================================'
      write(*,*)

      open(unit=44, file='output/ksfim_timeseries_35497.csv',
     &     status='unknown')
      write(44,'(A)') 'step,manom_deg,r_km,rank_c,rank_ks,'//
     &   'det_cartesian,reduced_det_ks,actual_ratio,'//
     &   'predicted_64r3,reldiff'

      rmin = 1.d30
      rmax = -1.d30
      reldiff_max = 0.d0

      do istep = 0, nstep-1
         mnow = m0 + mstep*dble(istep)
         mnow = dmod(mnow, 360.d0)
         a_kep(6) = mnow

         call oe2car(a_kep, x0, xd0, 1.d-12)
         w = dsqrt(-0.5d0*(amue/a_kep(1)))
         call car2ks(x0, xd0, u, us, w)

         rmag = dsqrt(x0(1)**2+x0(2)**2+x0(3)**2)
         if (rmag .lt. rmin) rmin = rmag
         if (rmag .gt. rmax) rmax = rmag

         call fim_build_cartesian(x0, nstat, stat_pos, sigma, Fc,
     &                             ierr)
         call fim_reduced_det_cartesian(Fc, det_c, rank_c, eigval_c,
     &                                    ierr)

         call fim_build_ks(u, w, nstat, stat_pos, sigma, Fk, ierr)
         call fim_reduced_det_ks(Fk, reduced_det, rank_est, eigval,
     &                            ierr2)

         actual_ratio = reduced_det/det_c
         predicted_ratio = 64.d0*rmag**3
         reldiff = dabs(actual_ratio-predicted_ratio)/predicted_ratio
         if (reldiff .gt. reldiff_max) reldiff_max = reldiff

         write(44,'(I4,A,F10.4,A,F12.3,A,I1,A,I1,A,E16.8,A,E16.8,'//
     &      'A,E16.8,A,E16.8,A,E12.4)')
     &      istep, ',', mnow, ',', rmag, ',', rank_c, ',', rank_est,
     &      ',', det_c, ',', reduced_det, ',', actual_ratio, ',',
     &      predicted_ratio, ',', reldiff
      end do

      close(44)

      write(*,'(A,F10.2,A,F10.2,A)') '  r range over 2 revolutions: ',
     &   rmin, ' - ', rmax, ' km'
      write(*,'(A,E12.4)') '  max |actual - predicted 64r^3| / '//
     &   'predicted, over all samples = ', reldiff_max
      if (reldiff_max .lt. 1.d-4) then
         write(*,*) '  [PASS] closed-form ratio=64*r^3 holds at'//
     &      ' every sampled point across both revolutions'
      else
         write(*,*) '  [FAIL] closed-form ratio broke down '//
     &      'somewhere in the time series'
      end if
      write(*,*)
      write(*,*) ' full time series -> output/'//
     &   'ksfim_timeseries_35497.csv'

      end
