c     ============================================================
c     ksfim_case_study.F -- reproduces the 2021 COSPAR presentation's
c     4 highly-elliptical-orbit case studies (objects 35497, 37151,
c     39615, 42928 -- the same orbital elements published in
c     "Regularized Orbit Observability for Resident Space Objects",
c     PEDAS.1-0023-21) and reports numbers side by side per case: the
c     raw 4x4 det(F_ks) (what the original presentation reported,
c     structurally rank-deficient per this repo's ALGORITHM.md), the
c     rank-corrected reduced determinant, det(F_cartesian) for the
c     same scenario, and (issue #4) the same reduced determinant
c     recomputed via the independent analytical KS gradient as a
c     cross-check -- so the "KS observability advantage" claim can be
c     read honestly instead of via the inflated raw 4x4 number.
c
c     Ground stations, orbital elements, AND epochs transcribed
c     directly from the slide deck (E:\Research\Conferences\c2021 -
c     COSPAR\Paper and Presentation\DROPBOX-31122020\...
c     FULLPRESENTATIONSLIDES.pdf, "Case Studies" slides) -- not
c     re-derived, so results are directly comparable to the original.
c
c     Issue #3 (resolved): stations are now placed at their published
c     geodetic lat/lon, rotated into each object's own ECI frame via
c     GMST at that object's actual observation epoch (src/timeconv.F)
c     -- a spherical, not WGS84-oblate, Earth model, consistent with
c     the rest of this repo's simplified station-geometry treatment.
c     This does not change the structural findings (rank-deficiency
c     and Jacobian-scaling are properties of the coordinate MAP, not
c     of tracking geometry) but makes the reported numbers a genuine
c     per-object, per-epoch snapshot rather than a fixed offset.
c     ============================================================
      program ksfim_case_study
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      double precision a_kep(7), x0(3), xd0(3)
      double precision u(4), us(4), w
      parameter(nstat_max = 4)
      double precision stat_lat(nstat_max), stat_lon(nstat_max)
      double precision stat_alt(nstat_max)
      double precision stat_pos(3,nstat_max), sigma(nstat_max)
      double precision Fc(3,3), Fk(4,4), Fk_an(4,4)
      double precision eigval(4), eigval_an(4)
      double precision det_c, det_k_raw, reduced_det, reduced_det_an
      integer rank_est, rank_an, ierr, ierr2, icase, k, nobj
      double precision epoch_cal(6), jd_epoch, gmst

      parameter(nobj = 4)
      integer norad(nobj)
      double precision perigee_alt(nobj), ecc(nobj), inc(nobj)
      double precision raan(nobj), aop(nobj), manom(nobj)
      double precision epoch_yr(nobj), epoch_mo(nobj), epoch_dy(nobj)
      double precision epoch_hr(nobj), epoch_mn(nobj), epoch_sc(nobj)
      character*30 objname(nobj)

      call init_constants()

      norad(1) = 35497
      objname(1) = 'Ariane 5 ESC-A'
      perigee_alt(1) = 243.13d0
      ecc(1) = 0.6049d0
      inc(1) = 5.7825d0
      raan(1) = 346.1476d0
      aop(1) = 154.0352d0
      manom(1) = 331.2056d0
      epoch_yr(1) = 2015.d0
      epoch_mo(1) = 6.d0
      epoch_dy(1) = 6.d0
      epoch_hr(1) = 14.d0
      epoch_mn(1) = 7.d0
      epoch_sc(1) = 22.d0

      norad(2) = 37151
      objname(2) = 'Long March 3B Upper Stage'
      perigee_alt(2) = 145.09d0
      ecc(2) = 0.4099d0
      inc(2) = 24.8672d0
      raan(2) = 204.061d0
      aop(2) = 163.8209d0
      manom(2) = 341.4979d0
      epoch_yr(2) = 2015.d0
      epoch_mo(2) = 3.d0
      epoch_dy(2) = 4.d0
      epoch_hr(2) = 16.d0
      epoch_mn(2) = 10.d0
      epoch_sc(2) = 37.d0

      norad(3) = 39615
      objname(3) = 'Proton-M Briz-M Aux Tank'
      perigee_alt(3) = 143.5d0
      ecc(3) = 0.3367d0
      inc(3) = 48.4571d0
      raan(3) = 118.9947d0
      aop(3) = 333.8816d0
      manom(3) = 12.4622d0
      epoch_yr(3) = 2017.d0
      epoch_mo(3) = 3.d0
      epoch_dy(3) = 2.d0
      epoch_hr(3) = 4.d0
      epoch_mn(3) = 20.d0
      epoch_sc(3) = 21.d0

      norad(4) = 42928
      objname(4) = 'IRNSS-1H/PSLV-PS4'
      perigee_alt(4) = 156.33d0
      ecc(4) = 0.3095d0
      inc(4) = 19.1455d0
      raan(4) = 315.4336d0
      aop(4) = 125.4322d0
      manom(4) = 330.6148d0
      epoch_yr(4) = 2017.d0
      epoch_mo(4) = 11.d0
      epoch_dy(4) = 1.d0
      epoch_hr(4) = 4.d0
      epoch_mn(4) = 15.d0
      epoch_sc(4) = 9.d0

c     Kwajalein, Ascension Island, SHAR, Tenerife. Sigma=10km per the
c     slide's own case-study settings box ("Std. deviation: 10 km") --
c     Phase 1 mistakenly reused the station ALTITUDE column (0.01,
c     0.01, 0.01, 2.39 km) for sigma; fixed here (Phase 3). Purely a
c     magnitude-normalization fix -- sigma enters both F_cartesian and
c     F_ks identically as a uniform 1/sigma^2 scalar per station, so
c     it doesn't change rank or any ratio-based finding, only the
c     absolute determinant magnitudes printed.
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

      write(*,*) '======================================'
      write(*,*) ' KS-FIM case-study reproduction'
      write(*,*) ' (COSPAR 2021 PEDAS.1-0023-21 objects)'
      write(*,*) ' (issue #3: GMST-accurate station placement)'
      write(*,*) '======================================'
      write(*,*)

      do icase = 1, nobj
         a_kep(1) = (R_earth+perigee_alt(icase))/(1.d0-ecc(icase))
         a_kep(2) = ecc(icase)
         a_kep(3) = inc(icase)
         a_kep(4) = raan(icase)
         a_kep(5) = aop(icase)
         a_kep(6) = manom(icase)
         a_kep(7) = 0.d0

         call oe2car(a_kep, x0, xd0, 1.d-12)
         w = dsqrt(-0.5d0*(amue/a_kep(1)))
         call car2ks(x0, xd0, u, us, w)

c        Object's own epoch -> JD -> GMST -> station ECI positions
         epoch_cal(1) = epoch_yr(icase)
         epoch_cal(2) = epoch_mo(icase)
         epoch_cal(3) = epoch_dy(icase)
         epoch_cal(4) = epoch_hr(icase)
         epoch_cal(5) = epoch_mn(icase)
         epoch_cal(6) = epoch_sc(icase)
         call cal2jd(epoch_cal, jd_epoch)
         call gmst_deg(jd_epoch, gmst)

         do k = 1, nstat_max
            call geodetic_to_eci(stat_lat(k), stat_lon(k),
     &         stat_alt(k), R_earth, gmst, stat_pos(1,k))
         end do

         call fim_build_cartesian(x0, nstat_max, stat_pos, sigma,
     &                             Fc, ierr)
         call mat3_det(Fc, det_c)

         call fim_build_ks(u, w, nstat_max, stat_pos, sigma, Fk,
     &                      ierr)
         call fim_reduced_det_ks(Fk, reduced_det, rank_est, eigval,
     &                            ierr2)

c        Issue #4 cross-check: same scenario, analytical gradient
         call fim_build_ks_analytical(u, w, nstat_max, stat_pos,
     &      sigma, Fk_an, ierr)
         call fim_reduced_det_ks(Fk_an, reduced_det_an, rank_an,
     &      eigval_an, ierr2)

c        Raw 4x4 determinant -- the quantity the original
c        presentation reported and compared directly to det_c
         det_k_raw = eigval(1)*eigval(2)*eigval(3)*eigval(4)

         write(*,'(A,I5,A,A)') ' Case ', norad(icase), ': ',
     &      objname(icase)
         write(*,'(A,F10.2,A,F8.4,A,F6.2,A)') '   perigee alt=',
     &      perigee_alt(icase), ' km  e=', ecc(icase),
     &      '  GMST=', gmst, ' deg'
         write(*,'(A,4E14.6)') '   KS eigenvalues (desc) = ', eigval
         write(*,'(A,I1,A)') '   rank_est = ', rank_est,
     &      ' (4 = presentation''s implicit assumption, 3 = this'//
     &      ' repo''s prediction)'
         write(*,'(A,E16.8)') '   det(F_cartesian)          = ',
     &      det_c
         write(*,'(A,E16.8)') '   det(F_ks), raw 4x4        = ',
     &      det_k_raw
         write(*,'(A,E16.8)') '   reduced det(F_ks), FD grad= ',
     &      reduced_det
         write(*,'(A,E16.8,A,I1)') '   reduced det(F_ks), '//
     &      'analytical  = ', reduced_det_an, '  rank=', rank_an
         write(*,'(A,E14.6)') '   raw ratio (KS/Cartesian, as'//
     &      ' presentation compared it)   = ', det_k_raw/det_c
         write(*,'(A,E14.6)') '   corrected ratio (reduced'//
     &      ' KS/Cartesian)                = ', reduced_det/det_c
         write(*,*)
      end do

c     ============================================================
c     Issue #1: GEO case (object 28868, ANIK F1R), 1- and 2-station
c     sub-cases, matching the presentation's own "Case Studies:
c     Geosynchronous Object" slide. That slide does not publish exact
c     orbital elements (unlike the HEO table) -- uses a representative
c     near-GEO state instead (a=42164.14 km geostationary radius,
c     small e/i typical of a station-kept GEO satellite) and a nominal
c     epoch, since none was published either. Documented here rather
c     than left implicit. Station data (Kermit, Fairbanks) and
c     sigma=10km ARE from the slide's own tables.
c     ============================================================
      call geo_case_study()

      end

c     ============================================================
      subroutine geo_case_study()
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      double precision a_kep(7), x0(3), xd0(3)
      double precision u(4), us(4), w
      double precision stat_lat(2), stat_lon(2), stat_alt(2)
      double precision stat_pos(3,2), sigma(2)
      double precision Fc(3,3), Fk(4,4), Fk_an(4,4)
      double precision eigval(4), eigval_an(4)
      double precision det_c, det_k_raw, reduced_det, reduced_det_an
      integer rank_est, rank_an, rank_c, ierr, ierr2, k, nstat
      double precision epoch_cal(6), jd_epoch, gmst

      write(*,*) '======================================'
      write(*,*) ' GEO case: 28868 ANIK F1R'
      write(*,*) ' (representative near-GEO state -- see'
      write(*,*) '  header comment, no orbital elements'
      write(*,*) '  published for this case in the source)'
      write(*,*) '======================================'
      write(*,*)

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

      epoch_cal(1) = 2020.d0
      epoch_cal(2) = 1.d0
      epoch_cal(3) = 1.d0
      epoch_cal(4) = 0.d0
      epoch_cal(5) = 0.d0
      epoch_cal(6) = 0.d0
      call cal2jd(epoch_cal, jd_epoch)
      call gmst_deg(jd_epoch, gmst)

c     Kermit, Texas; Fairbanks, Alaska
      stat_lat(1) = 31.964d0
      stat_lon(1) = -103.233d0
      stat_alt(1) = 0.855d0
      sigma(1) = 10.d0

      stat_lat(2) = 65.1299d0
      stat_lon(2) = -147.471d0
      stat_alt(2) = 0.213d0
      sigma(2) = 10.d0

      do k = 1, 2
         call geodetic_to_eci(stat_lat(k), stat_lon(k), stat_alt(k),
     &      R_earth, gmst, stat_pos(1,k))
      end do

      do nstat = 1, 2
         call fim_build_cartesian(x0, nstat, stat_pos, sigma, Fc,
     &                             ierr)
         call fim_reduced_det_cartesian(Fc, det_c, rank_c, eigval,
     &                                    ierr)

         call fim_build_ks(u, w, nstat, stat_pos, sigma, Fk, ierr)
         call fim_reduced_det_ks(Fk, reduced_det, rank_est, eigval,
     &                            ierr2)
         det_k_raw = eigval(1)*eigval(2)*eigval(3)*eigval(4)

         call fim_build_ks_analytical(u, w, nstat, stat_pos, sigma,
     &      Fk_an, ierr)
         call fim_reduced_det_ks(Fk_an, reduced_det_an, rank_an,
     &      eigval_an, ierr2)

         write(*,'(A,I1,A)') ' --- ', nstat, ' station(s) ---'
         write(*,'(A,4E14.6)') '   KS eigenvalues (desc) = ', eigval
         write(*,'(A,I1,A,I1)') '   rank(F_cartesian)=', rank_c,
     &      '  rank(F_ks)=', rank_est
         write(*,'(A,E16.8,A,E16.8)') '   det(F_cartesian) [reduced'//
     &      ' to rank_c] = ', det_c, '   det(F_ks), raw 4x4 = ',
     &      det_k_raw
         write(*,'(A,E16.8,A,I1)') '   reduced det(F_ks), FD/'//
     &      'analytical = ', reduced_det, ' / ', rank_an
         if (dabs(reduced_det_an) .gt. 1.d-300) then
            write(*,'(A,E12.4)') '   FD vs analytical relative'//
     &         ' diff = ',
     &         dabs(reduced_det-reduced_det_an)/dabs(reduced_det_an)
         end if
         if (dabs(det_c) .gt. 1.d-300) then
            write(*,'(A,E14.6,A,E14.6)') '   raw ratio = ',
     &         det_k_raw/det_c, '   corrected ratio = ',
     &         reduced_det/det_c
         else
            write(*,*) '   det(F_cartesian) ~ 0 (rank-deficient at'//
     &         ' this station count) -- ratio undefined'
         end if
         write(*,*)
      end do

      return
      end
