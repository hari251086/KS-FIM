c     ============================================================
c     ksfim_case_study.F -- reproduces the 2021 COSPAR presentation's
c     4 highly-elliptical-orbit case studies (objects 35497, 37151,
c     39615, 42928 -- the same orbital elements published in
c     "Regularized Orbit Observability for Resident Space Objects",
c     PEDAS.1-0023-21) and reports THREE numbers side by side per
c     case: the raw 4x4 det(F_ks) (what the original presentation
c     reported, structurally rank-deficient per this repo's ALGORITHM.
c     md), the rank-corrected reduced determinant (product of the top
c     3 eigenvalues only), and det(F_cartesian) for the same scenario
c     -- so the "KS observability advantage" claim can be read
c     honestly instead of via the inflated raw 4x4 number.
c
c     Ground stations and orbital elements transcribed directly from
c     the slide deck (E:\Research\Conferences\c2021 - COSPAR\Paper
c     and Presentation\DROPBOX-31122020\...FULLPRESENTATIONSLIDES.pdf,
c     "Case Studies" slides) -- not re-derived, so results are
c     directly comparable to the original numbers.
c
c     Simplification vs. the original: stations are placed at their
c     published geodetic lat/lon on a spherical Earth, in the SAME
c     inertial frame as the satellite state, at a nominal zero Earth-
c     rotation offset (no GMST epoch alignment). This does not affect
c     the structural question under test (rank-deficiency and Jacobian
c     -scaling are properties of the coordinate MAP, not of which
c     specific stations are used or their exact sidereal alignment) --
c     see ALGORITHM.md section 2 for why this simplification is valid
c     for this repo's purpose.
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
      double precision Fc(3,3), Fk(4,4), eigval(4)
      double precision det_c, det_k_raw, reduced_det
      integer rank_est, ierr, ierr2, icase, k, nobj

      parameter(nobj = 4)
      integer norad(nobj)
      double precision perigee_alt(nobj), ecc(nobj), inc(nobj)
      double precision raan(nobj), aop(nobj), manom(nobj)
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

      norad(2) = 37151
      objname(2) = 'Long March 3B Upper Stage'
      perigee_alt(2) = 145.09d0
      ecc(2) = 0.4099d0
      inc(2) = 24.8672d0
      raan(2) = 204.061d0
      aop(2) = 163.8209d0
      manom(2) = 341.4979d0

      norad(3) = 39615
      objname(3) = 'Proton-M Briz-M Aux Tank'
      perigee_alt(3) = 143.5d0
      ecc(3) = 0.3367d0
      inc(3) = 48.4571d0
      raan(3) = 118.9947d0
      aop(3) = 333.8816d0
      manom(3) = 12.4622d0

      norad(4) = 42928
      objname(4) = 'IRNSS-1H/PSLV-PS4'
      perigee_alt(4) = 156.33d0
      ecc(4) = 0.3095d0
      inc(4) = 19.1455d0
      raan(4) = 315.4336d0
      aop(4) = 125.4322d0
      manom(4) = 330.6148d0

c     Kwajalein, Ascension Island, SHAR, Tenerife
      stat_lat(1) = 8.717d0
      stat_lon(1) = 167.733d0
      stat_alt(1) = 0.01d0
      sigma(1) = 0.01d0

      stat_lat(2) = 7.972d0
      stat_lon(2) = -14.401d0
      stat_alt(2) = 0.01d0
      sigma(2) = 0.01d0

      stat_lat(3) = 13.719d0
      stat_lon(3) = 80.23d0
      stat_alt(3) = 0.01d0
      sigma(3) = 0.01d0

      stat_lat(4) = 28.3d0
      stat_lon(4) = -16.509d0
      stat_alt(4) = 2.39d0
      sigma(4) = 0.01d0

      do k = 1, nstat_max
         rs = R_earth + stat_alt(k)
         stat_pos(1,k) = rs*dcos(stat_lat(k)*d2r)*dcos(stat_lon(k)*
     &                    d2r)
         stat_pos(2,k) = rs*dcos(stat_lat(k)*d2r)*dsin(stat_lon(k)*
     &                    d2r)
         stat_pos(3,k) = rs*dsin(stat_lat(k)*d2r)
      end do

      write(*,*) '======================================'
      write(*,*) ' KS-FIM case-study reproduction'
      write(*,*) ' (COSPAR 2021 PEDAS.1-0023-21 objects)'
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

         call fim_build_cartesian(x0, nstat_max, stat_pos, sigma,
     &                             Fc, ierr)
         call mat3_det(Fc, det_c)

         call fim_build_ks(u, w, nstat_max, stat_pos, sigma, Fk,
     &                      ierr)
         call fim_reduced_det_ks(Fk, reduced_det, rank_est, eigval,
     &                            ierr2)

c        Raw 4x4 determinant -- the quantity the original
c        presentation reported and compared directly to det_c
         det_k_raw = eigval(1)*eigval(2)*eigval(3)*eigval(4)

         write(*,'(A,I5,A,A)') ' Case ', norad(icase), ': ',
     &      objname(icase)
         write(*,'(A,F10.2,A,F8.4)') '   perigee alt=',
     &      perigee_alt(icase), ' km  e=', ecc(icase)
         write(*,'(A,4E14.6)') '   KS eigenvalues (desc) = ', eigval
         write(*,'(A,I1,A)') '   rank_est = ', rank_est,
     &      ' (4 = presentation''s implicit assumption, 3 = this'//
     &      ' repo''s prediction)'
         write(*,'(A,E16.8)') '   det(F_cartesian)          = ',
     &      det_c
         write(*,'(A,E16.8)') '   det(F_ks), raw 4x4        = ',
     &      det_k_raw
         write(*,'(A,E16.8)') '   reduced det(F_ks), top-3  = ',
     &      reduced_det
         write(*,'(A,E14.6)') '   raw ratio (KS/Cartesian, as'//
     &      ' presentation compared it)   = ', det_k_raw/det_c
         write(*,'(A,E14.6)') '   corrected ratio (reduced'//
     &      ' KS/Cartesian)                = ', reduced_det/det_c
         write(*,*)
      end do

      end
