c     ============================================================
c     test_timeconv.F -- issue #3: sanity checks for gmst_deg and
c     geodetic_to_eci.
c     ============================================================
      program test_timeconv
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      double precision cal(6), jd, gmst, gmst2, jd2
      double precision r_eci(3), rmag
      integer npass, nfail

      call init_constants()
      npass = 0
      nfail = 0

      write(*,*) '=== test_timeconv ==='

c     Z1: well-known reference value -- GMST at J2000.0 epoch
c     (2000-01-01 12:00:00 UTC) is 280.46061837508 degrees (a widely
c     published constant, e.g. Vallado).
      cal(1) = 2000.d0
      cal(2) = 1.d0
      cal(3) = 1.d0
      cal(4) = 12.d0
      cal(5) = 0.d0
      cal(6) = 0.d0
      call cal2jd(cal, jd)
      call gmst_deg(jd, gmst)

      write(*,'(A,F16.10)') '  GMST at J2000.0 = ', gmst
      if (dabs(gmst - 280.46061837508d0) .lt. 1.d-4) then
         write(*,*) '[PASS] Z1 GMST matches published J2000.0 value'
         npass = npass + 1
      else
         write(*,*) '[FAIL] Z1 GMST matches published J2000.0 value'
         nfail = nfail + 1
      end if

c     Z2: GMST rate -- one sidereal day (~23h56m4s later) should
c     advance GMST by exactly 360 degrees; one solar (calendar) day
c     later should advance it by ~360.9856 degrees (Earth's extra
c     ~0.9856 deg/day orbital motion around the Sun).
      cal(3) = 2.d0
      call cal2jd(cal, jd2)
      call gmst_deg(jd2, gmst2)

      dgmst = gmst2 - gmst
      if (dgmst .lt. 0.d0) dgmst = dgmst + 360.d0

      write(*,'(A,F14.10)') '  GMST advance over 1 solar day '//
     &   '(mod 360) = ', dgmst
c     360.9856 deg/day mod 360 = 0.9856 deg/day (the ~4-minute-a-day
c     sidereal/solar mismatch, wrapped)
      if (dabs(dgmst - 0.9856d0) .lt. 1.d-2) then
         write(*,*) '[PASS] Z2 GMST solar-day rate ~360.9856 deg/day'
         npass = npass + 1
      else
         write(*,*) '[FAIL] Z2 GMST solar-day rate ~360.9856 deg/day'
         nfail = nfail + 1
      end if

c     Z3: geodetic_to_eci round-trip -- the ECI vector's magnitude
c     must equal R_earth+alt regardless of GMST (pure rotation).
      call geodetic_to_eci(28.3d0, -16.509d0, 2.39d0, R_earth,
     &   gmst, r_eci)
      rmag = dsqrt(r_eci(1)**2+r_eci(2)**2+r_eci(3)**2)

      if (dabs(rmag - (R_earth+2.39d0)) .lt. 1.d-8) then
         write(*,*) '[PASS] Z3 geodetic_to_eci preserves radius'
         npass = npass + 1
      else
         write(*,*) '[FAIL] Z3 geodetic_to_eci preserves radius'
         nfail = nfail + 1
      end if

c     Z4: at GMST=0, ECI and ECEF coincide -- a station at the
c     equator/prime-meridian should land exactly on the +x axis.
      call geodetic_to_eci(0.d0, 0.d0, 0.d0, R_earth, 0.d0, r_eci)
      if (dabs(r_eci(1)-R_earth) .lt. 1.d-8 .and.
     &    dabs(r_eci(2)) .lt. 1.d-8 .and.
     &    dabs(r_eci(3)) .lt. 1.d-8) then
         write(*,*) '[PASS] Z4 GMST=0 equator/meridian on +x axis'
         npass = npass + 1
      else
         write(*,*) '[FAIL] Z4 GMST=0 equator/meridian on +x axis'
         nfail = nfail + 1
      end if

      write(*,'(A,I2,A,I2)') ' Total pass=', npass, ' fail=', nfail
      if (nfail .gt. 0) then
         write(*,*) 'TESTS FAILED'
         stop 1
      end if
      write(*,*) 'ALL TESTS PASSED'

      end
