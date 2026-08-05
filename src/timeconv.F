c     ============================================================
c     timeconv.F -- GMST and geodetic-to-ECI conversion for KS-FIM
c
c     Issue #3: the original Phase 1 case study placed observing
c     stations at their published geodetic coordinates directly in
c     the satellite's inertial frame, with no Earth-rotation
c     (sidereal time) alignment to the object's actual epoch. That
c     was a deliberate simplification (doesn't affect the structural
c     rank-deficiency/Jacobian-scaling findings, which are properties
c     of the KS transformation itself, not of tracking geometry) but
c     means the case study wasn't placing stations where they'd
c     actually be at each object's real observation epoch. This
c     module fixes that with a standard GMST (IAU 1982) formula and
c     the standard ECEF->ECI rotation by GMST about the z-axis.
c
c     Author: Harishkumar Sellamuthu
c     ============================================================

c     ============================================================
      subroutine gmst_deg(jd_ut1, gmst)
c     Greenwich Mean Sidereal Time, IAU 1982 formula (e.g. Vallado,
c     "Fundamentals of Astrodynamics and Applications"). UT1~UTC is
c     an adequate approximation for this repo's purposes (station
c     geometry, not precision orbit determination).
c
c     Input
c       jd_ut1 : Julian date (UT1 ~ UTC)
c     Output
c       gmst   : Greenwich Mean Sidereal Time, degrees, in [0,360)
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision jd_ut1, gmst
      double precision t, gmst_sec

      t = (jd_ut1 - 2451545.0d0) / 36525.0d0

      gmst_sec = 67310.54841d0
     &   + (876600.d0*3600.d0 + 8640184.812866d0)*t
     &   + 0.093104d0*t*t
     &   - 6.2d-6*t*t*t

c     Seconds of sidereal time -> degrees (240 sidereal sec = 1 deg)
      gmst = gmst_sec/240.d0
      gmst = dmod(gmst, 360.d0)
      if (gmst .lt. 0.d0) gmst = gmst + 360.d0

      return
      end

c     ============================================================
      subroutine geodetic_to_eci(lat_deg, lon_deg, alt_km, r_earth,
     &                             gmst, r_eci)
c     Station position, spherical-Earth geodetic (lat/lon/altitude)
c     to ECI, via ECEF then rotating by GMST about the z-axis.
c     Spherical (not oblate/WGS84) Earth model, consistent with the
c     rest of this repo's simplified station-geometry treatment --
c     the structural findings under test don't depend on Earth's
c     oblateness at the station-position level.
c
c     Input
c       lat_deg, lon_deg : geodetic latitude/longitude, degrees
c       alt_km           : altitude above the spherical Earth, km
c       r_earth          : Earth radius, km (from KSROP's common /xy/)
c       gmst             : Greenwich Mean Sidereal Time, degrees
c     Output
c       r_eci(3) : station position in ECI, km
c     ============================================================
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth_dummy
      double precision lat_deg, lon_deg, alt_km, r_earth, gmst
      double precision r_eci(3)

      double precision rs, x_ecef, y_ecef, z_ecef, gmst_rad

      rs = r_earth + alt_km
      x_ecef = rs*dcos(lat_deg*d2r)*dcos(lon_deg*d2r)
      y_ecef = rs*dcos(lat_deg*d2r)*dsin(lon_deg*d2r)
      z_ecef = rs*dsin(lat_deg*d2r)

      gmst_rad = gmst*d2r
      r_eci(1) = x_ecef*dcos(gmst_rad) - y_ecef*dsin(gmst_rad)
      r_eci(2) = x_ecef*dsin(gmst_rad) + y_ecef*dcos(gmst_rad)
      r_eci(3) = z_ecef

      return
      end
