c     ============================================================
c     test_fim_ratio_formula.F -- closes out the analytical picture:
c     the rank-corrected KS/Cartesian ratio isn't just "constant
c     across station geometry" (test_fim_reduced_consistency) -- it
c     has an EXACT closed form, derivable from two standard KS
c     identities:
c
c       (1) L(u) is |u|-scaled orthogonal: L(u)^T L(u) = |u|^2 I_4
c           (hence also L(u)L(u)^T = |u|^2 I_4, same property for a
c           square matrix), so its 3-row physical submatrix L_3(u)
c           satisfies L_3(u)L_3(u)^T = |u|^2 I_3 -- L_3(u)'s 3
c           singular values are each exactly |u|.
c       (2) The standard KS regularizing identity |u|^2 = r (the
c           satellite's instantaneous radius).
c
c     Since F_ks = J^T F_cartesian J with J = dr/du = 2*L_3(u) (the
c     dr=2L(u)du identity used in range_grad_ks_analytical), F_ks's
c     nonzero eigenvalues equal 4*|u|^2 * (F_cartesian's eigenvalues)
c     exactly -- so:
c
c         reduced_det_ks / det_cartesian = (4*|u|^2)^3 = 64*r^3
c
c     exactly, independent of station geometry, count (>=3), or
c     sigma -- a single number determined entirely by the satellite's
c     instantaneous distance from the primary. Verified here across
c     many different orbits/anomalies/station geometries at once.
c     ============================================================
      program test_fim_ratio_formula
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      double precision a_kep(7), x0(3), xd0(3)
      double precision u(4), us(4), w
      double precision stat_pos(3,4), sigma(4)
      double precision Fc(3,3), Fk(4,4)
      double precision eigval(4), eigval_c(3)
      double precision det_c, reduced_det, rmag, umag2
      double precision predicted_ratio, actual_ratio, reldiff
      integer rank_est, rank_c, ierr, ierr2, npass, nfail
      integer icase, k

      call init_constants()
      npass = 0
      nfail = 0

      write(*,*) '=== test_fim_ratio_formula ==='

c     10 varied cases: different (a,e,i,RAAN,AOP,M) AND different
c     station geometries together, so this genuinely stresses both
c     the r-dependence and the geometry-independence at once.
      do icase = 1, 10

         a_kep(1) = 8000.d0 + 2000.d0*dble(icase)
         a_kep(2) = 0.05d0 + 0.08d0*dble(icase)
         if (a_kep(2) .gt. 0.85d0) a_kep(2) = 0.85d0
         a_kep(3) = 5.d0 + 7.d0*dble(icase)
         a_kep(4) = 20.d0*dble(icase)
         a_kep(5) = 15.d0*dble(icase)
         a_kep(6) = 37.d0*dble(icase)
         a_kep(7) = 0.d0

         call oe2car(a_kep, x0, xd0, 1.d-12)
         w = dsqrt(-0.5d0*(amue/a_kep(1)))
         call car2ks(x0, xd0, u, us, w)

         rmag = dsqrt(x0(1)**2+x0(2)**2+x0(3)**2)
         umag2 = u(1)*u(1)+u(2)*u(2)+u(3)*u(3)+u(4)*u(4)

c        Z_a: |u|^2 = r (the KS regularizing identity itself)
         if (dabs(umag2-rmag)/rmag .lt. 1.d-10) then
            npass = npass + 1
         else
            write(*,'(A,I2,A,E12.4)') ' [FAIL] |u|^2=r, case ',
     &         icase, ' reldiff=', dabs(umag2-rmag)/rmag
            nfail = nfail + 1
         end if

         do k = 1, 4
            ang = dble(k)*1.1d0 + 0.7d0*dble(icase)
            stat_pos(1,k) = R_earth*dcos(ang)*dcos(0.2d0*dble(icase))
            stat_pos(2,k) = R_earth*dsin(ang)
            stat_pos(3,k) = R_earth*dsin(0.2d0*dble(icase)+
     &                       0.5d0*dble(k))
            sigma(k) = 1.d0 + 0.3d0*dble(k)*dble(icase)
         end do

         call fim_build_cartesian(x0, 4, stat_pos, sigma, Fc, ierr)
         call fim_reduced_det_cartesian(Fc, det_c, rank_c, eigval_c,
     &                                    ierr)

         call fim_build_ks(u, w, 4, stat_pos, sigma, Fk, ierr)
         call fim_reduced_det_ks(Fk, reduced_det, rank_est, eigval,
     &                            ierr2)

         actual_ratio = reduced_det/det_c
         predicted_ratio = 64.d0*rmag**3

         reldiff = dabs(actual_ratio-predicted_ratio)/predicted_ratio

         write(*,'(A,I2,A,F10.2,A,E16.8,A,E16.8,A,E10.2)')
     &      ' case ', icase, ': r=', rmag, ' km  actual=',
     &      actual_ratio, '  predicted 64*r^3=', predicted_ratio,
     &      '  reldiff=', reldiff

c        Z_b: exact closed-form match
         if (rank_c .eq. 3 .and. rank_est .eq. 3 .and.
     &       reldiff .lt. 1.d-4) then
            npass = npass + 1
         else
            write(*,'(A,I2)') ' [FAIL] ratio=64*r^3, case ', icase
            nfail = nfail + 1
         end if
      end do

      write(*,'(A,I3,A,I3)') ' Total pass=', npass, ' fail=', nfail
      if (nfail .gt. 0) then
         write(*,*) 'TESTS FAILED'
         stop 1
      end if
      write(*,*) 'ALL TESTS PASSED'

      end
