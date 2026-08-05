c     ============================================================
c     test_fim_condition_number.F -- issue #5, first cut: is there
c     ANY rank-matched, non-dimensionalized comparison where KS
c     coordinates show a genuine observability advantage over
c     Cartesian, rather than just the trivial 64*r^3 scale factor
c     (Phase 4)?
c
c     Raw determinant fails as a metric because it isn't
c     dimensionless (Phase 1) or rank-matched (Phase 1 fix). The
c     natural NEXT candidate a nonlinear-filtering context would
c     actually care about is the FIM's condition number (ratio of
c     largest to smallest nonzero eigenvalue) -- it's dimensionless
c     by construction, and directly relevant to EKF/UKF numerical
c     conditioning and how "well-shaped" (vs. needle-like) the
c     information ellipsoid is, independent of its absolute size.
c
c     But this is already answered by the Phase 4 result before
c     writing a single test case: since EVERY nonzero eigenvalue of
c     F_ks equals EXACTLY 4*|u|^2 times the corresponding eigenvalue
c     of F_cartesian (same uniform scalar for all of them, at a given
c     instant), their RATIO to each other is invariant under the
c     KS map -- condition_number(F_ks) = condition_number(F_cartesian)
c     identically, not just approximately. Verified here across the
c     same 10 varied orbits/station geometries as
c     test_fim_ratio_formula, to confirm the analytical argument
c     rather than just assert it.
c
c     Disposition: this closes the "rank-matched, non-dimensionalized
c     LINEAR FIM comparison" reading of issue #5 -- there is no
c     conditioning advantage there either, exactly as the conformal-
c     map argument (ALGORITHM.md section 5) predicts. It does NOT
c     touch the deeper, genuinely different question of whether
c     nonlinear filter behavior (EKF linearization-point sensitivity
c     across a real estimation timeline, not a single-instant FIM)
c     shows a representation-dependent effect -- that would need an
c     actual filter-simulation study, out of scope here, and issue #5
c     stays open for that narrower remaining question.
c     ============================================================
      program test_fim_condition_number
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      double precision a_kep(7), x0(3), xd0(3)
      double precision u(4), us(4), w
      double precision stat_pos(3,4), sigma(4)
      double precision Fc(3,3), Fk(4,4)
      double precision eigval(4), eigval_c(3)
      double precision det_c, reduced_det, rmag
      double precision cond_c, cond_k, reldiff
      integer rank_est, rank_c, ierr, ierr2, npass, nfail
      integer icase, k

      call init_constants()
      npass = 0
      nfail = 0

      write(*,*) '=== test_fim_condition_number ==='

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
         w = dsqrt(0.25d0*(amue/a_kep(1)))
         call car2ks(x0, xd0, u, us, w)

         rmag = dsqrt(x0(1)**2+x0(2)**2+x0(3)**2)

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

         cond_c = eigval_c(1)/eigval_c(rank_c)
         cond_k = eigval(1)/eigval(rank_est)

         reldiff = dabs(cond_k-cond_c)/cond_c

         write(*,'(A,I2,A,F10.2,A,E16.8,A,E16.8,A,E10.2)')
     &      ' case ', icase, ': r=', rmag, ' km  cond_cart=',
     &      cond_c, '  cond_ks=', cond_k, '  reldiff=', reldiff

         if (rank_c .eq. 3 .and. rank_est .eq. 3 .and.
     &       reldiff .lt. 1.d-6) then
            npass = npass + 1
         else
            write(*,'(A,I2)') ' [FAIL] condition number invariance, '
     &         //'case ', icase
            nfail = nfail + 1
         end if
      end do

      write(*,'(A,I3,A,I3)') ' Total pass=', npass, ' fail=', nfail
      if (nfail .gt. 0) then
         write(*,*) 'TESTS FAILED'
         stop 1
      end if
      write(*,*) 'ALL TESTS PASSED'
      write(*,*)
      write(*,*) ' condition_number(F_ks) = condition_number('//
     &   'F_cartesian) EXACTLY, in every case -- KS offers no'
      write(*,*) ' conditioning advantage for the linear/FIM'//
     &   ' sense of observability either.'

      end
