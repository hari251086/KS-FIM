c     ============================================================
c     test_dynamics.F -- Phase 6 (issue #5 EKF study): validates the
c     new exact-flow-map + finite-difference-STM machinery
c     (src/dynamics.F) before it's trusted for the EKF study itself.
c
c     Three independent checks, all against already-established
c     ground truth (never a new hand-derived formula taken on
c     faith):
c       Z_a: round-trip -- propagating forward by dt then back by
c            -dt through the EXACT Kepler flow must recover the
c            original state to near machine precision.
c       Z_b: specific orbital energy is conserved by the unperturbed
c            flow (v^2/2 - mu/r unchanged by kepler_flow_cartesian).
c       Z_c: car2ks_jacobian's du/dv block is identically zero
c            (car2ks's u depends only on position, never velocity --
c            visible directly in Subrouts.F), and its du/dr block is
c            a genuine right-inverse of dr=2*L3(u)*du (Phase 4):
c            since r(u(r)) = r identically along whatever specific
c            branch car2ks picks of the KS map's gauge circle,
c            2*L3(u)*(du/dr) must equal the 3x3 identity exactly --
c            NOT the same thing as matching the Moore-Penrose
c            pseudo-inverse L3(u)^T/(2*|u|^2), which assumes a
c            particular (minimum-norm) branch that car2ks's actual
c            algebraic construction has no reason to follow.
c     Issue #5 (analytical KS STM, in progress): stm_ks_analytical in
c     dynamics.F implements the closed-form harmonic-oscillator
c     solution (user's own thesis eq 2.57/2.58/2.63) and is verified
c     correct in the sense that matters physically -- ks2car of its
c     predicted (u1,us1) reproduces the EXACT two-body (r1,v1) to
c     Kepler-solver precision (checked independently in Python, not a
c     Fortran test here yet). It does NOT, however, numerically match
c     stm_ks (FD) component-for-component: car2ks's inverse transform
c     fixes one of the 4 u-components to exactly zero based on the
c     sign of x at that instant (Table 2.2 of the thesis; section
c     2.3.4 calls this out as a known shortcoming), a discontinuous
c     gauge choice re-selected fresh every time car2ks is called.
c     stm_ks's FD differences straight through that re-selection at
c     t1, so its (u1,us1) is a different point on the same KS fiber
c     than the smoothly-evolved harmonic-oscillator state -- not a bug
c     in either routine, a genuine representation mismatch. Closing
c     this needs composing stm_ks_analytical with an analytical
c     ks2car Jacobian and the existing (FD) car2ks_jacobian via the
c     chain rule, not attempted yet -- see issue #5.
c     ============================================================
      program test_dynamics
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth

      double precision a_kep(7), r0(3), v0(3)
      double precision r1(3), v1(3), r2(3), v2(3)
      double precision u0(4), us0(4), w
      double precision Lmat(4,4), J(8,6)
      double precision Jpred, rmag2, reldiff, maxdiff
      double precision en0, en1
      integer icase, i, k, npass, nfail

      call init_constants()
      npass = 0
      nfail = 0

      write(*,*) '=== test_dynamics ==='

      do icase = 1, 3

         if (icase .eq. 1) then
            a_kep(1) = 7000.d0
            a_kep(2) = 0.01d0
            a_kep(3) = 28.5d0
         else if (icase .eq. 2) then
            a_kep(1) = 24500.d0
            a_kep(2) = 0.72d0
            a_kep(3) = 5.78d0
         else
            a_kep(1) = 42164.d0
            a_kep(2) = 0.001d0
            a_kep(3) = 0.05d0
         end if
         a_kep(4) = 30.d0*dble(icase)
         a_kep(5) = 50.d0*dble(icase)
         a_kep(6) = 80.d0*dble(icase)
         a_kep(7) = 0.d0

         call oe2car(a_kep, r0, v0, 1.d-12)
         w = dsqrt(0.25d0*(amue/a_kep(1)))
         call car2ks(r0, v0, u0, us0, w)

c        Z_a: round-trip through the exact flow
         call kepler_flow_cartesian(r0, v0, 1234.5d0, r1, v1)
         call kepler_flow_cartesian(r1, v1, -1234.5d0, r2, v2)

         maxdiff = 0.d0
         do i = 1, 3
            reldiff = dabs(r2(i)-r0(i))/dabs(r0(i))
            if (reldiff .gt. maxdiff) maxdiff = reldiff
            reldiff = dabs(v2(i)-v0(i))/dabs(v0(i))
            if (reldiff .gt. maxdiff) maxdiff = reldiff
         end do

         write(*,'(A,I2,A,E12.4)') ' case ', icase,
     &      ': round-trip maxreldiff=', maxdiff
         if (maxdiff .lt. 1.d-9) then
            npass = npass + 1
         else
            write(*,*) ' [FAIL] round-trip, case ', icase
            nfail = nfail + 1
         end if

c        Z_b: energy conservation across a forward propagation
         en0 = 0.5d0*(v0(1)**2+v0(2)**2+v0(3)**2)
     &       - amue/dsqrt(r0(1)**2+r0(2)**2+r0(3)**2)
         en1 = 0.5d0*(v1(1)**2+v1(2)**2+v1(3)**2)
     &       - amue/dsqrt(r1(1)**2+r1(2)**2+r1(3)**2)
         reldiff = dabs(en1-en0)/dabs(en0)

         write(*,'(A,I2,A,E12.4)') ' case ', icase,
     &      ': energy conservation reldiff=', reldiff
         if (reldiff .lt. 1.d-9) then
            npass = npass + 1
         else
            write(*,*) ' [FAIL] energy conservation, case ', icase
            nfail = nfail + 1
         end if

c        Z_c: car2ks_jacobian cross-checked two ways
         call car2ks_jacobian(r0, v0, w, J)
         call ks_lmatrix(u0, Lmat)
         rmag2 = u0(1)*u0(1)+u0(2)*u0(2)+u0(3)*u0(3)+u0(4)*u0(4)

c        du/dv block must be exactly zero (u depends only on r)
         maxdiff = 0.d0
         do i = 1, 4
            do k = 4, 6
               reldiff = dabs(J(i,k))
               if (reldiff .gt. maxdiff) maxdiff = reldiff
            end do
         end do

         write(*,'(A,I2,A,E12.4)') ' case ', icase,
     &      ': car2ks_jacobian du/dv block maxval=', maxdiff
         if (maxdiff .lt. 1.d-6) then
            npass = npass + 1
         else
            write(*,*) ' [FAIL] car2ks_jacobian du/dv block, case ',
     &         icase
            nfail = nfail + 1
         end if

c        right-inverse identity: 2*L3(u)*(du/dr) = I3 exactly, since
c        r(u(r))=r identically regardless of which gauge branch
c        car2ks picks
         maxdiff = 0.d0
         do i = 1, 3
            do k = 1, 3
               sum_ik = 0.d0
               do m = 1, 4
                  sum_ik = sum_ik + 2.d0*Lmat(i,m)*J(m,k)
               end do
               if (i .eq. k) then
                  reldiff = dabs(sum_ik-1.d0)
               else
                  reldiff = dabs(sum_ik)
               end if
               if (reldiff .gt. maxdiff) maxdiff = reldiff
            end do
         end do

         write(*,'(A,I2,A,E12.4)') ' case ', icase,
     &      ': 2*L3(u)*(du/dr) vs I3 maxdiff=', maxdiff
         if (maxdiff .lt. 1.d-6) then
            npass = npass + 1
         else
            write(*,*) ' [FAIL] right-inverse identity, case ',
     &         icase
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
