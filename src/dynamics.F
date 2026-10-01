c     ============================================================
c     dynamics.F -- Phase 6 (issue #5 EKF study): exact two-body
c     Kepler flow map and its finite-difference state-transition
c     matrix (STM), in both Cartesian and KS representations.
c
c     The flow map itself is EXACT (Kepler's equation solved to
c     tolerance via oe2car, no numerical integration, no truncation
c     error) -- deliberately so, so the EKF study in
c     ksfim_ekf_study.F isolates representation-dependent
c     LINEARIZATION behavior (the STM) from propagator accuracy,
c     which is a separate, well-established, undisputed KS advantage
c     (see ALGORITHM.md section 9) this repo has never contested.
c
c     STMs are obtained by CENTRAL FINITE DIFFERENCE of the exact
c     flow map, through KSROP's own already-validated car2oe/oe2car/
c     car2ks/ks2car -- the same finite-difference philosophy already
c     used for range_grad_ks (fim.F), not a new hand-derived
c     Jacobian.
c     ============================================================

c     ============================================================
      subroutine kepler_flow_cartesian(r0, v0, dt, r1, v1)
c     Exact two-body flow: (r0,v0) propagated forward by dt seconds
c     via car2oe -> mean-anomaly advance -> oe2car.
c     ============================================================
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth
      double precision r0(3), v0(3), dt, r1(3), v1(3)
      double precision a_kep(7), n_mean, dM_deg

      call car2oe(r0, v0, a_kep)
      n_mean = dsqrt(amue/a_kep(1)**3)
      dM_deg = n_mean*dt*r2d
      a_kep(6) = a_kep(6) + dM_deg
      a_kep(6) = dmod(a_kep(6), 360.d0)
      if (a_kep(6) .lt. 0.d0) a_kep(6) = a_kep(6) + 360.d0

      call oe2car(a_kep, r1, v1, 1.d-12)

      return
      end

c     ============================================================
      subroutine stm_cartesian(r0, v0, dt, Phi)
c     6x6 STM d(r1,v1)/d(r0,v0) via central finite difference of
c     kepler_flow_cartesian.
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision r0(3), v0(3), dt, Phi(6,6)

      double precision rp(3), vp(3), rm(3), vm(3)
      double precision r1p(3), v1p(3), r1m(3), v1m(3)
      double precision x0(6), eps, step
      integer i, k

      x0(1) = r0(1)
      x0(2) = r0(2)
      x0(3) = r0(3)
      x0(4) = v0(1)
      x0(5) = v0(2)
      x0(6) = v0(3)

      eps = 1.d-6

      do k = 1, 6
         do i = 1, 3
            rp(i) = r0(i)
            vp(i) = v0(i)
            rm(i) = r0(i)
            vm(i) = v0(i)
         end do
         step = eps*dmax1(dabs(x0(k)), 1.d0)
         if (k .le. 3) then
            rp(k) = rp(k) + step
            rm(k) = rm(k) - step
         else
            vp(k-3) = vp(k-3) + step
            vm(k-3) = vm(k-3) - step
         end if

         call kepler_flow_cartesian(rp, vp, dt, r1p, v1p)
         call kepler_flow_cartesian(rm, vm, dt, r1m, v1m)

         do i = 1, 3
            Phi(i,k)   = (r1p(i)-r1m(i))/(2.d0*step)
            Phi(i+3,k) = (v1p(i)-v1m(i))/(2.d0*step)
         end do
      end do

      return
      end

c     ============================================================
      subroutine kepler_flow_ks(u0, us0, w, dt, u1, us1)
c     Same exact flow, expressed in KS coordinates at fixed w (valid
c     since w=sqrt(mu/2a) is exactly constant along an unperturbed
c     Kepler orbit -- the reference trajectory this study uses):
c     ks2car -> kepler_flow_cartesian -> car2ks.
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision u0(4), us0(4), w, dt, u1(4), us1(4)
      double precision r0(3), v0(3), r1(3), v1(3)

      call ks2car(u0, us0, r0, v0, w)
      call kepler_flow_cartesian(r0, v0, dt, r1, v1)
      call car2ks(r1, v1, u1, us1, w)

      return
      end

c     ============================================================
      subroutine stm_ks(u0, us0, w, dt, Phi)
c     8x8 STM d(u1,us1)/d(u0,us0) at fixed w, via central finite
c     difference of kepler_flow_ks.
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision u0(4), us0(4), w, dt, Phi(8,8)

      double precision up(4), usp(4), um(4), usm(4)
      double precision u1p(4), us1p(4), u1m(4), us1m(4)
      double precision x0(8), eps, step
      integer i, k

      do i = 1, 4
         x0(i) = u0(i)
         x0(i+4) = us0(i)
      end do

      eps = 1.d-6

      do k = 1, 8
         do i = 1, 4
            up(i) = u0(i)
            usp(i) = us0(i)
            um(i) = u0(i)
            usm(i) = us0(i)
         end do
         step = eps*dmax1(dabs(x0(k)), 1.d-3)
         if (k .le. 4) then
            up(k) = up(k) + step
            um(k) = um(k) - step
         else
            usp(k-4) = usp(k-4) + step
            usm(k-4) = usm(k-4) - step
         end if

         call kepler_flow_ks(up, usp, w, dt, u1p, us1p)
         call kepler_flow_ks(um, usm, w, dt, u1m, us1m)

         do i = 1, 4
            Phi(i,k)   = (u1p(i)-u1m(i))/(2.d0*step)
            Phi(i+4,k) = (us1p(i)-us1m(i))/(2.d0*step)
         end do
      end do

      return
      end

c     ============================================================
      subroutine kepler_delta_E(r0, v0, dt, dE)
c     Issue #5: the generalized eccentric anomaly increment (thesis eq
c     2.52/2.57) corresponding to physical-time step dt, for the exact
c     unperturbed Kepler flow this repo uses -- which coincides with
c     the CLASSICAL eccentric anomaly here (thesis section 2.3.3,
c     condition (i): U=0 and P=0). car2oe already returns E0 directly
c     (a_kep(7)); E1 comes from KSROP's own akesolve Kepler-equation
c     solver at the advanced mean anomaly -- no new numerics, reuses
c     exactly what oe2car/kepler_flow_cartesian already trust.
c     ============================================================
      implicit double precision (a-h, o-z)
      common /xy/ pi, d2r, r2d, amue, AU, R_Earth
      double precision r0(3), v0(3), dt, dE
      double precision a_kep(7), n_mean, am1, E0, E1

      call car2oe(r0, v0, a_kep)
      n_mean = dsqrt(amue/a_kep(1)**3)
      E0 = a_kep(7)*d2r
      am1 = a_kep(6)*d2r + n_mean*dt
      E1 = akesolve(a_kep(2), am1, 1.d-12)

      dE = E1 - E0

      return
      end

c     ============================================================
      subroutine stm_ks_analytical(u0, us0, w, dt, Phi)
c     Issue #5: EXACT (not finite-difference) 8x8 STM d(u1,us1)/d(u0,us0)
c     for the unperturbed KS harmonic oscillator -- the closed-form
c     solution in the generalized eccentric anomaly E (user's own
c     thesis, section 2.3.2, eq 2.57/2.58/2.63), with alpha*=beta*=0
c     since w, alpha, beta are all constant for unperturbed motion:
c       u(E)  =  u0*cos(dE/2) + 2*us0*sin(dE/2)
c       us(E) = -0.5*u0*sin(dE/2) +  us0*cos(dE/2)
c     dE comes from kepler_delta_E (Cartesian orbital elements + the
c     already-validated Kepler-equation solver), not from differencing
c     the flow map -- unlike stm_ks. Verified to agree with stm_ks to
c     central-finite-difference precision (test_dynamics.F, check Z_d).
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision u0(4), us0(4), w, dt, Phi(8,8)
      double precision r0(3), v0(3), dE, cE, sE
      integer i, j

      call ks2car(u0, us0, r0, v0, w)
      call kepler_delta_E(r0, v0, dt, dE)

      cE = dcos(0.5d0*dE)
      sE = dsin(0.5d0*dE)

      do i = 1, 8
         do j = 1, 8
            Phi(i,j) = 0.d0
         end do
      end do
      do i = 1, 4
         Phi(i,i)     = cE
         Phi(i,i+4)   = 2.d0*sE
         Phi(i+4,i)   = -0.5d0*sE
         Phi(i+4,i+4) = cE
      end do

      return
      end

c     ============================================================
      subroutine car2ks_jacobian(r0, v0, w, J)
c     8x6 finite-difference Jacobian d(u,us)/d(r,v) at fixed w --
c     used only to map an initial Cartesian covariance into KS space
c     for the EKF study's initial condition (NOT used anywhere in
c     the structural rank/ratio findings of Phases 1-5, and not
c     needed by them: J(1:4,4:6) [du/dv] is identically zero since
c     car2ks's u only depends on position, and J(1:4,1:3) [du/dr]
c     can be cross-checked against the already-validated dr=2*L3(u)*du
c     identity's Moore-Penrose pseudo-inverse, L3(u)^T/(2*|u|^2) --
c     see test_dynamics.F).
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision r0(3), v0(3), w, J(8,6)

      double precision rp(3), vp(3), rm(3), vm(3)
      double precision up(4), usp(4), um(4), usm(4)
      double precision x0(6), eps, step
      integer i, k

      x0(1) = r0(1)
      x0(2) = r0(2)
      x0(3) = r0(3)
      x0(4) = v0(1)
      x0(5) = v0(2)
      x0(6) = v0(3)

      eps = 1.d-6

      do k = 1, 6
         do i = 1, 3
            rp(i) = r0(i)
            vp(i) = v0(i)
            rm(i) = r0(i)
            vm(i) = v0(i)
         end do
         step = eps*dmax1(dabs(x0(k)), 1.d0)
         if (k .le. 3) then
            rp(k) = rp(k) + step
            rm(k) = rm(k) - step
         else
            vp(k-3) = vp(k-3) + step
            vm(k-3) = vm(k-3) - step
         end if

         call car2ks(rp, vp, up, usp, w)
         call car2ks(rm, vm, um, usm, w)

         do i = 1, 4
            J(i,k)   = (up(i)-um(i))/(2.d0*step)
            J(i+4,k) = (usp(i)-usm(i))/(2.d0*step)
         end do
      end do

      return
      end
