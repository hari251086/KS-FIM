c     ============================================================
c     fim.F -- Fisher Information Matrix builders, Cartesian and
c     Kustaanheimo-Stiefel (KS) state representations
c
c     Motivation (see ALGORITHM.md for the full writeup): a 2021
c     COSPAR presentation claimed KS-space observability (measured as
c     det(FIM)) is "stronger" than Cartesian-space observability for
c     identical tracking scenarios. The claim is structurally
c     suspect: the physical position r=(x,y,z,0)=L(u)u depends on
c     only 3 true degrees of freedom, but the KS state u is 4
c     components -- range measurements (which only ever depend on
c     physical position, never on u directly) can therefore only
c     ever constrain a rank-<=3 subspace of the 4x4 KS-space FIM. A
c     4x4 FIM built this way should be rank-deficient (one
c     eigenvalue ~0) for ANY tracking geometry, not just some. This
c     module builds both FIMs so that claim can be checked directly
c     against KSROP's own real car2ks/ks2car transformation, rather
c     than argued from slides alone.
c
c     Author: Harishkumar Sellamuthu
c     ============================================================

c     ============================================================
      subroutine range_grad_cartesian(r, stat_r, grad)
c     Analytical gradient of range h(r) = |r - stat_r| w.r.t. the
c     3-vector r. This is exact (unit line-of-sight vector), no
c     finite-differencing needed -- used as the ground truth Fi
c     builder in Cartesian space.
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision r(3), stat_r(3), grad(3)
      double precision d(3), dmag
      integer i

      do i = 1, 3
         d(i) = r(i) - stat_r(i)
      end do
      dmag = dsqrt(d(1)*d(1)+d(2)*d(2)+d(3)*d(3))
      do i = 1, 3
         grad(i) = d(i)/dmag
      end do

      return
      end

c     ============================================================
      subroutine range_grad_ks(u, w, stat_r, grad, ierr)
c     Central-difference gradient of h(u) = |r(u) - stat_r| w.r.t.
c     the 4-vector u, where r(u) is KSROP's own ks2car forward map
c     (position only depends on u, not on the KS velocity us or the
c     oscillator frequency w -- both passed as dummies/zero to
c     ks2car since only x(1:3) is used here). Deliberately numerical
c     (not the analytical L(u)^-1 chain rule from the COSPAR slides)
c     so this doesn't silently inherit any error in that derivation
c     -- it independently re-derives the same gradient from KSROP's
c     actual, already-tested car2ks/ks2car code.
c
c     Output
c       grad(4) : d|r(u)-stat_r| / du_i, i=1..4
c       ierr    : 0 = ok, 1 = |u| too small for a stable KS map
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision u(4), w, stat_r(3), grad(4)
      integer ierr

      double precision us0(4), x(3), xd(3)
      double precision up(4), um(4), hp, hm
      double precision eps, ustep, umag
      double precision d(3)
      integer i, k

      ierr = 0
      umag = dsqrt(u(1)*u(1)+u(2)*u(2)+u(3)*u(3)+u(4)*u(4))
      if (umag .lt. 1.d-6) then
         ierr = 1
         return
      end if

      do i = 1, 4
         us0(i) = 0.d0
      end do

      eps = 1.d-6
      ustep = eps*umag

      do k = 1, 4
         do i = 1, 4
            up(i) = u(i)
            um(i) = u(i)
         end do
         up(k) = u(k) + ustep
         um(k) = u(k) - ustep

         call ks2car(up, us0, x, xd, w)
         do i = 1, 3
            d(i) = x(i) - stat_r(i)
         end do
         hp = dsqrt(d(1)*d(1)+d(2)*d(2)+d(3)*d(3))

         call ks2car(um, us0, x, xd, w)
         do i = 1, 3
            d(i) = x(i) - stat_r(i)
         end do
         hm = dsqrt(d(1)*d(1)+d(2)*d(2)+d(3)*d(3))

         grad(k) = (hp-hm)/(2.d0*ustep)
      end do

      return
      end

c     ============================================================
      subroutine ks_lmatrix(u, L)
c     Builds the 4x4 Levi-Civita/KS matrix L(u) from KSROP's own
c     convention (matches the presentation's equations.docx and
c     KSROP's internal ks2car construction: r=(x,y,z,0)=L(u)u).
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision u(4), L(4,4)

      L(1,1) = u(1)
      L(1,2) = -u(2)
      L(1,3) = -u(3)
      L(1,4) = u(4)
      L(2,1) = u(2)
      L(2,2) = u(1)
      L(2,3) = -u(4)
      L(2,4) = -u(3)
      L(3,1) = u(3)
      L(3,2) = u(4)
      L(3,3) = u(1)
      L(3,4) = u(2)
      L(4,1) = u(4)
      L(4,2) = -u(3)
      L(4,3) = u(2)
      L(4,4) = -u(1)

      return
      end

c     ============================================================
      subroutine range_grad_ks_analytical(u, w, stat_r, grad, ierr)
c     Issue #4 cross-check: analytical KS-space gradient, independent
c     of the finite-difference range_grad_ks above and independent of
c     the source presentation's own (possibly mis-transcribed/
c     ambiguous, see ALGORITHM.md) du/dr=L(u)^-1/2 chain-rule slide.
c     Uses the standard, textbook KS DIFFERENTIAL identity instead
c     (Stiefel & Scheifele): dr = 2*L(u)*du, i.e. dr_i/du_k = 2*L(u)_ik
c     for i=1,2,3 (r's 4th component is identically zero and doesn't
c     depend on the physical range h, so only L's first 3 rows enter
c     the chain rule). Then dh/du_k = sum_i (dh/dr_i)(dr_i/du_k)
c     = 2 * sum_{i=1}^3 grad_cartesian(i) * L(u)(i,k).
c
c     If this agrees with range_grad_ks's finite-difference result to
c     within numerical tolerance, Phase 1's findings are independently
c     confirmed via a second, analytically-derived method. If it
c     disagrees, that itself is a finding (see ALGORITHM.md).
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision u(4), w, stat_r(3), grad(4)
      integer ierr

      double precision us0(4), x(3), xd(3), d(3)
      double precision grad_cart(3), Lmat(4,4)
      double precision umag
      integer i, k

      ierr = 0
      umag = dsqrt(u(1)*u(1)+u(2)*u(2)+u(3)*u(3)+u(4)*u(4))
      if (umag .lt. 1.d-6) then
         ierr = 1
         return
      end if

      do i = 1, 4
         us0(i) = 0.d0
      end do
      call ks2car(u, us0, x, xd, w)

      do i = 1, 3
         d(i) = x(i) - stat_r(i)
      end do
      call range_grad_cartesian(x, stat_r, grad_cart)

      call ks_lmatrix(u, Lmat)

      do k = 1, 4
         grad(k) = 0.d0
         do i = 1, 3
            grad(k) = grad(k) + grad_cart(i)*Lmat(i,k)
         end do
         grad(k) = 2.d0*grad(k)
      end do

      return
      end

c     ============================================================
      subroutine fim_build_ks_analytical(u, w, nstat, stat_pos,
     &                                     sigma, F, ierr)
c     Same as fim_build_ks but using range_grad_ks_analytical instead
c     of the finite-difference gradient -- see that subroutine's
c     header for why (issue #4 cross-check).
c     ============================================================
      implicit double precision (a-h, o-z)
      integer nstat, ierr
      double precision u(4), w, stat_pos(3,nstat), sigma(nstat)
      double precision F(4,4)

      double precision grad(4), stat_r(3)
      integer i, j, k, ierr_g

      ierr = 0
      do i = 1, 4
         do j = 1, 4
            F(i,j) = 0.d0
         end do
      end do

      do k = 1, nstat
         stat_r(1) = stat_pos(1,k)
         stat_r(2) = stat_pos(2,k)
         stat_r(3) = stat_pos(3,k)
         call range_grad_ks_analytical(u, w, stat_r, grad, ierr_g)
         if (ierr_g .ne. 0) then
            ierr = 1
            goto 90
         end if
         do i = 1, 4
            do j = 1, 4
               F(i,j) = F(i,j) +
     &            grad(i)*grad(j)/(sigma(k)*sigma(k))
            end do
         end do
      end do

 90   continue
      return
      end

c     ============================================================
      subroutine fim_build_cartesian(r, nstat, stat_pos, sigma,
     &                                F, ierr)
c     3x3 FIM from nstat range measurements, Cartesian state.
c
c     Input
c       r(3)             : satellite ECI position (km)
c       nstat             : number of observing stations
c       stat_pos(3,nstat) : station ECI positions (km)
c       sigma(nstat)       : range-measurement 1-sigma (km)
c     Output
c       F(3,3) : Fisher Information Matrix
c       ierr   : 0 = ok
c     ============================================================
      implicit double precision (a-h, o-z)
      integer nstat, ierr
      double precision r(3), stat_pos(3,nstat), sigma(nstat)
      double precision F(3,3)

      double precision grad(3), stat_r(3)
      integer i, j, k

      ierr = 0
      do i = 1, 3
         do j = 1, 3
            F(i,j) = 0.d0
         end do
      end do

      do k = 1, nstat
         stat_r(1) = stat_pos(1,k)
         stat_r(2) = stat_pos(2,k)
         stat_r(3) = stat_pos(3,k)
         call range_grad_cartesian(r, stat_r, grad)
         do i = 1, 3
            do j = 1, 3
               F(i,j) = F(i,j) +
     &            grad(i)*grad(j)/(sigma(k)*sigma(k))
            end do
         end do
      end do

      return
      end

c     ============================================================
      subroutine fim_build_ks(u, w, nstat, stat_pos, sigma, F, ierr)
c     4x4 FIM from nstat range measurements, KS state.
c
c     Input
c       u(4), w            : KS state and oscillator frequency at
c                             the same physical instant as fim_build_
c                             cartesian's r (caller's responsibility
c                             to keep them consistent, e.g. via
c                             KSROP's car2ks on the same r)
c       nstat, stat_pos, sigma : as fim_build_cartesian
c     Output
c       F(4,4) : Fisher Information Matrix in KS-space
c       ierr   : 0 = ok, 1 = a station gradient failed (degenerate u)
c     ============================================================
      implicit double precision (a-h, o-z)
      integer nstat, ierr
      double precision u(4), w, stat_pos(3,nstat), sigma(nstat)
      double precision F(4,4)

      double precision grad(4), stat_r(3)
      integer i, j, k, ierr_g

      ierr = 0
      do i = 1, 4
         do j = 1, 4
            F(i,j) = 0.d0
         end do
      end do

      do k = 1, nstat
         stat_r(1) = stat_pos(1,k)
         stat_r(2) = stat_pos(2,k)
         stat_r(3) = stat_pos(3,k)
         call range_grad_ks(u, w, stat_r, grad, ierr_g)
         if (ierr_g .ne. 0) then
            ierr = 1
            goto 90
         end if
         do i = 1, 4
            do j = 1, 4
               F(i,j) = F(i,j) +
     &            grad(i)*grad(j)/(sigma(k)*sigma(k))
            end do
         end do
      end do

 90   continue
      return
      end

c     ============================================================
      subroutine fim_reduced_det_ks(F, reduced_det, rank_est,
     &                               eigval, ierr)
c     Eigen-decomposes a 4x4 KS-space FIM and reports:
c       - eigval(1:4), descending
c       - rank_est: count of eigenvalues above a relative threshold
c         (rel_tol * largest eigenvalue) -- expected to be 3, not 4,
c         if the rank-deficiency hypothesis is correct
c       - reduced_det: product of the top rank_est eigenvalues (the
c         "pseudo-determinant" restricted to the true physical
c         subspace actually constrained by the measurements -- NOT
c         always top-3: with fewer than 3 independent range
c         directions available, e.g. 2 stations, F_cartesian itself
c         is only rank-2, and the 3rd-largest KS eigenvalue is noise
c         around zero, not a genuine information direction). This is
c         the quantity that should be compared against
c         det(F_cartesian) restricted the same way, NOT det(F) of
c         the full 4x4 matrix
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision F(4,4), reduced_det, eigval(4)
      integer rank_est, ierr

      double precision a(4,4), rel_tol
      integer i, j, ierr_j

      parameter(rel_tol = 1.d-8)

      do i = 1, 4
         do j = 1, 4
            a(i,j) = F(i,j)
         end do
      end do

      call jacobi_eig(a, 4, 4, eigval, ierr_j)
      ierr = ierr_j
      call eig_sort_desc(eigval, 4)

      rank_est = 0
      do i = 1, 4
         if (eigval(i) .gt. rel_tol*eigval(1)) rank_est = rank_est+1
      end do

      reduced_det = 1.d0
      do i = 1, rank_est
         reduced_det = reduced_det*eigval(i)
      end do

      return
      end

c     ============================================================
      subroutine fim_reduced_det_cartesian(F, reduced_det, rank_est,
     &                                       eigval, ierr)
c     Same idea as fim_reduced_det_ks but for the 3x3 Cartesian FIM
c     -- with fewer than 3 independent range directions (e.g. 2
c     stations), F_cartesian itself is rank-deficient, and its raw
c     det() would be numerical noise around zero, not a meaningful
c     "observability" number. Used so KS-vs-Cartesian comparisons are
c     apples-to-apples at whatever the true rank actually is on both
c     sides, not just when both happen to be full rank.
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision F(3,3), reduced_det, eigval(3)
      integer rank_est, ierr

      double precision a(3,3), rel_tol
      integer i, j, ierr_j

      parameter(rel_tol = 1.d-8)

      do i = 1, 3
         do j = 1, 3
            a(i,j) = F(i,j)
         end do
      end do

      call jacobi_eig(a, 3, 3, eigval, ierr_j)
      ierr = ierr_j
      call eig_sort_desc(eigval, 3)

      rank_est = 0
      do i = 1, 3
         if (eigval(i) .gt. rel_tol*eigval(1)) rank_est = rank_est+1
      end do

      reduced_det = 1.d0
      do i = 1, rank_est
         reduced_det = reduced_det*eigval(i)
      end do

      return
      end
