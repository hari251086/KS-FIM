c     ============================================================
c     kalman.F -- Phase 6 (issue #5 EKF study): sequential scalar
c     Kalman covariance update.
c
c     Range measurements from different stations are independent
c     given the state, so a batch update can equivalently be done as
c     nstat sequential SCALAR updates (one station at a time) -- the
c     standard sequential/Bierman-style trick. This needs only a
c     scalar division, never an NxN matrix inverse, keeping this
c     routine generic (works for the 6-dim Cartesian state and the
c     8-dim KS state with the same code) and low-risk.
c     ============================================================

c     ============================================================
      subroutine kalman_scalar_update(P, n, H, sigma)
c     In-place scalar measurement update: P <- P - (P H^T)(H P) /
c     (H P H^T + sigma^2). Standard covariance-form Kalman update
c     for one scalar measurement with noise variance sigma^2 and
c     linearized measurement row H(1:n) (dh/dstate).
c     ============================================================
      implicit double precision (a-h, o-z)
      integer n
      double precision P(n,n), H(n), sigma

      double precision PH(n), s
      integer i, j

      do i = 1, n
         PH(i) = 0.d0
         do j = 1, n
            PH(i) = PH(i) + P(i,j)*H(j)
         end do
      end do

      s = 0.d0
      do i = 1, n
         s = s + H(i)*PH(i)
      end do
      s = s + sigma*sigma

      do j = 1, n
         do i = 1, n
            P(i,j) = P(i,j) - PH(i)*PH(j)/s
         end do
      end do

      return
      end
