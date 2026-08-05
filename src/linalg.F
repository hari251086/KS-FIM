c     ============================================================
c     linalg.F -- small dense linear-algebra helpers for KS-FIM
c
c     Fisher Information Matrices built from range measurements are
c     always real, symmetric, positive-semidefinite by construction
c     (sums of outer products grad*grad^T). Eigenvalues are the
c     quantity that actually answers "is this matrix rank-deficient"
c     -- a raw determinant can be a tiny nonzero float even when one
c     eigenvalue is genuinely zero to within numerical precision, so
c     KS-FIM's core question (is the 4x4 KS-space FIM structurally
c     rank<=3) is checked via eigenvalues, not determinant alone.
c
c     Author: Harishkumar Sellamuthu
c     ============================================================

c     ============================================================
      subroutine jacobi_eig(a, n, nmax, eigval, ierr)
c     Classic cyclic Jacobi eigenvalue algorithm for a real
c     symmetric n x n matrix (n <= nmax, nmax is the caller's fixed
c     array dimension -- this repo only ever calls with n=3 or n=4).
c     Destroys the input matrix `a` (used as scratch). Eigenvalues
c     returned UNSORTED in eigval(1:n) -- caller sorts if needed.
c
c     Input
c       a(nmax,nmax) : symmetric matrix, only a(1:n,1:n) used
c       n            : matrix dimension
c       nmax         : leading dimension of a
c     Output
c       eigval(n)    : eigenvalues (unsorted)
c       ierr         : 0 = converged, 1 = did not converge in the
c                      sweep budget (should not happen for n<=4)
c     ============================================================
      implicit double precision (a-h, o-z)
      integer n, nmax, ierr
      double precision a(nmax,nmax), eigval(n)

      parameter(maxsweep = 100, tol = 1.d-14)
      integer i, p, q, isweep
      double precision off, apq, app, aqq, phi, c, s, t, tau
      double precision aip, aiq

      ierr = 0

      do isweep = 1, maxsweep
c        Sum of squares of off-diagonal elements -- convergence test
         off = 0.d0
         do p = 1, n-1
            do q = p+1, n
               off = off + a(p,q)*a(p,q)
            end do
         end do
         if (off .lt. tol) goto 50

c        One cyclic sweep: zero every off-diagonal pair once
         do p = 1, n-1
            do q = p+1, n
               apq = a(p,q)
               if (dabs(apq) .lt. 1.d-300) goto 20

               app = a(p,p)
               aqq = a(q,q)
               phi = 0.5d0*(aqq-app)/apq
               if (phi .ge. 0.d0) then
                  t = 1.d0/(phi+dsqrt(1.d0+phi*phi))
               else
                  t = -1.d0/(-phi+dsqrt(1.d0+phi*phi))
               end if
               c = 1.d0/dsqrt(1.d0+t*t)
               s = t*c
               tau = s/(1.d0+c)

               a(p,p) = app - t*apq
               a(q,q) = aqq + t*apq
               a(p,q) = 0.d0
               a(q,p) = 0.d0

               do i = 1, n
                  if (i .eq. p .or. i .eq. q) goto 10
                  aip = a(i,p)
                  aiq = a(i,q)
                  a(i,p) = aip - s*(aiq + tau*aip)
                  a(p,i) = a(i,p)
                  a(i,q) = aiq + s*(aip - tau*aiq)
                  a(q,i) = a(i,q)
 10               continue
               end do
 20            continue
            end do
         end do
      end do

      ierr = 1

 50   continue
      do i = 1, n
         eigval(i) = a(i,i)
      end do

      return
      end

c     ============================================================
      subroutine mat3_det(a, det)
c     Determinant of a 3x3 matrix via cofactor expansion.
c     ============================================================
      implicit double precision (a-h, o-z)
      double precision a(3,3), det

      det = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))
     &    - a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1))
     &    + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

      return
      end

c     ============================================================
      subroutine eig_sort_desc(eigval, n)
c     In-place descending sort of a short eigenvalue list (n<=4,
c     insertion sort is plenty).
c     ============================================================
      implicit double precision (a-h, o-z)
      integer n
      double precision eigval(n)
      integer i, j
      double precision tmp

      do i = 2, n
         tmp = eigval(i)
         j = i - 1
         do while (j .ge. 1 .and. eigval(j) .lt. tmp)
            eigval(j+1) = eigval(j)
            j = j - 1
         end do
         eigval(j+1) = tmp
      end do

      return
      end
