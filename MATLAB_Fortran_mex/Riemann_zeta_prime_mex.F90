#include "fintrf.h"

!-----------------------------------------------------------------------
! [f1, f2] = Riemann_zeta_prime_mex(s)
!
! MEX gateway for the Riemann zeta function and its derivative.
! Input s can be a scalar or array (real or complex, double precision).
! Output f1 = zeta(s), f2 = zeta'(s), both complex double, same shape as s.
!-----------------------------------------------------------------------
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      use Riemann_zeta_module
      implicit none

!     MEX interface
      mwPointer :: plhs(*), prhs(*)
      integer   :: nlhs, nrhs

!     MEX API functions
      mwPointer :: mxGetPr, mxGetPi, mxCreateDoubleMatrix
      mwSize    :: mxGetM, mxGetN
      integer   :: mxIsComplex

!     Local variables
      mwPointer :: pr_in, pi_in, pr_out, pi_out
      mwSize    :: m, n, numel
      integer, parameter :: qp_local = selected_real_kind(33, 4931)
      real(8), allocatable :: sr(:), si(:)
      real(8), allocatable :: f1r(:), f1i(:), f2r(:), f2i(:)
      complex(kind=qp_local) :: s_qp, f_qp(1:2)
      integer :: k

!     Check arguments
      if (nrhs /= 1) then
          call mexErrMsgTxt('One input required.')
      end if
      if (nlhs > 2) then
          call mexErrMsgTxt('At most two outputs.')
      end if

!     Get input dimensions
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      numel = m * n

!     Allocate work arrays
      allocate(sr(numel), si(numel))
      allocate(f1r(numel), f1i(numel), f2r(numel), f2i(numel))

!     Copy input real part
      pr_in = mxGetPr(prhs(1))
      call mxCopyPtrToReal8(pr_in, sr, numel)

!     Copy input imaginary part (or set to zero if real)
      if (mxIsComplex(prhs(1)) == 1) then
          pi_in = mxGetPi(prhs(1))
          call mxCopyPtrToReal8(pi_in, si, numel)
      else
          si = 0.0d0
      end if

!     Evaluate zeta(s) and zeta'(s) for each element
      do k = 1, numel
          s_qp = cmplx(sr(k), si(k), kind=qp_local)
          f_qp = Riemann_zeta_prime(s_qp)
          f1r(k) = real(real(f_qp(1)), kind=8)
          f1i(k) = real(aimag(f_qp(1)), kind=8)
          f2r(k) = real(real(f_qp(2)), kind=8)
          f2i(k) = real(aimag(f_qp(2)), kind=8)
      end do

!     First output: f1 = zeta(s)
      plhs(1) = mxCreateDoubleMatrix(m, n, 1)
      pr_out = mxGetPr(plhs(1))
      pi_out = mxGetPi(plhs(1))
      call mxCopyReal8ToPtr(f1r, pr_out, numel)
      call mxCopyReal8ToPtr(f1i, pi_out, numel)

!     Second output: f2 = zeta'(s) (only if requested)
      if (nlhs >= 2) then
          plhs(2) = mxCreateDoubleMatrix(m, n, 1)
          pr_out = mxGetPr(plhs(2))
          pi_out = mxGetPi(plhs(2))
          call mxCopyReal8ToPtr(f2r, pr_out, numel)
          call mxCopyReal8ToPtr(f2i, pi_out, numel)
      end if

      deallocate(sr, si, f1r, f1i, f2r, f2i)
      return
      end subroutine mexFunction
