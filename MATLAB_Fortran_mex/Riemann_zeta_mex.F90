#include "fintrf.h"

!-----------------------------------------------------------------------
! f = Riemann_zeta_mex(s)
!
! MEX gateway for the Riemann zeta function.
! Input s can be a scalar or array (real or complex, double precision).
! Output f is complex double with the same shape as s.
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
      real(8), allocatable :: sr(:), si(:), fr(:), fi(:)
      complex(kind=qp_local) :: s_qp, f_qp
      integer :: k

!     Check arguments
      if (nrhs /= 1) then
          call mexErrMsgTxt('Riemann_zeta_mex requires one input.')
      end if
      if (nlhs > 1) then
          call mexErrMsgTxt('Riemann_zeta_mex returns one output.')
      end if

!     Get input dimensions
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      numel = m * n

!     Allocate work arrays
      allocate(sr(numel), si(numel), fr(numel), fi(numel))

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

!     Evaluate zeta(s) for each element
      do k = 1, numel
          s_qp = cmplx(sr(k), si(k), kind=qp_local)
          f_qp = Riemann_zeta(s_qp)
          fr(k) = real(real(f_qp), kind=8)
          fi(k) = real(aimag(f_qp), kind=8)
      end do

!     Create complex output matrix and copy results
      plhs(1) = mxCreateDoubleMatrix(m, n, 1)
      pr_out = mxGetPr(plhs(1))
      pi_out = mxGetPi(plhs(1))
      call mxCopyReal8ToPtr(fr, pr_out, numel)
      call mxCopyReal8ToPtr(fi, pi_out, numel)

      deallocate(sr, si, fr, fi)
      return
      end subroutine mexFunction
