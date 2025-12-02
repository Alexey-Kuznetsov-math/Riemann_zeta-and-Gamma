#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)
  use Riemann_zeta_module          
  implicit none

  integer*4 nlhs, nrhs
  mwPointer plhs(*), prhs(*)

  ! Fortran interface to MEX API functions
  integer*4 mxIsComplex, mxIsDouble
  mwPointer mxGetPr, mxGetPi, mxCreateDoubleMatrix
  mwSize    mxGetNumberOfElements, mxGetM, mxGetN

  ! MATLAB array handles and sizes
  mwPointer :: s_in_mx, f_out_mx
  mwPointer :: pr_re, pr_im, po_re, po_im
  mwSize   :: n_el, mrows, ncols

  ! Work arrays
  double precision, allocatable :: s_re(:), s_im(:), f_re(:), f_im(:)
  integer :: i
  complex(kind=16) :: s16, f16
  logical :: is_complex

  !------------------------------------------------------------
  ! Check number of inputs/outputs
  !------------------------------------------------------------
  if (nrhs .ne. 1) then
     call mexErrMsgTxt('One input (s) required.')
  end if
  if (nlhs .ne. 1) then
     call mexErrMsgTxt('One output (z) required.')
  end if

  s_in_mx = prhs(1)

  !------------------------------------------------------------
  ! We accept double input (real or complex)
  !------------------------------------------------------------
  if (mxIsDouble(s_in_mx) .eq. 0) then
     call mexErrMsgTxt('Input must be double (real or complex).')
  end if

  is_complex = (mxIsComplex(s_in_mx) .ne. 0)

  !------------------------------------------------------------
  ! Get sizes and data pointers
  !------------------------------------------------------------
  n_el  = mxGetNumberOfElements(s_in_mx)
  mrows = mxGetM(s_in_mx)
  ncols = mxGetN(s_in_mx)

  pr_re = mxGetPr(s_in_mx)
  if (is_complex) then
     pr_im = mxGetPi(s_in_mx)
  else
     pr_im = 0
  end if

  allocate(s_re(n_el), s_im(n_el), f_re(n_el), f_im(n_el))

  ! Copy real part
  call mxCopyPtrToReal8(pr_re, s_re, n_el)

  ! Copy imag part if complex; otherwise set to zero
  if (is_complex) then
     call mxCopyPtrToReal8(pr_im, s_im, n_el)
  else
     s_im = 0.0d0
  end if

  !------------------------------------------------------------
  ! Create complex output array
  !------------------------------------------------------------
  f_out_mx = mxCreateDoubleMatrix(mrows, ncols, .true.)  ! .true. => complex
  po_re    = mxGetPr(f_out_mx)
  po_im    = mxGetPi(f_out_mx)

  !------------------------------------------------------------
  ! Main loop: call Riemann_zeta in quad precision
  !------------------------------------------------------------
  do i = 1, n_el
     s16 = cmplx( real(s_re(i),kind=16), real(s_im(i),kind=16), kind=16 )
     f16 = Riemann_zeta(s16)
     f_re(i) = real(f16, kind=8)
     f_im(i) = aimag(f16)
  end do

  ! Copy results back to MATLAB
  call mxCopyReal8ToPtr(f_re, po_re, n_el)
  call mxCopyReal8ToPtr(f_im, po_im, n_el)

  plhs(1) = f_out_mx

  deallocate(s_re, s_im, f_re, f_im)

end subroutine mexFunction

