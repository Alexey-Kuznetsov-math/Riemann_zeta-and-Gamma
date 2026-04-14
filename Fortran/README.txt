Author: Alexey Kuznetsov
York University, Toronto, Canada
Website: https://kuznetsovmath.ca/
Email: akuznets@yorku.ca

Created: 28-Nov-2025
Last updated: 14-April-2026

License: BSD 3-Clause (https://opensource.org/licenses/BSD-3-Clause)

======================================================================
PACKAGE CONTENTS
======================================================================

  gamma_zeta_module.f90  -- Fortran module with four public functions
  test.f90               -- test program (12 tests)
  Makefile               -- build script for gfortran on Linux/macOS
  README.txt             -- this file

======================================================================
OVERVIEW
======================================================================

The module "gamma_zeta_module.f90" provides four public functions:
Riemann_zeta, Riemann_zeta_prime, ln_gamma, and psi. All inputs must
be scalars of type complex(kind=qp) (these functions are not
vectorized).

IMPORTANT: The kind parameter qp is not exported by the module. Users
must define it in their own program as:

    integer, parameter :: qp = selected_real_kind(33, 4931)

======================================================================
QUICK START
======================================================================

    use gamma_zeta_module
    integer, parameter    :: qp = selected_real_kind(33, 4931)
    complex(kind=qp)      :: s, z, h(1:2)

    ! Riemann zeta function at s = 0.5 + 14.134i
    s = cmplx(0.5_qp, 14.134_qp)
    print *, Riemann_zeta(s)

    ! Gamma(1) = 1
    z = cmplx(1.0_qp, 0.0_qp, kind=qp)
    print *, exp(ln_gamma(z))

    ! zeta(s) and zeta'(s) simultaneously
    h = Riemann_zeta_prime(s)   ! h(1)=zeta(s), h(2)=zeta'(s)

======================================================================
FUNCTION DESCRIPTIONS
======================================================================

----------------------------------------------------------------------
Riemann_zeta(s)
----------------------------------------------------------------------
Computes the Riemann zeta function zeta(s) for a complex argument s.


Precision:
    |Im(s)| < 100    -- approximately 31 correct decimal digits
    |Im(s)| < 1000   -- approximately 30 correct decimal digits
    |Im(s)| < 10000  -- approximately 29 correct decimal digits

For larger values of |Im(s)|, accuracy continues to decrease at a
rate of approximately one decimal digit per increase in |Im(s)| by 
a factor of ten. More details on this precision loss can be found in [1].

Computational complexity in the half-plane Re(s) >= 0.5: O(sqrt(|Im(s)|))

The computational complexity is O(sqrt(|Im(s)|)) in the critical 
strip 0<=Re(s)<=1. In other regions of the complex plane it is either 
O(sqrt(|Im(s)|))  (if we are close to the critical strip) or O(1). 

Depending on the value of s, one of three methods is used:
    i)   an approximation zeta_30(s) developed in [1],
    ii)  the Euler-Maclaurin summation formula, or
    iii) direct summation zeta(s) = sum_{n=1}^{inf} n^{-s}, truncated
         at N chosen to achieve the target precision.

----------------------------------------------------------------------
Riemann_zeta_prime(s)
----------------------------------------------------------------------
Returns an array h of type complex(kind=qp), dimension(2), where:
    h(1) = zeta(s)
    h(2) = zeta'(s)   (the derivative of zeta with respect to s)

The precision and computational complexity are the same as for
Riemann_zeta(s). Values of zeta'(s) are computed by differentiating
the approximations described in items (i)-(iii) above.

----------------------------------------------------------------------
ln_gamma(z)
----------------------------------------------------------------------
Computes the logarithm of the Gamma function, log(Gamma(z)), in the
entire complex plane to a relative accuracy of 10^{-33}.

The algorithm is based on Binet's formula, approximating the function

    f(x) = (exp(-x)/x) * (1/(exp(x)-1) + 1/2 - 1/x)

by an exponential sum with 30 terms using techniques from [3]. See
[2] for further details.

----------------------------------------------------------------------
psi(z)
----------------------------------------------------------------------
Computes the digamma function psi(z) = Gamma'(z)/Gamma(z) in the
entire complex plane to a relative accuracy of 10^{-32}, except when
z is close to the negative real axis, where accuracy degrades.

psi(z) is obtained by differentiating the approximation developed for
ln_gamma(z).

======================================================================
BUILDING AND RUNNING
======================================================================

Requirements:
  - gfortran (GCC Fortran compiler), version 5.0 or later
  - Quadruple precision (REAL128) support.

To compile and run on Linux or macOS:

    make && make run

To clean up compiled files:

    make clean

======================================================================
TEST DESCRIPTIONS
======================================================================

The program test.f90 contains 12 tests:

  Test 1:  Gamma(k) for integer k = 1, ..., 10
  Test 2:  Gamma(k - 1/2) for integer k = 1, ..., 5
  Test 3:  Testing the functional equation log(Gamma(z+1))=log(Gamma(z))+log(z)
  Test 4:  psi(1) = -gamma (Euler-Mascheroni constant)
  Test 5:  Testing the functional equation psi(z+1)=psi(z)+1/z 
  Test 6:  zeta(s) at s = 2, 4, 6, 8, 10 
  Test 7:  zeta(s) at s = 0, -1, -3, -5, -7, -9 
  Test 8:  zeta(1/2 + i*t_k) at the first five non-trivial zeros
  Test 9:  zeta(1/2 + i*t_k) at five zeros near height t = 200
  Test 10: zeta(1/2 + i*t_k) at five zeros near height t = 10^6
  Test 11: Newton's method to find the first non-trivial zero of zeta(s)
  Test 12: Newton's method for a zero near 0.5 + 401.8i

======================================================================
REFERENCES
======================================================================

[1] A. Kuznetsov, "Simple and accurate approximations to the Riemann
    zeta function", 2025, https://arxiv.org/abs/2503.09519

[2] A. Kuznetsov, "Computing the Barnes G-function and the gamma
    function in the entire complex plane", Journal of Computational
    and Applied Mathematics, Vol. 411, 2022, 114270.
    https://doi.org/10.1016/j.cam.2022.114270

[3] A. Kuznetsov, A. Mohammadioroojeh, "Approximating functions on
    R^+ by exponential sums", 2025, https://arxiv.org/abs/2508.19095
