Author: Alexey Kuznetsov
York University, Toronto, Canada
Website: https://kuznetsovmath.ca/
Email: akuznets@yorku.ca

Created: 29-Nov-2025
Last updated: 2-April-2026

License: BSD 3-Clause (https://opensource.org/licenses/BSD-3-Clause)

======================================================================
PACKAGE CONTENTS
======================================================================

  Riemann_zeta.py   -- Python implementation of the Riemann zeta function
  ln_gamma.py       -- Python implementation of the log-Gamma function
  test.py           -- correctness tests and plots
  README.txt        -- this file

======================================================================
OVERVIEW
======================================================================

This package provides two Python/NumPy functions, Riemann_zeta
and ln_gamma. 

Riemann_zeta(s) computes the Riemann zeta function for scalar,
list, or NumPy array inputs. The implementation is fully vectorized.

Accuracy:
  |Im(s)| < 100      -- approximately 13 correct decimal digits
  |Im(s)| < 1,000    -- approximately 12 correct decimal digits
  |Im(s)| < 10,000   -- approximately 11 correct decimal digits

For larger values of |Im(s)|, accuracy continues to decrease at a
rate of approximately one decimal digit per increase in |Im(s)| by 
a factor of ten. More details can be found in Section 1 of [1].

A higher-precision implementation is available in the companion package
Python_Fortran, which uses a compiled Fortran extension to achieve
full double precision. That implementation is more accurate but 
requires a build step.

======================================================================
QUICK START
======================================================================

Place Riemann_zeta.py and ln_gamma.py in the same directory (or on
the Python path), then:

    from Riemann_zeta import Riemann_zeta
    from ln_gamma import ln_gamma
    import numpy as np

    f = Riemann_zeta(2.0)                              # scalar input
    f = Riemann_zeta(np.array([2, 4, 6], dtype=complex))  # array input

ln_gamma is used internally by Riemann_zeta and need not be imported
directly, but it is also available as a standalone function:

    f = ln_gamma(z)

======================================================================
FUNCTION DESCRIPTIONS
======================================================================

----------------------------------------------------------------------
f = Riemann_zeta(s)
----------------------------------------------------------------------
Computes the Riemann zeta function zeta(s).

Input:
  s  -- complex scalar, list, or array (converted internally to
        complex128)

Output:
  f  -- complex scalar or array of the same shape as s

Requires ln_gamma.py to be on the Python path.

----------------------------------------------------------------------
f = ln_gamma(z)
----------------------------------------------------------------------
Computes the logarithm of the Gamma function log(Gamma(z)) in the
entire complex plane.

Input:
  z  -- complex scalar, list, or array (converted internally to
        complex128)

Output:
  f  -- complex scalar or array of the same shape as z

======================================================================
IMPLEMENTATION
======================================================================

----------------------------------------------------------------------
Riemann_zeta
----------------------------------------------------------------------
For Re(s) < 0.5, the functional equation is used to reduce the
computation to the half-plane Re(s) >= 0.5. Near s = 0, a Taylor
series is used instead (O(1) cost).

In the half-plane Re(s) >= 0.5, one of three methods is used:

    i)   an approximation zeta_8(s) developed in [1], used when
         Re(s) < 5 and |Im(s)| > 200,
    ii)  the Euler-Maclaurin summation formula, used when
         Re(s) < 5 and |Im(s)| <= 200, or
    iii) direct summation zeta(s) = sum_{n=1}^{N} n^{-s}, used
         when Re(s) >= 5.

The computational complexity is O(sqrt(|Im(s)|)) in the strip
-4 < Re(s) < 5, and O(1) everywhere else.

----------------------------------------------------------------------
ln_gamma
----------------------------------------------------------------------
The algorithm is based on Binet's formula, approximating the
integrand by an exponential sum with 12 terms (4 complex conjugate
pairs and 4 real terms) using techniques from [3]. See [2] for
further details.

======================================================================
REQUIREMENTS
======================================================================

  - Python 3.8 or later
  - NumPy; install with:  pip install numpy

No compilation step is required. matplotlib is needed only for Test 8
in test.py; install with:  pip install matplotlib

======================================================================
TESTS
======================================================================

Run test.py to verify correctness:

    $ python3 test.py

The script contains nine tests. Tests 1-8 print numerical results;
Test 9 produces a plot and runs last (requires matplotlib).

  Test 1:  exp(ln_gamma(k)) for integer k = 1, ..., 10
  Test 2:  exp(ln_gamma(k - 1/2)) for integer k = 1, ..., 5
  Test 3:  Test the functional equation ln_gamma(z+1)=ln_gamma(z)+log(z)
           for 10,000 random z with |Re(z)| < 20, |Im(z)| < 20
  Test 4:  zeta(s) for s = 2, 4, 6, 8, 10 (exact values pi^{2k}/...)
  Test 5:  zeta(s) for s = 0, -1, -3, -5, -7, -9 (known rational values)
  Test 6:  zeta(1/2 + i*t_k) at the first five non-trivial zeros
  Test 7:  zeta(1/2 + i*t_k) at five zeros near height t = 200
  Test 8:  zeta(1/2 + i*t_k) at five zeros near height t = 10^6
  Test 9:  plots the Riemann-Siegel function Z(t) over three ranges:
           0 < t < 100,  1000 < t < 1100,  10000 < t < 10100

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
