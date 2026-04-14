Author: Alexey Kuznetsov
York University, Toronto, Canada
Website: https://kuznetsovmath.ca/
Email: akuznets@yorku.ca

Created: 26-March-2026
Last updated: 14-April-2026

License: BSD 3-Clause (https://opensource.org/licenses/BSD-3-Clause)

======================================================================
PACKAGE CONTENTS
======================================================================

  Riemann_zeta_module.f90    -- Fortran module (core implementation)
  Riemann_zeta_wrapper.f90   -- f2py wrapper subroutines
  build.sh                   -- build script
  Riemann_zeta_python.py     -- Python interface module
  ln_gamma.py                -- pure Python/NumPy log-Gamma function
  test.py                    -- correctness tests and plots
  README.txt                 -- this file

======================================================================
OVERVIEW
======================================================================

This package provides two Python functions, Riemann_zeta and
Riemann_zeta_prime, which compute the Riemann zeta function and its
derivative for scalar or NumPy arrays of any shape. The core computation is
performed by the Fortran module Riemann_zeta_module.f90 in quadruple
precision; the results are returned to Python as standard
double-precision complex NumPy arrays.

The package also includes ln_gamma.py, a standalone pure Python/NumPy
implementation of the logarithm of the Gamma function, used in test.py
to compute the Riemann-Siegel Z function.

Accuracy: the returned values are accurate to full double precision
for |Im(s)| < 10^12.

======================================================================
QUICK START
======================================================================

Step 1: build the extension module (one-time setup, see BUILDING below)

    $ ./build.sh

Step 2: import and use the functions

    from Riemann_zeta_python import Riemann_zeta, Riemann_zeta_prime
    import numpy as np

    f = Riemann_zeta(2.0)                              # scalar input
    f = Riemann_zeta(np.array([2, 4, 6], dtype=complex))  # array input

    f1, f2 = Riemann_zeta_prime(2.0)          # f1=zeta(s), f2=zeta'(s)

======================================================================
FUNCTION DESCRIPTIONS
======================================================================

----------------------------------------------------------------------
f = Riemann_zeta(s)
----------------------------------------------------------------------
Computes the Riemann zeta function zeta(s).

Input:
  s  -- complex scalar or array (converted internally to
        complex128)

Output:
  f  -- complex scalar or array of the same shape as s

----------------------------------------------------------------------
f1, f2 = Riemann_zeta_prime(s)
----------------------------------------------------------------------
Computes the Riemann zeta function and its derivative simultaneously.

Input:
  s   -- complex scalar or array

Outputs:
  f1  -- zeta(s), complex scalar or array of the same shape as s
  f2  -- zeta'(s), complex scalar or array of the same shape as s

----------------------------------------------------------------------
f = ln_gamma(z)   [from ln_gamma.py]
----------------------------------------------------------------------
Computes the logarithm of the Gamma function log(Gamma(z)) in the
entire complex plane using a pure Python/NumPy implementation.

Input:
  z  -- complex scalar or array

Output:
  f  -- complex scalar or array of the same shape as z

This function is implemented directly in Python using a short
exponential sum approximation (no compiled extension required). It is
provided primarily for use in computing the Riemann-Siegel Z function,
as demonstrated in Test 7 of test.py.

======================================================================
IMPLEMENTATION
======================================================================

The package has a two-layer architecture:

  Layer 1 -- compiled extension (_Riemann_zeta_python.cpython-*.so):
    Built by f2py from Riemann_zeta_wrapper.f90, which provides two
    subroutines (Riemann_zeta_wrap and Riemann_zeta_prime_wrap). Each
    subroutine accepts a double-precision complex array, converts each
    element to quadruple precision (Fortran kind =
    selected_real_kind(33, 4931)), calls Riemann_zeta or
    Riemann_zeta_prime from the core module, and rounds the result
    back to double precision.

  Layer 2 -- Python interface (Riemann_zeta_python.py):
    Wraps the compiled extension to handle scalar inputs and arrays
    of arbitrary shape, presenting a clean NumPy-compatible API.

The computational complexity is O(sqrt(|Im(s)|)) in the critical 
strip 0<=Re(s)<=1. In other regions of the complex plane it is either 
O(sqrt(|Im(s)|))  (if we are close to the critical strip) or O(1). 

Depending on the value of s, one of three methods is used:
    i)   an approximation zeta_12(s) developed in [1],
    ii)  the Euler-Maclaurin summation formula, or
    iii) direct summation zeta(s) = sum_{n=1}^{inf} n^{-s}, truncated
         at N chosen to achieve the target precision.

Values of zeta'(s) are computed by differentiating the approximations
in items (i)-(iii). 

======================================================================
BUILDING
======================================================================

Prerequisites:
  - Python 3.8 or later
  - NumPy (provides f2py); install with:  pip install numpy
  - gfortran (GCC Fortran compiler)

To build, run:

    $ ./build.sh

If the script is not executable, first run:

    $ chmod +x build.sh

The script performs two steps:
  1. Compiles Riemann_zeta_module.f90 with gfortran (-O3) to produce
     Riemann_zeta_module.o and Riemann_zeta_module.mod.
  2. Calls f2py to compile Riemann_zeta_wrapper.f90 and link it
     against Riemann_zeta_module.o, producing the shared library
     _Riemann_zeta_python.cpython-*.so in the current directory.

After a successful build, a brief self-test is run automatically.

Both Riemann_zeta_python.py and _Riemann_zeta_python.cpython-*.so
must be in the same directory (or on the Python path) to use the
package.

======================================================================
TESTS
======================================================================

Run test.py to verify correctness:

    $ python3 test.py

The script contains nine tests. Tests 1-8 print numerical results;
Test 9 produces a plot (requires matplotlib).

  Test 1:  zeta(s) for s = 2, 4, 6, 8, 10 (exact values pi^{2k}/...)
  Test 2:  zeta(s) for s = 0, -1, -3, -5, -7, -9 (known rational values)
  Test 3:  zeta(1/2 + i*t_k) at the first five non-trivial zeros
  Test 4:  zeta(1/2 + i*t_k) at five zeros near height t = 10^6
  Test 5:  zeta(1/2 + i*t_k) at five zeros near height t = 10^9
  Test 6:  comparison with pre-computed quadruple-precision reference
           values at ten points across the complex plane
  Test 7:  Newton's method for the first non-trivial zero of zeta(s),
           starting from s = 1/2 + 14.1i
  Test 8:  Newton's method for the zero near 1/2 + 401.8i
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
