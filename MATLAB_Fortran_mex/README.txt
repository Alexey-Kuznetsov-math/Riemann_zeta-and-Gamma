Author: Alexey Kuznetsov
York University, Toronto, Canada
Website: https://kuznetsovmath.ca/
Email: akuznets@yorku.ca

Created: 1-Dec-2025
Last updated: 21-April-2026

License: BSD 3-Clause (https://opensource.org/licenses/BSD-3-Clause)

======================================================================
PACKAGE CONTENTS
======================================================================

  Riemann_zeta_module.f90      -- Fortran module (core implementation)
  Riemann_zeta_mex.F90         -- MEX gateway for Riemann_zeta_mex
  Riemann_zeta_prime_mex.F90   -- MEX gateway for Riemann_zeta_prime_mex
  build_mex.m                  -- MATLAB build script
  test.m                       -- accuracy and performance tests
  README.txt                   -- this file

======================================================================
OVERVIEW
======================================================================

This package provides two MATLAB functions, Riemann_zeta_mex and
Riemann_zeta_prime_mex, which compute the Riemann zeta function and
its derivative for scalar or array inputs. They are implemented as
MEX files that call a Fortran 90 core module (Riemann_zeta_module.f90)
in which all computations are performed in quadruple precision. 
The results are then rounded and returned to MATLAB as standard 
double-precision complex values.

The inputs to both MEX functions may be real or complex, and may be
scalars or vectors (rows or columns). The output has the same shape as the
input.

Accuracy: the returned values are accurate to full double precision
for |Im(s)| < 10^12.

======================================================================
QUICK START
======================================================================

Step 1: build the MEX files (one-time setup, see BUILDING below)

    >> mex -setup Fortran
    >> build_mex

Step 2: evaluate the Riemann zeta function

    >> f = Riemann_zeta_mex(s)

where s may be a real or complex scalar or array.

Step 3: evaluate zeta(s) and its derivative simultaneously

    >> [f, f1] = Riemann_zeta_prime_mex(s)

returns f = zeta(s) and f1 = zeta'(s). 

======================================================================
FUNCTION DESCRIPTIONS
======================================================================

----------------------------------------------------------------------
f = Riemann_zeta_mex(s)
----------------------------------------------------------------------
Computes the Riemann zeta function zeta(s).

Input:
  s  -- real or complex scalar or vector (double precision)

Output:
  f  -- complex array of the same size as s (double precision);
        f(k) = zeta(s(k)) for each element

----------------------------------------------------------------------
[f, f1] = Riemann_zeta_prime_mex(s)
----------------------------------------------------------------------
Computes the Riemann zeta function and its derivative simultaneously.

Input:
  s   -- real or complex scalar or vector (double precision)

Outputs:
  f   -- zeta(s), complex, same size as s (double precision)
  f1  -- zeta'(s), complex, same size as s (double precision)

======================================================================
IMPLEMENTATION
======================================================================

The MEX gateways extract the real and imaginary parts of each element
of s as double-precision values, promote them to quadruple
precision (Fortran kind = selected_real_kind(33, 4931)), and call
the Fortran routine Riemann_zeta or Riemann_zeta_prime. The results
are then rounded back to double precision and returned to MATLAB.

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
PERFORMANCE
======================================================================

Riemann_zeta_mex is significantly faster than MATLAB's built-in
zeta(s). Typical speedup factors measured on a desktop machine
(Lenovo ThinkCentre M90q Gen 3, Intel Core i5-12500, 16 GB RAM):

    Range of |Im(s)|     Speedup vs. zeta(s)
    ----------------     --------------------
    |Im(s)| < 100              100x
    |Im(s)| < 1,000            120x
    |Im(s)| < 10,000           482x
    |Im(s)| < 100,000        5,820x

Similarly, Riemann_zeta_prime_mex computes zeta(s) and zeta'(s)
simultaneously and is substantially faster than calling the MATLAB
built-in functions zeta(s) and zeta(1,s) separately.

======================================================================
BUILDING
======================================================================

Prerequisites:
  - MATLAB R2017b or later
  - A Fortran compiler supported by MATLAB's MEX interface. gfortran
    (GCC) is recommended. To check or configure the compiler, run:

        >> mex -setup Fortran

    This step only needs to be done once per MATLAB installation.

To compile both MEX files, place all package files in the same
directory, set that directory as the MATLAB working directory, and
run:

    >> build_mex

This script compiles Riemann_zeta_mex and Riemann_zeta_prime_mex with
the -R2017b flag (required for the separated complex API used by
mxGetPr/mxGetPi in MATLAB R2018a and later) and the -O3 optimization
flag. The compiled MEX files are placed in the same directory.

Once compiled, the MEX files can be called from any MATLAB script,
function, or Command Window session in the same directory (or any
directory on the MATLAB path).

======================================================================
TESTS
======================================================================

Run test.m to verify correctness and measure performance:

    >> test

The script contains nine tests:

  Testing the Riemann_zeta_prime_mex function: computing the first 
  non-trivial zero with Newton's method (starting with the approximate 
  value of s=1/2+14.1i).
  
  Testing the Riemann_zeta_prime_mex function: computing the 
  non-trivial zero close to 1/2+401.8i with Newton's method. 
  
  Test (a): Riemann_zeta_mex vs. zeta(s), 1000 random s with
            |Im(s)| < 100, -5 < Re(s) < 10
  Test (b): same comparison, |Im(s)| < 1,000
  Test (c): same comparison, |Im(s)| < 10,000
  Test (d): same comparison, |Im(s)| < 100,000

  Test (e): Riemann_zeta_prime_mex vs. zeta(s) and zeta(1,s),
            1000 random s with |Im(s)| < 100, -5 < Re(s) < 10
  Test (f): same comparison, |Im(s)| < 1,000
  Test (g): same comparison, |Im(s)| < 10,000

Each test reports computation times for both implementations and the
maximum mixed absolute/relative error, defined as:
  error(k) = |f(k) - g(k)|       if |g(k)| <= 1,
  error(k) = |f(k)/g(k) - 1|     if |g(k)| > 1.

A sample output from the test machine described above is included at
the end of this file.

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

======================================================================
SAMPLE OUTPUT
======================================================================

The following output was produced on a desktop machine
(Lenovo ThinkCentre M90q Gen 3, Intel Core i5-12500, 16 GB RAM):

>> test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Testing the Riemann_zeta_prime_mex function: computing the first non-trivial zero with Newton's method

s = 0.500489975147131 +14.134820075706351i
s = 0.499999913991052 +14.134725085038156i
s = 0.499999999999999 +14.134725141734689i
s = 0.500000000000000 +14.134725141734695i
s = 0.500000000000000 +14.134725141734695i
error = 0

Testing the Riemann_zeta_prime_mex function: computing the non-trivial zero close to 1/2+401.8i with Newton's method

s = 5.029976298275870e-01 + 4.018378529900980e+02i
s = 4.999922287853202e-01 + 4.018392517543437e+02i
s = 5.000000012889951e-01 + 4.018392286008850e+02i
s = 5.000000000000000e-01 + 4.018392286005332e+02i
s = 5.000000000000000e-01 + 4.018392286005332e+02i
error = 0


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Testing the Riemann_zeta_mex function against the MATLAB built-in function zeta(s)
 
Test (a): comparing the values of Riemann_zeta_mex(s) and zeta(s)
for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<100, -5<Re(s)<10:
 
Computation time (in seconds) of Riemann_zeta_mex(s) for these 1000 random numbers s_i:
    0.0548
 
Computation time (in seconds) of MATLAB built-in function zeta(s) for these 1000 random numbers:
    5.4532

the maximum mixed absolute/relative error is:
   4.4411e-16

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Testing the Riemann_zeta_mex function against the MATLAB built-in function zeta(s)
 
Test (b): comparing the values of Riemann_zeta_mex(s) and zeta(s)
for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<1000, -5<Re(s)<10:
 
Computation time (in seconds) of Riemann_zeta_mex(s) for these 1000 random numbers s_i:
    0.1204
 
Computation time (in seconds) of MATLAB built-in function zeta(s) for these 1000 random numbers:
   14.4796

the maximum mixed absolute/relative error is:
   2.2217e-16

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Testing the Riemann_zeta_mex function against the MATLAB built-in function zeta(s)
 
Test (c): comparing the values of Riemann_zeta_mex(s) and zeta(s)
for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<10 000, -5<Re(s)<10:
 
Computation time (in seconds) of Riemann_zeta_mex(s) for these 1000 random numbers s_i:
    0.1527
 
Computation time (in seconds) of MATLAB built-in function zeta(s) for these 1000 random numbers:
   73.6357

the maximum mixed absolute/relative error is:
   5.2735e-17

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Testing the Riemann_zeta_mex function against the MATLAB built-in function zeta(s)
 
Test (d): comparing the values of Riemann_zeta_mex(s) and zeta(s)
for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<100 000, -5<Re(s)<10:
 
Computation time (in seconds) of Riemann_zeta_mex(s) for these 1000 random numbers s_i:
    0.1893
 
Computation time (in seconds) of MATLAB built-in function zeta(s) for these 1000 random numbers:
   1.1019e+03

the maximum mixed absolute/relative error is:
   4.4792e-16

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Testing the Riemann_zeta_prime_mex function against the MATLAB built-in functions zeta(s) and zeta(1,s)
 
Test (e): comparing the values of Riemann_zeta_prime_mex(s) and zeta(s), zeta(1,s)
for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<100, -5<Re(s)<10:
 
Computation time (in seconds) of Riemann_zeta_prime_mex(s) for these 1000 random numbers s_i:
    0.0650

Computation time (in seconds) of MATLAB built-in functions zeta(s) and zeta(1,s) for these 1000 random numbers:
   12.8844

the maximum mixed absolute/relative error for zeta(s) is:
   1.1329e-16

the maximum mixed absolute/relative error for derivative of zeta(s) is:
   1.1102e-16

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Testing the Riemann_zeta_prime_mex function against the MATLAB built-in functions zeta(s) and zeta(1,s)
 
Test (f): comparing the values of Riemann_zeta_prime_mex(s) and zeta(s), zeta(1,s)
for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<1000, -5<Re(s)<10:
 
Computation time (in seconds) of Riemann_zeta_prime_mex(s) for these 1000 random numbers s_i:
    0.1408
 
Computation time (in seconds) of MATLAB built-in functions zeta(s) and zeta(1,s) for these 1000 random numbers:
   35.0687

the maximum mixed absolute/relative error for zeta(s) is:
   1.2143e-17

the maximum mixed absolute/relative error for derivative of zeta(s) is:
   2.2209e-16

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Testing the Riemann_zeta_prime_mex function against the MATLAB built-in functions zeta(s) and zeta(1,s)
 
Test (g): comparing the values of Riemann_zeta_prime_mex(s) and zeta(s), zeta(1,s)
for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<10 000, -5<Re(s)<10:
 
Computation time (in seconds) of Riemann_zeta_prime_mex(s) for these 1000 random numbers s_i:
    0.1774

Computation time (in seconds) of MATLAB built-in functions zeta(s) and zeta(1,s) for these 1000 random numbers:
  187.3539

the maximum mixed absolute/relative error for zeta(s) is:
   4.4430e-16

the maximum mixed absolute/relative error for derivative of zeta(s) is:
   2.2256e-16
