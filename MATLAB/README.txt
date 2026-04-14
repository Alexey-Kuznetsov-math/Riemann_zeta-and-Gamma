Author: Alexey Kuznetsov
York University, Toronto, Canada
Website: https://kuznetsovmath.ca/
Email: akuznets@yorku.ca

Created: 18-Nov-2025
Last updated: 2-April-2026

License: BSD 3-Clause (https://opensource.org/licenses/BSD-3-Clause)

======================================================================
PACKAGE CONTENTS
======================================================================

  Riemann_zeta.m    -- MATLAB implementation of the Riemann zeta function
  ln_gamma.m        -- MATLAB implementation of the log-Gamma function
  test.m            -- correctness tests and plots
  README.txt        -- this file

======================================================================
OVERVIEW
======================================================================

This package provides two MATLAB functions, Riemann_zeta and
ln_gamma.

Riemann_zeta(s) computes the Riemann zeta function for scalar, vector,
or array inputs. The implementation is fully vectorized.

Accuracy:
  |Im(s)| < 100      -- approximately 13 correct decimal digits
  |Im(s)| < 1,000    -- approximately 12 correct decimal digits
  |Im(s)| < 10,000   -- approximately 11 correct decimal digits

For larger values of |Im(s)|, accuracy continues to decrease at a
rate of approximately one decimal digit per increase in |Im(s)| by
a factor of ten. More details can be found in Section 1 of [1].

A higher-precision implementation is available in the companion package
MATLAB_Fortran_mex, which uses a Fortran MEX file to achieve full
double precision. That implementation is more accurate but slower 
by a factor of roughly 10 to 20 compared to this one.

======================================================================
QUICK START
======================================================================

Place Riemann_zeta.m and ln_gamma.m in the same directory (or on the
MATLAB path), then call:

    f = Riemann_zeta(s)

where s may be a real or complex scalar, vector, or array. ln_gamma
is used internally by Riemann_zeta and need not be called directly,
but it is also available as a standalone function:

    f = ln_gamma(z)

======================================================================
FUNCTION DESCRIPTIONS
======================================================================

----------------------------------------------------------------------
f = Riemann_zeta(s)
----------------------------------------------------------------------
Computes the Riemann zeta function zeta(s).

Input:
  s  -- real or complex scalar, vector, or array (double precision)

Output:
  f  -- complex array of the same size as s

Requires ln_gamma.m to be on the MATLAB path.

----------------------------------------------------------------------
f = ln_gamma(z)
----------------------------------------------------------------------
Computes the logarithm of the Gamma function log(Gamma(z)) in the
entire complex plane.

Input:
  z  -- real or complex scalar, vector, or array (double precision)

Output:
  f  -- complex array of the same size as z

======================================================================
IMPLEMENTATION
======================================================================

----------------------------------------------------------------------
Riemann_zeta
----------------------------------------------------------------------
For Re(s) < 0.5, the functional equation is used to reduce the
computation to the half-plane Re(s) >= 0.5. Near s = 0, a Taylor
series is used instead.

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
PERFORMANCE
======================================================================

Riemann_zeta is significantly faster than MATLAB's built-in zeta(s),
while also being faster than the companion MEX implementation due to
full vectorization. Typical speedup factors vs. the built-in, measured
on a desktop machine
(Lenovo ThinkCentre M90q Gen 3, Intel Core i5-12500, 16 GB RAM):

    Range of |Im(s)|     Speedup vs. zeta(s)
    ----------------     --------------------
    |Im(s)| < 100             567x
    |Im(s)| < 1,000         1,100x
    |Im(s)| < 10,000       12,256x

======================================================================
TESTS
======================================================================

Run test.m to verify correctness and measure performance:

    >> test

The script contains four tests:

  Test 1:  ln_gamma(z) vs. log(gamma(z)) for z in the range 0.5 to 10

  Test 2:  Testing the functional equation ln_gamma(z+1)=ln_gamma(z)+log(z)
           for 10,000 random z with |Re(z)| < 20, |Im(z)| < 20

  Test 3:  plots the Riemann-Siegel function Z(t) over three ranges:
           0 < t < 100,  1000 < t < 1100,  10000 < t < 10100

  Test 4:  Riemann_zeta vs. the MATLAB built-in zeta(s), with three
           subtests:
             (a) 1000 random s with |Im(s)| < 100,    -5 < Re(s) < 10
             (b) 1000 random s with |Im(s)| < 1,000,  -5 < Re(s) < 10
             (c) 1000 random s with |Im(s)| < 10,000, -5 < Re(s) < 10

Each subtest of Test 4 reports computation times for both
implementations and the maximum relative error max(|f/g - 1|).

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
Test 1: the maximum absolute error of ln_gamma(z)-log(gamma(z)) in the range 1/2<z<10 is
   1.110223024625157e-16

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Test 2: the maximum absolute value of ln_gamma(z+1)-ln_gamma(z)-ln(z)
for 10 000 random points, uniformly distributed in the square |Im(z)|<20, |Re(z)|<20:
   3.552713678800501e-15

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Test 3: plot the Riemann-Siegel function Z(t) for 0<t<100, 1000<t<1100 and 10000<t<10100

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Test 4: Testing the Riemann_zeta function against MATLAB built-in function zeta(s)
 
Test 4(a): computing the relative error Riemann_zeta(s)/zeta(s)-1
for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<100, -5<Re(s)<10:
 
Computation time (in seconds) of Riemann_zeta(s) for these 1000 random numbers s_i:
   0.009198000000000
 
Computation time (in seconds) of MATLAB built-in function zeta(s) for these 1000 random numbers:
   5.218138000000000

the maximum relative error is:
   2.371406605774137e-13

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Test 4(b): computing the maximum relative error Riemann_zeta(s)/zeta(s)-1
for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<1000, -5<Re(s)<10:
 
Computation time (in seconds) of Riemann_zeta(s) for these 1000 random numbers s_i:
   0.012621000000000
 
Computation time (in seconds) of MATLAB built-in function zeta(s) for these 1000 random numbers:
  13.858682000000000

the maximum relative error is:
   3.322460290492961e-12

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Test 4(c): computing the maximum relative error Riemann_zeta(s)/zeta(s)-1
for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<10 000, -5<Re(s)<10:
 
Computation time (in seconds) of Riemann_zeta(s) for these 1000 random numbers s_i:
   0.005968000000000
 
Computation time (in seconds) of MATLAB built-in function zeta(s) for these 1000 random numbers:
  73.057969000000000

the maximum relative error is:
   2.726475074058203e-11
