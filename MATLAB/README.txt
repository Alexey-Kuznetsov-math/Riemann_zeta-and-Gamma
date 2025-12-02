Author: Alexey Kuznetsov
York University, Toronto, Canada
Website: https://kuznetsovmath.ca/
Email: akuznets@yorku.ca

Created: 1-Dec-2025
Last updated: 2-Dec-2025  

License: BSD 3-Clause (https://opensource.org/licenses/BSD-3-Clause)
**********************************************************************************************************************
    f = Riemann_zeta(s) returns an approximation to zeta(s) for complex input s
    The input s can be a scalar, a vector or an array

  For |Im(s)|<100 the approximation is correct to around 13 decimal digits;
  For |Im(s)|<1000 the approximation is correct to around 12 decimal digits;
  For |Im(s)|<10 000 the approximation is correct to around 11 decimal digits;

For larger values of |Im(s)| the accuracy will continue to decrease in a similar way: 
with every increase of Im(s) by a factor of ten we lose approximately one decimal digit of precision. 
More details can be found at the end of Section 1 in [1] (see the list of references below).

A full double-precision implementation of Riemann zeta function can be found in the folder MATLAB_Fortran_mex. 
It includes a MEX file as an interface layer between MATLAB and the Fortran 90 function `Riemann_zeta(s)`, which computes zeta(s) to precision 
of 17 decimal digits (or higher) using quadruple precision numbers. That implementation is more precise than 
the current one (achieving full double-precision), but it is also slower by a factor of 20 to 50. 

When |Im(s)|>200 and -4<Re(s)<5 we use the approximation `zeta_8(s)` developed in [1]. 
For other values of s we use either Euler-Maclaurin formula or direct summation zeta(s)=\sum_{n=1}^{\infty} n^{-s}.
The computational complexity is O(sqrt(|Im(s)|)) in the strip -4 < Re(s) < 5, and O(1) everywhere else in the complex plane.


**********************************************************************************************************************
    f = ln_gamma(z) returns log(Gamma(z)) for any complex input z.
    The input z can be a scalar, a vector, or an array.

This algorithm was developed in [2,3]. 
**********************************************************************************************************************

The fourth test in the program `test.m` compares the accuracy and performance of `Riemann_zeta(s)` with the 
MATLAB built-in function `zeta(s)`. A typical output is shown below. We see that `Riemann_zeta` achieves the accuracy stated above and is significantly faster than the built-in `zeta(s)`:

  458× faster in the range |Im(s)| < 100
  1089× faster in the range |Im(s)| < 1000
  8558× faster in the range |Im(s)| < 10 000

Example output of `test.m` on a desktop machine
(Lenovo ThinkCentre M90q Gen 3, Intel Core i5-12500, 16 GB RAM):

Test 4(a): computing the relative error Riemann_zeta(s)/zeta(s)-1
for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<100, -5<Re(s)<10:
 
Computation time (in seconds) of Riemann_zeta(s) for these 1000 random numbers s_i:
    0.0119

 
Computation time (in seconds) of the MATLAB built-in function zeta(s) for these 1000 random numbers:
    5.4575

the maximum relative error is:
   4.4021e-13

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Test 4(b): computing the maximum relative error Riemann_zeta(s)/zeta(s)-1
for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<1000, -5<Re(s)<10:
 
Computation time (in seconds) of Riemann_zeta(s) for these 1000 random numbers s_i:
    0.0130

 
Computation time (in seconds) of the MATLAB built-in function zeta(s) for these 1000 random numbers:
   14.1584

the maximum relative error is:
   2.4620e-12

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Test 4(c): computing the maximum relative error Riemann_zeta(s)/zeta(s)-1
for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<10 000, -5<Re(s)<10:
 
Computation time (in seconds) of Riemann_zeta(s) for these 1000 random numbers s_i:
    0.0087

 
Computation time (in seconds) of the MATLAB built-in function zeta(s) for these 1000 random numbers:
   74.4586

the maximum relative error is:
   2.4216e-11


**********************************************************************************************************************
References:

 [1] A. Kuznetsov, "Simple and accurate approximations to the Riemann zeta function",
     2025, https://arxiv.org/abs/2503.09519

 [2] A. Kuznetsov, "Computing the Barnes G-function and the gamma function
     in the entire complex plane", Journal of Computational and Applied Mathematics,
     Vol. 411, 2022, 114270. https://doi.org/10.1016/j.cam.2022.114270

 [3] A. Kuznetsov, A. Mohammadioroojeh, "Approximating functions on R^+
     by exponential sums", 2025, https://arxiv.org/abs/2508.19095 
