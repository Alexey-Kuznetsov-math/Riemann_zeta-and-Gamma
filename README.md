# Riemann_zeta-and-Gamma
Fortran, MATLAB and Python code for computing the Riemann zeta and log Gamma functions. 

**********************************************************************************************************************

Fortran90 function Riemann_zeta(s) computes the Riemann zeta function for complex values of s.

The input s must be a scalar of type complex(kind=16) (note that this function is not vectorized).

The precision is close to quadruple:
    For |Im(s)| < 100   the result is correct to about 31 decimal digits;
    For |Im(s)| < 1000  the result is correct to about 30 decimal digits;
    For |Im(s)| < 10000 the result is correct to about 29 decimal digits.

For larger values of |Im(s)|, the accuracy continues to decrease in a similar way: with every increase of |Im(s)|
by a factor of ten, we lose approximately one decimal digit of precision.
More details on this loss of precision for large values of |Im(s)| can be found in [1] (see the list of references below).

When |Im(s)|>150 and -9<Re(s)<30 we use the approximation zeta_30(s) developed in [1]. 
For other values of s we use either Euler-Maclaurin formula or direct summation zeta(s)=\sum_{n=1}^{\infty} n^{-s}.
The computational complexity is O(sqrt(|Im(s)|)) in the strip -9 < Re(s) < 10, and O(1) everywhere else in the complex plane.

**********************************************************************************************************************

Fortran90 function ln_gamma(z) computes the logarithm of the Gamma function in the entire complex plane to quadruple precision.

The input z must be a scalar of type complex(kind=16) (note that this function is not vectorized).

This function is based on the approximation to log(Gamma(z)) developed in [2], with improved coefficients obtained in [3].

**********************************************************************************************************************

To compile and run the test.f90 program on Linux with the gfortran compiler, use the following command:

    gfortran test.f90 -o test.out && ./test.out

**********************************************************************************************************************   

MATLAB and Python function Riemann_zeta(s) returns an approximation to zeta(s) for complex input s.
The input s can be a scalar or vector.

  For |Im(s)|<100   the approximation is correct to around 13 decimal digits;
  For |Im(s)|<1000  the approximation is correct to around 12 decimal digits;
  For |Im(s)|<10000 the approximation is correct to around 11 decimal digits.
  
For larger values of |Im(s)| the accuracy will continue to decrease in a similar way: 
with every increase of Im(s) by a factor of ten we lose approximately one decimal digit of precision. 
More details can be found at the end of Section 1 in [1].

When |Im(s)|>200 and -4<Re(s)<5 we use the approximation zeta_8(s) developed in [1]. 
For other values of s we use either Euler-Maclaurin formula or direct summation zeta(s)=\sum_{n=1}^{\infty} n^{-s}.
The computational complexity is O(sqrt(|Im(s)|)) in the strip -4 < Re(s) < 5, and O(1) everywhere else in the complex plane.

**********************************************************************************************************************  

MATLAB and Python function ln_gamma(z) computes the logarithm of the Gamma function in the entire complex plane to double precision.
The input z can be a scalar, vector, or array. 

********************************************************************************************************************** 

References:

 [1] A. Kuznetsov, "Simple and accurate approximations to the Riemann zeta function",
     2025, https://arxiv.org/abs/2503.09519

 [2] A. Kuznetsov, "Computing the Barnes G-function and the gamma function
     in the entire complex plane", Journal of Computational and Applied Mathematics,
     Vol. 411, 2022, 114270. https://doi.org/10.1016/j.cam.2022.114270

 [3] A. Kuznetsov, A. Mohammadioroojeh, "Approximating functions on R^+
     by exponential sums", 2025, https://arxiv.org/abs/2508.19095
     
**********************************************************************************************************************

Author: Alexey Kuznetsov
York University, Toronto, Canada
Website: https://kuznetsovmath.ca/
Email: akuznets@yorku.ca

Created: 28-Nov-2025
Last updated: 28-Nov-2025     
