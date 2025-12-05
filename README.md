# Riemann_zeta-and-Gamma
Fortran, MATLAB and Python code for computing the Riemann zeta and log Gamma functions. 

**********************************************************************************************************************

Fortran90 function Riemann_zeta(s) computes the Riemann zeta function for complex values of s.

The input s must be a scalar of type complex(kind=16) (note that this function is not vectorized).

The precision is close to quadruple:
    For |Im(s)| < 100   the approximation is correct to about 31 decimal digits;
    For |Im(s)| < 1000  the approximation is correct to about 30 decimal digits;
    For |Im(s)| < 10000 the approximation is correct to about 29 decimal digits.

For larger values of |Im(s)|, the accuracy continues to decrease in a similar way: with every increase of |Im(s)|
by a factor of ten, we lose approximately one decimal digit of precision.
More details on this loss of precision for large values of |Im(s)| can be found in [1] (see the list of references below).

The computational complexity is O(sqrt(|Im(s)|)) in the critical strip. In other regions of the complex plane it is either O(sqrt(|Im(s)|))  
(if we are close to the critical strip) or O(1). Depending on the value of s, we use one of the following methods to compute zeta(s):
    i) an approximation zeta_30(s) developed in [1],
    ii) Euler-Maclaurin summation formula, or
    iii) direct summation zeta(s)=\sum_{n=1}^{\infty} n^{-s}, truncated at an appropriate N

**********************************************************************************************************************

Fortran90 function ln_gamma(z) computes the logarithm of the Gamma function in the entire complex plane to quadruple precision.

The input z must be a scalar of type complex(kind=16) (note that this function is not vectorized).

This function is based on the approximation to log(Gamma(z)) developed in [2], with improved coefficients obtained in [3].

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

The function Riemann_zeta(s), while not as accurate as the built-in MATLAB function zeta(s), is significantly faster than the built-in zeta(s):

  567× faster in the range |Im(s)| < 100
  1100× faster in the range |Im(s)| < 1000
  12256× faster in the range |Im(s)| < 10 000

**********************************************************************************************************************  

MATLAB and Python function ln_gamma(z) computes the logarithm of the Gamma function in the entire complex plane to double precision.
The input z can be a scalar, vector, or array. 

********************************************************************************************************************** 
MATLAB_Forran_mex: The Fortran90 function Riemann_zeta(s) computes the Riemann zeta function for complex scalar values of s (complex(kind=16)). 
The output of this function is also complex(kind=16) and the accuracy is at least 17 decimal digits for |Im(s)|<10^12. 

This Fortran function uses the approximation zeta_12(s) derived in [1], as well as the quadruple precision log(Gamma(s)) 
implementation derived in [2,3] (see the references below).

The MEX file provides an interface layer between MATLAB and the Fortran 90 implementation of the Riemann zeta function.
The routine mexFunction receives MATLAB inputs as mxArray pointers, extracts the real and (optionally) imaginary
parts as double-precision arrays, promotes them to complex(kind=16), and calls the Fortran routine Riemann_zeta for each element.
The results are then converted back to double-precision complex values and stored in a new mxArray, which is returned to MATLAB
as a standard MATLAB variable.

The computational complexity is O(sqrt(|Im(s)|)) in the critical strip. In other regions of the complex plane it is either O(sqrt(|Im(s)|))  
(if we are close to the critical strip) or O(1). Depending on the value of s, we use one of the following methods to compute zeta(s):
    i) an approximation zeta_12(s) developed in [1],
    ii) Euler-Maclaurin summation formula, or
    iii) direct summation zeta(s)=\sum_{n=1}^{\infty} n^{-s}, truncated at an appropriate N

The function Riemann_zeta_mex achieves full double-precision accuracy and is significantly faster
than the built-in MATLAB function zeta(s):

	58x faster in the range |Im(s)|<100
	68x faster in the range |Im(s)|<1000
	315x faster in the range |Im(s)|<10 000
	3955x faster in the range |Im(s)|<100 000    

********************************************************************************************************************** 

References:

 [1] A. Kuznetsov, "Simple and accurate approximations to the Riemann zeta function",
     2025, https://arxiv.org/abs/2503.09519

 [2] A. Kuznetsov, "Computing the Barnes G-function and the gamma function
     in the entire complex plane", Journal of Computational and Applied Mathematics,
     Vol. 411, 2022, 114270. https://doi.org/10.1016/j.cam.2022.114270

 [3] A. Kuznetsov, A. Mohammadioroojeh, "Approximating functions on R^+
     by exponential sums", 2025, https://arxiv.org/abs/2508.19095   
