Author: Alexey Kuznetsov
York University, Toronto, Canada
Website: https://kuznetsovmath.ca/
Email: akuznets@yorku.ca

Created: 28-Nov-2025
Last updated: 5-Dec-2025  

License: BSD 3-Clause (https://opensource.org/licenses/BSD-3-Clause)
**********************************************************************************************************************

The function Riemann_zeta(s) computes the Riemann zeta function in the entire complex plane.

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

The function ln_gamma(z) computes the logarithm of the Gamma function in the entire complex plane to quadruple precision.

The input z must be a scalar of type complex(kind=16) (note that this function is not vectorized).

This function is based on the approximation to log(Gamma(z)) developed in [2], with improved coefficients obtained in [3].
**********************************************************************************************************************

To compile and run the test.f90 program on Linux with the gfortran compiler, use the following command:
    gfortran test.f90 -o test.out && ./test.out


References:

 [1] A. Kuznetsov, "Simple and accurate approximations to the Riemann zeta function",
     2025, https://arxiv.org/abs/2503.09519

 [2] A. Kuznetsov, "Computing the Barnes G-function and the gamma function
     in the entire complex plane", Journal of Computational and Applied Mathematics,
     Vol. 411, 2022, 114270. https://doi.org/10.1016/j.cam.2022.114270

 [3] A. Kuznetsov, A. Mohammadioroojeh, "Approximating functions on R^+
     by exponential sums", 2025, https://arxiv.org/abs/2508.19095
    
    
