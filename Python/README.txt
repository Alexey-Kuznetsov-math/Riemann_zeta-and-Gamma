Author: Alexey Kuznetsov
York University, Toronto, Canada
Website: https://kuznetsovmath.ca/
Email: akuznets@yorku.ca

Created: 29-Nov-2025
Last updated: 5-Dec-2025  

License: BSD 3-Clause (https://opensource.org/licenses/BSD-3-Clause)
**********************************************************************************************************************
    f = Riemann_zeta(s) returns an approximation to zeta(s) for complex input s
    The input s can be a scalar, a list, or a NumPy array.

  For |Im(s)|<100 the approximation is correct to around 13 decimal digits;
  For |Im(s)|<1000 the approximation is correct to around 12 decimal digits;
  For |Im(s)|<10 000 the approximation is correct to around 11 decimal digits;

For larger values of |Im(s)| the accuracy will continue to decrease in a similar way: 
with every increase of Im(s) by a factor of ten we lose approximately one decimal digit of precision. 
More details can be found at the end of Section 1 in [1] (see the list of references below).

When |Im(s)|>200 and -4<Re(s)<5 we use the approximation `zeta_8(s)` developed in [1]. 
For other values of s we use either Euler-Maclaurin formula or direct summation zeta(s)=\sum_{n=1}^{\infty} n^{-s}.
The computational complexity is O(sqrt(|Im(s)|)) in the strip -4 < Re(s) < 5, and O(1) everywhere else in the complex plane.
**********************************************************************************************************************
    f = ln_gamma(z) returns log(Gamma(z)) for any complex input z.
    The input z can be  a scalar, a list, or a NumPy array.
    
This algorithm was developed in [2,3]. 
**********************************************************************************************************************
References:

 [1] A. Kuznetsov, "Simple and accurate approximations to the Riemann zeta function",
     2025, https://arxiv.org/abs/2503.09519

 [2] A. Kuznetsov, "Computing the Barnes G-function and the gamma function
     in the entire complex plane", Journal of Computational and Applied Mathematics,
     Vol. 411, 2022, 114270. https://doi.org/10.1016/j.cam.2022.114270

 [3] A. Kuznetsov, A. Mohammadioroojeh, "Approximating functions on R^+
     by exponential sums", 2025, https://arxiv.org/abs/2508.19095 
