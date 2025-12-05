Author: Alexey Kuznetsov
York University, Toronto, Canada
Website: https://kuznetsovmath.ca/
Email: akuznets@yorku.ca

Created: 1-Dec-2025
Last updated: 5-Dec-2025  

License: BSD 3-Clause (https://opensource.org/licenses/BSD-3-Clause)
**********************************************************************************************************************
The Fortran90 function Riemann_zeta(s) computes the Riemann zeta function for complex scalar values of s (complex(kind=16)). 
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
**********************************************************************************************************************
To compile and run this code, type the following commands in the MATLAB Command Window:

    >> mex -setup Fortran
    >> mex -c Riemann_zeta_module.f90
    >> mex -R2017b Riemann_zeta_mex.F90 Riemann_zeta_module.o

After successful compilation, values of the Riemann zeta function can be computed in any MATLAB script/function or in the Command Window via

    f = Riemann_zeta_mex(s)

where s may be a real or complex scalar, vector, or array.
**********************************************************************************************************************
The program test.m compares the accuracy and performance of Riemann_zeta_mex(s) with the MATLAB built-in function zeta(s).
A typical output is shown below. We see that Riemann_zeta_mex achieves full double-precision accuracy and is significantly faster
than the built-in zeta(s):

	58x faster in the range |Im(s)|<100
	68x faster in the range |Im(s)|<1000
	315x faster in the range |Im(s)|<10 000
	3955x faster in the range |Im(s)|<100 000

Example output of test.m on a desktop machine
(Lenovo ThinkCentre M90q Gen 3, Intel Core i5-12500, 16 GB RAM):

>> test
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Testing the Riemann_zeta_mex function against the MATLAB built-in function zeta(s)
 
Test (a): computing the relative error Riemann_zeta_mex(s)/zeta(s)-1
for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<100, -5<Re(s)<10:
 
Computation time (in seconds) of Riemann_zeta_mex(s) for these 1000 random numbers s_i:
   0.090132000000000

Computation time (in seconds) of MATLAB built-in function zeta(s) for these 1000 random numbers:
   5.253863000000000

the maximum relative error is:
     2.294551376473652e-16

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Test (b): computing the maximum relative error Riemann_zeta_mex(s)/zeta(s)-1
for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<1000, -5<Re(s)<10:
 
Computation time (in seconds) of Riemann_zeta_mex(s) for these 1000 random numbers s_i:
   0.200024000000000
 
Computation time (in seconds) of MATLAB built-in function zeta(s) for these 1000 random numbers:
  13.493294000000001

the maximum relative error is:
     2.247161154328318e-16

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Test (c): computing the maximum relative error Riemann_zeta_mex(s)/zeta(s)-1
for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<10 000, -5<Re(s)<10:
 
Computation time (in seconds) of Riemann_zeta_mex(s) for these 1000 random numbers s_i:
   0.229147000000000
 
Computation time (in seconds) of MATLAB built-in function zeta(s) for these 1000 random numbers:
  69.356413000000003

the maximum relative error is:
     2.220462989698777e-16

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Test (d): computing the maximum relative error Riemann_zeta_mex(s)/zeta(s)-1
for 1000 random numbers s_i, uniformly distributed in the rectangle |Im(s)|<100 000, -5<Re(s)<10:
 
Computation time (in seconds) of Riemann_zeta_mex(s) for these 1000 random numbers s_i:
   0.296897000000000
 
Computation time (in seconds) of MATLAB built-in function zeta(s) for these 1000 random numbers:
     1.174104960000000e+03

the maximum relative error is:
     7.695268066174558e-17

**********************************************************************************************************************
References:

 [1] A. Kuznetsov, "Simple and accurate approximations to the Riemann zeta function",
     2025, https://arxiv.org/abs/2503.09519

 [2] A. Kuznetsov, "Computing the Barnes G-function and the gamma function
     in the entire complex plane", Journal of Computational and Applied Mathematics,
     Vol. 411, 2022, 114270. https://doi.org/10.1016/j.cam.2022.114270  

 [3] A. Kuznetsov, A. Mohammadioroojeh, "Approximating functions on R^+
     by exponential sums", 2025, https://arxiv.org/abs/2508.19095 
     
