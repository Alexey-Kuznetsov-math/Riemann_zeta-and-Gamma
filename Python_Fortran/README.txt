Author: Alexey Kuznetsov
York University, Toronto, Canada
Website: https://kuznetsovmath.ca/
Email: akuznets@yorku.ca

Created: 26-March-2026
Last updated: 2-April-2026

License: BSD 3-Clause (https://opensource.org/licenses/BSD-3-Clause)
**********************************************************************************************************************
This package provides a Python interface to a Fortran 90 implementation of the Riemann zeta function.

The core routine is the Fortran function

    Riemann_zeta(s)

which computes \zeta(s) for complex scalar values of s using quadruple precision arithmetic
(complex(kind=16)). The Python interface accepts standard complex inputs, calls the Fortran code,
and returns the result in double precision.

The main files are:

    Riemann_zeta_module.f90   main Fortran implementation of the Riemann zeta function
    zeta_py_wrapper.f90       Fortran wrapper routines exposed to Python via f2py
    pyzeta.py                 Python interface function Riemann_zeta(s)
    __init__.py               package-level import
    build.sh                  build script for compiling the extension module
**********************************************************************************************************************
Depending on the value of s, the Fortran code uses one of several methods, including

    i)   the approximation zeta_12(s),
    ii)  Euler-Maclaurin summation, or
    iii) direct summation of the Dirichlet series.

The computational complexity is O(sqrt(|Im(s)|)) in the critical strip. In other regions of the
complex plane it is either O(sqrt(|Im(s)|)) or O(1), depending on the location of s.
**********************************************************************************************************************
To build the extension module on Linux, run

    chmod +x build.sh
    ./build.sh

This creates a compiled extension module named something like

    _pyzeta_fortran.cpython-310-x86_64-linux-gnu.so
**********************************************************************************************************************
After successful compilation, the package can be used from Python as follows:

    from pyzeta import Riemann_zeta

    print(Riemann_zeta(2))
    print(Riemann_zeta(0.5 + 14j))

It also accepts NumPy arrays.
**********************************************************************************************************************
References:

 [1] A. Kuznetsov, "Simple and accurate approximations to the Riemann zeta function",
     2025, https://arxiv.org/abs/2503.09519

 [2] A. Kuznetsov, "Computing the Barnes G-function and the gamma function
     in the entire complex plane", Journal of Computational and Applied Mathematics,
     Vol. 411, 2022, 114270. https://doi.org/10.1016/j.cam.2022.114270  

 [3] A. Kuznetsov, A. Mohammadioroojeh, "Approximating functions on R^+
     by exponential sums", 2025, https://arxiv.org/abs/2508.19095 
     
