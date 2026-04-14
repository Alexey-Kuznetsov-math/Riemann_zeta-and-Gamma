# Riemann_zeta-and-Gamma

This repository provides Fortran, MATLAB, and Python implementations
for computing the Riemann zeta function zeta(s), its derivative zeta'(s), and the logarithm of
the Gamma function log(Gamma(z)) for complex values of the variable.

Five packages are included:

  Fortran/             -- Fortran 90, ~30-31 digits, scalar only
  
  MATLAB/              -- MATLAB, ~13 digits, vectorized
  
  MATLAB_Fortran_mex/  -- MATLAB + Fortran, full double precision,
                          vectorized, requires build step
                          
  Python/              -- Python/NumPy, ~13 digits, vectorized
  
  Python_Fortran/      -- Python + Fortran, full double precision,
                          vectorized, requires build step


Fortran
=======================================

The module gamma_zeta_module.f90 provides four public functions. All
inputs must be scalars of type complex(kind=qp), where
qp = selected_real_kind(33, 4931).

  Riemann_zeta(s)        -- Riemann zeta function zeta(s)
  Riemann_zeta_prime(s)  -- returns array [zeta(s), zeta'(s)]
  ln_gamma(z)            -- log(Gamma(z)), relative accuracy 10^{-33}
  psi(z)                 -- digamma function psi(z)=Gamma'(z)/Gamma(z),
                            relative accuracy 10^{-32}

Precision of Riemann_zeta:

  |Im(s)| < 100      -- approximately 31 correct decimal digits
  |Im(s)| < 1,000    -- approximately 30 correct decimal digits
  |Im(s)| < 10,000   -- approximately 29 correct decimal digits

Accuracy decreases by approximately one digit per increase in
|Im(s)| by a factor of ten. See Section 1 of [1] for details.

MATLAB
=======================================

Two MATLAB functions:

  Riemann_zeta(s)  -- zeta(s) for scalar, vector, or array input
  ln_gamma(z)      -- log(Gamma(z)) for scalar, vector, or array input

Precision of Riemann_zeta:

  |Im(s)| < 100      -- approximately 13 correct decimal digits
  |Im(s)| < 1,000    -- approximately 12 correct decimal digits
  |Im(s)| < 10,000   -- approximately 11 correct decimal digits

While less accurate than MATLAB's built-in zeta(s), this
implementation is significantly faster due to full vectorization:

  |Im(s)| < 100        567x faster than built-in zeta(s)
  |Im(s)| < 1,000    1,100x faster than built-in zeta(s)
  |Im(s)| < 10,000  12,256x faster than built-in zeta(s)

MATLAB_Fortran_mex
=======================================

Two MATLAB functions backed by a quadruple-precision Fortran core,
accessed via MEX. All computations are performed internally in
quadruple precision; results are returned as standard double-precision
complex values.

  Riemann_zeta_mex(s)        -- zeta(s) for scalar or array input
  Riemann_zeta_prime_mex(s)  -- returns [zeta(s), zeta'(s)]
                                for scalar or array input

Accuracy: full double precision for |Im(s)| < 10^12.

Speedup vs. MATLAB built-in zeta(s):

  |Im(s)| < 100          79x
  |Im(s)| < 1,000        99x
  |Im(s)| < 10,000      463x
  |Im(s)| < 100,000   5,728x

Python
=======================================

Two  Python/NumPy functions:

  Riemann_zeta(s)  -- zeta(s) for scalar, list, or NumPy array input
  ln_gamma(z)      -- log(Gamma(z)) for scalar, list, or NumPy array

Precision of Riemann_zeta:

  |Im(s)| < 100      -- approximately 13 correct decimal digits
  |Im(s)| < 1,000    -- approximately 12 correct decimal digits
  |Im(s)| < 10,000   -- approximately 11 correct decimal digits

Requirements: Python 3.8 or later, NumPy. matplotlib is required
only for the plotting test in test.py.


Python_Fortran
=======================================

Python functions backed by a quadruple-precision Fortran core,
accessed via f2py. All computations are performed internally in
quadruple precision; results are returned as standard double-precision
complex NumPy arrays.

  Riemann_zeta(s)        -- zeta(s) for scalar or NumPy array input
  Riemann_zeta_prime(s)  -- returns (zeta(s), zeta'(s))
                            for scalar or NumPy array input
  ln_gamma(z)            -- log(Gamma(z)) for scalar or NumPy array

Accuracy: full double precision for |Im(s)| < 10^12.

Requirements: Python 3.8 or later, NumPy (provides f2py), gfortran.


REFERENCES
=======================================

[1] A. Kuznetsov, "Simple and accurate approximations to the Riemann
    zeta function", 2025, https://arxiv.org/abs/2503.09519

[2] A. Kuznetsov, "Computing the Barnes G-function and the gamma
    function in the entire complex plane", Journal of Computational
    and Applied Mathematics, Vol. 411, 2022, 114270.
    https://doi.org/10.1016/j.cam.2022.114270

[3] A. Kuznetsov, A. Mohammadioroojeh, "Approximating functions on
    R^+ by exponential sums", 2025, https://arxiv.org/abs/2508.19095
