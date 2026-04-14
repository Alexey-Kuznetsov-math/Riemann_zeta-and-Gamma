"""
Riemann zeta function and its derivative via quad-precision Fortran.

Functions
---------
Riemann_zeta(s)       -> f         : zeta(s)
Riemann_zeta_prime(s) -> (f1, f2)  : zeta(s), zeta'(s)

Both accept a scalar or array-like input and return the same shape.
"""

import numpy as np
from _Riemann_zeta_python import riemann_zeta_wrap, riemann_zeta_prime_wrap


def Riemann_zeta(s):
    """Compute the Riemann zeta function zeta(s).

    Parameters
    ----------
    s : complex scalar or array_like

    Returns
    -------
    f : complex scalar or ndarray (same shape as s)
    """
    s_arr = np.atleast_1d(np.asarray(s, dtype=np.complex128))
    shape = s_arr.shape
    f = riemann_zeta_wrap(s_arr.ravel())
    f = f.reshape(shape)
    if np.ndim(s) == 0:
        return f.item()
    return f


def Riemann_zeta_prime(s):
    """Compute the Riemann zeta function and its derivative.

    Parameters
    ----------
    s : complex scalar or array_like

    Returns
    -------
    f1 : complex scalar or ndarray – zeta(s)
    f2 : complex scalar or ndarray – zeta'(s)
    """
    s_arr = np.atleast_1d(np.asarray(s, dtype=np.complex128))
    shape = s_arr.shape
    f1, f2 = riemann_zeta_prime_wrap(s_arr.ravel())
    f1 = f1.reshape(shape)
    f2 = f2.reshape(shape)
    if np.ndim(s) == 0:
        return f1.item(), f2.item()
    return f1, f2
