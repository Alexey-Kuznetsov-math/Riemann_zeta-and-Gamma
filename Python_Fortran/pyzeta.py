import numpy as np
import _pyzeta_fortran


def Riemann_zeta(s):
    if np.isscalar(s):
        return complex(_pyzeta_fortran.zeta_py_wrapper.riemann_zeta_scalar(np.complex128(s)))

    a = np.asarray(s, dtype=np.complex128)
    b = _pyzeta_fortran.zeta_py_wrapper.riemann_zeta_vec(np.ravel(a))
    return np.asarray(b, dtype=np.complex128).reshape(a.shape)
