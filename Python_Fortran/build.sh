#!/usr/bin/env bash
set -e

python3 -m pip install --user meson ninja
export PATH="$HOME/.local/bin:$PATH"

python3 -m numpy.f2py -c -m _pyzeta_fortran \
    Riemann_zeta_module.f90 zeta_py_wrapper.f90 \
    only: riemann_zeta_scalar riemann_zeta_vec :
