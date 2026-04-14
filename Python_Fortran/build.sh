#!/bin/bash
# build.sh - Compile the Riemann zeta Python extension via f2py
#
# Usage:   ./build.sh
# Requires: gfortran, numpy (which provides f2py)
#
# This produces _Riemann_zeta_python.cpython-*.so in the current directory.
# Then just:  from Riemann_zeta_python import Riemann_zeta, Riemann_zeta_prime

set -e

# Step 1: compile the Fortran module to .o (also produces the .mod file)
gfortran -O3 -c Riemann_zeta_module.f90

# Step 2: f2py wraps only the wrapper file, links against the .o
python3 -m numpy.f2py -c \
    --f90flags='-O3' \
    -m _Riemann_zeta_python \
    Riemann_zeta_wrapper.f90 \
    Riemann_zeta_module.o

echo ""
echo "Build successful. Quick test:"
python3 -c "
from Riemann_zeta_python import Riemann_zeta, Riemann_zeta_prime
import numpy as np

# scalar test
f = Riemann_zeta(2.0)
print(f'  zeta(2)   = {f:.15f}  (expected {np.pi**2/6:.15f})')

# derivative test
f1, f2 = Riemann_zeta_prime(2.0)
print(f'  zeta(2)   = {f1:.15f}')
print(f\"  zeta\'(2)  = {f2:.15f}\")

# vector test
s = np.array([2, 4, 6], dtype=complex)
f = Riemann_zeta(s)
print(f'  zeta([2,4,6]) = {f}')
"
