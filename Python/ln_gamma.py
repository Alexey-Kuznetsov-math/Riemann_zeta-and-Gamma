import numpy as np
###########################################################################
def ln_gamma(z):
# ln_gamma  Computes the logarithm of the Gamma function in the entire complex plane.
#
#   f = ln_gamma(z) returns log(Gamma(z)) for any complex input z.
#   The input z can be a scalar, a list, or a NumPy array.
# -------------------------------------------------------------------------
# Author: Alexey Kuznetsov
# York University, Toronto, Canada
# Website: https://kuznetsovmath.ca/
# Email: akuznets@yorku.ca
#
# Created: 29-Nov-2025
# Last updated: 29-Nov-2025
#
# License: BSD 3-Clause (https://opensource.org/licenses/BSD-3-Clause)
# -------------------------------------------------------------------------
    scalar_input = np.isscalar(z)
    z = np.asarray(z, dtype=np.complex128)
    f = np.zeros_like(z, dtype=np.complex128)

    i1=z.real>=1.5
    i2=(z.real>=0.5)&(z.real<1.5)
    i3=(z.real<0.5)&(z.imag>=0)
    i4=(z.real<0.5)&(z.imag<0)

    if np.any(i1):
        f[i1]=ln_gamma_half_plane(z[i1])

    # use functional equation log(Gamma(z))=log(Gamma(z+1))-log(z)
    if np.any(i2):
        f[i2]=ln_gamma_half_plane(z[i2]+1)-np.log(z[i2])

    # use reflection formula, Im(z) >= 0
    if np.any(i3):
        w3=z[i3]
        f[i3]=np.log(1-w3)-ln_gamma_half_plane(2-w3)+1.83787706640935+1j*np.pi*(w3-0.5)-np.log(1-np.exp(2j*np.pi*w3))

    # use reflection formula, Im(z) < 0
    if np.any(i4):
        w4=z[i4]
        f[i4]=np.log(1-w4)-ln_gamma_half_plane(2-w4)+1.83787706640935-1j*np.pi*(w4-0.5)-np.log(1-np.exp(-2j*np.pi*w4))
        
    if scalar_input:
        return f.item()
    return f
###########################################################################
def ln_gamma_half_plane(z):
#Approximates log(Gamma(z)) in the half-plane Re(z) >= 1.5.
    z=np.asarray(z, dtype=np.complex128)
    c=np.array([
        -5.3035486658210424e-9 + 8.47546314707099368e-10j,
         1.0093321787093660e-6 + 2.72761804689954947e-7j,
        -2.9946104278440048e-5 + 3.68867702224184848e-5j,
         5.1324321142829496e-4 + 1.31127415401525033e-3j
    ], dtype=np.complex128)
    l=np.array([
         3.036399609619316394 + 2.183538146270238962j,
         2.619619673091566707 + 1.285513733113457465j,
         2.236468161644125811 - 6.451943678012780725e-1j,
         1.822691064348029853 + 1.260996493971840621e-1j
    ], dtype=np.complex128)
    cr=np.array([
         1.91267259422448627e-12,
        -1.34922683153592188e-3,
        -8.28182554958193984e-4,
        -1.80081775867242537e-4
    ], dtype=np.float64)
    lr=np.array([
        7.862008471885223854,
        1.278757134063109466,
        1.129114403699685475,
        1.037917494906115128
    ], dtype=np.float64)
    f=np.zeros_like(z,dtype=np.complex128)
    for j in range(4):
        f=f+c[j]/(z-1+l[j])**3+np.conj(c[j])/(z-1+np.conj(l[j]))**3+cr[j]/(z-1+lr[j])**3

    f=2*f+0.918938533204673+(z-0.5)*np.log(z)-z+0.083333333333333/z
    return f



