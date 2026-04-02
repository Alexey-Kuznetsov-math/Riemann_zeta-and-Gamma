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
# Last updated: 2-April-2026
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
        f[i3]=np.log(1-w3)-ln_gamma_half_plane(2-w3)+1.83787706640934548+1j*np.pi*(w3-0.5)-np.log(1-np.exp(2j*np.pi*w3))

    # use reflection formula, Im(z) < 0
    if np.any(i4):
        w4=z[i4]
        f[i4]=np.log(1-w4)-ln_gamma_half_plane(2-w4)+1.83787706640934548-1j*np.pi*(w4-0.5)-np.log(1-np.exp(-2j*np.pi*w4))
        
    if scalar_input:
        return f.item()
    return f
###########################################################################
def ln_gamma_half_plane(z):
#Approximates log(Gamma(z)) in the half-plane Re(z) >= 1.5.
    z=np.asarray(z, dtype=np.complex128)
    c=np.array([
        7.8489806116985789e-8+1.69933072157308014e-8j,
       -2.0915898156004773e-5-1.57315657837253062e-5j,
        8.0382728578667684e-4+1.46123728240593489e-3j,
       -1.7781861010073780e-2-3.91857202978460444e-2j
    ], dtype=np.complex128)
    l=np.array([
       3.089224411647167808+2.346905763041727244j,
       2.676547234010318606+1.370743226374250769j,
       2.266206347366603493+6.752681709263032345e-1j,
       1.829497420906493754+1.564512477809790312e-1j
    ], dtype=np.complex128)
    cr=np.array([
        -3.32926108478200017e-6,
         4.79008236266177865e-2,
         4.62832185218743242e-2,
         2.31503627111999874e-2
    ], dtype=np.float64)
    lr=np.array([
        3.754580109023743872,
        1.253693011087479497,  
        1.099995202351846400, 
        1.018671901578566989
    ], dtype=np.float64)
    f=np.zeros_like(z,dtype=np.complex128)
    for j in range(4):
        f=f+c[j]/(z-1+l[j])+np.conj(c[j])/(z-1+np.conj(l[j]))+cr[j]/(z-1+lr[j])

    f=f+0.91893853320467274+(z-0.5)*np.log(z)-z
    return f

