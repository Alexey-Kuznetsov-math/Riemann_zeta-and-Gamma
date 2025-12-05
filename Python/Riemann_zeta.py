import numpy as np
from ln_gamma import ln_gamma
###########################################################################
def Riemann_zeta(s):
# Riemann_zeta approximates the Riemann zeta function in the entire complex plane
#
#   f = Riemann_zeta(s) returns an approximation to zeta(s) for complex input s
#   The input s can be a scalar, a list, or a NumPy array.
#
# This function requires ln_gamma function (which computes the logarithm of the Gamma function 
# in the entire complex plane)
# -------------------------------------------------------------------------
# Author: Alexey Kuznetsov
# York University, Toronto, Canada
# Website: https://kuznetsovmath.ca/
# Email: akuznets@yorku.ca
#
# Created: 29-Nov-2025
# Last updated: 5-Dec-2025
#
# License: BSD 3-Clause (https://opensource.org/licenses/BSD-3-Clause)
#--------------------------------------------------------------------------
    scalar_input=np.isscalar(s)
    s=np.asarray(s,dtype=np.complex128)
    f=np.zeros_like(s,dtype=np.complex128)

    i1=s.real>=0.5
    i2=(s.real<0.5)&(s.imag>=0)
    i3=(s.real<0.5)&(s.imag<0)

    if np.any(i1):
        f[i1]=Riemann_zeta_half_plane(s[i1])

    if np.any(i2):
        w2=s[i2]
        f[i2]=1j*(2*np.pi)**(w2-1)*(1-np.exp(1j*np.pi*w2))*np.exp(-0.5j*np.pi*w2+ln_gamma(1-w2))*Riemann_zeta_half_plane(1-w2)

    if np.any(i3):
        w3=s[i3]
        f[i3]=-1j*(2*np.pi)**(w3-1)*(1-np.exp(-1j*np.pi*w3))*np.exp(0.5j*np.pi*w3+ln_gamma(1-w3))*Riemann_zeta_half_plane(1-w3)

    if scalar_input:
        return f.item()
    return f
###########################################################################
def Riemann_zeta_half_plane(s):
# Approximate zeta(s) in the half-plane Re(s)>0.5
    s=np.asarray(s,dtype=np.complex128)
    f=np.zeros_like(s,dtype=np.complex128)

    i1=s.real>=5
    i2=(s.real<5)&(np.abs(s.imag)<=200)
    i3=(s.real<5)&(s.imag>200)
    i4=(s.real<5)&(s.imag<-200)

    if np.any(i1):
        f[i1]=zeta_summation(s[i1])

    if np.any(i2):
        f[i2]=zeta_Euler_Maclaurin(s[i2])

    if np.any(i3):
        f[i3]=zeta_8(s[i3])

    if np.any(i4):
        f[i4]=np.conj(zeta_8(np.conj(s[i4])))

    return f
###########################################################################
def zeta_8(s):
# zeta_8 computes the approximation to the Riemann zeta function described in 
# [1] A. Kuznetsov, "Simple and accurate approximations to the Riemann zeta function", 2025, https://arxiv.org/abs/2503.09519
#--------------------------------------------------------------------------
# the cofficients lambda_j for j=1,2,...,8
    lam=np.array([
        0.152845417613666702426-0.119440685603870510384j,
        0.302346225128945757427-0.243989695504400621268j,
        0.451119584531782942888-0.378479770209444563858j,
        0.604563710297226464637-0.523486888629095259770j,
        0.765965706759629396959-0.678405572413543444272j,
        0.938371150977889047740-0.845332361280975174880j,
        1.128148837845288402558-1.030737947568157685685j,
        1.353030558654668162533-1.252503278108132307164j
    ],dtype=np.complex128)
#  the cofficients omega_j for j=0,1,2,...,8    
    omega0=1.926019633029103199063e-1+2.472986965795651842299e-2j
    omega=np.array([
        1.582954327321094104502e-1+4.149113569204600502105e-2j,
        7.826728293587305110862e-2+5.215518667623989653254e-2j,
        1.940595049247490540621e-2+2.977286598777633378610e-2j,
        1.691184771902755036966e-3+8.938933548999206800196e-3j,
       -2.994777986686168319731e-4+1.567541981830224487301e-3j,
       -9.837202592542590210980e-5+1.502108057352792742070e-4j,
       -9.346989286415688998740e-6+5.793852209955845432028e-6j,
       -2.451577304299235983015e-7+6.134784898751456953524e-9j
    ],dtype=np.complex128)
#--------------------------------------------------------------------------
    s=np.asarray(s,dtype=np.complex128)
# compute chi(s)=(2*pi)^s/(2*cos(pi*s/2)*gamma(s))
    chi=np.exp(s*np.log(2*np.pi)+0.5j*np.pi*s-np.log(1+np.exp(np.pi*1j*s))-ln_gamma(s))

    N=np.floor(np.sqrt(s.imag/(2*np.pi)))
    N_int=N.astype(int)
    f=1+chi
    for n in range(2,np.max(N_int)+1):
        u=n**(-s)
        f+=(n<=N_int)*(u+chi/(n*u))

    M=N_int+0.5
    s1=1-np.conj(s)
    I1=np.full_like(s,omega0,dtype=np.complex128)
    I2=I1.copy()
    for j in range(8):
        I1+=omega[j]*(np.exp(-2*np.pi*M*lam[j]-s*np.log(1+1j*lam[j]/M))+np.exp(2*np.pi*M*lam[j]-s*np.log(1-1j*lam[j]/M)))
        I2+=omega[j]*(np.exp(-2*np.pi*M*lam[j]-s1*np.log(1+1j*lam[j]/M))+np.exp(2*np.pi*M*lam[j]-s1*np.log(1-1j*lam[j]/M)))

    I1=np.exp(-s*np.log(M))*I1
    I2=np.conj(np.exp(-s1*np.log(M))*I2)
    f=f-0.5*(-1)**N_int*(I1+chi*I2)
    return f
###########################################################################
def zeta_summation(s):
# computes zeta(s)=\sum_{n=1}^{\infty} n^{-s}  
# We remove from this sum all numbers divisible by 2,3,5,7,11,13,17,19 and truncate the resulting sum at n=500
# This results in a more efficient algorithm (fewer terms in the main sum)
# This function gives a good approximation to zeta(s) in the half-plane Re(s)>5
#--------------------------------------------------------------------------
    s=np.asarray(s,dtype=np.complex128)
    n=np.array([
       23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131,137,139,
       149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263,
       269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,
       401,409,419,421,431,433,439,443,449,457,461,463,467,479,487,491,499
    ],dtype=float)
    f=np.ones_like(s,dtype=np.complex128)
    for k in n:
        f+=k**(-s)
        
    f=f/((1-2**(-s))*(1-3**(-s))*(1-5**(-s))*(1-7**(-s))*(1-11**(-s))*(1-13**(-s))*(1-17**(-s))*(1-19**(-s)))
    return f
###########################################################################
def zeta_Euler_Maclaurin(s):
# computes zeta(s) via Euler-Maclaurin summation 
# we use formula (25.2.10) at https://dlmf.nist.gov/25.2 with N=100 and n=15
#--------------------------------------------------------------------------
    s=np.asarray(s,dtype=np.complex128)
# precompute b(k)=B_{2k}/(2k)!, where B_{k} are Bernoulli numbers    
    b=np.array([
        8.33333333333333e-2,-1.38888888888889e-3,3.306878306878306e-5,-8.267195767195768e-7,2.087675698786810e-8,
        -5.284190138687493e-10,1.338253653068468e-11,-3.389680296322583e-13,8.586062056277845e-15,-2.174868698558061e-16, 
        5.509002828360229e-18, -1.395446468581252e-19, 3.534707039629467e-21,-8.953517427037546e-23,2.267952452337683e-24
    ],dtype=float)
    N=100
    m=s/N
    f=np.zeros_like(s,dtype=np.complex128)
    for k in range(1,16):
        f+=b[k-1]*m
        m=m*(s+2*k-1)*(s+2*k)/N**2
        
    f=1+N**(-s)*(f+0.5+N/(s-1))
    for k in range(2,N):
        f+=k**(-s)
        
    return f




