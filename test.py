import numpy as np
from Riemann_zeta import Riemann_zeta
from ln_gamma import ln_gamma
import matplotlib.pyplot as plt

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Test 1 : compute Gamma(k) for integer values of k')
f=np.zeros(10,dtype=np.complex128)
g=np.zeros(10,dtype=np.complex128)
for k in range(1,11):
    z=complex(k,0)
    f[k-1]=np.exp(ln_gamma(z))
    print(f[k-1])
g_vals=np.array([1,1,2,6,24,120,720,5040,40320,362880],dtype=np.float64)
g=g_vals.astype(np.complex128)
print('max relative error=',np.max(np.abs(f/g-1)))

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Test 2 : compute Gamma(1/2+k) for integer values of k')
f=np.zeros(10,dtype=np.complex128)
for k in range(1,6):
    z=k-0.5
    f[k-1]=np.exp(ln_gamma(z))
    print(f[k-1])
g=np.zeros(10,dtype=np.complex128)
g[:5]=np.sqrt(np.pi)*np.array([1.0,0.5,0.75,1.875,6.5625],dtype=np.float64)
print('max relative error=',np.max(np.abs(f[:5]/g[:5]-1)))

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Test 3 : compute |log(Gamma(z+1))-log(Gamma(z))-log(z)| for 10000 random values of z, uniformly distributed in the square |Im(z)|<20, |Re(z)|<20')
error=0.0
rng=np.random.default_rng(123456)
for _ in range(10000):
    U=rng.random(2)        # uniform in [0,1)
    z=20*complex(2*U[0]-1,2*U[1]-1)
    err_abs=abs(ln_gamma(z+1)-ln_gamma(z)-np.log(z))
    if err_abs>error:
        error=err_abs
print('maximum error=',error)

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Test 4 : compute zeta(s) for s=2,4,6,8,10')
f=np.zeros(10,dtype=np.complex128)
for k in range(1,6):
    s=2*k
    f[k-1]=Riemann_zeta(s)
    print(f[k-1])
g=np.zeros(10,dtype=np.complex128)
g[0]=np.pi**2/6
g[1]=np.pi**4/90
g[2]=np.pi**6/945
g[3]=np.pi**8/9450
g[4]=np.pi**10/93555
print('maximum error=',np.max(np.abs(f[:4]-g[:4])))

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Test 5 : compute zeta(s) for s=-1,-3,-5,-7,-9')
f=np.zeros(10,dtype=np.complex128)
for k in range(1,6):
    s=1-2*k
    f[k-1]=Riemann_zeta(s)
    print(f[k-1])
g=np.zeros(10,dtype=np.complex128)
g[0]=-1.0/12
g[1]=1.0/120
g[2]=-1.0/252
g[3]=1.0/240
g[4]=-1.0/132
print('maximum error=',np.max(np.abs(f[:5]-g[:5])))

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Test 5 : compute zeta(1/2+i*t_k) for the first five non-trivial zeros 1/2+i*t_k of zeta(s)')
t=np.zeros(5,dtype=np.float64)
t[0]=14.13472514173469379045725198356247027
t[1]=21.02203963877155499262847959389690277
t[2]=25.01085758014568876321379099256282181
t[3]=30.42487612585951321031189753058409132
t[4]=32.93506158773918969066236896407490348
for k in range(5):
    s=0.5+1j*t[k]
    print(Riemann_zeta(s))

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Test 6 : compute zeta(1/2+i*t_k) for the first five non-trivial zeros 1/2+i*t_k of zeta(s) at height t=200')
t[0]=201.2647519437037887330161334275481732
t[1]=202.4935945141405342776866606378643158
t[2]=204.1896718031045543307164383863136851
t[3]=205.3946972021632860252123793906930909
t[4]=207.9062588878062098615019679077536442
for k in range(5):
    s=0.5+1j*t[k]
    print(Riemann_zeta(s))

print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Test 7 : plot the Riemann-Siegel function Z(t) for 0<t<100, 1000<t<1100 and 10000<t<10100')
t=np.linspace(0,100,10000)
g=ln_gamma(0.25+0.5j*t)
f0=Riemann_zeta(0.5+1j*t)*np.exp(1j*g.imag-0.5j*t*np.log(np.pi))

t1=t+1000
g1=ln_gamma(0.25+0.5j*t1)
f1=Riemann_zeta(0.5+1j*t1)*np.exp(1j*g1.imag-0.5j*t1*np.log(np.pi))

t2=t+10000
g2=ln_gamma(0.25+0.5j*t2)
f2=Riemann_zeta(0.5+1j*t2)*np.exp(1j*g2.imag-0.5j*t2*np.log(np.pi))

fig,axes=plt.subplots(3,1,figsize=(8,12))  
axes[0].plot(t,f0.real)
axes[0].set_title('Z(t) for 0<t<100')
axes[0].grid(True)
axes[1].plot(t1,f1.real)
axes[1].set_title('Z(t) for 1000<t<1100')
axes[1].grid(True)
axes[2].plot(t2,f2.real)
axes[2].set_title('Z(t) for 10000<t<10100')
axes[2].grid(True)
plt.tight_layout()
plt.show()

