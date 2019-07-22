import numpy as np
import os
import matplotlib.pyplot as plt

f = os.getcwd() + '/init_struc_amrvac'

i, r, v, rho, Erad = np.loadtxt(f,unpack=True,skiprows=2)

Rsun = 6.96e10
Rstar = r[0]


arad = 7.5646e-15
kb = 1.380658e-16
mu = 0.6
mp = 1.6726231e-24
gamma = 1.66666667
c = 2.99792458e10
kappa = 0.33369889563715022
sigma = 5.67051e-5

mom = rho*v
T = (Erad/arad)**0.25
p = kb/(mp*mu)*T*rho
e = p/(gamma - 1) + rho*v*v/2
dEdr = np.zeros(len(r))
dEdr[1:-1] = (Erad[2:]-Erad[:-2])/(r[2:]-r[:-2])
dEdr[0] = (Erad[1]-Erad[0])/(r[1]-r[0])
dEdr[-1] = (Erad[-1]-Erad[-2])/(r[-1]-r[-2])
# dEdr.append(dEdr[-1])
F = -c/(3*kappa*rho)*dEdr




#> HEATING
heating = c*kappa*rho*Erad
cooling = 4*sigma*kappa*rho*T**4

f, (ax1,ax2,ax3,ax4) = plt.subplots(4,1,sharex = True)
ax1.semilogy(r/Rstar,rho,label='$\\rho$')
ax1.set_ylabel('$\\rho$')
ax1.legend(loc = 1)

ax2.semilogy(r/Rstar,mom,label='$\\rho \\vec{v}$')
ax2.semilogy(r/Rstar,mom*(r/Rstar)**2,label='$\\rho \\vec{v}r^2$')
ax2.set_ylabel('$\\vec{v} \\rho$')
ax2.legend(loc = 1)

ax3.plot(r/Rstar,F,'.',label='$\\vec{F}$')
ax3.plot(r/Rstar,F*(r/Rstar)**2,'.',label='$\\vec{F}r^2$')
ax3.set_ylabel('$\\vec{F}$')
ax3.legend(loc = 1)

ax4.semilogy(r/Rstar,Erad,'.',label='$E_r$')
ax4.semilogy(r/Rstar,e,'.',label='$e_g$')
ax4.set_ylabel('$E_r$')
ax4.legend(loc = 1)

ax4.set_xlabel('r/Rstar')

plt.figure()
plt.plot(r,heating,'x')
plt.plot(r,cooling,'+')

plt.show()
