import numpy as np
import os
import matplotlib.pyplot as plt

f_relax1 = os.getcwd() + '/Amrvac_profile1.okc'
A_relax1 = np.loadtxt(f_relax1,skiprows = 21)

x_rel, y_rel, z_rel, rho1, mom1, e1, er1, F1, kappa1, Gamma1 = np.transpose(A_relax1[0::256])


f_relax2 = os.getcwd() + '/Amrvac_profile2.okc'
A_relax2 = np.loadtxt(f_relax2,skiprows = 21)

x_rel, y_rel, z_rel, rho2, mom2, e2, er2, F2, kappa2, Gamma2 = np.transpose(A_relax2[0::256])

unit_length = 5075047950.3534832
unit_density = 8.8837999999999995e-9
unit_pressure = 128253.75816140750
unit_temperature = 106459.89999999999
unit_radflux = 487310493625.60626
unit_opacity = 2.217e-2
unit_time = 1335.68
unit_velocity = unit_length/unit_time

ny = len(y_rel)


rho1 = rho1*unit_density
mom1 = mom1*unit_density*unit_velocity
e1 = e1*unit_pressure
er1 = er1*unit_pressure
F1 = F1*unit_radflux
kappa1 = kappa1*unit_opacity

rho2 = rho2*unit_density
mom2 = mom2*unit_density*unit_velocity
e2 = e2*unit_pressure
er2 = er2*unit_pressure
F2 = F2*unit_radflux
kappa2 = kappa2*unit_opacity

v1 = mom1/rho1
p1 = (5./3.-1.)*(e1-0.5*v1*v1*rho1)

v2 = mom2/rho2
p2 = (5./3.-1.)*(e2-0.5*v2*v2*rho2)

mu = 0.6
mp = 1.6e-24
kb = 1.38e-16
Tg1 = p1/rho1 * (mu*mp)/kb
Tg2 = p2/rho2 * (mu*mp)/kb

R_star = 1252700000000.0
y_rel = y_rel*unit_length + R_star

dy = y_rel[1] - y_rel[0]

M_sun = 1.99e33
M_star = 66.8*M_sun

c = 2.99e10
G = 6.67e-8

grav = G*M_star/R_star**2
Gamma_e1 = 0.34*F1/(grav*c)
Gamma_e2 = 0.34*F2/(grav*c)

Gamma1 = kappa1*F1/(c*grav)
Gamma2 = kappa2*F2/(c*grav)

f,(ax1,ax2) = plt.subplots(2,1,sharex=True,figsize=(10,5))

rsun = 6.96e10
y_relsun = y_rel/rsun

ax1.set_title('1D profile stellar atmosphere $70M_\odot$', fontsize = 15)

ax1.plot(y_relsun,rho1, label = 'initial')
ax1.plot(y_relsun,rho2, label = 'final')
ax1.set_ylabel('$\\rho [g cm^{-3}]$', fontsize = 15)
ax1.legend()

ax2.plot(y_relsun,er1)
ax2.plot(y_relsun,er2)
ax2.set_ylabel('$E_r [erg cm^{-3}]$', fontsize = 15)
ax2.set_xlabel('r [$R_\odot$]', fontsize = 15)

plt.savefig('1Dprof')
plt.show()
