import numpy as np
import os
import matplotlib.pyplot as plt


f = os.getcwd() + '/FastWind_profile.txt'

r_fw,rho_fw,v_fw,pg_fw,T_fw,mu_fw,a2_fw = np.loadtxt(f,unpack=True,skiprows = 3)


f_relax1 = os.getcwd() + '/Amrvac_profile.okc'
A_relax1 = np.loadtxt(f_relax1,skiprows = 21)

x_rel, y_rel, z_rel, rho, mom, e, er, F, kappa, Gamma = np.transpose(A_relax1[0::256])

unit_length = 5075047950.3534832
unit_density = 8.8837999999999995e-9
unit_pressure = 128253.75816140750
unit_temperature = 106459.89999999999
unit_radflux = 487310493625.60626
unit_opacity = 2.217e-2
unit_time = 1335.68
unit_velocity = unit_length/unit_time

ny = len(y_rel)


rho = rho*unit_density
mom = mom*unit_density*unit_velocity
e = e*unit_pressure
er = er*unit_pressure
F = F*unit_radflux
kappa = kappa*unit_opacity

v = mom/rho
p = (5./3.-1.)*(e-0.5*v*v*rho)

mu = 0.6
mp = 1.6e-24
kb = 1.38e-16
Tg = p/rho * (mu*mp)/kb

R_star = 1252700000000.0
y_rel = y_rel*unit_length + R_star

dy = y_rel[1] - y_rel[0]

M_sun = 1.99e33
M_star = 66.8*M_sun

c = 2.99e10
G = 6.67e-8

grav = G*M_star/R_star**2
Gamma_e = 0.34*F/(grav*c)


Gamma = kappa*F/(c*grav)

tau = np.zeros(ny+1)
tau_0 = 0.34*p/abs(grav*(Gamma_e-1))



tau[-1] = tau_0[-1]
for i in range(ny-1,-1,-1):
     tau[i] = tau[i+1] + kappa[i]*rho[i]*dy
tau = tau[:-1]

print(tau_0[-1])

arad = 7.5646e-15
Trad = (er/arad)**0.25

delta_T = (Tg-Trad)/Tg


plt.figure()
plt.plot(y_rel,delta_T)
plt.ylabel('Detla_T')


plt.figure()
plt.plot(y_rel,Gamma,'o')
plt.plot(y_rel,Gamma_e,'o')
plt.ylabel('Gamma')
# plt.ylim([rho[0], rho[-1]])

plt.figure()
plt.semilogy(y_rel,rho,'o')
plt.semilogy(r_fw,rho_fw,'o')
plt.ylabel('density')
# plt.ylim([rho[0], rho[-1]])
plt.xlim([y_rel[0], y_rel[-1]])

plt.figure()
plt.semilogy(y_rel,p,'o')
plt.semilogy(r_fw,pg_fw,'o')
plt.ylabel('gas pressure')
# plt.ylim([p[0], p[-1]])
plt.xlim([y_rel[0], y_rel[-1]])


plt.figure()
plt.semilogy(y_rel,Tg,'o')
plt.semilogy(y_rel,Trad,'o')
plt.semilogy(r_fw,T_fw,'o')
plt.ylabel('gas temperature')
# plt.ylim([Tg[0], Tg[-1]])
plt.xlim([y_rel[0], y_rel[-1]])


plt.figure()
plt.plot(y_rel,tau,'o')
plt.ylabel('optical depth')

plt.figure()
plt.plot(y_rel,kappa,'o')
plt.ylabel('opacity')

plt.show()
