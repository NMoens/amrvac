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

R_star = 1252700000000.0
y_rel = y_rel*unit_length + R_star

dy = y_rel[1] - y_rel[0]

M_sun = 1.99e33
M_star = 66.8*M_sun

c = 2.99e10
G = 6.67e-8

grav = G*M_star/R_star**2
Gamma_e = 0.34*F/(grav*c)

tau = np.zeros(ny+1)
tau_0 = p/abs(grav*(Gamma-1))

tau[-1] = tau_0[-1]
for i in range(ny-1,-1,-1):
     tau[i] = tau[i+1] + kappa[i]*rho[i]*dy
tau = tau[:-1]

# plt.semilogy(r_fw,rho_fw,'o')
plt.plot(y_rel,tau,'o')
plt.show()
