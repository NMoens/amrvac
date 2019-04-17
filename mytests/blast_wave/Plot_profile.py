import numpy as np
import matplotlib.pyplot as plt
import os

f = os.getcwd() + '/blast_wave_output.okc'
F = np.loadtxt(f,skiprows=21)
x,y,z,rho,m1,m2,e,E_rad,F1,F2 = np.transpose(F)

unit_length = 1.0000000000000000
unit_velocity = 11645.084295622544
unit_density = 2.3416704877999998E-024
unit_pressure = 3.1754922399999997E-016
unit_flux = 3.6978874814895254E-012

kb = 1.380658e-16
arad = 7.5646e-15
sigma = 5.67051e-5
mp = 1.6726231e-24
mu = 0.6
gamma = 1.6666

x = x*unit_length
y = y*unit_length
rho = rho*unit_density
m1 = m1*unit_density*unit_velocity
m2 = m2*unit_density*unit_velocity
e = e*unit_pressure
E_rad = E_rad*unit_pressure
F1 = F1*unit_flux
F2 = F2*unit_flux

n_cells = len(x)
r = np.sqrt(x**2+y**2)
phi = np.arctan2(y,x)
mr = np.sqrt(m1**2+m2**2)
Fr = np.sqrt(F1**2+F2**2)
p = (gamma -1)*(e-0.5*((m1**2+m2**2)/rho))
Tgas = p/rho*mp*mu/kb
Trad = (E_rad/arad)**0.25

# for i in range(n_cells):
#     if (abs(np.tan(phi[i]))<0.3 or 1./abs(np.tan(phi[i]))<0.3):
#         plt.plot(r[i],rho[i],'r.')
#     else:
#         plt.plot(r[i],rho[i],'b.')

f,(ax1,ax2,ax3,ax4,ax5) = plt.subplots(5,1,sharex=True)

ax1.plot(r,E_rad,'.')
ax2.plot(r,Fr,'.')
ax3.plot(r,p,'.')
ax4.plot(r,Tgas,'.')
ax5.plot(r,rho,'.')


ax1.set_ylabel('$E_{rad}$')
ax2.set_ylabel('$F_{rad}$')
ax3.set_ylabel('$p_{gas}$')
ax4.set_ylabel('$T_{rad}$')
ax5.set_ylabel('$\\rho_{rad}$')


plt.show()
