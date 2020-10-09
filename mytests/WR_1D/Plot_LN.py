import numpy as np
import copy
import matplotlib.pyplot as plt

R_sun = 6.96e10
M_sun = 1.99e33
L_sun = 3.9e33
G = 6.67e-8
year = 365*24*60*60

# READING IN LUKA 1D HD

A = np.loadtxt('Output_Luka.blk',skiprows=3)

r_L = np.transpose(A)[0]
rho_L = np.transpose(A)[1]
v_L = np.transpose(A)[-4]/27
kappa_L = np.transpose(A)[-6]
T_L = np.transpose(A)[11]*1e5

x_L = 1 - r_L[0]/r_L

# plt.show()

A = np.loadtxt('test/2w_e1000.blk',skiprows=3)

r_SCD = np.transpose(A)[0]
rho_SCD = np.transpose(A)[1]
v_SCD = np.transpose(A)[2]
E_SCD = np.transpose(A)[3]
kappa_SCD = np.transpose(A)[5] + np.transpose(A)[6]
Gamma_SCD = np.transpose(A)[7]

unit_pressure = 1000000000
a_rad = 7.5646e-15
T_SCD = (E_SCD*unit_pressure/a_rad)**0.25

Mdot_SCD = 4*np.pi*rho_SCD*3.e-7*v_SCD*1e8*(r_SCD*R_sun)**2
Mdot_SCD = Mdot_SCD/M_sun*year

x_SCD = 1 - r_SCD[0]/r_SCD

#Compare Luka w FLD stable

plt.figure()
plt.title('velocity')
plt.plot(x_L,v_L,label='Luka')
plt.plot(x_SCD,v_SCD,label='Nico, CAK in D')
plt.ylabel('v [1d8 cm/s]')
plt.xlabel('x = 1-R/r')
plt.legend()

plt.figure()
plt.title('density')
plt.semilogy(x_L,rho_L,label='Luka')
plt.semilogy(x_SCD,rho_SCD,label='Nico, CAK in D')
plt.ylabel('rho [3d-7 cm/s]')
plt.xlabel('x = 1-R/r')
plt.legend()

plt.show()
