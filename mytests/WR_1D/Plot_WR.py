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

# A = np.loadtxt('stable/2test0070.blk',skiprows=3)
A = np.loadtxt('stable/2Cak_in_diff0050.blk',skiprows=3)

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

################################################################################

A = np.loadtxt('stable/2Cak_notin_diff0050.blk',skiprows=3)

r_S = np.transpose(A)[0]
rho_S = np.transpose(A)[1]
v_S = np.transpose(A)[2]
E_S = np.transpose(A)[3]
kappa_S = np.transpose(A)[5] + np.transpose(A)[6]
Gamma_S = np.transpose(A)[7]

unit_pressure = 1000000000
a_rad = 7.5646e-15
T_S = (E_S*unit_pressure/a_rad)**0.25

Mdot_S = 4*np.pi*rho_S*3.e-7*v_S*1e8*(r_S*R_sun)**2
Mdot_S = Mdot_S/M_sun*year

x_S = 1 - r_S[0]/r_S

################################################################################

A = np.loadtxt('unstable/2Cak_notin_diff0050.blk',skiprows=3)

r_U = np.transpose(A)[0]
rho_U = np.transpose(A)[1]
v_U = np.transpose(A)[2]
E_U = np.transpose(A)[3]
kappa_U = np.transpose(A)[5] + np.transpose(A)[6]
Gamma_U = np.transpose(A)[7]

unit_pressure = 1000000000
a_rad = 7.5646e-15
T_U = (E_U*unit_pressure/a_rad)**0.25

Mdot_U = 4*np.pi*rho_U*3.e-7*v_U*1e8*(r_U*R_sun)**2
Mdot_U = Mdot_U/M_sun*year

x_U = 1 - r_U[0]/r_U

################################################################################

A = np.loadtxt('unstable/2Cak_in_diff0050.blk',skiprows=3)

r_UCD = np.transpose(A)[0]
rho_UCD = np.transpose(A)[1]
v_UCD = np.transpose(A)[2]
E_UCD = np.transpose(A)[3]
kappa_UCD = np.transpose(A)[5] + np.transpose(A)[6]
Gamma_UCD = np.transpose(A)[7]

unit_pressure = 1000000000
a_rad = 7.5646e-15
T_UCD = (E_UCD*unit_pressure/a_rad)**0.25

Mdot_UCD = 4*np.pi*rho_UCD*3.e-7*v_UCD*1e8*(r_UCD*R_sun)**2
Mdot_UCD = Mdot_UCD/M_sun*year

x_UCD = 1 - r_UCD[0]/r_UCD


#Compare Luka w FLD stable

plt.figure()
plt.title('velocity')
plt.plot(x_L,v_L,label='Luka')
plt.plot(x_S,v_S,label='Nico')
plt.plot(x_SCD,v_SCD,label='Nico, CAK in D')
plt.ylabel('v [1d8 cm/s]')
plt.xlabel('x = 1-R/r')
plt.legend()

plt.figure()
plt.title('density')
plt.semilogy(x_L,rho_L,label='Luka')
plt.semilogy(x_S,rho_S,label='Nico')
plt.semilogy(x_SCD,rho_SCD,label='Nico, CAK in D')
plt.ylabel('rho [3d-7 cm/s]')
plt.xlabel('x = 1-R/r')
plt.legend()

plt.figure()
plt.title('Trad')
plt.plot(x_L,T_L,label='Luka')
plt.plot(x_S,T_S,label='Nico')
plt.plot(x_SCD,T_SCD,label='Nico, CAK in D')
plt.ylabel('T [K]')
plt.xlabel('x = 1-R/r')
plt.legend()

plt.figure()
plt.title('Mdot')
plt.hlines(2.06e-5,xmin=0,xmax=1,label='Luka')
plt.semilogy(x_S,Mdot_S,label='Nico')
plt.semilogy(x_SCD,Mdot_SCD,label='Nico, CAK in D')
plt.ylabel('Mdot [Msun/year]')
plt.xlabel('x = 1-R/r')
plt.legend()

#Compare stable w unstable fld

plt.figure()
plt.title('velocity')
plt.plot(x_U,v_U,label='unstable')
plt.plot(x_S,v_S,label='stable')
plt.plot(x_UCD,v_UCD,label='unstable, CAK in D')
plt.plot(x_SCD,v_SCD,label='stable, CAK in D')
plt.ylabel('v [1d8 cm/s]')
plt.xlabel('x = 1-R/r')
plt.legend()

plt.figure()
plt.title('density')
plt.semilogy(x_U,rho_U,label='unstable')
plt.semilogy(x_S,rho_S,label='stable')
plt.semilogy(x_UCD,rho_UCD,label='unstable, CAK in D')
plt.semilogy(x_SCD,rho_SCD,label='stable, CAK in D')
plt.ylabel('rho [3d-7 cm/s]')
plt.xlabel('x = 1-R/r')
plt.legend()

plt.figure()
plt.title('Trad')
plt.plot(x_U,T_U,label='unstable')
plt.plot(x_S,T_S,label='stable')
plt.plot(x_UCD,T_UCD,label='unstable, CAK in D')
plt.plot(x_SCD,T_SCD,label='stable, CAK in D')
plt.ylabel('T [K]')
plt.xlabel('x = 1-R/r')
plt.legend()

plt.figure()
plt.title('Mdot')
# plt.hlines(2.06e-5,xmin=0,xmax=1,label='Luka')
plt.semilogy(x_S,Mdot_S,label='stable')
plt.semilogy(x_U,Mdot_U,label='unstable')
plt.semilogy(x_UCD,Mdot_UCD,label='unstable, CAK in D')
plt.semilogy(x_SCD,Mdot_SCD,label='stable, CAK in D')
plt.ylabel('Mdot [Msun/year]')
plt.xlabel('x = 1-R/r')
plt.legend()

plt.figure()
plt.title('Wolf-Rayet wind')
plt.plot(r_SCD,v_SCD,'r-',label='stable')
plt.xlabel('$r/R_{sun}$',fontsize=15)
plt.ylabel('$v [10^8 cm/s]$',fontsize=15)


plt.show()
