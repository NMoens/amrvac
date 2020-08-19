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

A = np.loadtxt('stable/2test0050.blk',skiprows=3)

r_S = np.transpose(A)[0]
rho_S = np.transpose(A)[1]
v_S = np.transpose(A)[2]
E_S = np.transpose(A)[3]
kappa_S = np.transpose(A)[5] + np.transpose(A)[6]
Gamma_S = np.transpose(A)[7]

unit_pressure = 3000000000
a_rad = 7.5646e-15
T_S = (E_S*unit_pressure/a_rad)**0.25

Mdot_S = 4*np.pi*rho_S*3.e-7*v_S*1e8*(r_S*R_sun)**2
Mdot_S = Mdot_S/M_sun*year

x_S = 1 - r_S[0]/r_S

A = np.loadtxt('unstable/2test0050.blk',skiprows=3)

r_U = np.transpose(A)[0]
rho_U = np.transpose(A)[1]
v_U = np.transpose(A)[2]
E_U = np.transpose(A)[3]
kappa_U = np.transpose(A)[5] + np.transpose(A)[6]
Gamma_U = np.transpose(A)[7]

unit_pressure = 3000000000
a_rad = 7.5646e-15
T_U = (E_U*unit_pressure/a_rad)**0.25

Mdot_U = 4*np.pi*rho_U*3.e-7*v_U*1e8*(r_S*R_sun)**2
Mdot_U = Mdot_U/M_sun*year

x_U = 1 - r_U[0]/r_U

#Compare Luka w FLD stable

plt.figure()
plt.title('velocity')
plt.plot(x_L,v_L,label='Luka')
plt.plot(x_S,v_S,label='Nico')
plt.ylabel('v [1d8 cm/s]')
plt.xlabel('x = 1-R/r')
plt.legend()

plt.figure()
plt.title('density')
plt.semilogy(x_L,rho_L,label='Luka')
plt.semilogy(x_S,rho_S,label='Nico')
plt.ylabel('rho [3d-7 cm/s]')
plt.xlabel('x = 1-R/r')
plt.legend()

plt.figure()
plt.title('Trad')
plt.plot(x_L,T_L,label='Luka')
plt.plot(x_S,T_S,label='Nico')
plt.ylabel('T [K]')
plt.xlabel('x = 1-R/r')
plt.legend()

plt.figure()
plt.title('Mdot')
plt.hlines(2.06e-5,xmin=0,xmax=1,label='Luka')
plt.semilogy(x_S,Mdot_S,label='Nico')
plt.ylabel('Mdot [Msun/year]')
plt.xlabel('x = 1-R/r')
plt.legend()

#Compare stable w unstable fld

plt.figure()
plt.title('velocity')
plt.plot(x_U,v_U,label='unstable')
plt.plot(x_S,v_S,label='stable')
plt.ylabel('v [1d8 cm/s]')
plt.xlabel('x = 1-R/r')
plt.legend()

plt.figure()
plt.title('density')
plt.semilogy(x_U,rho_U,label='unstable')
plt.semilogy(x_S,rho_S,label='stable')
plt.ylabel('rho [3d-7 cm/s]')
plt.xlabel('x = 1-R/r')
plt.legend()

plt.figure()
plt.title('Trad')
plt.plot(x_U,T_U,label='unstable')
plt.plot(x_S,T_S,label='stable')
plt.ylabel('T [K]')
plt.xlabel('x = 1-R/r')
plt.legend()

plt.figure()
plt.title('Mdot')
# plt.hlines(2.06e-5,xmin=0,xmax=1,label='Luka')
plt.semilogy(x_S,Mdot_S,label='stable')
plt.semilogy(x_U,Mdot_U,label='unstable')
plt.ylabel('Mdot [Msun/year]')
plt.xlabel('x = 1-R/r')
plt.legend()



ind = 1
file = 'unstable/2test'

All_E = []
for i in range(10):
    A = np.loadtxt(file + '000'+str(i)+'.blk',skiprows=3)
    E_S = np.transpose(A)[ind]
    All_E.append(E_S)
    print('read ', i)
for i in range(10,100):
    A = np.loadtxt(file + '00'+str(i)+'.blk',skiprows=3)
    E_S = np.transpose(A)[ind]
    All_E.append(E_S)
    print('read ', i)

All_E = np.transpose(All_E)
mean_E = np.mean(All_E,axis=1)
nx,nt = np.shape(All_E)

Relldiff_E = copy.copy(All_E)
for x in range(nx):
    Relldiff_E[x] = (All_E[x]-mean_E[x])/mean_E[x]

plt.figure()
plt.pcolor(Relldiff_E,cmap='seismic',vmin=-0.5,vmax=0.5)
cbar = plt.colorbar()
cbar.set_label('Relative difference')


plt.show()
