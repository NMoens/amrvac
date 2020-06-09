import numpy as np
import matplotlib.pyplot as plt

R_sun = 6.96e10
M_sun = 1.99e33
L_sun = 3.9e33
G = 6.67e-8
year = 365*24*60*60

#READING IN LUKA 1D HD

A = np.loadtxt('out0300.blk',skiprows=3)

r_1D = np.transpose(A)[0]
# v = np.transpose(A)[2]*1e8
rho_1D = np.transpose(A)[-4]/np.transpose(A)[-3]
Mdot_1D = 4*np.pi*np.transpose(A)[0]**2*np.transpose(A)[-4]
Mdot_1D = Mdot_1D*R_sun**2*year/M_sun
v_1D = np.transpose(A)[-3]
T_1D = np.transpose(A)[-2]

x_1D = 1 - r_1D[0]/r_1D

# plt.plot(r,v_a)
# plt.plot(r,v)


#READING IN NICO 2D FLD
B = np.loadtxt('Output.okc',skiprows=23)
r_2D = np.transpose(B)[0]
rho_2D = np.transpose(B)[3]*1.e-7
v_2D = np.transpose(B)[4]*1.e8
T_2D = np.transpose(B)[5]
Mdot_2D = np.transpose(B)[6]
L_2D = np.transpose(B)[7]
kappa_2D = np.transpose(B)[8]
OPAL_2D = np.transpose(B)[9]
CAK_2D = np.transpose(B)[10]

x_2D = 1 - r_2D[0]/r_2D

L_w_2D = Mdot_2D*(v_2D**2/2 - G*M_sun/(r_2D*R_sun) + G*M_sun/(10*R_sun))/L_sun

#PLOTTING

f,(ax1,ax2) = plt.subplots(1,2,sharey=True)
f.suptitle('Density')
ax1.semilogy(r_1D,rho_1D,'b-',label='1D by LUKA')
ax1.semilogy(r_2D,rho_2D,'r.',label='2D by NICO')
ax2.semilogy(x_1D,rho_1D,'b-')
ax2.semilogy(x_2D,rho_2D,'r.')
ax1.set_ylabel('$\\rho [g/cm^3]$')
ax1.set_xlabel('$r/R_*$')
ax2.set_xlabel('$x = 1-R_*/r$')
ax1.legend()

f,(ax1,ax2) = plt.subplots(1,2,sharey=True)
f.suptitle('Velocity')
ax1.plot(r_1D,v_1D,'b-',label='1D by LUKA')
ax1.plot(r_2D,v_2D,'r.',label='2D by NICO')
ax2.plot(x_1D,v_1D,'b-')
ax2.plot(x_2D,v_2D,'r.')
ax1.set_ylabel('$v [cm/s]$')
ax1.set_xlabel('$r/R_*$')
ax2.set_xlabel('$x = 1-R_*/r$')
ax1.legend()


f,(ax1,ax2) = plt.subplots(1,2,sharey=True)
f.suptitle('Temperature')
ax1.plot(r_1D,T_1D,'b-',label='1D by LUKA')
ax1.plot(r_2D,T_2D,'r.',label='2D by NICO')
ax2.plot(x_1D,T_1D,'b-')
ax2.plot(x_2D,T_2D,'r.')
ax1.set_ylabel('$T_r,T_g [K]$')
ax1.set_xlabel('$r/R_*$')
ax2.set_xlabel('$x = 1-R_*/r$')
ax1.legend()


f,(ax1,ax2) = plt.subplots(1,2,sharey=True)
f.suptitle('Mass loss rate')
ax1.plot(r_1D,Mdot_1D,'b-',label='1D by LUKA')
ax1.plot(r_2D,Mdot_2D,'r.',label='2D by NICO')
ax2.plot(x_1D,Mdot_1D,'b-')
ax2.plot(x_2D,Mdot_2D,'r.')
ax1.set_ylabel('$\dot{M} [M_\odot/yr]$')
ax1.set_xlabel('$r/R_*$')
ax2.set_xlabel('$x = 1-R_*/r$')
ax1.legend()


f,(ax1,ax2) = plt.subplots(1,2,sharey=True)
f.suptitle('Luminosity')
ax1.hlines(10**5.28,1,70,'b',label='1D by LUKA')
ax1.plot(r_2D,L_2D,'r.',label='2D by Nico')
ax1.plot(r_2D,L_w_2D,'k.',label='Wind Luminosity')
ax2.hlines(10**5.28,0,1,'b')
ax2.plot(x_2D,L_2D,'r.')
ax2.plot(x_2D,L_w_2D,'k.')
ax1.set_ylim([0,1e6])
ax1.set_ylabel('$L/L_\odot$')
ax1.set_xlabel('$r/R_*$')
ax2.set_xlabel('$x = 1-R_*/r$')
ax1.legend()

f,(ax1,ax2) = plt.subplots(1,2,sharey=True)
f.suptitle('Opacity')
ax1.plot(r_2D,kappa_2D,'k.')
ax1.plot(r_2D,OPAL_2D,'b.')
ax1.plot(r_2D,CAK_2D,'r.')
ax2.plot(x_2D,kappa_2D,'k.',label='Total')
ax2.plot(x_2D,OPAL_2D,'b.',label='OPAL')
ax2.plot(x_2D,CAK_2D,'r.',label='CAK')
ax1.set_ylim([-0.1, 6])
ax1.set_ylabel('$\\kappa/\\kappa_e$')
ax1.set_xlabel('$r/R_*$')
ax2.set_xlabel('$x = 1-R_*/r$')
ax2.legend()

# v_1D = v_1D*(1+0.05*np.random.randn(len(v_1D)))

#GRADV RESOLUTION PROBLEM:
gradv1 = (v_1D[2:]-v_1D[:-2])/(r_1D[2:]-r_1D[:-2])
r1 = r_1D[1:-1]
gradv2 = (v_1D[4:]-v_1D[:-4])/(r_1D[4:]-r_1D[:-4])
r2 = r_1D[2:-2]
gradv3 = (v_1D[6:]-v_1D[:-6])/(r_1D[6:]-r_1D[:-6])
r3 = r_1D[3:-3]

gradv_combo = np.append(gradv1[:400],gradv2[401:800])
gradv_combo = np.append(gradv_combo,gradv3[801:])
r_combo = np.append(r1[:400],r2[401:800])
r_combo = np.append(r_combo,r3[801:])

plt.figure()
plt.plot(r1, gradv1, label='$v_{i+1}-v_{i-1}$')
plt.plot(r2, gradv2, label='$v_{i+2}-v_{i-2}$')
plt.plot(r3, gradv3, label='$v_{i+3}-v_{i-3}$')
plt.plot(r_combo,gradv_combo,'k-',label='amr situation')
plt.legend()


# https://en.wikipedia.org/wiki/Finite_difference_coefficient
gradv_o1 = gradv1[2:-2]
gradv_o2 = 2*(2/3*gradv1[2:-2] - 1/12*gradv2[1:-1])
gradv_o3 = 2*(3/4*gradv1[2:-2] - 3/20*gradv2[1:-1] + 1/60*gradv3)

plt.figure()
plt.plot(r3, gradv_o1, label='$O(\Delta x^2)$')
plt.plot(r3, gradv_o2, label='$O(\Delta x^4)$')
plt.plot(r3, gradv_o3, label='$O(\Delta x^6)$')
plt.legend()

plt.show()
