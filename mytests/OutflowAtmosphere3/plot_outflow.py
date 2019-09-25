import numpy as np
import os
import matplotlib.pyplot as plt

f_relax1 = os.getcwd() + '/FirstSnap.okc'
A_relax1 = np.loadtxt(f_relax1,skiprows = 17)

f_relax2 = os.getcwd() + '/LastSnap.okc'
A_relax2 = np.loadtxt(f_relax2,skiprows = 17)


def convert_data(A):
    nlines,ncols = np.shape(A)

    x, y, z, rho, v, p, er, gamma = np.sort(np.transpose(A),axis=1)


    print(rho)
    print(v)
    print(gamma)
    return x, y, z, rho, v, p, er, gamma

x_rel, y_rel, z_rel, rho1, v1, p1, er1, Gamma1 = convert_data(A_relax1)
x_rel, y_rel, z_rel, rho2, v2, p2, er2, Gamma2 = convert_data(A_relax2)

unit_length = 69599000000.000000
unit_density = 1.7098979425032789e-009
unit_pressure = 43950.111561334961
unit_temperature = 189666.48683552662
unit_radflux = 222820363223.95514
unit_opacity = 8.4028538615729537e-003
unit_time = 13728.026336097881
unit_velocity = 5069847.5000000000

ny = len(y_rel)


rho1 = rho1*unit_density
v1 = v1*unit_velocity
p1 = p1*unit_pressure
er1 = er1*unit_pressure


rho2 = rho2*unit_density
v2 = v2*unit_velocity
p2 = p2*unit_pressure
er2 = er2*unit_pressure

f,(ax1,ax2) = plt.subplots(2,1,sharex=True,figsize=(10,5))
ax1.set_title('Density and Velocity', fontsize = 15)

ax1.plot(y_rel,rho1, label = 'initial')
for i in range(128):
    ax1.plot(y_rel,rho2, 'r',label = 'final')
    # ax1.semilogy(y_rel[i*256::(i+1)*256],rho2[i*256:(i+1)*256], 'r',label = 'final')
ax1.set_ylabel('$\\rho [g cm^{-3}]$', fontsize = 15)
ax1.legend()

ax2.plot(y_rel,v1)
ax2.plot(y_rel,v2)
ax2.set_ylabel('$radial Velocity [erg cm^{-3}]$', fontsize = 15)
ax2.set_xlabel('r [$R_\odot$]', fontsize = 15)

plt.savefig('1Dprof')
plt.show()
