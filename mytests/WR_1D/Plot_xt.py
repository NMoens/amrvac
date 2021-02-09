import numpy as np
import copy
import matplotlib.pyplot as plt


R_sun = 6.96e10
M_sun = 1.99e33
L_sun = 3.9e33
G = 6.67e-8
year = 365*24*60*60



ind = 2
# file = 'unstable/2test'
# file = 'stable/Cak_in_diff'
# file = 'unstable/2Cak_in_diff'
# file = 'unstable/2Cak_notin_diff'
# file = 'test/w_e'
file = 'test/u_w_e'

unit_time = 695.9
All_E = []
t_axis = []
A = np.loadtxt(file + '0000.blk',skiprows=3)
r_axis = np.transpose(A)[0]
x_axis = 1. - 1./r_axis


for i in range(100,3000):
    t = str(i).zfill(4)
    try:
        A = np.loadtxt(file + t + '.blk',skiprows=3)
        E_S = np.transpose(A)[ind]
        All_E.append(E_S)

        t = open(file + t +'.blk', "r").readlines()[2]
        t_axis.append(float(t))
        print('read ', i)
    except:
        print('Failed ', i)

t_axis = np.multiply(t_axis,unit_time)

t_mesh,r_mesh = np.meshgrid(t_axis,r_axis)
x_mesh = 1.-1./r_mesh

first_E = All_E[0]
All_E = np.transpose(All_E)
mean_E = np.mean(All_E,axis=1)
nx,nt = np.shape(All_E)

try:
    Logdiff_E = copy.copy(All_E)
    Relldiff_E = copy.copy(All_E)
    for x in range(nx):
        # Logdiff_E[x] = np.log10(All_E[x]/mean_E[x])
        # Logdiff_E[x] = np.log10(All_E[x]/first_E[x])
        # Relldiff_E[x] = (All_E[x]-mean_E[x])/mean_E[x]
        Relldiff_E[x] = (All_E[x]-first_E[x])/first_E[x]
except:
    print("Logdiff or Relldiff not defined")

plt.figure()
plt.pcolor(t_mesh,x_mesh,Relldiff_E,cmap='seismic',vmin=-0.01,vmax=0.01)
# plt.pcolor(t_mesh,x_mesh,np.log10(All_E),cmap='viridis',vmin=0,vmax=0.6)

cbar = plt.colorbar()
cbar.set_label('$Log(\\rho/<\\rho>)$')
# plt.contour(t_mesh,x_mesh,Logdiff_E, [0.0], colors = 'k')

plt.ylabel('x = 1-R/r')
plt.xlabel('t [s]')

plt.figure()
plt.plot(r_axis,mean_E)


plt.show()
