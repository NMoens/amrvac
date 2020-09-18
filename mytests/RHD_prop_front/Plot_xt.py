import numpy as np
import copy
import matplotlib.pyplot as plt
import matplotlib.colors as colors


ind = 4
# file = 'unstable/2test'
# file = 'stable/Cak_in_diff'
file = 'test/test'
# file = 'unstable/2Cak_notin_diff'


unit_time = 1.3370309040864490E-006
unit_pressure = 13984827615.086012
unit_length = 1.e0
c_light = 2.99e10
All_E = []
t_axis = []
A = np.loadtxt(file + '0000.blk',skiprows=3)
r_axis = np.transpose(A)[0]
x_axis = 1. - 1./r_axis

for i in range(0,10):
    try:
        A = np.loadtxt(file + '000'+str(i)+'.blk',skiprows=3)
        E_S = np.transpose(A)[ind]
        All_E.append(E_S)

        t = open(file + '000'+str(i)+'.blk', "r").readlines()[2]
        t_axis.append(float(t))
        print('read ', i)
    except:
        print('Failed ', i)

for i in range(10,100):
    try:
        A = np.loadtxt(file + '00'+str(i)+'.blk',skiprows=3)
        E_S = np.transpose(A)[ind]
        All_E.append(E_S)

        t = open(file + '00'+str(i)+'.blk', "r").readlines()[2]
        t_axis.append(float(t))
        print('read ', i)
    except:
        print('Failed ', i)

for i in range(100,300):
    try:
        A = np.loadtxt(file + '0'+str(i)+'.blk',skiprows=3)
        E_S = np.transpose(A)[ind]
        All_E.append(E_S)

        t = open(file + '0'+str(i)+'.blk', "r").readlines()[2]
        t_axis.append(float(t))
        print('read ', i)
    except:
        print('Failed ', i)

t_axis = np.multiply(t_axis,unit_time)
r_axis = np.multiply(r_axis,unit_length)
t_mesh,r_mesh = np.meshgrid(t_axis,r_axis)

All_E = np.transpose(All_E)*unit_pressure
mean_E = np.mean(All_E,axis=1)

plt.figure()
# PC = plt.pcolor(r_mesh,t_mesh,All_E,cmap='viridis')
PC = plt.pcolor(r_mesh,t_mesh,All_E,norm=colors.LogNorm(vmin=1.e-3, vmax=1.e0),cmap='viridis')


x_an = np.linspace(0,2.5,1000)
t_an = x_an/c_light
t_an_1 = x_an/(c_light*0.95)
t_an_2 = x_an/(c_light*0.90)
t_an_3 = x_an/(c_light*0.85)

CF = plt.contour(r_mesh,t_mesh,All_E, [1e-1], colors = 'r')
plt.annotate('E=0.5 contour',(0.3,0.6))

cbar = plt.colorbar(PC)
cbar.add_lines(CF)
cbar.set_label('$E_{rad}$ [erg/cm^3]')

plt.plot(x_an,t_an,c='1.0',label='speed of light')
plt.plot(x_an,t_an_1,c='0.5',label='95% c')
plt.plot(x_an,t_an_2,c='0.5',label='90% c')
plt.plot(x_an,t_an_3,c='0.25',label='85% c')

# plt.xlim([-0.1,0.5])
# plt.ylim([0.0,2.e-11])

plt.xlim([-0.1,1.4])
plt.ylim([0.0,4.e-11])

plt.xlabel('x [cm]')
plt.ylabel('t [s]')

plt.legend(loc='best')
plt.show()
