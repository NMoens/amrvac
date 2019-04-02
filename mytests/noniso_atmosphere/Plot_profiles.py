import numpy as np
import os
import matplotlib.pyplot as plt


f_tot = os.getcwd() + '/initial_conditions/init_struc_total'
f_vac = os.getcwd() + '/initial_conditions/init_struc_amrvac'

A_tot = np.loadtxt(f_tot)
A_vac = np.loadtxt(f_vac)

i_vac, y_vac, rho_vac, tg_vac, pg_vac, er_vac, tau_vac = np.transpose(A_vac)
rho_tot, tr_tot, pg_tot, tau_tot, mc_tot, gammar_tot, er_tot, y_tot = np.transpose(A_tot)

print('#cells total structure: ', len(y_tot))
print('#cells amrvac structure: ', len(y_vac))

r_bottom_cm = float(y_vac[0])
r_top_cm = float(y_vac[-1])

print('bottom bound = ' + str(r_bottom_cm) + ' cm')
print('upper bound = ' + str(r_top_cm) + ' cm')

hp_top_sr = 2.29e-3
hp_bottom_sr = 6.644e-3

hp_bottom_cm = hp_bottom_sr*r_bottom_cm

print('one scaleheight is ' + str(hp_bottom_sr) + ' stellar radii')
print('one scaleheight is ' + str(hp_bottom_cm) + ' cm')

r_bottom_hp = r_bottom_cm/ hp_bottom_cm
r_top_hp = r_top_cm/ hp_bottom_cm

print('In code units (hp), the star bottom starts at ' + str(r_bottom_hp) + ' hp')
print('In code units (hp), the star top ends at ' + str(r_top_hp) + ' hp')

print('The numerical domain entails [0, ' + str(r_top_hp - r_bottom_hp) + ' ]')
print('If there are ' + str(len(y_vac)) + ' cells, dx should be ' + str((r_top_hp - r_bottom_hp)/len(y_vac)))


################################################################################


f_relax1 = os.getcwd() + '/profile_init.okc'
A_relax1 = np.loadtxt(f_relax1,skiprows = 15)

x_rel, y_rel, z_rel, int_rho1, int_p1, int_r_e1, int_dt1 = np.transpose(A_relax1[0::128])

f_relax2 = os.getcwd() + '/relaxed_profile_2.okc'
A_relax2 = np.loadtxt(f_relax2,skiprows = 15)

x_rel, y_rel, z_rel, int_rho2, int_p2, int_r_e2, int_dt2 = np.transpose(A_relax2[0::128])

unit_length = 8963898879.7622433
unit_density = 3.2439805779365456E-008
unit_pressure = 530577.09957521548
unit_temperature = 120610.46245923516

rho_rel1 = int_rho1/int_dt1*unit_density
p_rel1 = int_p1/int_dt1*unit_pressure
er_rel1 = int_r_e1/int_dt1*unit_pressure

rho_rel2 = int_rho2/int_dt2*unit_density
p_rel2 = int_p2/int_dt2*unit_pressure
er_rel2 = int_r_e2/int_dt2*unit_pressure

y_rel = y_rel*unit_length + r_bottom_cm

f,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)


def rho_p_2_T(rho,p):
    mu = 0.6
    mp = 1.672e-24
    kb = 1.3806e-16
    return mp*mu/kb * p/rho

def E_r_2_T(E_r):
    a = 7.5646e-15
    return (E_r/a)**(1./4)

ax1.plot(y_tot, rho_tot, label = 'Initial conditions total')
ax1.plot(y_vac, rho_vac, label = 'Initial conditions amrvac')
ax1.plot(y_rel, rho_rel1,'o', label = 'time-averaged conditions')
ax1.plot(y_rel, rho_rel2,'o', label = 'time-averaged conditions')

ax1.legend()
ax1.set_ylabel('$\\rho$ [g/cm3]')

ax2.plot(y_tot, rho_p_2_T(rho_tot,pg_tot), label = 'Initial conditions total')
ax2.plot(y_vac, rho_p_2_T(rho_vac,pg_vac), label = 'Initial conditions amrvac')
ax2.plot(y_rel, rho_p_2_T(rho_rel1,p_rel1),'o', label = 'Initial conditions')
ax2.plot(y_rel, rho_p_2_T(rho_rel2,p_rel2),'o', label = 'time-averaged conditions')

ax2.legend()
ax2.set_ylabel('$T_{gas}$ [g/cm3]')

ax3.plot(y_tot, E_r_2_T(er_tot), label = 'Initial conditions total')
ax3.plot(y_vac, E_r_2_T(er_vac), label = 'Initial conditions amrvac')
ax3.plot(y_rel, E_r_2_T(er_rel1),'o', label = 'Initial conditions')
ax3.plot(y_rel, E_r_2_T(er_rel2),'o', label = 'time-averaged conditions')

ax3.legend()
ax3.set_xlabel('r [cm]')
ax3.set_ylabel('$T_{rad}$ [g/cm3]')



plt.figure()
plt.plot(y_tot,tau_tot)
plt.xlabel('r [cm]')
plt.ylabel('$\\tau$')

plt.show()
