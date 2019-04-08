import numpy as np
import os
import matplotlib.pyplot as plt

c = 2.99e10
lam = 1./3.
kap = 0.34
G = 6.67e-8
R_sun = 6.96e10
M_sun = 1.99e33
R_star = 20.0*R_sun
M_star = 50.0*M_sun

f_tot = os.getcwd() + '/initial_conditions/init_struc_total'
f_vac = os.getcwd() + '/initial_conditions/init_struc_amrvac'
f_it = os.getcwd() + '/iterated_profile'

A_tot = np.loadtxt(f_tot)
A_vac = np.loadtxt(f_vac)
A_it = np.loadtxt(f_it)


i_vac, y_vac, rho_vac, tg_vac, pg_vac, er_vac, tau_vac = np.transpose(A_vac)
rho_tot, tr_tot, pg_tot, tau_tot, mc_tot, gammar_tot, er_tot, y_tot = np.transpose(A_tot)
i_it, y_it, rho_it, pg_it, tr_it, er_it = np.transpose(A_it)
tg_it = tr_it

print (len(tg_it), len(y_it))


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


f_relax1 = os.getcwd() + '/relaxed_profile_1.okc'
A_relax1 = np.loadtxt(f_relax1,skiprows = 17)

x_rel, y_rel, z_rel, int_rho1, int_p1, int_r_e1, int_dt1, F_rel1 = np.transpose(A_relax1[0::128])

# f_relax2 = os.getcwd() + '/relaxed_profile_2.okc'
# A_relax2 = np.loadtxt(f_relax2,skiprows = 17)
#
# x_rel, y_rel, z_rel, int_rho2, int_p2, int_r_e2, int_dt2, F_rel2 = np.transpose(A_relax2[0::128])

f_final1 = os.getcwd() + '/final_profile_1.okc'
A_final1 = np.loadtxt(f_final1,skiprows = 15)

x_fin, y_fin, z_fin, rho_fin, p_fin, er_fin, F_fin = np.transpose(A_final1[0::128])

unit_length = 8963898879.7622433
unit_density = 3.2439805779365456E-008
unit_pressure = 530577.09957521548
unit_temperature = 120610.46245923516
unit_radflux = 2024138883985.9893

rho_rel1 = int_rho1/int_dt1*unit_density
p_rel1 = int_p1/int_dt1*unit_pressure
er_rel1 = int_r_e1/int_dt1*unit_pressure
F_rel1 = F_rel1*unit_radflux

# rho_rel2 = int_rho2/int_dt2*unit_density
# p_rel2 = int_p2/int_dt2*unit_pressure
# er_rel2 = int_r_e2/int_dt2*unit_pressure
# F_rel2 = F_rel2*unit_radflux

rho_fin = rho_fin*unit_density
p_fin = p_fin*unit_pressure
er_fin = er_fin*unit_pressure
F_fin = F_fin*unit_radflux

y_rel = y_rel*unit_length + r_bottom_cm
y_fin = y_fin*unit_length + r_bottom_cm

f,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)


def rho_p_2_T(rho,p):
    mu = 0.6
    mp = 1.672e-24
    kb = 1.3806e-16
    return mp*mu/kb * p/rho

def E_r_2_T(E_r):
    a = 7.5646e-15
    return (E_r/a)**(1./4)

def F_2_Gamma(F,rho):
    g_rad = kap*rho[1:]/c*F[:]
    g_grav = G*M_star/R_star**2
    return g_rad/g_grav

def tem_ratio(Trad,Tgas):
    return Tgas/Trad


# ax1.semilogy(y_tot, rho_tot, label = 'Initial conditions total')
ax1.semilogy(y_vac, rho_vac, label = 'Initial conditions amrvac')
ax1.semilogy(y_rel, rho_rel1,'o', label = 'Average')
ax1.semilogy(y_fin, rho_fin,'o', label = 'Final')
ax1.semilogy(y_it, rho_it,'o', label = 'iterated')

ax1.legend()
ax1.set_ylabel('$\\rho$ [g/cm3]')

# ax2.plot(y_tot, rho_p_2_T(rho_tot,pg_tot))
# ax2.plot(y_vac, rho_p_2_T(rho_vac,pg_vac))
# ax2.plot(y_rel, rho_p_2_T(rho_rel1,p_rel1),'o')
# ax2.plot(y_rel, rho_p_2_T(rho_fin,p_fin),'o')
# ax2.plot(y_it, tg_it,'o')
# ax2.plot(y_tot, pg_tot)
ax2.plot(y_vac, pg_vac)
ax2.plot(y_rel, p_rel1,'o')
ax2.plot(y_rel, p_fin,'o')
ax2.plot(y_it, pg_it,'o')


ax2.set_ylabel('$p_{gas}$ [g/cm3]')

# ax3.plot(y_tot, E_r_2_T(er_tot))
# ax3.plot(y_vac, E_r_2_T(er_vac))
# ax3.plot(y_rel, E_r_2_T(er_rel1),'o')
# ax3.plot(y_rel, E_r_2_T(er_fin),'o')
# ax3.plot(y_it, tg_it,'o')
# ax3.plot(y_tot, er_tot)
ax3.plot(y_vac, er_vac)
ax3.plot(y_rel, er_rel1,'o')
ax3.plot(y_rel, er_fin,'o')
ax3.plot(y_it, er_it,'o')

ax3.set_xlabel('r [cm]')
ax3.set_ylabel('$E_{rad}$ [g/cm3]')


flux_tot = -lam*c/(kap*rho_tot[1:])*(er_tot[1:] - er_tot[:-1])/(y_tot[1:] - y_tot[:-1])
flux_vac = -lam*c/(kap*rho_vac[1:])*(er_vac[1:] - er_vac[:-1])/(y_vac[1:] - y_vac[:-1])
flux_it = -lam*c/(kap*rho_it[1:])*(er_it[1:] - er_it[:-1])/(y_it[1:] - y_it[:-1])


plt.figure()
# plt.plot(y_tot[1:],flux_tot,label='total profile')
plt.plot(y_vac[1:],flux_vac,label='initial amrvac profile')
plt.plot(y_it[1:],flux_it,label='iterated profile')
plt.plot(y_rel,F_rel1,label='Averaged')
plt.plot(y_fin,F_fin,label='Final')

# plt.plot(y_tot,rho_p_2_T(rho_tot,pg_tot)/E_r_2_T(er_tot),label='total profile')
# plt.plot(y_vac,rho_p_2_T(rho_vac,pg_vac)/E_r_2_T(er_vac),label='initial amrvac profile')
# plt.plot(y_it,rho_p_2_T(rho_it,pg_it)/E_r_2_T(er_it),label='iterated profile')
# plt.plot(y_rel,rho_p_2_T(rho_rel1,p_rel1)/E_r_2_T(er_rel1),label='Average')
# plt.plot(y_rel,rho_p_2_T(rho_fin,p_fin)/E_r_2_T(er_fin),label='Final')

# plt.plot(y_tot[1:],F_2_Gamma(flux_tot,rho_tot),label='total profile')
# plt.plot(y_vac[1:],F_2_Gamma(flux_vac,rho_vac),label='initial amrvac profile')
# plt.plot(y_it[1:],F_2_Gamma(flux_it,rho_it),label='iterated profile')
# plt.plot(y_rel[1:],F_2_Gamma(F_rel1[1:],rho_rel1),label='Pomraning')
# plt.plot(y_rel[1:],F_2_Gamma(F_rel2[1:],rho_rel2),label='Diffusion')
plt.legend()

plt.show()
