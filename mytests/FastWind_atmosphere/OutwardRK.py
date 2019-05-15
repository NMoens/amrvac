import numpy as np
import matplotlib.pyplot as plt
import os

f_relax1 = os.getcwd() + '/Initial_Amrvac.okc'
A_relax1 = np.loadtxt(f_relax1,skiprows = 21)

x_rel, y_rel, z_rel, rho_vac, mom_vac, e_vac, er_vac, F_vac, kappa_vac, Gamma_vac = np.transpose(A_relax1[0::256])

unit_length = 5075047950.3534832
unit_density = 8.8837999999999995e-9
unit_pressure = 128253.75816140750
unit_temperature = 106459.89999999999
unit_radflux = 487310493625.60626
unit_opacity = 2.217e-2
unit_time = 1335.68
unit_velocity = 3799580.6174531388

ny = len(y_rel)

rho_vac = rho_vac*unit_density
mom_vac = mom_vac*unit_density*unit_velocity
e_vac = e_vac*unit_pressure
er_vac = er_vac*unit_pressure
F_vac = F_vac*unit_radflux
kappa_vac = kappa_vac*unit_opacity

v_vac = mom_vac/rho_vac
p_vac = (5./3.-1.)*(e_vac-0.5*v_vac*v_vac*rho_vac)

mu = 0.6
mp = 1.6e-24
kb = 1.38e-16
arad = 7.5646e-15
Tg_vac = p_vac/rho_vac * (mu*mp)/kb

R_star = 1252700000000.0
y_rel = y_rel*unit_length + R_star

dy = y_rel[1] - y_rel[0]

M_sun = 1.99e33
M_star = 66.8*M_sun
c = 2.99e10
G = 6.67e-8
kappa_0 = 0.34
akram = 13.1351597305
bkram = -4.5182188206
grav = G*M_star/R_star**2


y_rk = y_rel

rho_rk = [rho_vac[0]]
p_rk = [p_vac[0]]
E_rk = [er_vac[0]]
T_rk = [Tg_vac[0]]
F_rk = np.mean(F_vac)

a2 = p_vac[0]/rho_vac[0]
kappa_rk = [kappa_0*(1.+10**akram*rho_rk[0]*(a2/1.e12)**bkram)]
Gamma_rk = [kappa_rk[0]*F_rk/(grav*c)]


def get_rho(p,E):
    temp = (E/arad)**0.25
    rho = p/temp*mp*mu/kb
    return rho

def get_kappa(p,E):
    rho = get_rho(p,E)

    speedofsound = p/rho
    kappa = kappa_0*(1.+10**akram*rho*(speedofsound/1.e12)**bkram)
    return kappa

def get_Gamma(p,E):
    kappa = get_kappa(p,E)
    Gamma = kappa*F_rk/(grav*c)
    return Gamma


def f_pressure(p,E):
    return -get_rho(p,E)*grav*(1-get_Gamma(p,E))

def f_radiation(p,E):
    return -3*F_rk/c*get_kappa(p,E)*get_rho(p,E)


def Runge_Kutta():
    for i in range(1,len(y_rk)):
        k1 = dy*f_pressure(p_rk[-1], E_rk[-1])
        k2 = dy*f_pressure(p_rk[-1]+k1/2, E_rk[-1])
        k3 = dy*f_pressure(p_rk[-1]+k2/2, E_rk[-1])
        k4 = dy*f_pressure(p_rk[-1]+k3, E_rk[-1])

        p_rk.append(p_rk[-1] + 1/6*(k1+2*k2+2*k3+k4))

        l1 = dy*f_radiation(p_rk[-1],E_rk[-1])
        l2 = dy*f_radiation(p_rk[-1],E_rk[-1]+l1/2)
        l3 = dy*f_radiation(p_rk[-1],E_rk[-1]+l2/2)
        l4 = dy*f_radiation(p_rk[-1],E_rk[-1]+l3)

        E_rk.append(E_rk[-1] + 1/6*(l1+2*l2+2*l3+l4))


        rho_rk.append(get_rho(p_rk[-1],E_rk[-1]))
        kappa_rk.append(get_kappa(p_rk[-1],E_rk[-1]))
        Gamma_rk.append(get_Gamma(p_rk[-1],E_rk[-1]))
        T_rk.append((E_rk[-1]/arad)**0.25)

        print(y_rk[i], rho_rk[i])

    return

# for i in range(1,len(y_rk)):
#
#     dp = -rho_rk[i-1]*grav*(1-Gamma_rk[i-1])
#     dE = -3*F_rk*kappa_rk[i-1]*rho_rk[i-1]/c
#
#     p_rk.append(p_rk[i-1]+dy*dp)
#     E_rk.append(E_rk[i-1]+dy*dE)
#
#     T_rk.append((E_rk[i]/arad)**0.25)
#     rho_rk.append(p_rk[i]/T_rk[i]*mp*mu/kb)
#
#     a2 = p[i]/rho[i]
#     kappa_rk.append(kappa_0*(1.+10.e0**akram*rho_rk[i]*(a2/1.e12)**bkram))
#     Gamma_rk.append(kappa_rk[i]*F_rk/(grav*c))
#
#     print(y_rk[i], rho_rk[i], T_rk[i])


Runge_Kutta()

plt.figure()
plt.plot(y_rk,E_rk,'-',label='outward RK')
plt.plot(y_rel,er_vac,'-',label='initial conditions')
plt.legend()
plt.show()
