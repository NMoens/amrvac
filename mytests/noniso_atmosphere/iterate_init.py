import numpy as np
import os
import matplotlib.pyplot as plt


f_vac = os.getcwd() + '/initial_conditions/init_struc_amrvac'

A_vac = np.loadtxt(f_vac)

i_vac, y_vac, rho_vac, tg_vac, pg_vac, er_vac, tau_vac = np.transpose(A_vac)

mu = 0.6
mp = 1.672e-24
arad = 7.564e-15
kb = 1.38e-16
G_newt = 6.672e-8

M_sun = 1.9e33
R_sun = 6.96e10

M_star = 50.0*M_sun
R_star = 20.0*R_sun

g_star = G_newt*M_star/R_star**2
gamma = 0.5

rho_0 = rho_vac[0]
p_0 = pg_vac[0]

print(rho_vac[0], pg_vac[0], er_vac[0],g_star,gamma)

y_new = (y_vac-y_vac[0])

print y_new

def IterateConditions(y_array):
    #Initial iteration:
    g0 = g_star*(gamma-1)

    a2 = (p_0/rho_0)**0.5*np.ones((len(y_array)))
    scaleheight = a2/g0*np.ones((len(y_array)))
    p = p_0*np.exp(-y_array/scaleheight)
    rho = rho_0*np.exp(-y_array/scaleheight)


    #plt.plot(y_array,p,label = str(0))

    for i in range(100):

        rho = rho_0*np.exp(-y_array/scaleheight)
        temp = (3.*gamma/(1.-gamma)*rho*kb/mu/mp/arad)**(1./3)
        p = kb*rho*temp/mu/mp

        a2 = p/rho
        scaleheight = a2/g0

        if (i > 10):
            plt.plot(y_array,p,label = str(i+1))


    return rho, p

plt.figure()

rho, p = IterateConditions(rho_vac[0],pg_vac[0],y_new)


plt.legend()
plt.show()
