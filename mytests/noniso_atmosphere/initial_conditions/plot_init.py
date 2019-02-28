import numpy as np
import os
import matplotlib.pyplot as plt


f = os.getcwd() + '/init_struc_total'
g = os.getcwd() + '/init_struc_amrvac'

rho_tot,tr_tot,pg_tot,tau_tot,mc_tot,gammar_tot,er_tot,yarr_tot = np.loadtxt(f,unpack=True)
rho_amr,tr_amr,pg_amr,tau_amr,mc_amr,gammar_amr,er_amr,yarr_amr = np.loadtxt(g,unpack=True)

plt.plot(y_arr_tot,er_tot,'b-')
plt.plot(y_arr_amr,er_amr,'ro')

plt.show()
