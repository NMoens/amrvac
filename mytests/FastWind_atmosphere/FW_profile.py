import numpy as np
import os
import matplotlib.pyplot as plt


f = os.getcwd() + '/FastWind_profile.txt'

r_fw,rho_fw,v_fw,pg_fw,T_fw,mu_fw,a2_fw = np.loadtxt(f,unpack=True,skiprows = 3)


f_relax1 = os.getcwd() + '/Amrvac_profile.okc'
A_relax1 = np.loadtxt(f_relax1,skiprows = 17)

x_rel, y_rel, z_rel, int_rho1, int_p1, int_r_e1, int_dt1, F_rel1 = np.transpose(A_relax1[0::128])

unit_length = 5075047950.3534832
unit_density = 8.8837999999999995e-9
unit_pressure = 128253.75816140750
unit_temperature = 106459.89999999999
unit_radflux = 487310493625.60626
unit_opacity = 2.217e-2
unit_time = 1335.68

rho_rel1 = int_rho1/int_dt1*unit_density
p_rel1 = int_p1/int_dt1*unit_pressure
er_rel1 = int_r_e1/int_dt1*unit_pressure
F_rel1 = F_rel1*unit_radflux

R_star = 1252700000000.0
y_rel = y_rel*unit_length + R_star

plt.semilogy(r_fw,rho_fw,'o')
plt.semilogy(y_rel,rho_rel1,'o')
plt.show()
