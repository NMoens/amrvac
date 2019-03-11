import numpy as np
import os
import matplotlib.pyplot as plt


f_tot = os.getcwd() + '/init_struc_total'
f_vac = os.getcwd() + '/init_struc_amrvac'

A_tot = np.loadtxt(f_tot,dtype = 'string')
A_vac = np.loadtxt(f_vac,dtype = 'string')

def stringtodouble(M):
    I,J = np.shape(M)
    for i in range(I):
        for j in range(J):
            m = M[i][j]
            m = float(m.replace('D','e'))
            M[i][j] = m
    return M

A_tot = stringtodouble(A_tot)
A_vac = stringtodouble(A_vac)

i_vac, y_vac, rho_vac, tg_vac, pg_vac, er_vac, tau_vac = np.transpose(A_vac)
rho_tot, tr_tot, pg_tot, tau_tot, mc_tot, gammar_tot, er_tot, y_tot = np.transpose(A_tot)

print '#cells total structure: ', len(y_tot)
print '#cells amrvac structure: ', len(y_vac)

r_bottom_cm = float(y_vac[0])
r_top_cm = float(y_vac[-1])

print 'bottom bound = ' + str(r_bottom_cm) + ' cm'
print 'upper bound = ' + str(r_top_cm) + ' cm'

hp_top_sr = 2.29e-3
hp_bottom_sr = 6.644e-3

hp_bottom_cm = hp_bottom_sr*r_bottom_cm

print 'one scaleheight is ' + str(hp_bottom_sr) + ' stellar radii'
print 'one scaleheight is ' + str(hp_bottom_cm) + ' cm'

r_bottom_hp = r_bottom_cm/ hp_bottom_cm
r_top_hp = r_top_cm/ hp_bottom_cm

print 'In code units (hp), the star bottom starts at ' + str(r_bottom_hp) + ' hp'
print 'In code units (hp), the star top ends at ' + str(r_top_hp) + ' hp'

print 'The numerical domain entails [0, ' + str(r_top_hp - r_bottom_hp) + ' ]'
print 'If there are ' + str(len(y_vac)) + ' cells, dx should be ' + str((r_top_hp - r_bottom_hp)/len(y_vac))

plt.plot(y_tot, rho_tot)
plt.plot(y_vac, rho_vac)
plt.show()
