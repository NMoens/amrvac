import numpy as np
import matplotlib.pyplot as plt

f = 'Halley1_1.d2'

A = np.loadtxt(f)

t = A.transpose()[0]
e_g = A.transpose()[1]
E_r = A.transpose()[2]
e_eq = A.transpose()[3]
T_g = A.transpose()[4]
T_r = A.transpose()[5]


fig, ax1 = plt.subplots()

ax1.set_xlabel('time (s)')
ax1.set_ylabel('Energy density [equilibrium]', color='r')
ax1.loglog(t, e_g, 'rx', label='gas energy ')
ax1.loglog(t, e_eq, 'r--', label='radiation energy ')
ax1.loglog(t, E_r, 'r+',label='Radiation energy ')
ax1.tick_params(axis='y', labelcolor='r')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('$Log(T)$', color='b')  # we already handled the x-label with ax1
ax2.loglog(t, T_g, 'bx', label = 'gas temperature')
ax2.loglog(t,T_r, 'b+', label ='radiation temperature')
ax2.tick_params(axis='y', labelcolor='b')

ax1.legend(loc='left')
ax2.legend(loc='right')

fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
