import numpy as np
import matplotlib.pyplot as plt
import os

def ReadE(filename):
    f = os.getcwd() + filename

    A = np.loadtxt(f)

    t = A.transpose()[0]
    e_eq = A.transpose()[1]
    e_r = A.transpose()[2]
    e_g = A.transpose()[3]

    return t, e_eq, e_r, e_g

plt.title('Energy interaction schemes')

t, e_eq, e_r, e_g = ReadE('/newton_11')
plt.loglog(t,e_g,'r1',label='dt = 1e-10')

t, e_eq, e_r, e_g = ReadE('/newton_12')
plt.loglog(t,e_g,'b2',label='dt = 1e-8')

t, e_eq, e_r, e_g = ReadE('/newton_13')
plt.loglog(t,e_g,'y3',label='dt = 1e-6')

t, e_eq, e_r, e_g = ReadE('/newton_21')
plt.loglog(t,e_g,'r1')

t, e_eq, e_r, e_g = ReadE('/newton_22')
plt.loglog(t,e_g,'b2')

t, e_eq, e_r, e_g = ReadE('/newton_23')
plt.loglog(t,e_g,'y3')

plt.loglog(t,e_eq,'k--', label='equilibrium level')

plt.xlabel('time')
plt.ylabel('gas energy')

plt.legend()
plt.show()
