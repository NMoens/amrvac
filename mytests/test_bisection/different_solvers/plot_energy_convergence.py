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

t, e_eq, e_r, e_g = ReadE('/bisect_1')
plt.loglog(t,e_g,'r1',label='bisect 1')

t, e_eq, e_r, e_g = ReadE('/bisect_2')
plt.loglog(t,e_g,'b1',label='bisect 2')

t, e_eq, e_r, e_g = ReadE('/newton_1')
plt.loglog(t,e_g,'r2',label='Newton 1')

t, e_eq, e_r, e_g = ReadE('/newton_2')
plt.loglog(t,e_g,'b2',label='Newton 2')

t, e_eq, e_r, e_g = ReadE('/halley_1')
plt.loglog(t,e_g,'r3',label='Halley 1')

t, e_eq, e_r, e_g = ReadE('/Halley_updated')
plt.loglog(t,e_g,'b3',label='Halley New')

plt.loglog(t,e_eq,'k--', label='equilibrium level')

plt.xlabel('time')
plt.ylabel('gas energy')

plt.legend()
plt.show()
