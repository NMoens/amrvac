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


t, e_eq, e_r, e_g = ReadE('/energy_1')
plt.loglog(t,e_g,'r--')
plt.loglog(t,e_r,'r-')

t, e_eq, e_r, e_g = ReadE('/energy_2')
plt.loglog(t,e_g,'b--')
plt.loglog(t,e_r,'b-')

t, e_eq, e_r, e_g = ReadE('/energy_3')
plt.loglog(t,e_g,'y--')
plt.loglog(t,e_r,'y-')

plt.loglog(t,e_eq,'k-')
plt.show()
