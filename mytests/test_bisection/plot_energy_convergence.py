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
plt.loglog(t,e_g,'ro')
plt.loglog(t,e_r,'bo')

t, e_eq, e_r, e_g = ReadE('/energy_2')
plt.loglog(t,e_g,'ro')
plt.loglog(t,e_r,'bo')

plt.loglog(t,e_eq,'r--')
plt.show()
