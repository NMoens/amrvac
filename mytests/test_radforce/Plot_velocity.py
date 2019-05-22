import numpy as np
import matplotlib.pyplot as plt

import os

f = os.getcwd() + '/f1_1'
t,momentum,theoretical,flux = np.loadtxt(f,unpack = True)

plt.figure()
plt.plot(t,theoretical,'b-',label = 'theoretical')
plt.plot(t,momentum,'r.',label = 'simulation')

plt.title('Testcase radiation force',fontsize=20)
plt.xlabel('$t$',fontsize=15)
plt.ylabel('$\\vec{v}\\rho$',fontsize=15)
plt.legend()
plt.show()
