import numpy as np
import matplotlib.pyplot as plt

import os

f = os.getcwd() + '/f1_1'
t,momentum,theoretical,flux = np.loadtxt(f,unpack = True)

plt.figure()
plt.plot(t,theoretical,'b-',label = 'theoretical')
plt.plot(t,momentum,'r.',label = 'simulation')

plt.figure()
plt.semilogy(flux)

plt.show()
