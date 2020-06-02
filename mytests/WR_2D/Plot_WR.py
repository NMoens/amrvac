import numpy as np
import matplotlib.pyplot as plt

A = np.loadtxt('out0300.blk',skiprows=3)

r = np.transpose(A)[0]
# v = np.transpose(A)[2]*1e8
v_a = np.transpose(A)[-3]
T = np.transpose(A)[-2]
kappa = np.transpose(A)[5]

x = 1 - r[0]/r

# plt.plot(r,v_a)
# plt.plot(r,v)
plt.plot(r,T)
plt.show()
