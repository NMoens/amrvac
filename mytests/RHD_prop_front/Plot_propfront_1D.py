from amrvac_tools.datfiles.reading import amrvac_reader
from scipy import interpolate as itp
import numpy as np
import os
import matplotlib.pyplot as plt

unit_length = 1.0000000000000000
unit_numberdensity = 1.4936755387807165e22
unit_temperature = 3390.6871563953396

c = 2.99792458e10
year = 3.1536e7
arad = 7.5657e-15
kb = 1.380648520e-16
mp = 1.6737236e-24
mu = 0.60869565217391308
# hd_gamma = 1.6667

def AMRVAC_single_profile(file,variable):
    A = np.loadtxt(file,skiprows=3)
    x = np.transpose(A)[0]

    re = np.transpose(A)[4]
    return x, re

# folder = 'output'



x0,Er_Euler1 = AMRVAC_single_profile('test/Euler0010.blk','re')
x0,Er_SP1 = AMRVAC_single_profile('test/SP0010.blk','re')
x0,Er_Midpoint1 = AMRVAC_single_profile('test/SP0010.blk','re')
x0,Er_ARS31 = AMRVAC_single_profile('test/SP0010.blk','re')


t0 = 0.e0
t1 = 10e-6
t2 = 20e-6

f0 = 0.1 + t0*c
f1 = 0.1 + t1*c
f2 = 0.1 + t2*c

f,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)

# ax1.plot(r0,Er0,'bx',label='Pomraning \n t = ' + str(t0))
# ax1.plot(r1,Er1,'rx',label='Minerbo \n t = ' + str(t0))
ax1.plot(x0,Er_Euler1,'bx',label='Pomraning \n t = ' + str(t0))
ax1.vlines(f0,0,1)
ax1.set_ylabel('$E_r [erg cm^{-3}]$')
ax1.legend()

# ax2.plot(r2,Er2,'bx',label='t = ' + str(t2))
# ax2.plot(r3,Er3,'rx',label='t = ' + str(t3))
ax2.vlines(f1,0,1)
ax2.set_ylabel('$E_r [erg cm^{-3}]$')
ax2.legend()

# ax3.plot(r4,Er5,'bx',label='t = ' + str(t4))
# ax3.plot(r5,Er5,'rx',label='t = ' + str(t5))
ax3.vlines(f2,0,1)
ax3.set_ylabel('$E_r [erg cm^{-3}]$')
ax3.legend()

ax3.set_xlabel('$x/cm$')

plt.show()
