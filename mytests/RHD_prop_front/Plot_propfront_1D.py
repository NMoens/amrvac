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


plt.figure()
i_min = 10
i_max = 300
n_i = 20
for ii in np.linspace(i_min,i_max,n_i):
    t = int(ii)
    t = str(t).zfill(4)
    color = (ii - i_min)/(i_max - i_min)

    x0,Er = AMRVAC_single_profile('test/test' + t +'.blk','re')
    plt.plot(x0, Er, c= str(color))

plt.title('Optically thin propagation')
plt.xlabel('x [cm]')
plt.ylabel('$E_{rad}$ erg/cm3')

plt.show()
