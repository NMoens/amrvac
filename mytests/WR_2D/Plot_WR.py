from amrvac_tools.datfiles.reading import amrvac_reader
from scipy import interpolate as itp
import numpy as np
import os
import matplotlib.pyplot as plt

unit_length = 100000.00000000000
unit_numberdensity = 5.9747021551228651e21
unit_temperature = 6061367450.7107716

c = 2.99792458e10
year = 3.1536e7
arad = 7.5657e-15
kb = 1.380648520e-16
mp = 1.6737236e-24
mu = 0.60869565217391308
hd_gamma = 1.6667

wvl = 777363184079.60193
x = np.linspace(0,25*wvl,10000)
