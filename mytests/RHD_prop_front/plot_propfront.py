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
    ds = amrvac_reader.load_file(file)
    ds.get_info()
    ds.units.set_units(unit_length=unit_length, unit_numberdensity=unit_numberdensity, unit_temperature=unit_temperature,He_abundance=0.0)
    if variable == 'time':
        return ds.get_time() * ds.units.unit_time

    ad = ds.load_all_data()

    data = [[0]]

    x,y = ds.get_coordinate_arrays()

    if variable == 'rho':
        data = ad['rho']
        data = data*ds.units.unit_density

    if variable == 'v':
        data = ad['m1']/ad['rho']
        data = data*ds.units.unit_velocity

    if variable == 'e':
        data = ad['e']
        data = data*ds.units.unit_pressure

    if variable == 'p':
        rho = ad['rho']
        v = ad['m1']/ad['rho']
        e = ad['e']

        data = (hd_gamma - 1)*(e - 0.5*rho*v**2)
        data = data*ds.units.unit_pressure

    if variable == 're':
        data = ad['r_e']
        data = data*ds.units.unit_pressure

    if variable == 'D':
        data = ad['D']
        data = data*ds.units.unit_velocity*ds.units.unit_length

    if variable == 'a':
        rho = ad['rho']
        v = ad['m2']/ad['rho']
        e = ad['e']

        data = (hd_gamma - 1)*(e - 0.5*rho*v**2)
        p = data*ds.units.unit_pressure

        data = np.sqrt(p/(ad['rho']*ds.units.unit_density))

    data = np.mean(data,axis=1)
    x = x*ds.units.unit_length
    return x,data,ds.get_time() * ds.units.unit_time


# folder = 'output'
amrvac_outfile0 = '/lhome/nicolasm/amrvac/mytests/RHD_prop_front/output/propfront_pom0000.dat'
amrvac_outfile1 = '/lhome/nicolasm/amrvac/mytests/RHD_prop_front/output/propfront_min0000.dat'
amrvac_outfile2 = '/lhome/nicolasm/amrvac/mytests/RHD_prop_front/output/propfront_pom0010.dat'
amrvac_outfile3 = '/lhome/nicolasm/amrvac/mytests/RHD_prop_front/output/propfront_min0010.dat'
amrvac_outfile4 = '/lhome/nicolasm/amrvac/mytests/RHD_prop_front/output/propfront_pom0020.dat'
amrvac_outfile5 = '/lhome/nicolasm/amrvac/mytests/RHD_prop_front/output/propfront_min0020.dat'

r0,Er0,t0 = AMRVAC_single_profile(amrvac_outfile0,'re')
r1,Er1,t1 = AMRVAC_single_profile(amrvac_outfile1,'re')
r2,Er2,t2 = AMRVAC_single_profile(amrvac_outfile2,'re')
r3,Er3,t3 = AMRVAC_single_profile(amrvac_outfile3,'re')
r4,Er4,t4 = AMRVAC_single_profile(amrvac_outfile4,'re')
r5,Er5,t5 = AMRVAC_single_profile(amrvac_outfile5,'re')

f0 = 0.1 + t0*c
f2 = 0.1 + t2*c
f4 = 0.1 + t4*c

f,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)

ax1.plot(r0,Er0,'bx',label='Pomraning \n t = ' + str(t0))
ax1.plot(r1,Er1,'rx',label='Minerbo \n t = ' + str(t1))
ax1.vlines(f0,0,1)
ax1.set_ylabel('$E_r [erg cm^{-3}]$')
ax1.legend()

ax2.plot(r2,Er2,'bx',label='t = ' + str(t2))
ax2.plot(r3,Er3,'rx',label='t = ' + str(t3))
ax2.vlines(f2,0,1)
ax2.set_ylabel('$E_r [erg cm^{-3}]$')
ax2.legend()

ax3.plot(r4,Er5,'bx',label='t = ' + str(t4))
ax3.plot(r5,Er5,'rx',label='t = ' + str(t5))
ax3.vlines(f4,0,1)
ax3.set_ylabel('$E_r [erg cm^{-3}]$')
ax3.legend()

ax3.set_xlabel('$x/cm$')

plt.show()
