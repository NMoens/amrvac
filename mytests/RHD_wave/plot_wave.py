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

def AMRVAC_single_profile(file,variable):
    ds = amrvac_reader.load_file(file)
    ds.get_info()
    ds.units.set_units(unit_length=unit_length, unit_numberdensity=unit_numberdensity, unit_temperature=unit_temperature, He_abundance=0.0)
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
    return x,data


# folder = 'output'
amrvac_outfile0 = '/lhome/nicolasm/amrvac/mytests/RHD_shock/output/wave_tau_thick0010.dat'


r0,rho0 = AMRVAC_single_profile(amrvac_outfile0,'rho')
r0,v0 = AMRVAC_single_profile(amrvac_outfile0,'v')
r0,p0 = AMRVAC_single_profile(amrvac_outfile0,'p')
r0,e0 = AMRVAC_single_profile(amrvac_outfile0,'e')
r0,Er0 = AMRVAC_single_profile(amrvac_outfile0,'re')
r0,a0 = AMRVAC_single_profile(amrvac_outfile0,'a')
r0,Diff0 = AMRVAC_single_profile(amrvac_outfile0,'D')

f,(ax1,ax2,ax3,ax4) = plt.subplots(4,1,sharex=True)

ax1.plot(r0,rho0,'rx')
ax1.set_ylabel('$\\rho$')

e_sc = (e0-0.5*rho0*v0**2)/10**14
ax2.plot(r0,e_sc,'rx')
ax2.set_ylabel('$e_g$')

ax3.plot(r0,Er0,'rx')
ax3.set_ylabel('$E_r$')

ax4.plot(r0,v0,'rx')
ax4.set_ylabel('$v$')
ax4.set_xlabel('$x/\\lambda$')


plt.show()
