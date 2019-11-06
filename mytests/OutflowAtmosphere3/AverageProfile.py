from amrvac_tools.datfiles.reading import amrvac_reader
import numpy as np
import os, glob
import matplotlib.pyplot as plt

unit_length=69599000000.0
unit_numberdensity=729723637293892.50
unit_temperature=189666.48683552662

dir_path = os.getcwd() + '/output/'
files = list(glob.glob(os.path.join(dir_path, 'const*.dat')))
files.sort()

ds = amrvac_reader.load_file(files[0])
ds.units.set_units(unit_length=unit_length, unit_numberdensity=unit_numberdensity, unit_temperature=unit_temperature)
ad = ds.load_all_data()
orig_variables = ds.get_varnames()
x,y = ds.get_coordinate_arrays()
nx = len(x)
ny = len(y)

def calc_qtt(ad,qtt):

    r = []
    for i in range(nx):
        r.append(y)
    r = np.array(r)

    if qtt in orig_variables:
        data = ad[qtt]
    elif qtt == 'v':
        data = ad['m2']/ad['rho']
    elif qtt == 'mdot':
        data = 4*np.pi*r**2*ad['m2']
    else:
        print('Not implemented yet')
    return data

def average_quantity(qtt):

    # Initialise time iteration counter
    nt = 0

    # Initialise average snap array
    average_snap = np.zeros((nx,ny))

    #Loop over timesnaps
    for f in files:
        nt = nt+1
        ds = amrvac_reader.load_file(f)
        ds.units.set_units(unit_length=unit_length, unit_numberdensity=unit_numberdensity, unit_temperature=unit_temperature)
        ad = ds.load_all_data()
        average_snap = average_snap + calc_qtt(ad,qtt)
    average_snap = average_snap/nt

    # Initialize average profile array
    average_profile = np.mean(average_snap,axis=0)
    return average_profile

# rho_mean = average_quantity('rho')
mdot_mean = average_quantity('mdot')
M_sun = 1.9891000e33

unit_mdot = ds.units.unit_density*ds.units.unit_length**3/ds.units.unit_time*(365.25*24*60*60)/M_sun

mdot_mean = mdot_mean*unit_mdot

# plt.plot(y,rho_mean)

plt.plot(y[1:],np.log10(mdot_mean[1:]))

plt.title('Mass Loss rate')
plt.xlabel('r/r_star')
plt.ylabel('log Mdot in Msun/yr')
plt.show()
