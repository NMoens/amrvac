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
x_a = np.linspace(0,28*wvl,10000)

def get_data(myfile):
    A = np.loadtxt(myfile,skiprows=3)

    x0 = np.transpose(A)[0]
    rho0 = np.transpose(A)[1]
    v0 = np.transpose(A)[2]
    e0 = np.transpose(A)[3]
    E0 = np.transpose(A)[4]

    delta_rho = 100*(rho0-np.mean(rho0))
    delta_v = v0
    delta_e = e0-np.mean(e0)
    delta_E = E0-np.mean(E0)

    return x0, delta_rho, delta_v, delta_e, delta_E

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

    data_1D = np.zeros(len(x))
    for i in range(len(x)):
        data_1D[i] = data[i][i]
    x = np.sqrt(x**2+y**2)*ds.units.unit_length
    return x,data_1D

def analytic(ii,jj,cheat):
    a = 3.e6
    omega = 1.8769477562401875e-005
    wvl = 777363184079.60193

    tau_a = 1.e3/(2*np.pi)
    tau_c = 1.e4*tau_a
    Bo = 1.e-3
    r = 0.1
    hd_gamma = 5./3.

    z4 = (1. - 1j*16*tau_a/Bo)
    z3 = 0.
    z2 = (3*tau_a**2*(1+1j/tau_c)**2 - 1. + 1j*(16*hd_gamma*tau_a/Bo) + (16*r)*tau_a**2*(1+1j/tau_c)*(5+1j/(3*(hd_gamma-1)*tau_c) +16./3*r*hd_gamma/(hd_gamma-1)))
    z1 = 0.
    z0 = -3*tau_a**2*((1+1j/tau_c)**2 + 16*hd_gamma*r*(1 + 1j/tau_c))
    coeffs = [z4, z3, z2, z1, z0]
    z = np.roots(coeffs)
    k = omega*z/a

    L_damp = 1./np.imag(k[ii])
    L_oscillate = 1./np.real(k[jj])

    if cheat:
        L_damp = L_damp/(2*np.pi)

    osc = np.sin(x_a/L_oscillate)
    damp = np.exp((x_a-wvl)/L_damp)
    sol = osc
    for i in range(len(sol)):
        if x_a[i] > wvl:
            sol[i] = sol[i]*damp[i]
    return sol

# folder = 'output'
amrvac_outfile0 = '/lhome/nicolasm/amrvac/mytests/RHD_wave_diag/output/2wave_thick0039.dat'


r0,rho0 = AMRVAC_single_profile(amrvac_outfile0,'rho')
# r0,v0 = AMRVAC_single_profile(amrvac_outfile0,'v')
# r0,p0 = AMRVAC_single_profile(amrvac_outfile0,'p')
# r0,e0 = AMRVAC_single_profile(amrvac_outfile0,'e')
# r0,Er0 = AMRVAC_single_profile(amrvac_outfile0,'re')
# r0,a0 = AMRVAC_single_profile(amrvac_outfile0,'a')
# r0,Diff0 = AMRVAC_single_profile(amrvac_outfile0,'D')

r_1D , rho_1D, v_1D, e_1D, Er_1D = get_data('/lhome/nicolasm/amrvac/mytests/RHD_wave/convergence/Midpoint_mp5_4000039.blk')


analytic = analytic(0,0,False)


plt.figure()

# plt.title('Density perturbation')
plt.grid()
plt.plot(x_a/wvl,analytic,'k--',label='analytic')
plt.plot((r0)/unit_length+0.25,1e4*(rho0-np.mean(rho0)),'r-',label='Diagonal')
plt.plot(r_1D,(rho_1D-np.mean(rho_1D)),'b-',label='1D')
plt.ylabel('$(\\rho - \\rho_0) / A_\\rho$',fontsize=10)
plt.xlabel('$x/\\lambda$',fontsize=10)
plt.legend()


# f,(ax1,ax2,ax3,ax4) = plt.subplots(4,1,sharex=True)
#
#
# ax1.plot(r0,rho0,'rx')
# ax1.set_ylabel('$\\rho$')
#
# e_sc = (e0-0.5*rho0*v0**2)
# ax2.plot(r0,e_sc,'rx')
# ax2.set_ylabel('$e_g$')
#
# ax3.plot(r0,Er0,'rx')
# ax3.set_ylabel('$E_r$')
#
# ax4.plot(r0,v0,'rx')
# ax4.set_ylabel('$v$')
# ax4.set_xlabel('$x/\\lambda$')


plt.show()
