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

    osc = np.sin(x/L_oscillate)
    damp = np.exp((x-wvl)/L_damp)
    sol = osc
    for i in range(len(sol)):
        if x[i] > wvl:
            sol[i] = sol[i]*damp[i]
    return sol

# folder = 'output'
amrvac_outfile0 = '/lhome/nicolasm/amrvac/mytests/RHD_wave/test/fs_hll_lim_mp50100.dat'

r0 , rho0, v0, e0, Er0 = get_data('test/fs_hll_lim_mp50100.blk')

analytic_cheat = analytic(0,0,True)
analytic_mix = analytic(2,0,False)




# plt.figure()
#
# plt.title('Density perturbation')
# plt.plot(r0/unit_length,1e4*(rho0-np.mean(rho0)),'rx',label='simulation')
# plt.plot(x/wvl,analytic_cheat,'k--',label='cheaty')
# plt.plot(x/wvl,analytic_mix,'b--',label='mixed')
# plt.ylabel('perturbation')
# plt.xlabel('$x/\\lambda$')
# plt.legend()



# amrvac_outfile1 = '/lhome/nicolasm/amrvac/mytests/RHD_wave/test/fs_hll_lim_mp50100.dat'
# r1,rho1 = AMRVAC_single_profile(amrvac_outfile1,'rho')
#
# plt.figure()
#
# plt.title('Density perturbation')
# plt.plot(r0/unit_length,1e4*(rho0-np.mean(rho0)),'r-',label='sim, cfl = 0.5')
# plt.plot(r1/unit_length,1e4*(rho0-np.mean(rho1)),'kx',label='sim, cfl = 0.005')
# plt.ylabel('perturbation')
# plt.xlabel('$x/\\lambda$')
# plt.legend()
#


f,(ax1,ax2,ax3,ax4) = plt.subplots(4,1,sharex=True)


ax1.plot(r0,rho0,'rx')
ax1.set_ylabel('$\\delta \\rho$',fontsize=15)

ax2.plot(r0,v0,'rx')
ax2.set_ylabel('$v$',fontsize=15)

e_sc = (e0-0.5*rho0*v0**2)
ax3.plot(r0,e_sc,'rx')
ax3.set_ylabel('$\\delta e_g$',fontsize=15)

Er_cheat = Er0*0.e0
for i in range(len(r0)):
    if r0[i] > 1.0:
        Er_cheat[i] = 0.002/0.025*e_sc[i]

ax4.plot(r0,Er_cheat,'rx')
ax4.set_ylabel('$\\delta E_r$',fontsize=15)
ax4.set_xlabel('$x/\\lambda$',fontsize=15)

plt.tight_layout()

plt.show()
