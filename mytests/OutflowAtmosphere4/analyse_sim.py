from amrvac_tools.datfiles.reading import amrvac_reader
from scipy import interpolate as itp
import numpy as np
import os
import matplotlib.pyplot as plt

unit_length=69599000000.0
unit_numberdensity=2218242342924238.2
unit_temperature=296199.82122247218

R = 1
M = 10
L = 1.9e5

Rsun = 6.96e10
Msun = 1.99e33
Lsun = 3.9e33
G = 6.67259e-8
c = 2.99792458e10
year = 356.25*24*60*60
arad = 7.5646e-15
kb = 1.380658e-16
mp = 1.6726231e-24
mu = 0.6
hd_gamma = 5.0/3.0

def AMRVAC_single_profile(file,variable):
    ds = amrvac_reader.load_file(file)
    ds.get_info()
    ds.units.set_units(unit_length=unit_length, unit_numberdensity=unit_numberdensity, unit_temperature=unit_temperature)
    if variable == 'time':
        return ds.get_time() * ds.units.unit_time

    ad = ds.load_all_data()

    data = [[0]]

    x,y = ds.get_coordinate_arrays()

    if variable == 'rho':
        data = ad['rho']
        data = data*ds.units.unit_density

    if variable == 'v':
        data = ad['m2']/ad['rho']
        data = data*ds.units.unit_velocity

    if variable == 'e':
        data = ad['e']
        data = data*ds.units.unit_pressure

    if variable == 'p':
        rho = ad['rho']
        v = ad['m2']/ad['rho']
        e = ad['e']

        data = (hd_gamma - 1)*(e - 0.5*rho*v**2)
        data = data*ds.units.unit_pressure

    if variable == 're':
        data = ad['r_e']
        data = data*ds.units.unit_pressure

    if variable == 'a':
        rho = ad['rho']
        v = ad['m2']/ad['rho']
        e = ad['e']

        data = (hd_gamma - 1)*(e - 0.5*rho*v**2)
        p = data*ds.units.unit_pressure

        data = np.sqrt(p/(ad['rho']*ds.units.unit_density))

    if variable == 'Mdot':
        unit_Mdot = ds.units.unit_length**3*ds.units.unit_density/ds.units.unit_time
        data = 4*np.pi*ad['m2']
        data = np.mean(data,axis=0)
        data = data*y**2*unit_Mdot
    else:
        data = np.mean(data,axis=0)

    if variable == 'Gamma':
        E = ad['r_e']
        gradE = np.mean(E,axis=0)

        gradE = np.gradient(gradE,y)

        grad = -gradE/(3*np.mean(ad['rho'],axis=0))

        ggrav = G*M*Msun/(y*ds.units.unit_length)**2\
        *(ds.units.unit_time**2/ds.units.unit_length)
        data = grad/ggrav

    if variable == 'grada':
        data = ad['r_e']
        data = np.mean(data,axis=0)
        data = np.gradient(data,y)

    y = y*ds.units.unit_length
    return y,data


def Get_sonic_point(r,v,a):

    tck_v = itp.splrep(r,v)
    tck_a = itp.splrep(r,a)

    rnew = np.linspace(0,3*r[0],1000)
    vnew = itp.splev(rnew,tck_v)
    anew = itp.splev(rnew,tck_a)

    fitness = abs(vnew - anew)

    i_sp = np.argmin(fitness)
    r_sp = rnew[i_sp]
    v_sp = itp.splev(r_sp,tck_v)
    print(r_sp, v_sp)

    return r_sp, v_sp


folder = 'grid4_output'
amrvac_outfile0 = '/lhome/nicolasm/amrvac/mytests/OutflowAtmosphere4/'+ folder +'/G2m02_d20.00'+str(100)+'.dat'

r0,Mdot0 = AMRVAC_single_profile(amrvac_outfile0,'Mdot')
r0,rho0 = AMRVAC_single_profile(amrvac_outfile0,'rho')
r0,v0 = AMRVAC_single_profile(amrvac_outfile0,'v')
r0,p0 = AMRVAC_single_profile(amrvac_outfile0,'p')
r0,Er0 = AMRVAC_single_profile(amrvac_outfile0,'re')
r0,a0 = AMRVAC_single_profile(amrvac_outfile0,'a')
r0,Gamma0 = AMRVAC_single_profile(amrvac_outfile0,'Gamma')
r_sp0, v_sp0 = Get_sonic_point(r0,v0,a0)


#Energy=energy density*volume
#volume element ~ r**2
U_pot = -rho0*G*M*Msun/r0 *r0**2
# U_pot = U_pot - min(U_pot)
E_k = 1./2.*rho0*v0**2 *r0**2
E_int = p0/(hd_gamma-1) *r0**2
E_rad = Er0 *r0**2
E_tot = U_pot + E_k + E_int + E_rad


#photon-tiring term:
gradV = np.gradient(v0,r0)
Ptt = gradV*Er0/3
Ptt_pp = 2./3.*v0*Er0/r0
Ptt_ppc = Ptt + Ptt_pp


#Heating and cooling terms:
error_b = 10.0
kappa_b = 0.61885728378588867
kappa_0 = 1.3752384084130858
kappa = kappa_b + (1+np.erf((r0-R)/R*error_b-error_b/2))*(kappa_0-kappa_b)/2
Tg = p0*mp*mu/(kb*rho0)
cool = c*kappa*rho0*arad*Tg**4
heat = c*kappa*rho0*Er0
q_dot = heat-cool

# plt.figure()
# plt.semilogy(r0/Rsun,U_pot,'r--',label='Potential energy')
# plt.semilogy(r0/Rsun,E_k,'g--',label='Kinetic energy')
# plt.semilogy(r0/Rsun,E_int,'b--',label='Internal energy')
# plt.semilogy(r0/Rsun,E_rad,'c--',label='Radiation energy')
# plt.semilogy(r0/Rsun,E_tot,'k-',label='Total energy')
# plt.legend()

plt.figure()
plt.plot(r0/Rsun,Ptt,'r-',label='photon-tiring term')
plt.plot(r0/Rsun,Ptt_pp,'b-',label='pseudo-planar correction')
plt.plot(r0/Rsun,Ptt_ppc,'k-',label='pseudo-planar corrected \n photon-tiring term')
plt.legend()

plt.figure()
plt.semilogy(r0/Rsun,cool,'r-',label='cooling term')
plt.semilogy(r0/Rsun,heat,'b-',label='heating term')
plt.semilogy(r0/Rsun,q_dot,'k-',label='Net energy exch for gas')
plt.legend()


plt.show()
