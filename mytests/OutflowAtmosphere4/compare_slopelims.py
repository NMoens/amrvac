from amrvac_tools.datfiles.reading import amrvac_reader
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

    # if variable == 'grada':
    #     rho = ad['rho']
    #     v = ad['m2']/ad['rho']
    #     e = ad['e']
    #
    #     data = (hd_gamma - 1)*(e - 0.5*rho*v**2)
    #     p = data*ds.units.unit_pressure
    #
    #     data = np.sqrt(p/(ad['rho']*ds.units.unit_density))
    #     data = np.mean(data,axis=0)
    #     data = np.gradient(data,y)

    if variable == 'grada':
        data = ad['r_e']
        data = np.mean(data,axis=0)
        data = np.gradient(data,y)

    y = y*ds.units.unit_length
    return y,data

#=============================================================================================

# koren, mp5, woodward crash

d1 = 50
d2 = 60
d3 = 70
d4 = 80
d5 = 90

it = 40

amrvac_outfile0 = '/lhome/nicolasm/amrvac/mytests/OutflowAtmosphere4/test_slopelims/G2m02_minmod00'+str(it)+'.dat'
amrvac_outfile1 = '/lhome/nicolasm/amrvac/mytests/OutflowAtmosphere4/test_slopelims/G2m02_cada00'+str(it)+'.dat'
amrvac_outfile2 = '/lhome/nicolasm/amrvac/mytests/OutflowAtmosphere4/test_slopelims/G2m02_cada300'+str(it)+'.dat'
amrvac_outfile3 = '/lhome/nicolasm/amrvac/mytests/OutflowAtmosphere4/test_slopelims/G2m02_weno300'+str(it)+'.dat'
amrvac_outfile4 = '/lhome/nicolasm/amrvac/mytests/OutflowAtmosphere4/test_slopelims/G2m02_ppm00'+str(it+10)+'.dat'
amrvac_outfile5 = '/lhome/nicolasm/amrvac/mytests/OutflowAtmosphere4/test_slopelims/G2m02_vanleer00'+str(it)+'.dat'


r0,Mdot0 = AMRVAC_single_profile(amrvac_outfile0,'Mdot')
r1,Mdot1 = AMRVAC_single_profile(amrvac_outfile1,'Mdot')
r2,Mdot2 = AMRVAC_single_profile(amrvac_outfile2,'Mdot')
r3,Mdot3 = AMRVAC_single_profile(amrvac_outfile3,'Mdot')
r4,Mdot4 = AMRVAC_single_profile(amrvac_outfile4,'Mdot')
r5,Mdot5 = AMRVAC_single_profile(amrvac_outfile5,'Mdot')

r0,rho0 = AMRVAC_single_profile(amrvac_outfile0,'rho')
r1,rho1 = AMRVAC_single_profile(amrvac_outfile1,'rho')
r2,rho2 = AMRVAC_single_profile(amrvac_outfile2,'rho')
r3,rho3 = AMRVAC_single_profile(amrvac_outfile3,'rho')
r4,rho4 = AMRVAC_single_profile(amrvac_outfile4,'rho')
r5,rho5 = AMRVAC_single_profile(amrvac_outfile5,'rho')

r0,v0 = AMRVAC_single_profile(amrvac_outfile0,'v')
r1,v1 = AMRVAC_single_profile(amrvac_outfile1,'v')
r2,v2 = AMRVAC_single_profile(amrvac_outfile2,'v')
r3,v3 = AMRVAC_single_profile(amrvac_outfile3,'v')
r4,v4 = AMRVAC_single_profile(amrvac_outfile4,'v')
r5,v5 = AMRVAC_single_profile(amrvac_outfile5,'v')

r0,a0 = AMRVAC_single_profile(amrvac_outfile0,'a')
r1,a1 = AMRVAC_single_profile(amrvac_outfile1,'a')
r2,a2 = AMRVAC_single_profile(amrvac_outfile2,'a')
r3,a3 = AMRVAC_single_profile(amrvac_outfile3,'a')
r4,a4 = AMRVAC_single_profile(amrvac_outfile4,'a')
r5,a5 = AMRVAC_single_profile(amrvac_outfile5,'a')

r0,Gamma0 = AMRVAC_single_profile(amrvac_outfile0,'Gamma')
r1,Gamma1 = AMRVAC_single_profile(amrvac_outfile1,'Gamma')
r2,Gamma2 = AMRVAC_single_profile(amrvac_outfile2,'Gamma')
r3,Gamma3 = AMRVAC_single_profile(amrvac_outfile3,'Gamma')
r4,Gamma4 = AMRVAC_single_profile(amrvac_outfile4,'Gamma')
r5,Gamma5 = AMRVAC_single_profile(amrvac_outfile5,'Gamma')

plt.figure(1)
plt.title('Mdot(dinflo)')
plt.grid()
plt.plot(r0/r0[0],Mdot0*year/Msun,'kx',label='minmod')
plt.plot(r1/r0[0],Mdot1*year/Msun,'rx',label='cada')
plt.plot(r2/r0[0],Mdot2*year/Msun,'bx',label='cada3')
plt.plot(r3/r0[0],Mdot3*year/Msun,'gx',label='weno3')
plt.plot(r4/r0[0],Mdot4*year/Msun,'mx',label='ppm')
plt.plot(r5/r0[0],Mdot5*year/Msun,'cx',label='vanleer')
plt.xlabel('r/R')
plt.ylabel('Mdot')
plt.ylim([-5e-5,5e-4])
plt.legend()

plt.savefig('figs/lim_mdot')

plt.figure(2)
plt.title('rho(dinflo)')
plt.grid()
plt.loglog(r0/r0[0],rho0,'kx',label='minmod')
plt.loglog(r1/r0[0],rho1,'rx',label='cada')
plt.loglog(r2/r0[0],rho2,'bx',label='cada3')
plt.loglog(r3/r0[0],rho3,'gx',label='weno3')
plt.loglog(r4/r0[0],rho4,'mx',label='ppm')
plt.loglog(r5/r0[0],rho5,'cx',label='vanleer')
plt.xlabel('r/R')
plt.ylabel('rho')
plt.legend()

plt.savefig('figs/lim_rho')

plt.figure(3)
plt.title('v(dinflo)')
plt.grid()
plt.plot(r0/r0[0],v0,'kx',label='minmod')
plt.plot(r1/r0[0],v1,'rx',label='cada')
plt.plot(r2/r0[0],v2,'bx',label='cada3')
plt.plot(r3/r0[0],v3,'gx',label='weno3')
plt.plot(r4/r0[0],v4,'mx',label='ppm')
plt.plot(r5/r0[0],v5,'cx',label='vanleer')
plt.plot(r0/r0[0],a0,'k--',label='speed of sound')
plt.plot(r1/r0[0],a1,'r--')
plt.plot(r2/r0[0],a2,'b--')
plt.plot(r3/r0[0],a3,'g--')
plt.plot(r4/r0[0],a4,'m--')
plt.plot(r5/r0[0],a5,'c--')
plt.xlabel('r/R')
plt.ylabel('v')
plt.legend()

plt.savefig('figs/lim_v')

plt.figure(4)
plt.title('Gamma(dinflo)')
plt.grid()
plt.plot(r0/r0[0],Gamma0,'kx',label='minmod')
plt.plot(r1/r0[0],Gamma1,'rx',label='cada')
plt.plot(r2/r0[0],Gamma2,'bx',label='cada3')
plt.plot(r3/r0[0],Gamma3,'gx',label='weno3')
plt.plot(r4/r0[0],Gamma4,'mx',label='ppm')
plt.plot(r5/r0[0],Gamma5,'cx',label='vanleer')

plt.plot(r0[0:3]/r0[0],Gamma0[0:3],'k-')
plt.plot(r1[0:3]/r0[0],Gamma1[0:3],'r-')
plt.plot(r2[0:3]/r0[0],Gamma2[0:3],'b-')
plt.plot(r3[0:3]/r0[0],Gamma3[0:3],'g-')
plt.plot(r4[0:3]/r0[0],Gamma4[0:3],'m-')
plt.plot(r5[0:3]/r0[0],Gamma5[0:3],'c-')

plt.xlabel('r/R')
plt.ylabel('Gamma')
plt.legend()

plt.savefig('figs/lim_gamma')

plt.show()
