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

#=============================================================================================



folder = 'grid4_output'


amrvac_outfile0 = '/lhome/nicolasm/amrvac/mytests/OutflowAtmosphere4/'+ folder +'/G2m02_d20.00000.dat'
amrvac_outfile1 = '/lhome/nicolasm/amrvac/mytests/OutflowAtmosphere4/'+ folder +'/G2m02_d20.000'+str(40)+'.dat'
amrvac_outfile2 = '/lhome/nicolasm/amrvac/mytests/OutflowAtmosphere4/'+ folder +'/G2m02_d20.000'+str(60)+'.dat'
amrvac_outfile3 = '/lhome/nicolasm/amrvac/mytests/OutflowAtmosphere4/'+ folder +'/G2m02_d20.000'+str(80)+'.dat'
amrvac_outfile4 = '/lhome/nicolasm/amrvac/mytests/OutflowAtmosphere4/'+ folder +'/G2m02_d20.000'+str(90)+'.dat'
amrvac_outfile5 = '/lhome/nicolasm/amrvac/mytests/OutflowAtmosphere4/'+ folder +'/G2m02_d20.00'+str(100)+'.dat'

d1 = AMRVAC_single_profile(amrvac_outfile1,'time')
d2 = AMRVAC_single_profile(amrvac_outfile2,'time')
d3 = AMRVAC_single_profile(amrvac_outfile3,'time')
d4 = AMRVAC_single_profile(amrvac_outfile4,'time')
d5 = AMRVAC_single_profile(amrvac_outfile5,'time')

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

r_sp0, v_sp0 = Get_sonic_point(r0,v0,a0)
r_sp1, v_sp1 = Get_sonic_point(r1,v1,a1)
r_sp2, v_sp2 = Get_sonic_point(r2,v2,a2)
r_sp3, v_sp3 = Get_sonic_point(r3,v3,a3)
r_sp4, v_sp4 = Get_sonic_point(r4,v4,a4)
r_sp5, v_sp5 = Get_sonic_point(r5,v5,a5)

plt.figure(1)
plt.title('Mdot(dinflo)')
plt.grid()
plt.plot(r0/r0[0],Mdot0*year/Msun,'k',label='initial conditions')
plt.plot(r1/r0[0],Mdot1*year/Msun,'rx',label='t = ' + str(d1))
plt.plot(r2/r0[0],Mdot2*year/Msun,'bx',label='t = ' + str(d2))
plt.plot(r3/r0[0],Mdot3*year/Msun,'gx',label='t = ' + str(d3))
plt.plot(r4/r0[0],Mdot4*year/Msun,'mx',label='t = ' + str(d4))
plt.plot(r5/r0[0],Mdot5*year/Msun,'cx',label='t = ' + str(d5))
plt.vlines(r_sp0/r0[0],max(Mdot5*year/Msun),0,'k')
plt.vlines(r_sp1/r0[0],max(Mdot5*year/Msun),0,'r')
plt.vlines(r_sp2/r0[0],max(Mdot5*year/Msun),0,'b')
plt.vlines(r_sp3/r0[0],max(Mdot5*year/Msun),0,'g')
plt.vlines(r_sp4/r0[0],max(Mdot5*year/Msun),0,'m')
plt.vlines(r_sp5/r0[0],max(Mdot5*year/Msun),0,'c')
plt.xlabel('r/R')
plt.ylabel('Mdot')
plt.ylim([-5e-5,5e-4])
plt.legend()

plt.savefig('figs/d_mdot')

plt.figure(2)
plt.title('rho(dinflo)')
plt.grid()
plt.loglog(r0/r0[0],rho0,'k',label='initial conditions')
plt.loglog(r1/r0[0],rho1,'rx',label='t = ' + str(d1))
plt.loglog(r2/r0[0],rho2,'bx',label='t = ' + str(d2))
plt.loglog(r3/r0[0],rho3,'gx',label='t = ' + str(d3))
plt.loglog(r4/r0[0],rho4,'mx',label='t = ' + str(d4))
plt.loglog(r5/r0[0],rho5,'cx',label='t = ' + str(d5))
plt.vlines(r_sp0/r0[0],max(rho5),min(rho0),'k')
plt.vlines(r_sp1/r0[0],max(rho5),min(rho0),'r')
plt.vlines(r_sp2/r0[0],max(rho5),min(rho0),'b')
plt.vlines(r_sp3/r0[0],max(rho5),min(rho0),'g')
plt.vlines(r_sp4/r0[0],max(rho5),min(rho0),'m')
plt.vlines(r_sp5/r0[0],max(rho5),min(rho0),'c')
plt.xlabel('r/R')
plt.ylabel('rho')
plt.legend()

# plt.savefig('figs/d_rho')

plt.figure(3)
plt.title('v(dinflo)')
plt.grid()
plt.plot(r0/r0[0],v0,'black',label='initial conditions')
# plt.plot(r0/r0[0] + 0.2,v0,'gray',label='initial conditions')
# plt.plot(r0/r0[0] + 0.4,v0,'dimgray',label='initial conditions')
# plt.plot(r0/r0[0] + 0.8,v0,'lightgray',label='initial conditions')
plt.plot(r1/r0[0],v1,'rx',label='t = ' + str(d1))
plt.plot(r2/r0[0],v2,'bx',label='t = ' + str(d2))
plt.plot(r3/r0[0],v3,'gx',label='t = ' + str(d3))
plt.plot(r4/r0[0],v4,'mx',label='t = ' + str(d4))
plt.plot(r5/r0[0],v5,'cx',label='t = ' + str(d5))
plt.plot(r0/r0[0],a0,'k--',label='speed of sound')
plt.plot(r1/r0[0],a1,'r--')
plt.plot(r2/r0[0],a2,'b--')
plt.plot(r3/r0[0],a3,'g--')
plt.plot(r4/r0[0],a4,'m--')
plt.plot(r5/r0[0],a5,'c--')
plt.plot(r_sp0/r0[0],v_sp0,'ko')
plt.plot(r_sp1/r0[0],v_sp1,'ro')
plt.plot(r_sp2/r0[0],v_sp2,'bo')
plt.plot(r_sp3/r0[0],v_sp3,'go')
plt.plot(r_sp4/r0[0],v_sp4,'mo')
plt.plot(r_sp5/r0[0],v_sp5,'co')
plt.xlabel('r/R')
plt.ylabel('v')
plt.legend()

# plt.savefig('figs/d_v')

plt.figure(4)
plt.title('Gamma(dinflo)')
plt.grid()
plt.plot(r0/r0[0],Gamma0,'k',label='initial conditions')
plt.plot(r1/r0[0],Gamma1,'rx',label='t = ' + str(d1))
plt.plot(r2/r0[0],Gamma2,'bx',label='t = ' + str(d2))
plt.plot(r3/r0[0],Gamma3,'gx',label='t = ' + str(d3))
plt.plot(r4/r0[0],Gamma4,'mx',label='t = ' + str(d4))
plt.plot(r5/r0[0],Gamma5,'cx',label='t = ' + str(d5))
plt.vlines(r_sp0/r0[0],max(Gamma1),min(Gamma5),'k')
plt.vlines(r_sp1/r0[0],max(Gamma1),min(Gamma5),'r')
plt.vlines(r_sp2/r0[0],max(Gamma1),min(Gamma5),'b')
plt.vlines(r_sp3/r0[0],max(Gamma1),min(Gamma5),'g')
plt.vlines(r_sp4/r0[0],max(Gamma1),min(Gamma5),'m')
plt.vlines(r_sp5/r0[0],max(Gamma1),min(Gamma5),'c')
plt.xlabel('r/R')
plt.ylabel('Gamma')
plt.legend()

# plt.savefig('figs/d_gamma')

plt.show()
#=============================================================================================
