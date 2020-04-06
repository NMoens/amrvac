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
G = 6.67191e-8
c = 2.99792458e10
year = 3.1536e7
arad = 7.5657e-15
kb = 1.380648520e-16
mp = 1.6737236e-24
mu = 0.60869565217391308
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


# folder = 'output'
amrvac_outfile0 = '/lhome/nicolasm/amrvac/mytests/OutflowAtmosphere4/grid4_output/G2m02_d20.00'+str(100)+'.dat'
amrvac_outfile0 = '/lhome/nicolasm/amrvac/mytests/OutflowAtmosphere4/output/G2m02_3L00'+str(36)+'.dat'

r0,Mdot0 = AMRVAC_single_profile(amrvac_outfile0,'Mdot')
r0,rho0 = AMRVAC_single_profile(amrvac_outfile0,'rho')
r0,v0 = AMRVAC_single_profile(amrvac_outfile0,'v')
r0,p0 = AMRVAC_single_profile(amrvac_outfile0,'p')
r0,e0 = AMRVAC_single_profile(amrvac_outfile0,'e')
r0,Er0 = AMRVAC_single_profile(amrvac_outfile0,'re')
r0,a0 = AMRVAC_single_profile(amrvac_outfile0,'a')
r0,Diff1 = AMRVAC_single_profile(amrvac_outfile0,'D')
r0,Gamma0 = AMRVAC_single_profile(amrvac_outfile0,'Gamma')
r_sp0, v_sp0 = Get_sonic_point(r0,v0,a0)

#photon-tiring term:
gradV = np.gradient(v0,r0)
Ptt = gradV*Er0/3/Er0
Ptt_pp = 2./3.*v0*Er0/r0/Er0
Ptt_ppc = Ptt + Ptt_pp

#Heating and cooling terms, relative to Er:
error_b = 10.0#*R*Rsun
kappa_b = 0.61885728378588867
kappa_0 = 1.3752384084130858
kappa = kappa_b + (1 + np.erf((r0/(R*Rsun) - 1)*error_b-error_b/2))*(kappa_0-kappa_b)/2
Tg = p0*mp*mu/(kb*rho0)
cool = c*kappa*rho0*arad*Tg**4/Er0
heat = c*kappa*rho0*Er0/Er0
q_dot = -heat+cool

#Advection flux Div(Ev + F), relative to Er
Diff2 = c/(3*kappa*rho0)
F = -Diff2*np.gradient(Er0,r0)
Div_sph_F = 1/r0**2*np.gradient(r0**2*F,r0)/Er0
Div_sph_Ev = 1/r0**2*np.gradient(r0**2*Er0*v0,r0)/Er0

#PseudoPlanar correction for advection flux
DivF = np.gradient(F,r0)/Er0
DivEv = np.gradient(Er0*v0,r0)/Er0
ppc_DivF = -2*F/r0/Er0
ppc_DivEv = -2*v0*Er0/r0/Er0
DivF_ppc = DivF - ppc_DivF
DivEv_ppc = DivEv - ppc_DivEv

#Steady state Continuity equation:
div_rhov = np.gradient(rho0*v0,r0)/rho0
div_sph_rhov = 1/r0**2*np.gradient(r0**2*rho0*v0,r0)/rho0
ppc_rhov = -2*rho0*v0/r0/rho0

#Steady state momentum equation:
div_vrhov_p = np.gradient(rho0*v0**2+p0,r0)/(rho0*v0)
ppc_vrhov_p = -2*(rho0*v0**2+p0)/r0/(rho0*v0)
g_rad_grav = rho0*(kappa*F/c-G*M*Msun/r0**2)/(rho0*v0)

#Steady state gas-energy equation:
div_ev_pv = np.gradient(e0*v0+p0*v0,r0)/e0
ppc_ev_pv = -2*(e0*v0+p0*v0)/r0/e0
vg_rad_grav = v0*rho0*(kappa*F/c-G*M*Msun/r0**2)/e0
e_qdot = -cool+heat

#Steady state radiation-energy equation:
div_F_Ev = np.gradient(F+Er0*v0,r0)/Er0
ppc_F_Ev = -2*(F+Er0*v0)/r0/Er0

#mass loss parameter
Mdot = 4*np.pi*r0**2*rho0*v0
L_adv = 4*np.pi*r0**2*(4./3.*Er0*v0)
L_cmf = 4*np.pi*r0**2*F
L_obs = L_adv+L_cmf

m = Mdot*G*M*Msun/(r0*L_obs)
m1 = Mdot[0]*G*M*Msun/(r0[0]*L_obs[0])
m2 = max(Mdot)*G*M*Msun/(r0[0]*L_obs[0])


plt.figure()
plt.title('Photon tiring term')
plt.plot(r0/Rsun,Ptt,'r-',label='photon-tiring term')
plt.plot(r0/Rsun,Ptt_pp,'b-',label='pseudo-planar correction')
plt.plot(r0/Rsun,Ptt_ppc,'k-',label='pseudo-planar corrected \n photon-tiring term')
plt.legend()

plt.figure()
plt.title('Heating and Cooling')
plt.plot(r0/Rsun,cool,'r-',label='cooling term')
plt.plot(r0/Rsun,heat,'b-',label='heating term')
plt.plot(r0/Rsun,abs(q_dot),'k-',label='|Net energy exch for gas|')
# plt.plot(r0/Rsun,abs(heat/cool),'k--',label='ratio')
plt.legend()

plt.figure()
plt.title('Radiation energy source terms')
plt.plot(r0/Rsun,q_dot,'r.',label='cooling-heating')
plt.plot(r0/Rsun,Ptt_ppc,'r--',label='photon tiring term')
plt.plot(r0/Rsun,Div_sph_F,'b.',label='$\\nabla \cdot \\vec{F}$')
plt.plot(r0/Rsun,Div_sph_Ev,'b--',label='$\\nabla \cdot E\\vec{v}$')
plt.plot(r0/Rsun,Div_sph_F+Div_sph_Ev-q_dot-Ptt_ppc,'k-',label='$\\nabla \cdot (Ev+F)$ \n $= \dot{q} - \\nabla v : P$')
plt.legend()

plt.figure()
plt.title('PseudoPlanar correction E_rad')
plt.plot(r0/Rsun,DivF,'r.',label='$\\nabla_c \cdot F$')
plt.plot(r0/Rsun,DivEv,'b.',label='$\\nabla_c \cdot \\vec{v}E$')
plt.plot(r0/Rsun,Div_sph_F,'r--',label='$\\nabla_{sph} \cdot F$')
plt.plot(r0/Rsun,Div_sph_Ev,'b--',label='$\\nabla_{sph} \cdot \\vec{v}E$')
plt.plot(r0/Rsun,ppc_DivF,'r-.',label='$S_{ppc}^F$')
plt.plot(r0/Rsun,ppc_DivEv,'b-.',label='$S_{ppc}^{Ev}$')
plt.plot(r0/Rsun,DivEv_ppc - Div_sph_F,'r-',label='$\\nabla_c \cdot F + S_{ppc}^F - \\nabla_{sph} \cdot \\vec{F}$')
plt.plot(r0/Rsun,DivF_ppc - Div_sph_Ev,'b-',label='$\\nabla_c \cdot Ev + S_{ppc}^{vE} - \\nabla_{sph} \cdot \\vec{v}E$')
plt.plot(r0/Rsun,DivEv_ppc - Div_sph_F+DivF_ppc - Div_sph_Ev,'k-.',label='total')
plt.legend()

plt.figure()
plt.title('Steady state continuity equation')
plt.plot(r0/Rsun,div_rhov,'r',label='$\\nabla \cdot (\\rho v)$')
plt.plot(r0/Rsun,div_sph_rhov,'k-.',label='$\\nabla_s \cdot (\\rho v)$')
plt.plot(r0/Rsun,ppc_rhov,'b',label='$S_{ppc}^{\\rho}$')
plt.plot(r0/Rsun,div_rhov-ppc_rhov,'k',label='total')
plt.legend()

plt.figure()
plt.title('Steady state momentum equation')
plt.plot(r0/Rsun,div_vrhov_p,'r',label='$\\nabla(v\\rho v + p)$')
plt.plot(r0/Rsun,ppc_vrhov_p,'b',label='$S_{ppc}^{v \\rho}$')
plt.plot(r0/Rsun,g_rad_grav,'g',label='$g_{rad}-g_{grav}$')
plt.plot(r0/Rsun,div_vrhov_p-ppc_vrhov_p-g_rad_grav,'k',label='total')
plt.legend()

plt.figure()
plt.title('Steady state gas-energy equation')
plt.plot(r0/Rsun,div_ev_pv,'r',label='$\\nabla(ev+pv)$')
plt.plot(r0/Rsun,ppc_ev_pv,'b',label='$S_{ppc}^e$')
plt.plot(r0/Rsun,vg_rad_grav,'g',label='$vg_{rad}-vg_{grav}$')
plt.plot(r0/Rsun,e_qdot,'c',label='$\dot{q}$')
plt.plot(r0/Rsun,div_ev_pv-ppc_ev_pv-vg_rad_grav-e_qdot,'k',label='total')
plt.legend()

plt.figure()
plt.title('Steady state radiation-energy equation')
plt.plot(r0/Rsun,div_F_Ev,'r',label='$\\nabla(Ev+F)$')
plt.plot(r0/Rsun,ppc_F_Ev,'b',label='$S_{ppc}^E$')
plt.plot(r0/Rsun,Ptt_ppc,'g',label='$\\nabla_s v: P$')
plt.plot(r0/Rsun,q_dot,'c',label='$-\dot{q}$')
plt.plot(r0/Rsun,div_F_Ev-ppc_F_Ev+Ptt_ppc-q_dot,'k',label='total')
plt.legend()

plt.figure()
plt.title('Opacity')
plt.plot(r0/Rsun,kappa,'r.')

plt.figure()
plt.title('check diffusion coefficient')
plt.plot(r0/Rsun,Diff1/1e20,'r-',label='amrvac readout: $D$')
plt.plot(r0/Rsun,Diff2/1e20,'b-',label='recalculated: $\\frac{c\\lambda}{\\kappa \\rho}$')
plt.plot(r0/Rsun,Diff1/Diff2-1,'k--', label='relative difference')
plt.legend()

plt.figure()
plt.title('mass loss parameter')
plt.plot(r0/Rsun,m,'r',label= 'Local value' )
plt.hlines(m1,r0[0]/Rsun,r0[-1]/Rsun,'r','--',label= 'Value at first cell')
plt.hlines(m2,r0[0]/Rsun,r0[-1]/Rsun,'b','--',label= 'Value using max(Mdot)')
plt.legend()




plt.show()
