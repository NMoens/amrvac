from amrvac_tools.datfiles.reading import amrvac_reader
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import special as sp

R = 1
M = 10
L = 1.9e5
Mdot = 0.00012015046769275339 #8.3e-6

kappa_0 = 306.57376321732482
kappa_b = 223.87947156182562

m = 0.2
Gamma = 2.0


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

unit_length=69599000000.0
unit_numberdensity=2218242342924238.2
unit_temperature=296199.82122247218


print(Mdot*G*M*Msun*Msun/year/(R*Rsun*L*Lsun))

v_esc = (2*G*M*Msun/(R*Rsun))**0.5
print('v_esc',v_esc)
print('v_inf',np.sqrt(2.20441348 -1)*v_esc)

####################################################
### Stuff for Amrvac profile
####################################################

def AMRVAC_profile(file,variable):
    ds = amrvac_reader.load_file(file)
    ds.get_info()
    ds.units.set_units(unit_length=unit_length, unit_numberdensity=unit_numberdensity, unit_temperature=unit_temperature)
    ad = ds.load_all_data()

    if variable == 'rho':
        data = ad['int_r']/ad['int_dt']
        data = data*ds.units.unit_density

    if variable == 'v':
        data = ad['int_v']/ad['int_dt']
        data = data*ds.units.unit_velocity

    if variable == 'e':
        data = ad['int_e']/ad['int_dt']
        data = data*ds.units.unit_pressure

    if variable == 'p':
        rho = ad['int_r']/ad['int_dt']
        v = ad['int_v']/ad['int_dt']
        e = ad['int_e']/ad['int_dt']

        data = (hd_gamma - 1)*(e - 0.5*rho*v**2)
        data = data*ds.units.unit_pressure

    if variable == 're':
        data = ad['int_re']/ad['int_dt']
        data = data*ds.units.unit_pressure

    x,y = ds.get_coordinate_arrays()

    if variable == 'kappa':
        unit_kappa = 1/(ds.units.unit_length*ds.units.unit_density)

        print(unit_kappa)
        data = kappa_b + sp.erf((y-1)*10)*(kappa_0-kappa_b)
        data = data*unit_kappa
    else:
        data = np.mean(data,axis=0)


    y = y*ds.units.unit_length

    return y,data

def AMRVAC_single_profile(file,variable):
    ds = amrvac_reader.load_file(file)
    ds.get_info()
    ds.units.set_units(unit_length=unit_length, unit_numberdensity=unit_numberdensity, unit_temperature=unit_temperature)
    ad = ds.load_all_data()

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

    x,y = ds.get_coordinate_arrays()

    if variable == 'kappa':
        unit_kappa = 1/(ds.units.unit_length*ds.units.unit_density)
        kappa_0 = 306.57376321732482
        kappa_b = 128.78925554545540
        print(unit_kappa)
        data = kappa_b + sp.erf((y-1)*100)*(kappa_0-kappa_b)
        data = data*unit_kappa
    else:
        data = np.mean(data,axis=0)


    y = y*ds.units.unit_length

    return y,data

####################################################
### Stuff for Gamma,m - profile
####################################################

class MyWindModel():
    def __init__(self,fname,Gamma,m):
        self.fname = fname
        self.wg,self.p,self.x = np.loadtxt(self.fname,unpack=True)
        self.Gamma = Gamma
        self.w = self.wg*(self.Gamma-1)
        self.m = m
        self.eta = 4*self.p*np.sqrt(self.w)/(self.m*(1-self.x)**2)

    def SetStarParams(self,R,M,L):
        self.R = R*Rsun
        self.M = M*Msun
        self.L = L*Lsun
        self.Mdot = self.m*self.R*self.L/(G*self.M) #8.3e-6*Msun/year

        self.v_esc = (2*G*self.M/self.R)**0.5

        self.r = self.R/(1-self.x)
        self.v = np.sqrt(self.w)*self.v_esc
        self.Prad = self.p*self.L/(4*np.pi*self.R**2*self.v_esc)
        self.rho = self.Mdot/(4*np.pi*self.r**2*self.v)
        self.Erad = 3*self.Prad
        self.T = abs(self.Erad/arad)**(1./4)
        self.pg = kb*self.T/(mp*mu) *self.rho
        self.e = self.pg/(hd_gamma-1)+0.5*self.v**2*self.rho

        self.Ltot = (1-self.m*(self.w+self.x))*self.L
        self.Ladv = (self.m*self.eta)*self.L
        self.Ldiff = (1-self.m*(self.w+self.x+self.eta))*self.L

def SE_UNIFIED_profile(file,variable):

    Mwm = MyWindModel(file,Gamma,m)
    Mwm.SetStarParams(R,M,L)

    if variable == 'rho':
        data = Mwm.rho
    if variable == 'v':
        data = Mwm.v
    if variable == 'e':
        data = Mwm.e
    if variable == 'p':
        data = Mwm.pg
    if variable == 're':
        data = Mwm.Erad

    return Mwm.r, data


def Get_SE_Luminosities(file):
    m = 0.1
    Gamma = 2.0

    Mwm = MyWindModel(file,Gamma,m)
    Mwm.SetStarParams(R,M,L)

    return Mwm.r, Mwm.Ltot, Mwm.Ladv, Mwm.Ldiff

####################################################
### General stuff
####################################################

def GetGamma(r,rho,re):
    grad_re = np.gradient(re,r)
    g_rad = grad_re/(3*rho)
    g_gra = G*M*Msun/r**2

    Gamma = -g_rad/g_gra
    return Gamma

def GetRadflux(r,rho,re,kappa):
    grad_re = np.gradient(re,r)

    F = -c*1/3/(kappa*rho)*grad_re
    return F

def GetMdot(r,rho,v):
    return 4*np.pi*r**2*rho*v


def OptDepth(r,rho,v,kappa):
    #Outer boundary condition
    Mdot = GetMdot(r,rho,v)
    vinf,beta = FitBetaLaw(r,v)
    Outer_tau = kappa[-1]*np.mean(Mdot)/(4*np.pi*vinf)*(1-R*Rsun/r[-1])**(1-b)/((1-b)*r[-1])

    #Integrate inwards:
    tau_arr = [Outer_tau]
    for i in range(len(r)-1,-1,-1):
        # print(i,tau_arr[-1])
        tau_arr.append(tau_arr[-1] + kappa[i]*rho[i]*(r[i]-r[i-1]))

    return np.flip(np.array(tau_arr[:-1]))


def beta_law(r, vinf, b):
	return vinf * (1.- r[0]/r)**b

def FitBetaLaw(r,v):
    popt, pcov = curve_fit(beta_law, r, v, bounds=(0,[100*np.max(v),3]))
    vinf = popt[0]
    beta = popt[1]
    return vinf,beta

def GetSoundSpeed2(re):
    T = abs(re)**(1/4)/arad
    a2 = kb*T/(mp*mu)
    return a2

def GetHeff(r,Gamma,a2):
    g_eff = G*M*Msun/r**2*abs(Gamma-1)
    H_eff = a2/g_eff
    return H_eff

def GetL_tot(r,rho,v,re,kappa):
    #L = 4pr2F_obs
    #F_obs = F_cmf + 4/3 v Er
    F_cmf = GetRadflux(r,rho,re,kappa)
    F_obs = F_cmf + 4/3*v*re
    L = 4*np.pi*r**2*F_obs
    return L

def GetL_diff(r,rho,re,kappa):
    #L = 4pr2F_obs
    #F_obs = F_cmf + 4/3 v Er
    F_cmf = GetRadflux(r,rho,re,kappa)
    L = 4*np.pi*r**2*F_cmf
    return L


def GetL_adv(r,v,re):
    #L = 4pr2F_obs
    #F_obs = F_cmf + 4/3 v Er
    F_adv = 4/3*v*re
    L = 4*np.pi*r**2*F_adv
    return L

def GetSEL(r,rho,v):
    #L = L0 - Mdot (v2/2 - GM/r + GM/R)
    Mdot = GetMdot(r,rho,v)
    LSE = L*Lsun - Mdot*(v**2/2 - G*M*Msun/r +G*M*Msun/(R*Rsun))
    return LSE

def Getm(r,rho,v,L):
    Mdot = GetMdot(r,rho,v)
    m_arr = Mdot*G*M*Msun/(r*L)
    return m_arr

# amrvac_outfile = '../output/G2m020094.dat'
amrvac_outfile = '../output/G2m02_L0080.dat'
SE_infile = 'model_G2_m0.2'

r_A,rho_A = AMRVAC_profile(amrvac_outfile,'rho')
r_A,v_A = AMRVAC_profile(amrvac_outfile,'v')
r_A,p_A = AMRVAC_profile(amrvac_outfile,'p')
r_A,e_A = AMRVAC_profile(amrvac_outfile,'e')
r_A,re_A = AMRVAC_profile(amrvac_outfile,'re')
r_A, kappa_A = AMRVAC_profile(amrvac_outfile,'kappa')
Gamma_A = GetGamma(r_A,rho_A,re_A)
v,b = FitBetaLaw(r_A,v_A)
v_fit_A = beta_law(r_A,v,b)
tau_A = OptDepth(r_A,rho_A,v_A,kappa_A)
a2_A = GetSoundSpeed2(re_A)
Heff_A = GetHeff(r_A,Gamma_A,a2_A)
Ltot_A = GetL_tot(r_A,rho_A,v_A,re_A,kappa_A)
Ldiff_A = GetL_diff(r_A,rho_A,re_A,kappa_A)
Ladv_A = GetL_adv(r_A,v_A,re_A)
F_A = GetRadflux(r_A,rho_A,re_A,kappa_A)
m_A = Getm(r_A,rho_A,v_A,Ltot_A)


r_As,rho_As = AMRVAC_single_profile(amrvac_outfile,'rho')
r_As,v_As = AMRVAC_single_profile(amrvac_outfile,'v')
r_As,p_As = AMRVAC_single_profile(amrvac_outfile,'p')
r_As,e_As = AMRVAC_single_profile(amrvac_outfile,'e')
r_As,re_As = AMRVAC_single_profile(amrvac_outfile,'re')
r_As, kappa_As = AMRVAC_single_profile(amrvac_outfile,'kappa')
Gamma_As = GetGamma(r_As,rho_As,re_As)
v,b = FitBetaLaw(r_As,v_As)
v_fit_As = beta_law(r_As,v,b)
tau_As = OptDepth(r_As,rho_As,v_As,kappa_As)
a2_As = GetSoundSpeed2(re_As)
Heff_As = GetHeff(r_As,Gamma_As,a2_As)
Ltot_As = GetL_tot(r_As,rho_As,v_As,re_As,kappa_As)
Ldiff_As = GetL_diff(r_As,rho_As,re_As,kappa_As)
Ladv_As = GetL_adv(r_As,v_As,re_As)
F_As = GetRadflux(r_As,rho_As,re_As,kappa_As)
m_As = Getm(r_As,rho_As,v_As,Ltot_As)

print('amrvac', v,b)

r_S,rho_S = SE_UNIFIED_profile(SE_infile,'rho')
r_S,v_S = SE_UNIFIED_profile(SE_infile,'v')
r_S,p_S = SE_UNIFIED_profile(SE_infile,'p')
r_S,e_S = SE_UNIFIED_profile(SE_infile,'e')
r_S,re_S = SE_UNIFIED_profile(SE_infile,'re')
r_S,L_tot,L_adv,L_diff = Get_SE_Luminosities(SE_infile)
Gamma_S = GetGamma(r_S,rho_S,re_S)
v,b = FitBetaLaw(r_S,v_S)
v_fit_S = beta_law(r_S,v,b)
a2_S = GetSoundSpeed2(re_S)
Heff_S = GetHeff(r_S,Gamma_S,a2_S)
L_S = GetSEL(r_S,rho_S,v_S)
m_S = Getm(r_S,rho_S,v_S,L_S)
Mdot_S = GetMdot(r_S,rho_S,v_S)

print('SE', v,b)

fig,(ax1,ax2,ax3,ax4) = plt.subplots(4,1,sharex=True)
ax1.grid()
ax2.grid()
ax3.grid()
ax4.grid()

ax1.loglog(r_A/(R*Rsun),rho_A,'r',label='AMRVAC')
ax1.loglog(r_As/(R*Rsun),rho_As,'r.',label='final')
ax1.loglog(r_S/(R*Rsun),rho_S,'b',label='SE_Unified')
ax2.semilogx(r_A/(R*Rsun),v_A/1e5,'r')
ax2.semilogx(r_As/(R*Rsun),v_As/1e5,'r.')
ax2.semilogx(r_S/(R*Rsun),v_S/1e5,'b')
ax3.semilogx(r_A/(R*Rsun),e_A,'r')
ax3.semilogx(r_As/(R*Rsun),e_As,'r.')
ax3.semilogx(r_S/(R*Rsun),e_S,'b')
ax4.loglog(r_A/(R*Rsun),re_A,'r')
ax4.loglog(r_As/(R*Rsun),re_As,'r.')
ax4.loglog(r_S/(R*Rsun),re_S,'b')

ax1.legend()
ax1.set_title('Conserved quantities')
ax1.set_ylabel('rho [g/cm3]')
ax2.set_ylabel('v [km/s]')
ax3.set_ylabel('e [erg/cm3]')
ax4.set_ylabel('Erad [erg/cm3]')
ax4.set_xlabel('r [R*]')
ax1.set_xlim([1,11])

plt.figure(2)

plt.title('Velocity profile, beta law')
plt.plot(r_A/(R*Rsun),v_A,'r-',label='AMRVAC')
plt.plot(r_S/(R*Rsun),v_S,'b-',label='SE_unified')

plt.plot(r_A/(R*Rsun),np.sqrt(a2_A)/1e5,'r--',label='soundspeed')
plt.plot(r_S/(R*Rsun),np.sqrt(a2_S)/1e5,'b--')

plt.plot(r_A/(R*Rsun),v_fit_A,'r.',label='betalaw')
plt.plot(r_S/(R*Rsun),v_fit_S,'b.')

plt.xlabel('r [R*]')
plt.ylabel('v [km/s]')

plt.xlim(min(r_A/(R*Rsun)),max(r_A/(R*Rsun)))
plt.legend()
# plt.ylim(min(v_A),max(v_A))

plt.figure(3)

plt.title('Eddington factor')
plt.plot(r_A/(R*Rsun),Gamma_A,'r-',label='AMRVAC')
plt.plot(r_As/(R*Rsun),Gamma_As,'r.',label='final')
plt.plot(r_S/(R*Rsun),Gamma_S,'b-',label='SE_unified')

plt.xlabel('r [R*]')
plt.ylabel('Gamma')

plt.xlim(min(r_A/(R*Rsun)),max(r_A/(R*Rsun)))
plt.ylim(0,3)
plt.legend()


plt.figure(4)

plt.title('Optical depth')
plt.semilogy(r_A/(R*Rsun),tau_A,'r-',label='AMRVAC')
plt.semilogy(r_As/(R*Rsun),tau_As,'r.',label='final')
# plt.loglog(r_S/(R*Rsun),tau_S,'b-',label='SE_unified')

plt.xlabel('r [R*]')
plt.ylabel('tau')

plt.xlim(min(r_A/(R*Rsun)),max(r_A/(R*Rsun)))
plt.legend()

plt.figure(5)
plt.title('speed of sound')
plt.plot(r_A/(R*Rsun),np.sqrt(a2_A)/1e5,'r-',label='AMRVAC')
plt.plot(r_As/(R*Rsun),np.sqrt(a2_As)/1e5,'r.',label='final')
plt.plot(r_S/(R*Rsun),np.sqrt(a2_S)/1e5,'b-',label='SE_unified')

plt.xlabel('r [R*]')
plt.ylabel('a [km/s]')

plt.xlim(min(r_A/(R*Rsun)),max(r_A/(R*Rsun)))
# plt.ylim(0,3)
plt.legend()


plt.figure(6)
plt.title('Effective scaleheight')
plt.semilogy(r_A/(R*Rsun),Heff_A/(R*Rsun)**2,'r-',label='AMRVAC')
plt.semilogy(r_As/(R*Rsun),Heff_As/(R*Rsun)**2,'r.',label='final')
plt.semilogy(r_S/(R*Rsun),Heff_S/(R*Rsun)**2,'b-',label='SE_unified')
dx = (r_A[2]-r_A[1])/(R*Rsun)
plt.semilogy(r_A/(R*Rsun),dx*np.ones(len(r_A)),'r--',label='dx in amrvac')


plt.xlabel('r [R*]')
plt.ylabel('H_eff [R*]')

plt.xlim(min(r_A/(R*Rsun)),max(r_A/(R*Rsun)))
# plt.ylim(0,3)
plt.legend()

plt.figure(7)
plt.title('Luminosity')
plt.plot(r_A/(R*Rsun),Ltot_A/Lsun,'r-',label='AMRVAC')
plt.plot(r_A/(R*Rsun),Ladv_A/Lsun,'r.',label='AMRVAC')
plt.plot(r_A/(R*Rsun),Ldiff_A/Lsun,'r--',label='AMRVAC')
# plt.plot(r_As/(R*Rsun),Ltot_As/Lsun,'r.',label='AMRVAC')
# plt.plot(r_S/(R*Rsun),L_S/Lsun,'b-',label='SE_unified')
plt.plot(r_S/(R*Rsun),L_tot/Lsun,'k-',label='SE_L_tot')
plt.plot(r_S/(R*Rsun),L_adv/Lsun,'k.',label='SE_L_adv')
plt.plot(r_S/(R*Rsun),L_diff/Lsun,'k--',label='SE_L_diff')


plt.xlabel('r [R*]')
plt.ylabel('L [L*]')

plt.xlim(min(r_A/(R*Rsun)),max(r_A/(R*Rsun)))
# plt.ylim(0,3)
plt.legend()


plt.figure(8)
plt.title('Flux')
plt.plot(r_A/(R*Rsun),F_A)

plt.figure(9)
plt.title('m')
plt.plot(r_A/(R*Rsun),m_A,'r-',label='AMRVAC')
plt.plot(r_As/(R*Rsun),m_As,'r.',label='final')
plt.plot(r_S/(R*Rsun),m_S,'b-',label='SE_unified')
plt.plot(r_S/(R*Rsun),Mdot_S*G*M*Msun/(R*Rsun*L*Lsun),'b--',label='Input value')
plt.xlim(min(r_A/(R*Rsun)),max(r_A/(R*Rsun)))


plt.figure(10)
plt.title('Mdot')
plt.plot(r_A/(R*Rsun),GetMdot(r_A,rho_A,v_A),'r-',label='AMRVAC')
plt.plot(r_As/(R*Rsun),GetMdot(r_As,rho_As,v_As),'r.',label='final')
plt.plot(r_S/(R*Rsun),GetMdot(r_S,rho_S,v_S),'b-',label='SE_unified')
plt.plot(r_S/(R*Rsun),Mdot_S,'b--',label='Input value')
# plt.xlim(min(r_A/(R*Rsun)),max(r_A/(R*Rsun)))


plt.legend()

plt.show()
