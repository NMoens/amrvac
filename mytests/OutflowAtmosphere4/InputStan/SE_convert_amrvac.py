import numpy as np
import matplotlib.pyplot as plt

import os
from scipy import interpolate

Rsun = 6.96e10
Msun = 1.99e33
Lsun = 3.9e33
G = 6.67259e-8
year = 356.25*24*60*60
arad = 7.5646e-15
kb = 1.380658e-16
mp = 1.6726231e-24
mu = 0.6
clight = 2.99792458e10
hd_gamma = 5.0/3.0

m = 0.2
Gamma = 2
error_b = 10.0
# error_b = 5 #this is used in the grids


R = 1*Rsun
M = 10*Msun
L = 1.9e5*Lsun

n_cells = 4096
n_gh = 2
R_min = 1.*R
R_max = 21.*R

Mdot = m*R*L/(G*M) #8.3e-6*Msun/year

v_esc = (2*G*M/R)**0.5

def convert_all(wg,p,x):
    r = R/(1-x)
    w = wg*(Gamma-1)
    v = np.sqrt(w)*v_esc
    Prad = p*L/(4*np.pi*R**2*v_esc)
    rho = Mdot/(4*np.pi*r**2*v)
    Erad = 3*Prad
    T = (Erad/arad)**(1./4.)
    pg = kb*T/(mp*mu) *rho
    e = pg/(hd_gamma-1)+0.5*v**2*rho
    return r, rho, v, pg, Erad

def GetProfiles():
    fname = 'model_G' + str(Gamma) + '_m' + str(m)
    # fname = 'My_model2'
    print('reading ', fname)

    wg,p,x = np.loadtxt(fname,unpack=True)
    r,rho,v,pg,Erad = convert_all(wg,p,x)
    return r,rho,v,pg,Erad

# strc = os.getcwd() + '/struct_steep.txt'
r,rho,v,pg,er = GetProfiles()

# def GetSonicPoint(r,rho,v,er):
#     #Find the sonic point (where v == sqrt(kb T/mp mu))
#     T = (er/arad)**0.25
#     a_s = kb*T/(mp*mu)
#     sonic_fit = v - a_s
#     i_sp = 10 #np.argmin(abs(sonic_fit))
#
#     #Get values at Sonic point
#     r_sp = r[i_sp]
#     rho_sp = rho[i_sp]
#     v_sp = v[i_sp]
#     er_sp = er[i_sp]
#
#     #Integrate downwards into hydrostatic core
#     # dE = Gamma rho G M / (lambda r**2) dx
#     dx = r[i_sp+1] - r[i_sp]
#     nx_constv = int((r_sp-0.5*r[0])/dx)
#
#     print(0.1*r[0],r_sp, nx_constv)
#
#     r_constv = np.linspace(0.5*r[0],r_sp, nx_constv)
#     rho_constv = np.ones(len(r_constv))*rho_sp
#     v_constv = np.ones(len(r_constv))*v_sp
#     er_constv = er_sp - 3*rho_sp*G*M*(1/r_sp - 1/r_constv)
#
#     #concatenate constant velocity region with SE region
#     r = np.concatenate((r_constv[:-1],r[i_sp:]))
#     rho = np.concatenate((rho_constv[:-1],rho[i_sp:]))
#     v = np.concatenate((v_constv[:-1],v[i_sp:]))
#     er = np.concatenate((er_constv[:-1],er[i_sp:]))
#
#     return r, rho, v, er
#
# r, rho, v, er = GetSonicPoint(r,rho,v,er)

# plt.figure()
# plt.plot(r,er)
# plt.show()
# exit()

dx = (R_max - R_min)/(n_cells)
my_amrvac_grid = np.linspace(R_min - (n_gh-0.5)*dx,R_max + (n_gh-0.5)*dx,n_cells+2*n_gh)
# my_amrvac_grid = np.linspace(R_min - (n_gh-0.5)*dx + (n_gh-1)*dx,R_max + (n_gh-0.5)*dx + (n_gh-1)*dx,n_cells+2*n_gh)

spl_v = interpolate.splrep(r, v, s=0)
vnew = interpolate.splev(my_amrvac_grid, spl_v, der=0,  ext=3)

spl_rho = interpolate.splrep(r, rho, s=0)
rhonew = interpolate.splev(my_amrvac_grid, spl_rho, der=0,  ext=3)

# vnew[0:2] = v[0]
# rhonew[0:2] = Mdot/(4*np.pi*my_amrvac_grid[0:2]**2*vnew[0:2])
dinflo = rhonew[1]

spl_er = interpolate.splrep(r, er, s=0)
ernew = interpolate.splev(my_amrvac_grid, spl_er, der=0,  ext=3)

gradE = (er[1:] - er[:-1])/(r[1:] - r[:-1])

spl_gradE = interpolate.splrep(r[:-1],gradE, s=0)
gradEnew = interpolate.splev(my_amrvac_grid,spl_gradE, der=0, ext=3)

################################################################################
################################################################################



################################################################################
################################################################################

L_cmf = L - Mdot*(v**2/2 + M*G/R - M*G/r)
L_cmfnew = L - Mdot*(vnew**2/2 + M*G/R - M*G/my_amrvac_grid)

kappa_b = 4*np.pi*G*M*clight/L
kappa_0 = 4*np.pi*G*M*clight*Gamma/L

kappa_prof = kappa_b + (1+np.erf((r-R)/R*error_b-error_b/2))*(kappa_0-kappa_b)/2
kappa_profnew = kappa_b + (1+np.erf((my_amrvac_grid-R)/R*error_b-error_b/2))*(kappa_0-kappa_b)/2

F = -clight/(3*kappa_prof*rho)*np.gradient(er,r)
Fnew = -clight/(3*kappa_profnew*rhonew)*np.gradient(ernew,my_amrvac_grid)

Gamma_1 = np.ones(len(r))*Gamma
Gamma_2 = kappa_0*L_cmf/(4*np.pi*G*M*clight)
Gamma_3 = kappa_prof*L_cmf/(4*np.pi*G*M*clight)
Gamma_4 = kappa_profnew*L_cmfnew/(4*np.pi*G*M*clight)

Lum_1 = np.ones(len(r))*L
Lum_2 = L_cmf
Lum_3 = 4*np.pi*r**2*F
Lum_4 = 4*np.pi*my_amrvac_grid**2*Fnew
Lum_5 = 4*np.pi*r**2*F

m_1 = np.ones(len(r))*m
m_2 = Mdot*G*M/(R*Lum_1)
m_3 = Mdot*G*M/(r*Lum_2)

np.savetxt('structure_amrvac_G' + str(Gamma) + '_m' +str(m) + '_20R' + '.txt',np.transpose([my_amrvac_grid,vnew,rhonew,ernew]))

params = ['Gamma', 'm', 'M_Msun', 'L_Lsun', 'R_Rsun', 'huh', 'Mdot', 'rho0', 'error_b']
paramvals = [Gamma,m,M/Msun,L/Lsun,R/Rsun,R/Rsun,Mdot*year/Msun,dinflo,error_b]

import csv
with open('params_G' + str(Gamma) + '_m' +str(m) + '_20R'+ '.txt', 'w') as f:
    writer = csv.writer(f, delimiter='\t')
    writer.writerows(zip(params,paramvals))

plt.figure()
plt.title('velocity')
plt.plot(r/R_min,v,'r')
plt.plot(my_amrvac_grid/R_min,vnew,'bx')
plt.xlim([0,10])

plt.figure()
plt.title('radiation energy')
plt.plot(r/R_min,er,'r')
plt.plot(my_amrvac_grid/R_min,ernew,'bx')
plt.xlim([0,10])

plt.figure()
plt.semilogy(r/R_min,rho,'r')
plt.title('density')
plt.semilogy(my_amrvac_grid/R_min,rhonew,'bx')
plt.xlim([0,10])

plt.figure()
plt.title('Flux')
plt.plot(r/R_min,F,'r')
plt.plot(my_amrvac_grid/R_min,Fnew,'bx')
plt.plot(r/R_min,F*(r/R_min)**2,'r--')
plt.plot(my_amrvac_grid/R_min,Fnew*(my_amrvac_grid/R_min)**2,'b--')
plt.xlim([0,10])

plt.figure()
plt.title('opacity')
plt.plot(r/R_min,kappa_b*np.ones(len(r)),'b--')
plt.plot(r/R_min,kappa_0*np.ones(len(r)),'k--')
plt.plot(r/R_min,kappa_prof,'r')
plt.plot(my_amrvac_grid/R_min,kappa_profnew,'rx')

plt.figure()
plt.title('Gamma')
plt.grid()
plt.plot(r/R_min,Gamma_1,'b--',label='base model')
plt.plot(r/R_min,Gamma_2,'k',label='Stans model')
plt.plot(r/R_min,Gamma_3,'r',label='bound')
plt.plot(my_amrvac_grid/R_min,  Gamma_4  ,'rx',label='Amrvac')
plt.xlim([1,10])
plt.ylim([0.8,1.2*Gamma])
plt.xlabel('r/R')
plt.ylabel('Gamma')
plt.legend()

plt.figure()
plt.title('Luminosity')
plt.grid()
plt.plot(r/R_min,Lum_1/L,'b--',label='base model')
plt.plot(r/R_min,Lum_2/L,'k',label='Stans model')
plt.plot(r/R_min,Lum_3/L,'r',label='bound')
plt.plot(my_amrvac_grid/R_min,Lum_4/L,'rx',label='Amrvac')
plt.plot(r/R_min,Lum_5/L,'g',label='bound')
plt.xlim([1,10])
plt.ylim([0.8,1.05])
plt.xlabel('r/R')
plt.ylabel('L/L_0')
plt.legend()

plt.figure()
plt.title('m')
plt.grid()
plt.plot(r/R_min,m_1,'b--',label='base model')
plt.plot(r/R_min,m_2,'r--',label='L drops')
plt.plot(r/R_min,m_3,'k',label='L,r drops')

plt.show()
