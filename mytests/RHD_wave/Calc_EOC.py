from scipy.interpolate import interp1d
import numpy as np
import os
import matplotlib.pyplot as plt

c = 2.99792458e10
year = 3.1536e7
arad = 7.5657e-15
kb = 1.380648520e-16
mp = 1.6737236e-24
mu = 0.60869565217391308
hd_gamma = 1.6667

wvl = 777363184079.60193

def analytic(ii,jj,cheat,x):
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

def get_data(myfile):
    A = np.loadtxt(myfile,skiprows=3)

    x0 = np.transpose(A)[0]
    rho0 = np.transpose(A)[1]
    v0 = np.transpose(A)[2]
    e0 = np.transpose(A)[3]
    E0 = np.transpose(A)[4]

    # delta_rho = 100*(rho0-np.mean(rho0))
    # delta_v = v0
    # delta_e = e0-np.mean(e0)
    # delta_E = E0-np.mean(E0)

    return x0, rho0, v0, e0, E0

def get_error(x1,y1,x2,y2):
    f1 = interp1d(x1,y1)
    err = y2 - f1(x2)
    return abs(err)

def Calculate_EOC(Imex,limiter,it,norm,plot):

    # Imex = 'SP'
    # limiter = 'weno5'
    # it = '0100'
    ampl = 1.e-2

    x_100 , rho_100, v_100, e_100, E_100 = get_data('./convergence/' + Imex + '_' + limiter + '_100' + it + '.blk')
    x_200 , rho_200, v_200, e_200, E_200 = get_data('./convergence/' + Imex + '_' + limiter + '_200' + it + '.blk')
    x_400 , rho_400, v_400, e_400, E_400 = get_data('./convergence/' + Imex + '_' + limiter + '_400' + it + '.blk')
    x_800 , rho_800, v_800, e_800, E_800 = get_data('./convergence/' + Imex + '_' + limiter + '_800' + it + '.blk')
    x_1600 , rho_1600, v_1600, e_1600, E_1600 = get_data('./convergence/' + Imex + '_' + limiter + '_1600' + it + '.blk')
    x_3200 , rho_3200, v_3200, e_3200, E_3200 = get_data('./convergence/' + Imex + '_' + limiter + '_3200' + it + '.blk')

    # rho_100 = v_100
    # rho_200 = v_200
    # rho_400 = v_400
    # rho_800 = v_800
    # rho_1600 = v_1600
    # rho_3200 = v_3200

    x_max = x_3200
    rho_max = rho_3200

    analytic_100 = 1 + ampl*analytic(0,0,False,x_100*wvl)
    analytic_200 = 1 + ampl*analytic(0,0,False,x_200*wvl)
    analytic_400 = 1 + ampl*analytic(0,0,False,x_400*wvl)
    analytic_800 = 1 + ampl*analytic(0,0,False,x_800*wvl)
    analytic_1600 = 1 + ampl*analytic(0,0,False,x_1600*wvl)
    analytic_3200 = 1 + ampl*analytic(0,0,False,x_3200*wvl)

    # err_200 = max(abs(analytic_200 - rho_200))
    # err_400 = max(abs(analytic_400 - rho_400))
    # err_800 = max(abs(analytic_800 - rho_800))
    err_100 = abs(analytic_100 - rho_100)
    err_200 = abs(analytic_200 - rho_200)
    err_400 = abs(analytic_400 - rho_400)
    err_800 = abs(analytic_800 - rho_800)
    err_1600 = abs(analytic_1600 - rho_1600)
    err_3200 = abs(analytic_3200 - rho_3200)

    err2_100 = get_error(x_max,rho_max,x_100,rho_100)
    err2_200 = get_error(x_max,rho_max,x_200,rho_200)
    err2_400 = get_error(x_max,rho_max,x_400,rho_400)
    err2_800 = get_error(x_max,rho_max,x_800,rho_800)
    err2_1600 = get_error(x_max,rho_max,x_1600,rho_1600)

    if plot:
        fig,(ax1,ax2,ax3) = plt.subplots(3,1,sharex=True)

        ax1.set_title(Imex+' & '+limiter)

        ax1.plot(x_100,rho_100,'x',label='$N_c = 100$')
        ax1.plot(x_200,rho_200,'x',label='$N_c = 200$')
        ax1.plot(x_400,rho_400,'x',label='$N_c = 400$')
        ax1.plot(x_800,rho_800,'x',label='$N_c = 800$')
        ax1.plot(x_1600,rho_1600,'x',label='$N_c = 1600$')
        ax1.plot(x_3200,rho_3200,'x',label='$N_c = 3200$')
        ax1.plot(x_3200,analytic_3200,'k-',label='analytic')
        ax1.set_ylabel('$\\rho$')
        ax1.legend(loc='right')

        ax2.plot(x_100,err_100,'x')
        ax2.plot(x_200,err_200,'x')
        ax2.plot(x_400,err_400,'x')
        ax2.plot(x_800,err_800,'x')
        ax2.plot(x_1600,err_1600,'x')
        ax2.plot(x_3200,err_3200,'x')
        ax2.set_ylabel('error wrt analytic')

        ax3.plot(x_100,err2_100,'x')
        ax3.plot(x_200,err2_200,'x')
        ax3.plot(x_400,err2_400,'x')
        ax3.plot(x_800,err2_800,'x')
        ax3.plot(x_1600,err2_1600,'x')
        ax3.set_ylabel('error wrt $N_max$')
        ax3.set_xlabel('$x/\\lambda$')


    # EOC_200 = -np.log(np.linalg.norm(err_200,norm)/np.linalg.norm(err_100,norm))/np.log(200/100)
    # EOC_400 = -np.log(np.linalg.norm(err_400,norm)/np.linalg.norm(err_200,norm))/np.log(400/200)
    # EOC_800 = -np.log(np.linalg.norm(err_800,norm)/np.linalg.norm(err_400,norm))/np.log(800/400)
    # EOC_1600 = -np.log(np.linalg.norm(err_1600,norm)/np.linalg.norm(err_800,norm))/np.log(1600/800)

    EOC2_200 = np.log(np.linalg.norm(err2_200,norm)/np.linalg.norm(err2_100,norm))/np.log(0.5)
    EOC2_400 = np.log(np.linalg.norm(err2_400,norm)/np.linalg.norm(err2_200,norm))/np.log(0.5)
    EOC2_800 = np.log(np.linalg.norm(err2_800,norm)/np.linalg.norm(err2_400,norm))/np.log(0.5)
    EOC2_1600 = np.log(np.linalg.norm(err2_1600,norm)/np.linalg.norm(err2_800,norm))/np.log(0.5)

    print(Imex+' & '+limiter)
    print('EOC_200', round(EOC2_200,2))
    print('EOC_400', round(EOC2_400,2))
    print('EOC_800', round(EOC2_800,2))
    print('EOC_1600', round(EOC2_1600,2))

    return [EOC2_200, EOC2_400, EOC2_800, EOC2_1600]

res = [200, 400, 800, 1600]

nrm = 2
it = '0100'

sp_minmod = Calculate_EOC('SP','minmod',it,nrm,True)
# sp_weno5 = Calculate_EOC('SP','weno5',it,nrm,True)
Midpoint_mp5 = Calculate_EOC('Midpoint','mp5',it,nrm,True)
Midpoint_weno5 = Calculate_EOC('Midpoint','weno5',it,nrm,True)
ARS3_weno5 = Calculate_EOC('ARS3','weno5',it,nrm,True)

plt.figure()
plt.title('RHD wave EOC at it ' + it)
plt.loglog(res,sp_minmod,'b^',label='SP & minmod')
plt.loglog(res,sp_minmod,'b-')
# plt.plot(res,sp_weno5,'bs',label='SP & Weno5')
# plt.plot(res,sp_weno5,'b')
plt.loglog(res,Midpoint_mp5,'rs',label='Midpoint & mp5')
plt.loglog(res,Midpoint_mp5,'r-')
plt.loglog(res,Midpoint_weno5,'r^',label='Midpoint & Weno5')
plt.loglog(res,Midpoint_weno5,'r-')
plt.loglog(res,ARS3_weno5,'g^',label='ARS3 & Weno5')
plt.loglog(res,ARS3_weno5,'g-')
plt.xlabel('$N_{cells}$')
plt.xticks(res,res)
plt.ylabel('$EOC_{N}$')
plt.hlines(1,200,2000,linestyles='--',colors='blue')
plt.hlines(2,200,2000,linestyles='--',colors='red')
plt.hlines(3,200,2000,linestyles='--',colors='green')
plt.legend()

imex = 'Midpoint'
limiter = 'mp5'
nrm = 2#np.inf
# its = ['0010','0020','0030','0040','0050']
# its = ['0050','0060','0070','0080','0090','0100']
its = ['0100']
cls = 1- np.linspace(0.1,0.9,len(its))

plt.figure()
for i in range(len(its)):
    EOC_arr = Calculate_EOC('SP','minmod',its[i],nrm,False)
    plt.plot(res,EOC_arr,'bs')
    # plt.plot(res,EOC_arr,'-',c=str(cls[i]),label= its[i])
    plt.plot(res,EOC_arr,'-',c='b',label= its[i])

    EOC_arr = Calculate_EOC('Midpoint','weno5',its[i],nrm,False)
    plt.plot(res,EOC_arr,'rs')
    # plt.plot(res,EOC_arr,'-',c=str(cls[i]))
    plt.plot(res,EOC_arr,'-',c='r')

    EOC_arr = Calculate_EOC('ARS3','weno5',its[i],nrm,False)
    plt.plot(res,EOC_arr,'gs')
    # plt.plot(res,EOC_arr,'-',c=str(cls[i]))
    plt.plot(res,EOC_arr,'-',c='g')


plt.title('RHD wave EOC using '+str(nrm)+'-norm')
plt.hlines(1,200,2000,linestyles='--',colors='blue',label='SP')
plt.hlines(2,200,2000,linestyles='--',colors='red',label='Midoint')
plt.hlines(3,200,2000,linestyles='--',colors='green',label='ARS3')
plt.xlabel('$N_{cells}$')
plt.xticks(res,res)
plt.ylabel('$EOC_{N}$')

plt.legend()



plt.show()
