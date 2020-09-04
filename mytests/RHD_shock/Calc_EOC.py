from scipy.interpolate import interp1d
import numpy as np
import os
import matplotlib.pyplot as plt

def get_data(myfile):
    A = np.loadtxt(myfile,skiprows=3)

    x0 = np.transpose(A)[0]
    rho0 = np.transpose(A)[1]
    v0 = np.transpose(A)[2]
    e0 = np.transpose(A)[3]
    E0 = np.transpose(A)[4]

    return x0, rho0, v0, e0, E0

def get_error(x1,y1,x2,y2):
    f1 = interp1d(x1,y1)
    err = y2 - f1(x2)
    return abs(err)

def Calculate_EOC(Imex,limiter,it,norm,plot,N_c):

    # Imex = 'SP'
    # limiter = 'weno5'
    # it = '0100'
    ampl = 1.e-2

    x_init, rho_init, v_init, e_init, E_init = get_data('./convergence/' + Imex + '_' + limiter + '_' + N_c[-1]  + '0000' + '.blk')

    x_ = []
    rho_ = []

    for n in N_c:
        x , rho, v, e, E = get_data('./convergence/' + Imex + '_' + limiter + '_' + n + it + '.blk')
        x_.append(x)
        rho_.append(rho)

    x_max = x_[-1]
    rho_max = rho_[-1]

    err_ = []

    for i in range(len(N_c)):
        err_.append(get_error(x_max,rho_max,x_[i],rho_[i]))

    if plot:
        fig,(ax1,ax2) = plt.subplots(2,1,sharex=True)

        ax1.set_title(Imex+' & '+limiter)

        for i in range(len(N_c)):
            ax1.plot(x_[i],rho_[i],'x',label='$N_c =$' + N_c[i])

        ax1.plot(x_init,rho_init,'k-',label='initial')
        ax1.set_ylabel('$\\rho$')
        ax1.legend()

        for i in range(len(N_c)-1):
            ax2.plot(x_[i],err_[i],'x')

        ax2.set_ylabel('error wrt $N_max$')
        ax2.set_xlabel('$x/\\lambda$')

    EOC_ = [0]

    for i in range(1,len(N_c)-1):
        EOC_.append(np.log(np.linalg.norm(err_[i],norm)/np.linalg.norm(err_[i-1],norm))/np.log(0.5))

    print(Imex+' & '+limiter)
    for i in range(1,len(N_c)-1):
        print('EOC ',N_c[i], ': ' , round(EOC_[i],2))

    return EOC_[1:]

# N_c = ['50', '100', '200', '400', '800', '1600', '3200','6400']
# res = [100, 200, 400, 800, 1600, 3200]
#
# nrm = np.inf
# it = '0014'
#
#
# # Euler_weno5 = Calculate_EOC('Euler','weno5',it,nrm,True,N_c)
# # SP_weno5 = Calculate_EOC('SP','weno5',it,nrm,True,N_c)
# Midpoint_mp5 = Calculate_EOC('Midpoint','mp5',it,nrm,True,N_c)
# # ARS3_mp5 = Calculate_EOC('ARS3','mp5',it,nrm,True,N_c)
#
# plt.figure()
# plt.title('RHD Shock EOC at it ' + it)
# # plt.loglog(res,Euler_weno5,'b^',label='Euler & weno5')
# # plt.loglog(res,Euler_weno5,'b-')
# # plt.loglog(res,SP_weno5,'bs',label='SP & weno5')
# # plt.loglog(res,SP_weno5,'b-')
# plt.loglog(res,Midpoint_mp5,'r^',label='Midpoint & mp5')
# plt.loglog(res,Midpoint_mp5,'r-')
# # plt.loglog(res,ARS3_mp5,'g^',label='ARS3 & mp5')
# # plt.loglog(res,ARS3_mp5,'g-')
# plt.xlabel('$N_{cells}$')
# plt.xticks(res,res)
# plt.ylabel('$EOC_{N}$')
# plt.hlines(1,100,2000,linestyles='--',colors='blue')
# plt.hlines(2,100,2000,linestyles='--',colors='red')
# plt.hlines(3,100,2000,linestyles='--',colors='green')
# plt.legend()



# limiter = 'weno5'
nrm = np.inf
# its = ['0021','0022','0023','0024','0025','0026','0027','0028','0029','0030','0031','0032','0033','0034','0035','0036','0037','0038','0039','0040']
its = ['0001','0002','0003','0004','0005','0006','0007','0008','0009','0010','0011','0012','0013','0014','0015','0016']

cls = 1- np.linspace(0.01,0.99,len(its))

N_c = ['50', '100', '200', '400', '800', '1600', '3200', '6400']
res = [100, 200, 400, 800, 1600, 3200]

# N_c = ['50', '100', '200', '400', '800', '1600', '3200']
# res = [100, 200, 400, 800, 1600]

plt.figure()
for i in range(len(its)):
    # EOC_arr = Calculate_EOC('Euler','weno5',its[i],nrm,False, N_c)
    # plt.semilogx(res,EOC_arr,'-',c=str(cls[i]),label= its[i])
    # plt.semilogx(res,EOC_arr,'bs')

    # EOC_arr = Calculate_EOC('SP','weno5',its[i],nrm,False, N_c)
    # plt.semilogx(res,EOC_arr,'-',c=str(cls[i]))
    # plt.semilogx(res,EOC_arr,'b^')

    EOC_arr = Calculate_EOC('Midpoint','mp5',its[i],nrm,False, N_c)
    plt.semilogx(res,EOC_arr,'-',c=str(cls[i]))
    plt.semilogx(res,EOC_arr,'rs')

    # EOC_arr = Calculate_EOC('ARS3','mp5',its[i],nrm,False, N_c)
    # plt.semilogx(res,EOC_arr,'-',c=str(cls[i]))
    # plt.semilogx(res,EOC_arr,'gs')

plt.title('RHD Shock EOC using '+str(nrm)+'-norm')
plt.hlines(1,100,4000,linestyles='--',colors='blue',label='SP')
plt.hlines(2,100,4000,linestyles='--',colors='red',label='Midoint')
plt.hlines(3,100,4000,linestyles='--',colors='green',label='ARS3')
plt.xlabel('$N_{cells}$')
plt.xticks(res,res)
plt.ylabel('$EOC_{N}$')
plt.legend()





plt.show()
