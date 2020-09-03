from scipy import interpolate as itp
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
x = np.linspace(0,25*wvl,10000)

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

analytic1 = analytic(0,0,False)

x0 , rho_minmod, v_minmod, e_minmod, E_minmod = get_data('test/fs_hll_lim_minmod0100.blk')
x0 , rho_koren, v_koren, e_koren, E_koren = get_data('test/fs_hll_lim_koren0100.blk')
x0 , rho_mp5, v_mp5, e_mp5, E_mp5 = get_data('test/fs_hll_lim_mp50100.blk')
x0 , rho_weno5, v_weno5, e_weno5, E_weno5 = get_data('test/fs_hll_lim_weno50100.blk')

# x0 , rho_Euler, v_Euler, e_Euler, E_Euler = get_data('1D_output/Euler_thick0100.blk')
# x0 , rho_SP, v_SP, e_SP, E_SP = get_data('1D_output/SP_thick0100.blk')
# x0 , rho_Midpoint, v_Midpoint, e_Midpoint, E_Midpoint = get_data('1D_output/Midpoint_thick0100.blk')
# x0 , rho_ARS3, v_ARS3, e_ARS3, E_ARS3 = get_data('1D_output/ARS3_thick0100.blk')

x0 , rho_Euler, v_Euler, e_Euler, E_Euler = get_data('1D_output/Euler_weno50100.blk')
x0 , rho_SP, v_SP, e_SP, E_SP = get_data('1D_output/SP_weno50100.blk')
x0 , rho_Midpoint, v_Midpoint, e_Midpoint, E_Midpoint = get_data('1D_output/Midpoint_weno50100.blk')
x0 , rho_ARS3, v_ARS3, e_ARS3, E_ARS3 = get_data('1D_output/ARS3_weno50100.blk')

plt.figure()
plt.plot(x/wvl,analytic1,'k-',label='analytic')
plt.plot(x0, rho_minmod, '--',label='hll + minmod')
plt.plot(x0, rho_koren, '--',label='hll + koren')
plt.plot(x0, rho_mp5, '--',label='hll + mp5')
plt.plot(x0, rho_weno5, '--',label='hll + weno5')
plt.ylabel('$\\delta \\rho$')
plt.xlabel('$x/\\lambda$')
plt.legend()

plt.figure()
plt.plot(x/wvl,analytic1,'k-',label='analytic')
plt.plot(x0, rho_Euler, '--',label='Euler')
plt.plot(x0, rho_SP, '--',label='SP')
plt.plot(x0, rho_Midpoint, '--',label='Midpoint')
plt.plot(x0, rho_ARS3, '--',label='ARS3')
plt.ylabel('$\\delta \\rho$')
plt.xlabel('$x/ \\lambda$')
plt.legend()


plt.show()
