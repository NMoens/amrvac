import numpy as np
import matplotlib.pyplot as plt

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

L_damp = 1./np.imag(k)
L_oscillate = 1./np.real(k)

L_G = -L_damp/L_oscillate


print(L_damp/wvl)
print(L_oscillate/wvl)
print(L_G)
print(L_G/(2*np.pi), '---> Cheaty')

x = np.linspace(0,25*wvl,10000)

def Mode(L_d,L_o,x):
    osc = np.sin(x/L_o)
    damp = np.exp((x-wvl)/L_d)
    sol = osc
    for i in range(len(sol)):
        if x[i] > wvl:
            sol[i] = sol[i]*damp[i]
    return sol

def Mode2(k,x):
    osc = np.exp(-1j*np.real(k)*x)
    damp = np.exp((x-wvl)*np.imag(k))
    sol = osc
    for i in range(len(sol)):
        if x[i] > wvl:
            sol[i] = sol[i]*damp[i]
    return sol



plt.figure(1)
plt.plot(x/wvl,Mode(L_damp[0],L_oscillate[0],x),'b-',label='Material dominated')
plt.plot(x/wvl,Mode(L_damp[2],L_oscillate[2],x),'r-',label='Radiation dominated')
# plt.plot(x/wvl,Mode(L_damp[1],L_oscillate[1],x),'b--',label='Acoustic dominated')
# plt.plot(x/wvl,Mode(L_damp[3],L_oscillate[3],x),'k--',label='Radiation dominated')
# plt.plot(x/wvl,Mode(L_damp[2],L_oscillate[2]/(2*np.pi),x),'k.-',label='cheaty Radiation dominated')
plt.plot(x/wvl,Mode(L_damp[0]/(2*np.pi),L_oscillate[0],x),'k--',label='cheaty Acoustic dominated')
plt.plot(x/wvl,Mode(L_damp[2],L_oscillate[0],x),'k-',label='Material wvl, Radiation dampening')
plt.ylabel('perturbation')
plt.xlabel('$x/\\lambda$')
plt.legend()

plt.figure(2)
plt.plot(x/wvl,Mode2(k[0],x),'b-',label='Acoustic dominated')
plt.plot(x/wvl,Mode2(k[2],x),'r-',label='Radiation dominated')
# plt.plot(x/wvl,Mode2(k[1],x),'b--',label='Acoustic dominated')
# plt.plot(x/wvl,Mode2(k[3],x),'k--',label='Radiation dominated')
plt.ylabel('perturbation')
plt.xlabel('$x/\\lambda$')
plt.legend()
plt.show()
