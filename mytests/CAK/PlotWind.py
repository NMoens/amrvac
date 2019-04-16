import numpy as np
import matplotlib.pyplot as plt
import os


def plot_profile(sn):
    if sn < 10:
        n = '000' + str(sn)
    elif sn < 100:
        n = '00' + str(sn)
    elif sn < 1000:
        n = '0' + str(sn)

    f = os.getcwd() + '/output/CAK' + n + '.blk'

    r,rho,v, g_cak, g_eff, f_fd = np.loadtxt(f,skiprows=3, unpack=True)

    ax1.semilogy(r,rho,'.',label=n)
    ax2.plot(r,f_fd,'.',label=n)

    return n

f,(ax1,ax2) = plt.subplots(2,1,sharex=True)

n = plot_profile(0)
n = plot_profile(18)
n = plot_profile(9)

ax1.grid(which='both')
ax2.grid(which='both')

ax1.set_title('snapshot ' + n, fontsize=15)
ax1.set_ylabel('$\\frac{\\rho(r)}{\\rho_{bound}}$', fontsize=15)
ax2.set_ylabel('$\\frac{v(r)}{v_{sound}}$', fontsize=15)
ax2.set_xlabel('$\\frac{r}{R_*}$', fontsize=15)

ax2.legend(loc='best')

plt.show()
