import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

def Get_data(still,adv,ind,Deltax):

    # u_t = 4.7999999999999996e-7
    # t = 10
    # v = 5.e7

    A = np.loadtxt(still,skiprows=3)
    B = np.loadtxt(adv,skiprows=3)

    x = np.transpose(A)[0]

    still = np.transpose(A)[ind]
    adv = np.transpose(B)[ind]

    if ind == 4:
        still = np.transpose(A)[4] + 

    spl_still = interpolate.splrep(x, still, s=0)
    spl_adv = interpolate.splrep(x, adv, s=0)

    x_adv = x - Deltax

    x_new = np.linspace(-10,10,500)

    still_new = interpolate.splev(x_new, spl_still, der=0)
    adv_new = interpolate.splev(x_new + Deltax, spl_adv, der=0)

    reldiff = (still_new - adv_new)/still_new

    return x_new, x, x_adv,  still, adv, reldiff

ind = 4
dx = 10.08

still_file = '1D_output/still_Euler0010.blk'
adv_file = '1D_output/adv_Euler0010.blk'
x_new, x_still, x_adv, still, adv, delta_Euler = Get_data(still_file, adv_file, ind, dx)

# still_file = '1D_output/still_SP0010.blk'
# adv_file = '1D_output/adv_SP0010.blk'
# x_new, x_still, x_adv, still, adv, delta_SP = Get_data(still_file, adv_file, ind, dx)
#
# still_file = '1D_output/still_Midpoint0010.blk'
# adv_file = '1D_output/adv_Midpoint0010.blk'
# x_new, x_still, x_adv, still, adv, delta_Midpoint = Get_data(still_file, adv_file, ind, dx)
#
# still_file = '1D_output/still_ARS30010.blk'
# adv_file = '1D_output/adv_ARS30010.blk'
# x_new, x_still, x_adv, still, adv, delta_ARS3 = Get_data(still_file, adv_file, ind, dx)


# #FIND OPTIMAL dx
# dx_arr = []
# err_arr = []
# for dx in np.linspace(10,10.1,500):
#     print(dx)
#     x_new, x_still, x_adv, still, adv, delta_ARS3 = Get_data(still_file, adv_file, 1, dx)
#
#     dx_arr.append(dx)
#     err_arr.append(max(abs(delta_ARS3)))
#
# plt.plot(dx_arr,err_arr)

plt.figure()
plt.plot(x_still,still)
plt.plot(x_still,adv)

# plt.figure()
# plt.plot(x_still,still)
# plt.plot(x_adv,adv)

plt.figure()
plt.title('advection pulse: relative difference')
plt.plot(x_new,delta_Euler,label='Euler')
# plt.plot(x_new,delta_SP,label='SP')
# plt.plot(x_new,delta_Midpoint,label='Midpoint')
# plt.plot(x_new,delta_ARS3,label='ARS3')

plt.ylabel('$(\\rho_s - \\rho_a)/\\rho_s$',fontsize = 10)
plt.xlabel('$x$',fontsize = 10)

plt.legend()

plt.show()
