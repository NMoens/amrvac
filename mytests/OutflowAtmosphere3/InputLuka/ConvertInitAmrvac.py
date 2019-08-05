import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import interpolate

rsun = 6.96e10
msun = 1.99e33

strc = os.getcwd() + '/struct.txt'
i, r, v, rho, er = np.loadtxt(strc,skiprows = 2, unpack=True)

Mdot = 4*np.pi*r**2*rho*v

n_cells = 256
n_gh = 2
R_min = 1.*rsun
R_max = 7.*rsun

inew = range(n_cells+2*n_gh)

dx = (R_max - R_min)/(n_cells)
my_amrvac_grid = np.linspace(R_min - n_gh*dx,R_max + n_gh*dx,n_cells+2*n_gh)

spl_v = interpolate.splrep(r, v, s=0)
vnew = interpolate.splev(my_amrvac_grid, spl_v, der=0)

spl_rho = interpolate.splrep(r, rho, s=0)
rhonew = interpolate.splev(my_amrvac_grid, spl_rho, der=0)

dinflo = 4.26875892E-09
for i in range(n_gh):
    rhonew[i] = max(dinflo,rhonew[i])

vnew[0:2] = v[0]
rhonew[0:2] = Mdot[2]/(4*np.pi*my_amrvac_grid[0:2]**2*vnew[0:2])

spl_er = interpolate.splrep(r, er, s=0)
ernew = interpolate.splev(my_amrvac_grid, spl_er, der=0)


np.savetxt('structure_amrvac.txt',np.transpose([my_amrvac_grid,vnew,rhonew,ernew]))

plt.figure()
plt.plot(r/R_min,v,'r')
plt.plot(my_amrvac_grid/R_min,vnew,'bx')
plt.xlim([0,10])

plt.figure()
plt.plot(r/R_min,er,'r')
plt.plot(my_amrvac_grid/R_min,ernew,'bx')
plt.xlim([0,10])

plt.figure()
plt.semilogy(r/R_min,rho,'r')
plt.semilogy(my_amrvac_grid/R_min,rhonew,'bx')
plt.xlim([0,10])

plt.show()
