import numpy as np
import os
import matplotlib.pyplot as plt

rho0 = 3.1198e-8
pg0 = 503463.97
a20 = pg0/rho0
er0 = 1499628.401102

pr0 = er0/3.

kb = 1.3806e-16
cc = 2.99e10
mp = 1.672e-24
grav = 6.67e-8
arad = 7.5646e-15

tr0 = (er0/arad)**0.25
mu = 0.6
tg0 = pg0*mu*mp/rho0/kb

rsun = 6.96e10
msun = 1.99e33

rstar = 19.370383317782409*rsun
mstar = 50.*msun
gstar = grav*mstar/rstar**2.
gamma = 0.5
geff = gstar*(1.-gamma)

delta = pr0 - pg0*gamma*(1-gamma)

dy = 138087438.2

scale0 = a20/geff

y_arr = [rstar]
y_new = rstar

while y_new < rstar+2*scale0:
    y_new = y_new+dy
    y_arr.append(y_new)

ny = len(y_arr)


################################################################################

scale = scale0*np.ones((ny))
temp = a20*mu*mp/kb
rho = rho0*np.exp(-(np.array(y_arr)-rstar)/scale)
pg = rho*kb*temp/mu/mp
pr = pg*gamma/(1.-gamma) + delta
er = pr*3.
trad = (er/arad)**0.25

max_err = 100000


while max_err > 1e-10:
    temp = trad
    a2 = kb*temp/mu/mp
    scale = a2/geff
    rho = rho0*np.exp(-(np.array(y_arr)-rstar)/scale)
    pg = rho*kb*temp/mu/mp
    #pg = pg0*np.exp(-(np.array(y_arr)-rstar)/scale)

    # print(pr0,pg0)
    pg0 = rho[0]*kb*temp[0]/mu/mp
    # pr0 = 1/3*(trad[0]*arad)**4

    delta = pr0 - pg0*gamma*(1-gamma)

    pr = pg*gamma/(1.-gamma) + delta
    er = pr*3.
    trad = (er/arad)**0.25

    max_err = max(abs(trad-temp)/temp)
    print(delta)

file = open('iterated_profile','w')
for i in range(ny):
    # print(i,y_arr[i],rho[i],pg[i],trad[i],er[i])
    file.write(str(i)+'\t'+str(y_arr[i])+'\t'+str(rho[i])+'\t'+str(pg[i])+'\t'+str(trad[i])+'\t'+str(er[i]) + '\n')
file.close()

# plt.semilogy(y_arr,rho)
# plt.show()
