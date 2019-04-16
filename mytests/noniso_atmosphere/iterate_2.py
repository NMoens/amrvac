import numpy as np
import matplotlib.pyplot as plt

rho0 = 0.30543206e-07
er0 = 0.14784754e+07

pr0 = 1/3*er0
gamma = 0.5
pg0 = 0.48944484e+06 # (1-gamma)*pr0/gamma #

delta = pr0 - pg0*gamma/(1 - gamma)

print(delta, delta/pr0, delta/pg0)

a20 = pg0/rho0
print(a20)

rsun = 6.96e10
msun = 1.99e33
grav = 6.67e-8
mu = 0.6
kb = 1.3806e-16
cc = 2.99e10
mp = 1.672e-24
arad = 7.5646e-15


rstar = 19.370383317782409*rsun
rmax = 20.0*rsun
mstar = 50.*msun
gstar = grav*mstar/rstar**2.
geff = gstar*(1.-gamma)

Heff0 = a20/geff
print(Heff0)


dy = 138087438.2

y_arr = np.arange(rstar,rmax,dy)

ny = len(y_arr)

#Make everything in arrays
rho = rho0*np.exp(-(y_arr-rstar)/Heff0)
pg = rho*a20
pr = gamma*pg/(1-gamma) + delta
er = 3.0*pr
trad = (er/arad)**0.25
tgas = mp*mu/kb*pg/rho

print(pr)

print(tgas)
print(trad)

plt.figure(1)
plt.plot(y_arr,tgas,'r')
plt.plot(y_arr,trad,'b')

plt.figure(2)
plt.semilogy(y_arr,rho,'k')


max_err = 1000
while max_err > 1.e-15:
    tgas = trad
    a2 = kb*tgas/(mu*mp)
    Heff = a2/geff
    rho = rho0*np.exp(-(np.array(y_arr)-rstar)/Heff)
    pg = rho*a2
    pr = gamma*pg/(1-gamma) + delta
    er = 3.0*pr
    trad = (er/arad)**0.25
    tgas = mp*mu/kb*pg/rho
    plt.figure(1)
    plt.plot(y_arr,tgas,'r')
    plt.plot(y_arr,trad,'b')
    plt.figure(2)
    plt.semilogy(y_arr,rho,'k')

    max_err = max(abs(tgas-trad)/tgas)

# plt.figure(1)
# plt.plot(y_arr,tgas,'ko')
# plt.show()


file = open('iterated_gamma05','w')
for i in range(ny):
    file.write(str(i)+'\t'+str(y_arr[i])+'\t'+str(rho[i])+'\t'+str(pg[i])+'\t'+str(trad[i])+'\t'+str(er[i]) + '\n')
file.close()
