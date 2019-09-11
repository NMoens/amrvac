import numpy as np
import matplotlib.pyplot as plt
import os

def Open_file(file):
    file = os.getcwd() + '/' + file
    x,y,z,rho,m,e,Er = np.loadtxt(file, unpack = True, skiprows = 15)
    return x, rho, m, e, Er

x0, rho0, m0, e0, Er0 = Open_file('tau_1d0.okc')
x3, rho3, m3, e3, Er3 = Open_file('tau_1d3.okc')
x9, rho9, m9, e9, Er9 = Open_file('tau_1d9.okc')

x = np.linspace(0,20,1000)

def Signal(L,x):
    ampl = 5.9999998800000023e-003*6.2831853071795862
    print(ampl)
    s=[]
    for xx in x:
        if xx >= 1:
            s.append(1+ampl*np.sin(xx*2*np.pi)*np.exp(-(xx-1)/L))
        else:
            s.append(1+ampl*np.sin(xx*2*np.pi))
    return s


def DampingLength(tau):
    L = 1/(np.sqrt(3)*tau)
    return L

# plt.plot(x0, rho0, '.', label= '$\\tau_\\Lambda = 1d0$')
# plt.plot(x3, rho3, '.', label= '$\\tau_\\Lambda = 1d3$')

plt.plot(x9, rho9, '.', label= '$\\tau_\\Lambda = 1d9$')
L = DampingLength(1.e0)
plt.plot(x,Signal(L,x))


plt.axvspan(0., 1., facecolor='b', alpha=0.5)

plt.xlim([0,20])

plt.legend()

plt.show()
