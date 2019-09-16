import numpy as np
import matplotlib.pyplot as plt
import os

def Open_file(file):
    file = os.getcwd() + '/' + file
    x,y,z,rho,m,e,Er = np.loadtxt(file, unpack = True, skiprows = 15)
    return x, rho, m, e, Er

# x0, rho0, m0, e0, Er0 = Open_file('tau_1d0.okc')
x3, rho3, m3, e3, Er3 = Open_file('tau_1d3.okc')
# x9, rho9, m9, e9, Er9 = Open_file('tau_1d9.okc')

x = np.linspace(0,20,1000)

def Signal(L,x):
    base = np.mean(Er3)
    ampl = np.max(Er3) - np.mean(Er3)

    s=[]
    for xx in x:
        if xx >= 1:
            s.append(base+ampl*np.sin(xx*2*np.pi)*np.exp(-(xx-1)/L))
        else:
            s.append(base+ampl*np.sin(xx*2*np.pi))
    return s


def DampingLength(tau):
    Bo = 1e-3
    gamma = 1.666667
    omega = 6.28
    c= 2.99e10/2998295.0

    L = 2*c*tau/omega

    # L = Bo/(8*tau*(gamma-1))
    # L = 1/(np.sqrt(3)*tau)

    print(L)
    return L

plt.title('Radiation energy density')
# plt.plot(x0, Er0, 'b.', label= '$\\tau_\\Lambda = 1d0$')
plt.plot(x3, Er3, 'r.', label= '$\\tau_\\Lambda = 1d3$')
# plt.plot(x9, Er9, 'g.', label= '$\\tau_\\Lambda = 1d9$')

# plt.plot(x,Signal(DampingLength(1.e0),x),'b')
plt.plot(x,Signal(DampingLength(1.e3),x),'r')
# plt.plot(x,Signal(DampingLength(1.e9),x),'g')


plt.axvspan(0., 1., facecolor='b', alpha=0.5)

plt.xlim([0,20])

plt.legend()

plt.show()
