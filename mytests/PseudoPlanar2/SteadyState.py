import numpy as np
import matplotlib.pyplot as plt

M_dot = 1.
P_dot = 1.
e_dot = 1.
E_dot = 1.

gamma = 5./3.

ny = 100
y = np.linspace(1,15,ny)

def Get_rho(v):
    rho = M_dot/(4*np.pi*y**2*v)
    return rho

def Get_p(rho,v):
    p = P_dot/(4*np.pi*y**2) -v**2*rho
    return p

def Get_e(v,p):
    e = e_dot/(4*np.pi*y**2*v) - p
    return e

def Get_E(v):
    E = E_dot/(4*np.pi*y**2*v)
    return E

def Get_v(p):
    a = M_dot/(8*np.pi*y**2)
    b = p*gamma/(gamma -1)
    c = -e_dot/(4*np.pi*y**2)
    D = b**2 - 4*a*c
    v = (-b - np.sqrt(D))/(2*a)
    return v

rho_arr = np.ones([ny])
v_arr = np.ones([ny])
e_arr = np.ones([ny])
p_arr = np.ones([ny])
E_arr = np.ones([ny])

for i in range(3):
    v_arr = Get_v(p_arr)
    rho_arr = Get_rho(v_arr)
    p_arr = Get_p(rho_arr,v_arr)

    plt.plot(y,p_arr)


e_arr = Get_e(v_arr,p_arr)
E_arr = Get_E(v_arr)

plt.show()
