import numpy as np
import matplotlib.pyplot as plt
import copy

m = 0.2 #1.3694919e-2
Gamma = 2 #2.2044134

epsilon_in = 1.e-7
epsilon_out0 = 0.5
dx0 = 1e-5
it1_toll = 1e-10
itot = 0

def System(w,q,x):
    # dwdx = Gamma*(1-m*(w+x)+4*q*np.sqrt(abs(w)))-1
    dwdx = Gamma*(1-m*(w+x)-4*q*np.sqrt(abs(w)))-1
    dqdx = -m*Gamma/np.sqrt(abs(w))*(1-m*(w+x))+2*q/(1-x)+4*m*Gamma*q
    return dwdx,-1*dqdx

q_in = 1

#Solving the system:
def Solve_system_out(w_zpe,q_zpe, epsilon_out):
    w = copy.copy(w_zpe)
    q = copy.copy(q_zpe)
    x = copy.copy(epsilon_in)
    w_array = [w]
    q_array = [q]
    x_array = [x]

    while x < 1-epsilon_out:
        w_min = copy.copy(w)
        q_min = copy.copy(q)
        err = 1
        tol = 1e-10
        dx = dx0

        while err > tol:
            w = abs(w)
            w_min = abs(w_min)
            # w_min = min(w_min,(x+dx)*(Gamma-1))
            # w_new = (w + dx*Gamma -m*dx*Gamma*(x+dx)+4*q_min*np.sqrt(abs(w_min))*dx*Gamma - dx)/(1 + dx*Gamma*m)
            w_new = (w + dx*Gamma -m*dx*Gamma*(x+dx)-4*q_min*np.sqrt(abs(w_min))*dx*Gamma - dx)/(1 + dx*Gamma*m)
            q_new = (q + dx*(-m*Gamma)/np.sqrt(abs(w_min))*(1-m*(w_min+x+dx)))/(1 - 2*dx/(1-x-dx) - 4*dx*m*Gamma)

            w_min = w_new
            q_min = q_new

            err = (((w_min - w_new)/w_new)**2 + ((q_min - q_new)/q_new)**2)**0.5

        # w_new = min(w_new,(x+dx)*(Gamma-1))

        w = w_new
        q = q_new
        x = x+dx

        w_array.append(w)
        q_array.append(q)
        x_array.append(x)

    return w_array, q_array, x_array

def Solve_system_in(w_arr,q_ome,x_arr):
    q = copy.copy(q_ome)
    q_array = [q]

    dx = -dx0
    for i in range(len(x_arr)-1,-1,-1):
        err = 1
        tol = 1e-7
        w = w_arr[i]
        x = x_arr[i]

        q_new = (q + dx*(-m*Gamma)/np.sqrt(abs(w))*(1-m*(w+x+dx)))/(1 - 2*dx/(1-x-dx) - 4*dx*m*Gamma)


        q = q_new

        q_array.append(q)
    return q


def IterativeScheme(w_zpe, q_zpe, epsilon_out, data):
    #choose initial value for w and q:
    #w and q at epsilon_in (zero plus epsilon_in = zpe)
    global itot
    err = 1
    it = 0

    while abs(err) > it1_toll and  it < 30:
        #integrate both w and q outward from w_zpe and q_zpe:
        w_arr,q_arr,x_arr = Solve_system_out(w_zpe,q_zpe,epsilon_out)

        plt.figure(1)
        plt.loglog(x_arr,w_arr)

        #for this outer value of w, set 'correct' outer value of q
        w_ome = w_arr[-1]
        q_ome = epsilon_out*(m*Gamma/np.sqrt(w_ome)*(1-m*w_ome-m+m*epsilon_out))/(1+2+4*m*Gamma*epsilon_out)

        #integrate only q inward from 'correct' outer value for q_zpe
        #Hereby setting the new value for q_zpe
        q_zpe_old = q_zpe
        q_zpe = Solve_system_in(w_arr,q_ome,x_arr)

        old_err = err
        err = (q_zpe - q_zpe_old)/q_zpe_old

        if abs(err) > abs(old_err):
            print('diverging')
            exit()

        it = it+1
        itot = itot+1

        print(it, err, q_zpe)

        plt.figure(6)
        plt.plot(itot,q_zpe,'bo')

    if data:
        return w_arr,q_arr,x_arr
    else:
        return q_zpe

def xtoxpo(epsilon):
    x = 1-epsilon
    r = 1/(1-x)
    r = r + r
    x = 1 - 1/r
    epsilon = 1 - x
    return epsilon

def Second_iteration():

    p_in = q_in*(1-epsilon_in)**2
    q_zpe = p_in -2*m*Gamma*np.sqrt(epsilon_in/(Gamma-1))+2*p_in*(1+2*m*Gamma)*epsilon_in
    w_zpe = (Gamma-1)*epsilon_in

    epsilon_out = epsilon_out0
    it2 = 0

    while epsilon_out > 1e-3 and it2 < 100:
        q_zpe = IterativeScheme(w_zpe, q_zpe, epsilon_out, False)
        w_arr,q_arr,x_arr = Solve_system_out(w_zpe,q_zpe,epsilon_out)

        plt.figure(2)
        plt.loglog(x_arr,w_arr)

        # epsilon_out = 0.9*epsilon_out
        epsilon_out = xtoxpo(epsilon_out)

        it2 = it2 + 1

        x_out = 1-epsilon_out
        r_out = 1/(1-x_out)

        print('##########')
        print('# ', it2, epsilon_out, r_out)
        print('##########')
    return w_arr,q_arr,x_arr


ww,qq,xx = Second_iteration()

dwdx = copy.copy(ww)
min_dqdx = copy.copy(qq)
for i in range(len(xx)):
    dwdx[i],min_dqdx[i] = System(ww[i],qq[i],xx[i])

dw = copy.copy(ww)
dq = copy.copy(qq)
for i in range(1,len(xx)-1):
    dw[i] = (ww[i+1]-ww[i-1])/(xx[i+1]-xx[i-1])
    dq[i] = -(qq[i+1]-qq[i-1])/(xx[i+1]-xx[i-1])

dw[0] = dw[1]
dw[-1] = dw[-2]
dq[0] = dq[1]
dq[-1] = dq[-2]

print('converting ......')
wg = [w/(Gamma-1) for w in ww]
p = copy.copy(qq)
for i in range(len(qq)):
    p[i] = qq[i]*(1-xx[i])**2
print('done')

wg = np.array(wg)
p = np.array(p)
xx = np.array(xx)

fname = 'model_G' + str(Gamma) + '_m' + str(m)
np.savetxt(fname,np.transpose([wg,p,xx]))

plt.figure(3)
plt.title('w/(gamma -1)')
plt.semilogy(xx,wg, 'r-')
plt.loglog(xx,xx,'k')
plt.xlim([1e-5,1])
plt.ylim([1e-5,1])

plt.figure(4)
plt.title('q(1-x)^2')
plt.loglog(xx,p, 'b-')
plt.xlim([1e-5,1])
plt.ylim([1e-5,1])

plt.figure(5)
plt.title('Derivatives')
plt.loglog(xx,dwdx,label='dwdx')
plt.loglog(xx,min_dqdx,label='dqdx')
plt.loglog(xx,dw,'rx',label='dwdx numerical')
plt.loglog(xx,dq,'bx',label='dqdx numerical')
plt.xlim([1e-5,1])
plt.legend()

print('Gamma = ',Gamma)
print('m =', m)
print('epsilon_in =', epsilon_in)

plt.show()
