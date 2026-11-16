import numpy as np
S=np.load('solv/longrun.npy'); T=np.load('solv/longrun_t.npy'); R=0.1
r=S[0,:,0]; dr=r[1]-r[0]; dV=4*np.pi*((r+dr/2)**3-(r-dr/2)**3)/3; A=4*np.pi*r**2
G=2.1898172642416400E-003; M=24524.996128252151; gam=1.667
Ls=np.median(A*S[0,:,5]/S[0,:,18])
rho,v,p,E=S[...,1],S[...,2],S[...,3],S[...,4]
eg=p/(gam-1)+0.5*rho*v**2; ep=-rho*G*M/r
Ltot=S[...,18]+S[...,19]+S[...,20]+A*rho*v*(-G*M/r)/Ls
print("Budget  t=0.3-3.3, shell [r_a, 4.89]:  in - out  vs  d/dt(content), content with E and E/R")
m=T>=0.3; T_=T[m]; dt=T_[-1]-T_[0]
ib=np.argmin(abs(r-4.89))
for ra in [1.02,1.10,1.15,1.20]:
    ia=np.argmin(abs(r-ra)); w=np.zeros(len(r)); w[ia+1:ib]=1; w[ia]=w[ib]=0.5
    C=lambda x: np.sum(x*dV*w,axis=-1)/Ls
    Lin=np.trapz(Ltot[m,ia],T_)/dt; Lout=np.trapz(Ltot[m,ib],T_)/dt
    dg=(C(eg+ep)[m][-1]-C(eg+ep)[m][0])/dt; dE=(C(E)[m][-1]-C(E)[m][0])/dt
    print(f"  r_a={r[ia]:.3f}: in-out = {Lin-Lout:+.3f} L*   d(gas+pot)/dt = {dg:+.3f}   dE/dt = {dE:+.3f}  -> dE/dt/R = {dE/R:+.3f}   residual(E) = {Lin-Lout-dg-dE:+.3f}  residual(E/R) = {Lin-Lout-dg-dE/R:+.3f}")
# radiative-only view: time-averaged L_F profile and E/R storage per region
print("\nRadiative energy stored (E/R) per region, rate over t=0.3-3.3 in L*:")
i0=np.where(m)[0][0]
for a,b in [(1.0,1.10),(1.10,1.20),(1.20,1.30),(1.30,4.9)]:
    mm=(r>=a)&(r<b); d=np.sum((E[-1]-E[i0])[mm]*dV[mm])/Ls/dt
    print(f"  {a:.2f}-{b:.2f}: dE/dt = {d:+.4f} -> /R = {d/R:+.3f}")
# effective thermal timescale (E/R + e_int) / L*
for ra in [1.0,1.10,1.15]:
    mm=r>=ra
    print(f"r>{ra}: E/R content = {np.sum(E[0][mm]*dV[mm])/Ls/R:.2f} L*.t, gas internal = {np.sum((p[0]/(gam-1))[mm]*dV[mm])/Ls:.2f} L*.t  (pulsation period 0.23)")
# time-averaged radiative L profile
print("\n<L_F>(r), t=0.3-3.3:", ' '.join(f"{x:.2f}:{np.trapz(S[m,np.argmin(abs(r-x)),18],T_)/dt:.3f}" for x in [1.01,1.05,1.10,1.15,1.20,1.25,1.3,2,4.9]))
# trend of photospheric E content: is E above 1.1 still rising?
mm=r>=1.10
Ec=np.sum(E[:,mm]*dV[mm],axis=1)/Ls/R
print("E/R content r>1.10 vs t:", ' '.join(f"{T[k]:.1f}:{Ec[k]:.2f}" for k in range(0,len(T),30)))
mm=r>=1.0
Ec=np.sum(E[:,mm]*dV[mm],axis=1)/Ls/R
print("E/R content r>1.00 vs t:", ' '.join(f"{T[k]:.1f}:{Ec[k]:.2f}" for k in range(0,len(T),30)))
