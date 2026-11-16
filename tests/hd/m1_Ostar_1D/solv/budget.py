import numpy as np,re,glob
def rd(f):
    fh=open(f); fh.readline(); n=int(fh.readline()); t=float(fh.readline())
    return t,np.array([[float(re.sub(r"(\d)([+-]\d{3})$",r"\1E\2",x)) for x in fh.readline().split()] for _ in range(n)])
fs=sorted(glob.glob('output/Star8_1d0*.blk'))
T=[];S=[]
for f in fs: t,d=rd(f); T.append(t); S.append(d)
T=np.array(T); S=np.array(S); np.save('solv/longrun.npy',S); np.save('solv/longrun_t.npy',T)
r=S[0,:,0]; dr=r[1]-r[0]; fl=r-dr/2; fr=r+dr/2; dV=4*np.pi*(fr**3-fl**3)/3; A=4*np.pi*r**2
G=2.1898172642416400E-003; M=24524.996128252151; gam=1.667
# L* in code units (energy/time) from lum column: lum = A F / Lstar_code
Ls=np.median(A*S[0,:,5]/S[0,:,18])
rho,v,p,E=S[...,1],S[...,2],S[...,3],S[...,4]
eg=p/(gam-1); ek=0.5*rho*v**2; ep=-rho*G*M/r
# total energy flux (code), in L*: radiative + advective (gas,rad) from diagnostics + kinetic... lum_gas uses e (total gas) + p -> includes kinetic
Lgrav=A*rho*v*(-G*M/r)/Ls
Ltot=S[...,18]+S[...,19]+S[...,20]+Lgrav
print('L* code =',Ls, ' t_end',T[-1], ' nsnap',len(T))
ra,rb=1.02,4.9; ia=np.argmin(abs(r-ra)); ib=np.argmin(abs(r-rb)); sl=slice(ia+1,ib)  # cells strictly between centers ia and ib
# energy content between cell centres ia and ib (approx: whole cells ia+1..ib-1 plus half cells)
w=np.zeros(len(r)); w[ia+1:ib]=1; w[ia]=0.5; w[ib]=0.5
def content(x): return np.sum(x*dV*w,axis=-1)/Ls
Cg,Ck,Cr,Cp=content(eg),content(ek),content(E),content(ep); Ct=Cg+Ck+Cr+Cp
print(f"budget for shell r={r[ia]:.4f}..{r[ib]:.4f}")
for a,b in [(0,0.5),(0.5,1),(1,2),(2,T[-1]+1e-9),(0.3,T[-1]+1e-9)]:
    m=(T>=a)&(T<=b); i0,i1=np.where(m)[0][[0,-1]]; dt=T[i1]-T[i0]
    Lin=np.trapz(Ltot[m,ia],T[m])/dt; Lout=np.trapz(Ltot[m,ib],T[m])/dt
    LFout=np.trapz(S[m,ib,18],T[m])/dt; LFin=np.trapz(S[m,ia,18],T[m])/dt
    dC=lambda C:(C[i1]-C[i0])/dt
    print(f" t {T[i0]:.2f}-{T[i1]:.2f}: <L_in,tot>={Lin:.3f} (rad {LFin:.3f}) <L_out,tot>={Lout:.3f} (rad {LFout:.3f}, grav+adv {Lout-LFout:+.4f}) | in-out={Lin-Lout:+.3f}  dE/dt={dC(Ct):+.3f} [gas {dC(Cg):+.3f} kin {dC(Ck):+.3f} rad {dC(Cr):+.3f} pot {dC(Cp):+.3f}]  residual={Lin-Lout-dC(Ct):+.3f}")
# where is energy stored: regions
print('\nStored energy change (L* x time units) by region, t=0 -> end:')
for a,b in [(1.02,1.10),(1.10,1.20),(1.20,1.30),(1.30,2.0),(2.0,4.9)]:
    mm=(r>=a)&(r<b)
    d=lambda x: (np.sum((x[-1]-x[0])[mm]*dV[mm]))/Ls
    print(f"  {a:.2f}-{b:.2f}: gas int {d(eg):+.3f} kin {d(ek):+.3f} rad {d(E):+.3f} pot {d(ep):+.3f}  total {d(eg+ek+E+ep):+.3f}   mass change {100*np.sum((rho[-1]-rho[0])[mm]*dV[mm])/np.sum(rho[0]*dV):+.2f}% of total")
# compare energy timescale: total stored energy above 1.1 vs L*
mm=r>1.10; print('\nThermal+rad energy above r=1.10 at t=0 (in L* x time):', np.sum((eg[0]+E[0])[mm]*dV[mm])/Ls, ' gravitational binding:',np.sum(ep[0][mm]*dV[mm])/Ls)
# Mass flux / wind mechanical luminosity at 4.9
md=A*rho*v; Lw=md*(0.5*v**2+G*M/r)/Ls
print('wind: <Mdot v^2/2 + GM/r ... > at 4.9 (L*):', np.mean(md[:,ib]*0.5*v[:,ib]**2/Ls), ' Mdot*GM/R*(=lift from 1):', np.mean(md[:,ib]*G*M/1.0/Ls))
print('mass flux at base r=1.02 mean (code):',np.mean(md[T>=0.3,ia]),' at 4.9:',np.mean(md[T>=0.3,ib]))
print('Time series (t, L_rad out, L_tot in, stored E_tot rel t0, R_ph(theta=.5)):')
for k in range(0,len(T),20):
    th=S[k,:,22]; j=np.where(th<0.5)[0][0]
    print(f"  {T[k]:.2f} {S[k,ib,18]:.3f} {Ltot[k,ia]:+.3f} {Ct[k]-Ct[0]:+.4f} {r[j]:.3f}")
