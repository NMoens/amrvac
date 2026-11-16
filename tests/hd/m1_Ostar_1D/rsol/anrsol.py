import numpy as np,re,glob,matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
D0='/home/nicolasm/Codes/amrvac3/amrvac/tests/hd/m1_Ostar_1D/rsol'
def rd(f):
    fh=open(f); fh.readline(); n=int(fh.readline()); t=float(fh.readline())
    return t,np.array([[float(re.sub(r"(\d)([+-]\d{3})$",r"\1E\2",x)) for x in fh.readline().split()] for _ in range(n)])
Rs=[0.1,0.3,1.0]; col={0.1:'k',0.3:'#1f77b4',1.0:'#d62728'}
D={}
for R in Rs:
    T=[];S=[]
    for f in sorted(glob.glob(f'{D0}/R{R}0*.blk')): t,d=rd(f); T.append(t); S.append(d)
    D[R]=(np.array(T),np.array(S))
r=D[0.1][1][0,:,0]; dr=r[1]-r[0]; dV=4*np.pi*((r+dr/2)**3-(r-dr/2)**3)/3; A=4*np.pi*r**2
G=2.1898172642416400E-003; M=24524.996128252151; gam=1.667
I=lambda s,c,x: np.interp(x,r,s[:,c]); m=(r>1.03)&(r<1.25)
ib=np.argmin(abs(r-4.89))
def rph(s):
    j=np.where(s[:,22]<0.5)[0]; return r[j[0]] if len(j) else np.nan
print("R_hat  nsnap t_end | <L_out> t=.1-.5  t=.25-.5 | min L_out (t) | L_out at .5 | <L_F> at r=1.10/1.15/1.20/1.25/1.5 (t=.25-.5) | R_ph t=0/.25/.5 | Gmax range t>.1 | v1.15 rms | E/R content r>1.1 change")
for R in Rs:
    T,S=D[R]; Ls=np.median(A*S[0,:,5]/S[0,:,18])
    a=(T>=0.1); b=(T>=0.25)
    Lo=S[:,ib,18]; k=np.argmin(Lo)
    LF=[np.mean(S[b,np.argmin(abs(r-x)),18]) for x in [1.10,1.15,1.20,1.25,1.5]]
    rp=[rph(S[np.argmin(abs(T-x))]) for x in [0,0.25,0.5]]
    G_=[z[m,9].max() for z in S[a]]
    v15=np.std([I(z,2,1.15) for z in S[a]])
    mm=r>1.10; Ec=np.sum(S[:,mm,4]*dV[mm],axis=1)/Ls/R
    print(f"{R:4.1f}  {len(T):3d} {T[-1]:.3f} | {Lo[a].mean():.3f}  {Lo[b].mean():.3f} | {Lo[k]:.3f} ({T[k]:.2f}) | {Lo[-1]:.3f} | "+'/'.join(f"{x:.2f}" for x in LF)+f" | {rp[0]:.3f}/{rp[1]:.3f}/{rp[2]:.3f} | {min(G_):.2f}-{max(G_):.2f} | {v15:.3f} | {Ec[0]:.2f}->{Ec[-1]:.2f}")
# budget above r_a with E/R, t=0.1-end
print("\nBudget shell [r_a,4.89], t=0.1-end (L*): in-out, d(gas+pot)/dt, d(E/R)/dt, residual")
for R in Rs:
    T,S=D[R]; Ls=np.median(A*S[0,:,5]/S[0,:,18]); a=T>=0.1; Ta=T[a]; dt=Ta[-1]-Ta[0]
    rho,v,p,E=S[...,1],S[...,2],S[...,3],S[...,4]; eg=p/(gam-1)+0.5*rho*v**2-rho*G*M/r
    Lt=S[...,18]+S[...,19]+S[...,20]+A*rho*v*(-G*M/r)/Ls
    out=[]
    for ra in [1.15,1.20]:
        ia=np.argmin(abs(r-ra)); w=np.zeros(len(r)); w[ia+1:ib]=1; w[ia]=w[ib]=0.5
        C=lambda x: np.sum(x[a]*dV*w,axis=-1)/Ls
        io=(np.trapz(Lt[a,ia],Ta)-np.trapz(Lt[a,ib],Ta))/dt; dg=(C(eg)[-1]-C(eg)[0])/dt; dE=(C(E)[-1]-C(E)[0])/dt/R
        out.append(f"r_a={ra}: {io:+.3f} {dg:+.3f} {dE:+.3f} res {io-dg-dE:+.3f}")
    print(f"  R={R}: "+' | '.join(out))
# plots
fig,ax=plt.subplots(2,3,figsize=(16,9))
for R in Rs:
    T,S=D[R]; c=col[R]; lab=f'RSOL={R}'
    ax[0,0].plot(T,S[:,ib,18],c=c,label=lab)
    ax[0,0].plot(T,S[:,ib,18]+S[:,ib,19]+S[:,ib,20],c=c,ls=':',lw=1.6)
    ax[0,1].plot(T,[I(z,2,1.15) for z in S],c=c,label=lab)
    ax[0,2].plot(T,[rph(z) for z in S],c=c,label=lab)
    b=T>=0.25
    ax[1,0].plot(r,S[b,:,18].mean(axis=0),c=c,label=lab)
    ax[1,0].plot(r,(S[b,:,18]+S[b,:,19]+S[b,:,20]).mean(axis=0),c=c,ls=':',lw=1.6)
    ax[1,1].semilogy(r,S[-1,:,1],c=c,label=lab)
    ax[1,2].plot(r,S[-1,:,2],c=c,label=lab)
ax[0,0].set(xlabel='t',ylabel=r'$L(r=4.9)/L_*$',title='emergent luminosity (solid $L_F$, dotted total)'); ax[0,0].axhline(1,c='gray',ls=':')
ax[0,1].set(xlabel='t',ylabel=r'$v_r(1.15)$',title='iron-bump velocity')
ax[0,2].set(xlabel='t',ylabel=r'$R(\theta=0.5)$',title='photosphere radius')
ax[1,0].set(xlabel='r',ylabel=r'$\langle L\rangle/L_*$',title='mean L over t=0.25-0.5 (solid $L_F$, dotted total)',xlim=(1,1.6),ylim=(0,2.5)); ax[1,0].axhline(1,c='gray',ls=':')
ax[1,1].set(xlabel='r',ylabel=r'$\rho$',title='density at t=0.5',xlim=(1,1.6))
ax[1,2].set(xlabel='r',ylabel=r'$v_r$',title='velocity at t=0.5',xlim=(1,1.6),ylim=(-3,6))
for a in ax.flat: a.grid(alpha=.3); a.legend(fontsize=8)
fig.suptitle('M1 Star8 1D: RSOL sensitivity'); fig.tight_layout(); fig.savefig(f'{D0}/rsol_comparison.png',dpi=130); print('saved')
