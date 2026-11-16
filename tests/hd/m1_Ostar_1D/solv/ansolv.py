import numpy as np,re,glob
D0='/home/nicolasm/Codes/amrvac3/amrvac/tests/hd/m1_Ostar_1D/solv'
def rd(f):
    fh=open(f); fh.readline(); n=int(fh.readline()); t=float(fh.readline())
    return t,np.array([[float(re.sub(r"(\d)([+-]\d{3})$",r"\1E\2",x)) for x in fh.readline().split()] for _ in range(n)])
runs=['Euler_tvdlf_minmod','Midpoint_tvdlf_minmod','Euler_hll_minmod','Euler_tvdlf_koren','Euler_tvdlf_weno3']
D={}
for nm in runs:
    T=[];S=[]
    for f in sorted(glob.glob(f'{D0}/{nm}0*.blk')): t,d=rd(f); T.append(t); S.append(d)
    D[nm]=(np.array(T),np.array(S))
T0,S0=D[runs[0]]; r=S0[0,:,0]; I=lambda s,c,x: np.interp(x,r,s[:,c])
f=np.zeros(len(r)+1); f[0]=1.0
for i in range(len(r)): f[i+1]=2*r[i]-f[i]
dV=(f[1:]**3-f[:-1]**3)/3; m=(r>1.03)&(r<1.25)
# columns: 0 r 1 rho 2 v 3 p 4 E 5 F 9 Gamma 18 lum 19 lum_gas 20 lum_rad 22 theta
print(f"{'run':22s} nsnap tend  rms v1.05 rms v1.15 | v1.15 t=.05/.10/.15  | Gmax  | Lout mean/rms | L_F base mean/rms | dM%   | <|drho|/rho> vs ref r<1.3 | NaN")
for nm in runs:
    T,S=D[nm]; k=min(len(T),len(T0))
    v15=np.array([I(z,2,1.15) for z in S]); v05=np.array([I(z,2,1.05) for z in S])
    at=lambda tt: v15[np.argmin(abs(T-tt))]
    dm=100*(np.sum(S[-1,:,1]*dV)-np.sum(S[0,:,1]*dV))/np.sum(S[0,:,1]*dV)
    mm=r<1.3; drho=np.mean(abs(S[:k,mm,1]-S0[:k,mm,1])/S0[:k,mm,1])
    print(f"{nm:22s} {len(T):3d}  {T[-1]:.3f}  {np.std(v05):.3f}     {np.std(v15):.3f}   | {at(0.05):+.3f}/{at(0.10):+.3f}/{at(0.15):+.3f} | {max(z[m,9].max() for z in S):.2f}  | {S[:,-1,18].mean():.3f}/{S[:,-1,18].std():.3f}  | {S[:,0,18].mean():.3f}/{S[:,0,18].std():.3f}       | {dm:+.2f} | {drho:.4f} | {np.isnan(S).sum()}")
print('\nAt t_end: rho, v, Gamma, L_F, theta')
for x in [1.02,1.08,1.12,1.15,1.20,1.25,1.5,3.0,4.9]:
    print(f"r={x:4.2f}: "+' | '.join(f"{I(D[nm][1][-1],1,x):.3e} {I(D[nm][1][-1],2,x):+.3f} {I(D[nm][1][-1],9,x):.2f} {I(D[nm][1][-1],18,x):.2f} {I(D[nm][1][-1],22,x):.2f}" for nm in runs))
print('\nv(1.15) time series (every 0.01):')
for nm in runs:
    T,S=D[nm]; print(f"{nm:22s} "+' '.join(f"{I(z,2,1.15):+.2f}" for z in S))
print('\nGamma_max (1.03<r<1.25) time series:')
for nm in runs:
    T,S=D[nm]; print(f"{nm:22s} "+' '.join(f"{z[m,9].max():.2f}" for z in S))
print('\nEmergent L time series:')
for nm in runs:
    T,S=D[nm]; print(f"{nm:22s} "+' '.join(f"{z[-1,18]:.2f}" for z in S))
