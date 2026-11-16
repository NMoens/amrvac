"""Fit the measured early drift dF/dt (code - exact) to the drift predicted by each
suspected source-term error in mod_m1.t. Usage: attribute.py <blk> <iprob>"""
import numpy as np,re,sys
f,ip=sys.argv[1],int(sys.argv[2]); V=0.1; F0=0.5; c=1.0
fh=open(f); hdr=fh.readline().split(); dims=[int(x) for x in fh.readline().split()]; t=float(fh.readline())
a=np.array([[float(re.sub(r"(\d)([+-]\d{3})$",r"\1E\2",x)) for x in l.split()] for l in fh if l.strip()])
C={n:i for i,n in enumerate(hdr)}
r=a[:,C['r']]; th=a[:,C['Theta']]; ph=a[:,C['Phi']] if 'Phi' in C else 0*r
def sph(ax,ay,az): return (ax*np.sin(th)*np.cos(ph)+ay*np.sin(th)*np.sin(ph)+az*np.cos(th),
                            ax*np.cos(th)*np.cos(ph)+ay*np.cos(th)*np.sin(ph)-az*np.sin(th), -ax*np.sin(ph)+ay*np.cos(ph))
Fr,Ft,Fp=sph(F0,0,0) if ip in (5,6) else sph(0,0,F0)
vr,vt,vp={1:sph(V,0,0),2:sph(0,0,V),5:sph(0,0,0),6:sph(0,0,V)}[ip]
cot=1/np.tan(th)
meas=np.array([(a[:,C['F1']]-Fr),(a[:,C['F2']]-Ft),(a[:,C['F3']]-Fp)])/t if 'F3' in C else None
z=0*r
# suspect 1: geom_corr added with the wrong sign -> code - correct = 2*geom_corr
gc=np.array([-(Ft*vt+Fp*vp)/r, Ft*vr/r-Fp*vp*cot/r, Fp*vr/r+Fp*vt*cot/r])
if dims and len(dims)==2: gc[2]=0*r   # 2.5D: no idim=3 loop, phi correction never applied (code - correct = -(-gc3)... handled below)
S1=2*gc
# suspect 2: T_r_theta used instead of T_theta_r, T_r_phi instead of T_phi_r in m1_add_source_geom
S2=np.array([z, -(vr*Ft-vt*Fr)/r, -(vr*Fp-vp*Fr)/r])
# suspect 3: 2*ppt with ppt = v_theta F_phi instead of 1*(v_phi F_theta)
S3=np.array([z, z, -cot*(2*vt*Fp-vp*Ft)/r])
# suspect 4 (2.5D only): phi-component of (F.grad)v correction missing entirely
S4=np.array([z,z,-(Fp*vr/r+Fp*vt*cot/r)])
# suspect 5: c^2 P_phi_theta counted with the wrong weight in the phi source.
# Code now uses 1*P_phi_theta; S5 = drift if the correct weight were 2 (code - correct = +cot c^2 P_pt / r).
ff=(F0/(c*1.0))**2; chi=(3+4*ff)/(5+2*np.sqrt(4-3*ff)); nt=Ft/F0; npp=Fp/F0
Ppt=(3*chi-1)/2*nt*npp*1.0
S5=np.array([z,z,cot*c**2*Ppt/r])
# interior cells only (2 layers from each boundary)
idx=np.arange(len(r)); n1=dims[0]; i1=idx%n1; i2=(idx//n1)%dims[1]; m=(i1>=2)&(i1<n1-2)&(i2>=2)&(i2<dims[1]-2)
if len(dims)==3: i3=idx//(n1*dims[1]); m&=(i3>=2)&(i3<dims[2]-2)
names=['S1 geom_corr sign','S2 T_r.th/T_th.r order','S3 2*ppt (advective)','S4 2.5D phi corr missing','S5 c2P_pt weight 2 vs 1']
Ss=[S1,S2,S3,S4,S5]
A=np.stack([s[:,m].ravel() for s in Ss],1); b=meas[:,m].ravel()
keep=[k for k in range(5) if np.abs(A[:,k]).max()>1e-12]
coef,res,*_=np.linalg.lstsq(A[:,keep],b,rcond=None)
print(f"{f.split('/')[-1]}  t={t:.3f}  |measured dF/dt| rms = {np.sqrt(np.mean(b**2)):.3e}")
for k,cf in zip(keep,coef): print(f"   {names[k]:28s} coefficient {cf:+.3f}   (rms of term {np.sqrt(np.mean(A[:,k]**2)):.3e})")
fit=A[:,keep]@coef; print(f"   residual rms after fit = {np.sqrt(np.mean((b-fit)**2)):.3e}  ({100*np.sqrt(np.mean((b-fit)**2))/np.sqrt(np.mean(b**2)):.1f}% of measured)")
for k in keep:
    print(f"   {names[k]:28s} alone with coef 1: residual {100*np.sqrt(np.mean((b-A[:,k])**2))/np.sqrt(np.mean(b**2)):.1f}%")
