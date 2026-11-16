import numpy as np,re,glob,matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
D0='/home/nicolasm/Codes/amrvac3/amrvac/tests/hd/m1_Ostar_1D/solv'
def rd(f):
    fh=open(f); fh.readline(); n=int(fh.readline()); t=float(fh.readline())
    return t,np.array([[float(re.sub(r"(\d)([+-]\d{3})$",r"\1E\2",x)) for x in fh.readline().split()] for _ in range(n)])
runs=[('Euler_tvdlf_minmod','reference (IMEX_Euler, TVDLF, minmod)','k','-'),
      ('Midpoint_tvdlf_minmod','IMEX_Midpoint','#1f77b4','--'),
      ('Euler_hll_minmod','HLL','#d62728','-.'),
      ('Euler_tvdlf_koren','koren','#2ca02c',':'),
      ('Euler_tvdlf_weno3','weno3','#ff7f0e',(0,(5,1,1,1,1,1)))]
fig,ax=plt.subplots(2,2,figsize=(12,8.5),sharex=True)
# columns: 0 r 1 rho 2 v 14 Tgas 15 Trad 18 lum(L_F) 19 lum_gas 20 lum_rad
for nm,lab,c,ls in runs:
    f=sorted(glob.glob(f'{D0}/{nm}0*.blk'))[-1]; t,d=rd(f); r=d[:,0]
    kw=dict(color=c,ls=ls,lw=1.4,label=f'{lab}  (t={t:.3f})')
    ax[0,0].semilogy(r,d[:,1],**kw); ax[0,1].plot(r,d[:,2],**kw)
    ax[1,0].semilogy(r,d[:,14],**kw)
    ax[1,1].plot(r,d[:,18],**kw)
    ax[1,1].plot(r,d[:,18]+d[:,19]+d[:,20],color=c,ls=ls,lw=0.7,alpha=0.5)
ax[0,0].set_ylabel(r'$\rho$ [code units]'); ax[0,1].set_ylabel(r'$v_r$ [code units]')
ax[1,0].set_ylabel(r'$T_{\rm gas}$ [K]'); ax[1,1].set_ylabel(r'$L/L_*$  (bold: radiative flux $L_F$; thin: total incl. advection)')
ax[1,1].set_ylim(0,2)
for a in ax.flat: a.set_xscale('log'); a.grid(alpha=0.3); a.set_xlim(1,5)
for a in ax[1]: a.set_xlabel(r'$r\,[R_*]$')
ax[0,0].legend(fontsize=8); fig.suptitle('M1 Star8 1D: solver comparison, final profiles')
fig.tight_layout(); fig.savefig(f'{D0}/profiles_final.png',dpi=140)
# zoom on the deep envelope / iron bump
for a in ax.flat: a.set_xlim(1,1.4); a.set_xscale('linear')
ax[1,1].set_ylim(-1,2.5); ax[0,1].set_ylim(-2.5,2.5); ax[0,0].set_ylim(1e-10,2)
fig.savefig(f'{D0}/profiles_final_zoom.png',dpi=140)
print('saved')
