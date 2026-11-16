import numpy as np,re,glob,os,sys
d=sys.argv[1] if len(sys.argv)>1 else '.'
order={'static':0,'para':1,'perp':2,'rot':3,'xstatic':4,'xpara':5}
fs=sorted([g for g in glob.glob(f'{d}/output/*.log') if os.path.basename(g).split('_')[0] in order],key=lambda f:(order.get(os.path.basename(f).split('_')[0],9),int(os.path.basename(f)[:-4].split('_')[1])))
print(f"{'run':11s}  final: {'eE':>8s} {'eF':>8s} {'eFx':>8s} {'eFy':>8s} {'eFz':>8s} {'ev':>7s} | eF at t=0.1/0.5/1/2/3")
for f in fs:
    L=open(f).read().splitlines(); h=L[0].split('|')[0].split(); c={n:i for i,n in enumerate(h)}
    a=np.array([[float(re.sub(r"(\d)([+-]\d{3})$",r"\1E\2",x)) for x in l.split('|')[0].split()] for l in L[1:]]); t=a[:,1]
    g=lambda n:a[-1,c[n]]
    print(f"{os.path.basename(f)[:-4]:11s}  t={t[-1]:.2f} {g('err_E'):8.2e} {g('err_F'):8.2e} {g('err_Fx'):8.2e} {g('err_Fy'):8.2e} {g('err_Fz'):8.2e} {g('err_v'):7.1e} | "+' '.join(f"{a[np.argmin(abs(t-x)),c['err_F']]:.2e}" for x in [0.1,0.5,1,2,3]))
