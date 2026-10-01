"""Measured convergence orders, full window and t<=8 kyr (c2), for a set of runs vs a reference."""
import sys, numpy as np, compare as C
CUT=1e-12
def errs(ref, run, tmax):
    nr,Tr,Yr=C.read(ref); ns,Ts,Ys=C.read(run)
    jco=[nr.index('J%02dCO'%k) for k in range(1,21)]; co=nr.index('CO'); _,_,_,chem,_=C.classes(nr)
    e=dict(jco=0.,co=0.,all=0.)
    for i,t in enumerate(Ts):
        if t<=0 or t>tmax: continue
        j=np.argmin(abs(Tr-t)); yr,ys=Yr[j],Ys[i]
        for k in jco:
            if yr[k]>CUT: e['jco']=max(e['jco'],abs(ys[k]-yr[k])/yr[k])
        e['co']=max(e['co'],abs(ys[co]-yr[co])/yr[co])
        sel=[k for k in chem if yr[k]>CUT]; e['all']=max(e['all'],(np.abs(ys[sel]-yr[sel])/yr[sel]).max())
    return e
ref=sys.argv[1]; pre=sys.argv[2]; dts=[1000,500,250,100,50,25]
for tmax,lab in ((1e4,'full window'),(8000,'t<=8 kyr (c2)')):
    E=[errs(ref,f'{pre}{dt}',tmax) for dt in dts]
    print(f"--- {pre}*  {lab}")
    for key in ('jco','co','all'):
        v=[e[key] for e in E]; p=[np.log(v[i]/v[i+1])/np.log(dts[i]/dts[i+1]) for i in range(len(v)-1)]
        print(f"  {key:4s} err: "+" ".join(f"{x:.1e}" for x in v)+"   local order: "+" ".join(f"{x:4.2f}" for x in p))
