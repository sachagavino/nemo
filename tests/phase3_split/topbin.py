"""Item 4: top-bin ice sweep-up -- CO-ice share & monolayers per bin over time (reference), and the
per-bin CO-ice error profile of a split run (which bins carry the headline error)."""
import sys, numpy as np, compare as C
ref, run = sys.argv[1], sys.argv[2]
n,T,Y = C.read(ref); ns,Ts,Ys = C.read(run); ng=20
rows=[l.split() for l in open(ref+'/dust_grid_table.in') if l.strip() and not l.startswith('!')]
a=np.array([float(r[0]) for r in rows]); sites=4*np.pi*a**2*8e14
jco=[n.index('J%02dCO'%k) for k in range(1,ng+1)]
gr=[(n.index('GRAIN%02d'%k),n.index('GRAIN%02d-'%k)) for k in range(1,ng+1)]
Jall=[[i for i,s in enumerate(n) if s.startswith('J%02d'%k)] for k in range(1,ng+1)]
print("t[yr]  CO-ice share bins18-20  | monolayers bins 1,4,8,12,16,20 | split rel.err of JkCO, bins 1..20 (x1e-3)")
for i,t in enumerate(T):
    if t<=0: continue
    y=Y[i]; ys=Ys[np.argmin(abs(Ts-t))]
    co=y[jco]; share=co[17:].sum()/co.sum()
    ngr=np.array([y[g0]+y[g1] for g0,g1 in gr]); ml=np.array([y[Jall[k]].sum() for k in range(ng)])/(ngr*sites)
    err=np.abs(ys[jco]-co)/co
    print(f"{t:6.0f}  {share:6.3f}   | "+" ".join(f"{ml[k]:6.1f}" for k in (0,3,7,11,15,19))+" | "+" ".join(f"{1e3*e:4.1f}" for e in err))
