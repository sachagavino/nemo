"""Gate/money-plot comparison of a split run against a coupled reference (NEMO outputs)."""
import sys, re, struct, numpy as np
import os
NET=os.path.join(os.environ.get('NEMO_ROOT', os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..')),'networks','reduced_CHO')+'/'
YR=3.15576e7; CUT=1e-12
def read(d):
    names=dict((int(a),b) for a,b in re.findall(r'(\d+)\)\s+(\S+)',open(d+'/species.out').read()))
    nb=max(names); raw=open(d+'/abundances.out','rb').read(); p=0; T=[]; Y=[]
    def rec(p):
        (n,)=struct.unpack_from('<i',raw,p); return raw[p+4:p+4+n],p+8+n
    while p<len(raw):
        t,p=rec(p); _,p=rec(p); a,p=rec(p)
        T.append(struct.unpack('<d',t)[0]/YR); Y.append(struct.unpack('<%dd'%nb,a))
    return [names[i] for i in range(1,nb+1)], np.array(T), np.array(Y)
def comp():
    c={}
    for f in ('gas_species.in','grain_species.in'):
        for l in open(NET+f):
            s=l.split()
            if s and not s[0].startswith('!'): c[s[0]]=np.array([int(x) for x in s[1:]])
    return c
COMP=comp()
def species_vec(name):
    if name=='e-': v=np.zeros(14,int); v[0]=-1; return v
    m=re.fullmatch(r'GRAIN(\d\d)(-?)',name)
    if m: v=np.zeros(14,int); v[0]=-1 if m.group(2) else 0; return v
    m=re.fullmatch(r'J(\d\d)(.+)',name)
    if m: return COMP['J'+m.group(2)]
    return COMP[name]
def classes(names):
    ng=sum(1 for n in names if re.fullmatch(r'GRAIN\d\d',n))
    g0=[names.index('GRAIN%02d'%k) for k in range(1,ng+1)]; gm=[names.index('GRAIN%02d-'%k) for k in range(1,ng+1)]
    gi=set(g0+gm); chem=[i for i,n in enumerate(names) if i not in gi]
    jco=[names.index('J%02dCO'%k) for k in range(1,ng+1)]
    return ng,g0,gm,chem,jco
def invariants(names,Y,mass):
    ng,g0,gm,chem,jco=classes(names)
    V=np.array([species_vec(n) for n in names])        # col0 charge, col1 H, col2 He, col3 C, col5 O
    tot=Y@V
    dust=(Y[:,g0]+Y[:,gm])@mass
    return dict(charge=tot[:,0], H=tot[:,1], C=tot[:,3], O=tot[:,5], dust=dust)
def gate(ref, run, mass, verbose=True):
    nr,Tr,Yr=read(ref); ns,Ts,Ys=read(run); assert nr==ns
    ng,g0,gm,chem,jco=classes(nr); co=nr.index('CO')
    rows=[]
    for i,t in enumerate(Ts):
        j=np.argmin(abs(Tr-t))
        if abs(Tr[j]-t)>1e-6*max(t,1) or t<=0: continue
        yr,ys=Yr[j],Ys[i]
        sel=[k for k in chem if yr[k]>CUT]
        rel_chem=np.abs(ys[sel]-yr[sel])/yr[sel]
        gtr=yr[g0]+yr[gm]; gts=ys[g0]+ys[gm]
        rel_gr=np.abs(gts-gtr)/gtr
        rel_jco=np.abs(ys[jco]-yr[jco])/yr[jco]
        L1=np.abs(ys[sel]-yr[sel]).sum()/yr[sel].sum(); L2=np.sqrt(((ys[sel]-yr[sel])**2).sum()/(yr[sel]**2).sum())
        worst=nr[sel[int(np.argmax(rel_chem))]]
        rows.append(dict(t=t,max_chem=rel_chem.max(),worst=worst,max_grain=rel_gr.max(),max_jco=rel_jco.max(),
                         jco_bin=int(np.argmax(rel_jco))+1,co_gas=abs(ys[co]-yr[co])/yr[co],
                         jco_tot=abs(ys[jco].sum()-yr[jco].sum())/yr[jco].sum(),L1=L1,L2=L2,nsel=len(sel)))
    if verbose:
        print(f"{'t[yr]':>8} {'max_rel(gas+ice>1e-12)':>22} {'worst':>10} {'max_rel(grain bins)':>19} {'max_rel(JkCO)':>13} {'bin':>3} {'CO gas':>9} {'sum JkCO':>9} {'L1':>9} {'L2':>9}")
        for r in rows:
            print(f"{r['t']:8.0f} {r['max_chem']:22.3e} {r['worst']:>10} {r['max_grain']:19.3e} {r['max_jco']:13.3e} {r['jco_bin']:3d} {r['co_gas']:9.2e} {r['jco_tot']:9.2e} {r['L1']:9.2e} {r['L2']:9.2e}")
    return rows
def conservation(d, mass):
    n,T,Y=read(d); inv=invariants(n,Y,mass)
    out={}
    for k,v in inv.items():
        scale = inv['H'][0] if k=='charge' else abs(v[0])   # charge: relative to H nuclei
        out[k]=np.max(np.abs(v-v[0]))/scale
    return out
def mass_grid(d):
    rows=[l.split() for l in open(d+'/dust_grid_table.in') if l.strip() and not l.startswith('!')]
    a=np.array([float(r[0]) for r in rows]); return 4/3*np.pi*a**3*3.0
if __name__=='__main__':
    ref,run=sys.argv[1],sys.argv[2]; m=mass_grid(ref)
    gate(ref,run,m)
    for d in (ref,run): print(d, 'max drift:', {k:f"{v:.1e}" for k,v in conservation(d,m).items()})
