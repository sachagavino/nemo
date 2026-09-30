"""Step-0 per-bin timescales from a coupled NEMO run (reads NEMO outputs; kernel,
weights and desorption re-evaluated with NEMO's own formulas/constants)."""
import sys, re, struct, numpy as np
KB=1.3806488e-16; PI=3.1415926535898; AMU=1.66053892e-24; YR=3.15576e7
E_CO=1150.; M_CO=28.; NS=8e14; ZETA=1.3e-17; FE=3e-14; AREF=1e-5; TPK=70.; DUR=1e-5
def read_run(d):
    names=dict((int(a),b) for a,b in re.findall(r'(\d+)\)\s+(\S+)',open(d+'/species.out').read()))
    nb=max(names); ng=sum(1 for n in names.values() if re.fullmatch(r'GRAIN\d\d',n))
    raw=open(d+'/abundances.out','rb').read(); p=0; out=[]
    def rec(p):
        (n,)=struct.unpack_from('<i',raw,p); pl=raw[p+4:p+4+n]; return pl,p+8+n
    while p<len(raw):
        t,p=rec(p); phys,p=rec(p); ab,p=rec(p)
        t=struct.unpack('<d',t)[0]; ph=struct.unpack('<%dd'%(len(phys)//8),phys)
        out.append((t/YR, ph[0], np.array(ph[1:1+ng]), ph[1+ng], np.array(struct.unpack('<%dd'%nb,ab))))
    return names, ng, out
def grid(d):
    rows=[l.split() for l in open(d+'/dust_grid_active.out') if l.strip() and not l.lstrip().startswith('!')]
    a=np.array([float(r[0]) for r in rows]); m=4/3*PI*a**3*3.0
    return a,m
def weights(m):
    N=len(m); W=np.zeros((N,N,N))  # W[i,j,k]: fraction of (i,j) product into bin k; overflow -> all zero
    for i in range(N):
        for j in range(N):
            mp=m[i]+m[j]
            if mp>m[-1]*(1+1e-9): continue
            k=np.searchsorted(m,mp*(1-1e-12),side='right')-1
            if k>=N-1: W[i,j,N-1]=1; continue
            e=(m[k+1]-mp)/(m[k+1]-m[k]); W[i,j,k]+=e; W[i,j,k+1]+=1-e
    return W
def analyse(d, t_want):
    names,ng,out=read_run(d); a,m=grid(d); W=weights(m); idx={v:k-1 for k,v in names.items()}
    t=np.array([o[0] for o in out]); it=int(np.argmin(abs(np.log(t/t_want)))); tt,Tg,Td,nH,Y=out[it]
    n=np.array([Y[idx['GRAIN%02d'%k]]+Y[idx['GRAIN%02d-'%k]] for k in range(1,ng+1)])
    mu=np.outer(m,m)/np.add.outer(m,m); K=np.add.outer(a,a)**2*np.sqrt(8*PI*KB*Tg/mu)
    ok=W.sum(2)>0                                   # non-overflow pairs (only these are generated)
    Kg=np.where(ok,K,0.)
    gross=(Kg*n[None,:]).sum(1)*nH                  # per-grain (and per-ice-molecule) removal rate [s^-1]
    stay=np.array([[W[k,j,k] for j in range(ng)] for k in range(ng)])
    ice_net=(Kg*(1-stay)*n[None,:]).sum(1)*nH       # ice actually leaving bin k
    # grain net d n_k/dt (Smoluchowski with NEMO weights; self pair counted with 1/2)
    R=0.5*Kg*np.outer(n,n)*nH                       # unordered-pair-symmetric rate/2 per ordered pair
    gain=np.einsum('ij,ijk->k',R,W); loss=(Kg*np.outer(n,n)).sum(1)*nH
    dndt=gain-loss
    nu=np.sqrt(2*KB/PI/PI/AMU*NS*E_CO/M_CO)
    kth=nu*np.exp(-E_CO/Td); kcr=(ZETA/1.3e-17)*nu*FE*(a/AREF)**2*DUR*np.exp(-E_CO/TPK)*(a<=1e-4)
    vth=np.sqrt(8*KB*Tg/(PI*M_CO*AMU)); kacc=PI*a**2*vth*n*nH   # per gas-CO molecule, S=1
    jco=np.array([Y[idx['J%02dCO'%k]] for k in range(1,ng+1)]); gco=Y[idx['CO']]
    res=dict(t=tt,nH=nH,Tg=Tg,Td=Td,a=a,n=n,tau_coag_gross=1/gross/YR,tau_coag_net=n/np.maximum(abs(dndt),1e-300)/YR,
             tau_ice_net=1/np.maximum(ice_net,1e-300)/YR,tau_des=1/(kth+kcr)/YR,tau_acc_bin=1/kacc/YR,
             tau_fo=1/kacc.sum()/YR,jco=jco,gco=gco,ice_frac_bin=jco/(jco.sum()+gco))
    return res
def report(d,t_want,dt_head=1e3):
    r=analyse(d,t_want)
    print(f"\n== {d}  t={r['t']:.3g} yr  nH={r['nH']:.3g}  Tgas={r['Tg']:.3g} K  tau_fo(CO)={r['tau_fo']:.3g} yr"
          f"  CO gas={r['gco']:.2e}  sum JkCO={r['jco'].sum():.2e}")
    print(" bin  a[nm]   Td    n_k/H     tau_coag_gross tau_coag_net tau_ice_net  tau_des(CO)  tau_acc,k  JkCO/H   dt/tau_ice_net")
    for k in range(len(r['a'])):
        print(f" {k+1:3d} {r['a'][k]*1e7:6.1f} {r['Td'][k]:5.2f}  {r['n'][k]:.2e}   {r['tau_coag_gross'][k]:.2e}     {r['tau_coag_net'][k]:.2e}   "
              f"{r['tau_ice_net'][k]:.2e}   {r['tau_des'][k]:.2e}   {r['tau_acc_bin'][k]:.2e}  {r['jco'][k]:.2e}  {dt_head/r['tau_ice_net'][k]:.2e}")
    return r
if __name__=='__main__':
    for tw in sys.argv[2:]: report(sys.argv[1], float(tw))
