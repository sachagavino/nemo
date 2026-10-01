"""Phase III correctness gate as defined in the design-thread decision note (sec. 1)."""
import sys, re, numpy as np, compare as C
TOL=1e-4; CUT=1e-12
def gate(ref, run, verbose=True):
    nr,Tr,Yr=C.read(ref); ns,Ts,Ys=C.read(run); ng=20
    jco=[nr.index('J%02dCO'%k) for k in range(1,ng+1)]; co=nr.index('CO')
    g=[(nr.index('GRAIN%02d'%k),nr.index('GRAIN%02d-'%k)) for k in range(1,ng+1)]
    ng_,g0,gm,chem,_=C.classes(nr)
    e_jco=e_co=e_tot=e_gr=e_all=0.; worst=''
    for i,t in enumerate(Ts):
        if t<=0: continue
        j=np.argmin(abs(Tr-t)); yr,ys=Yr[j],Ys[i]
        for k in jco:
            if yr[k]>CUT: e_jco=max(e_jco,abs(ys[k]-yr[k])/yr[k])
        if yr[co]>CUT: e_co=max(e_co,abs(ys[co]-yr[co])/yr[co])
        e_tot=max(e_tot,abs(ys[jco].sum()-yr[jco].sum())/yr[jco].sum())
        for a,b in g: e_gr=max(e_gr,abs(ys[a]+ys[b]-yr[a]-yr[b])/(yr[a]+yr[b]))
        sel=[k for k in chem if yr[k]>CUT]; r=np.abs(ys[sel]-yr[sel])/yr[sel]
        if r.max()>e_all: e_all=r.max(); worst=f"{nr[sel[int(np.argmax(r))]]} @ {t:.0f} yr"
    m=C.mass_grid(ref); cr=C.conservation(ref,m); cs=C.conservation(run,m)
    res=dict(jco=e_jco,co=e_co,tot=e_tot,grain=e_gr)
    ok_i=max(e_jco,e_co,e_tot)<TOL; ok_ii=e_gr<TOL
    # (iii) 'conserved to the coupled tolerance': executor-proposed explicit bound 1e-10 relative
    # (100x below rtol 1e-8); the ratio to the coupled drift is reported alongside.
    ok_iii=all(cs[k]<1e-10 for k in ('charge','H','C','O','dust'))
    log=open(run+'/run.log').read()
    m_aud=re.search(r'chem block slots\s+=\s+\[(\d+),(\d+)\].*?dust block slots\s+=\s+\[(\d+),(\d+)\].*?total = (\d+)',log,re.S)
    ok_iv=bool(m_aud) and int(m_aud[1])==1 and int(m_aud[3])==int(m_aud[2])+1 and int(m_aud[4])==int(m_aud[5])
    if verbose:
        print(f"== GATE {run} vs {ref}")
        print(f" (i)  per-bin CO ice {e_jco:.2e}  gas CO {e_co:.2e}  total CO ice {e_tot:.2e}   -> {'PASS' if ok_i else 'FAIL'} (<{TOL:g})")
        print(f" (ii) grain bins (totals) {e_gr:.2e}  -> {'PASS' if ok_ii else 'FAIL'}")
        print(f" (iii) drift run/ref: "+"  ".join(f"{k} {cs[k]:.1e}/{cr[k]:.1e}" for k in ('charge','H','C','O','dust'))+f"  -> {'PASS' if ok_iii else 'FAIL'}")
        print(f" (iv) mask audit: chem [{m_aud[1]},{m_aud[2]}] dust [{m_aud[3]},{m_aud[4]}] of {m_aud[5]}  -> {'PASS' if ok_iv else 'FAIL'}")
        print(f" diagnostic (not gated): all-species max (gas+ice>1e-12) {e_all:.2e}  [{worst}]")
    return ok_i and ok_ii and ok_iii and ok_iv, res
if __name__=='__main__':
    for run in sys.argv[2:]: gate(sys.argv[1],run)
