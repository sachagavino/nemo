"""Collect split-vs-reference error metrics and solver cost per run into a TSV."""
import sys, re, glob, os, numpy as np
import compare as C
def stats(d):
    s={}
    for l in open(d+'/run.log'):
        for key,lab in (('NST','nst'),('NFE','nfe'),('NJE','nje'),('NLU','nlu'),('cold starts','cold'),
                        ('macro-steps','macro'),('restarts','fail'),('CPU time','cpu')):
            if key in l and '=' in l: s[lab]=float(l.split('=')[-1])
    return s
def metrics(ref,run):
    rows=C.gate(ref,run,None,verbose=False)
    g=lambda k: max(r[k] for r in rows)
    at=lambda k,t: [r[k] for r in rows if abs(r['t']-t)<1e-6][0]
    return dict(max_chem=g('max_chem'), max_grain=g('max_grain'), max_jco=g('max_jco'), co_gas=g('co_gas'),
                jco_tot=g('jco_tot'), L1=g('L1'), L2=g('L2'), jco_t1e3=at('max_jco',1000), co_t1e3=at('co_gas',1000),
                worst=max(rows,key=lambda r:r['max_chem'])['worst'])
if __name__=='__main__':
    ref=sys.argv[1]; out=sys.argv[2]; dirs=sys.argv[3:]
    cols=['run','variant','order','dt_yr','rtol','max_chem','worst','max_grain','max_jco','jco_t1e3','co_gas','co_t1e3','jco_tot','L1','L2',
          'nst','nfe','nje','nlu','cold','macro','fail','cpu']
    with open(out,'w') as f:
        f.write('\t'.join(cols)+'\n')
        for d in dirs:
            if not os.path.exists(d+'/STATUS') or open(d+'/STATUS').read().strip()!='ok':
                if d.startswith('ref') and os.path.exists(d+'/abundances.out'): pass
                else: print('skip',d); continue
            p=open(d+'/parameters.in').read()
            g=lambda k,dflt='': (re.search(r'^%s\s*=\s*(\S+)'%k,p,re.M) or [None,dflt])[1]
            m=metrics(ref,d); s=stats(d)
            row=dict(run=d,variant=g('split_variant','coupled'),order=g('split_order','-'),dt_yr=float(g('split_dt','0')),
                     rtol=float(g('relative_tolerance')),**m,**s)
            f.write('\t'.join(str(row.get(c,'')) if not isinstance(row.get(c),float) else f"{row[c]:.4e}" for c in cols)+'\n')
    os.system(f"column -t -s $'\\t' {out} | cut -c1-230")
