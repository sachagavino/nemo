import csv, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt; 
from matplotlib.ticker import  FixedLocator, FixedFormatter, NullFormatter;
from matplotlib.lines import Line2D
def load(f): return list(csv.DictReader(open(f),delimiter='\t'))
R=load('sweep_all.tsv'); W=load('wp_all.tsv')
fl=lambda r,k: float(r[k])
# ---- accuracy panel: diagnostic error vs dt, all four (variant, order) series
fig,ax=plt.subplots(1,2,figsize=(11,4.3))
sty={('A','CDC'):('C0','o','-'),('B','CDC'):('C1','s','-'),('A','DCD'):('C0','o','--'),('B','DCD'):('C1','s','--')}
for (v,o),(c,m,ls) in sty.items():
    rr=sorted([r for r in R if r['variant']==v and r['order']==o],key=lambda r:fl(r,'dt_yr'))
    if not rr: continue
    dt=[fl(r,'dt_yr') for r in rr]
    ax[0].loglog(dt,[fl(r,'max_jco') for r in rr],color=c,marker=m,ls=ls,label=f'{v}, {o}: per-bin CO ice')
    ax[1].loglog(dt,[fl(r,'co_gas') for r in rr],color=c,marker=m,ls=ls,label=f'{v}, {o}: gas CO')
for a in ax:
    x=np.array([10,1000]); a.loglog(x,1e-3*(x/100)**1,'k:',lw=1,label='slope 1'); a.loglog(x,1e-4*(x/100)**2,'k-.',lw=1,label='slope 2')
    a.axvspan(100,1000,color='0.9',zorder=0); a.axvline(1000,color='0.5',lw=.8)
    a.axhline(1e-4,color='r',lw=1.); a.set_xlabel(r'$\Delta t_{\rm split}$ [yr]'); #a.grid(alpha=.3,which='both')
ax[0].set_ylabel('max relative error vs coupled (rtol 1e-8)'); ax[0].legend(fontsize=7, handlelength=3.0); ax[1].legend(fontsize=7, handlelength=3.0)
ax[0].set_title('per-bin CO ice (max over bins)'); ax[1].set_title('gas-phase CO')
plt.tight_layout(); plt.savefig('fig_accuracy.pdf',dpi=130)
#fig.suptitle(r'$T_{\rm gas}=21.23$ K, $n_{\rm H}=2.2\times10^8$, $T_d$ 25$\to$13 K',fontsize=9);
# ---- work-precision: error (vs rtol 1e-10 reference) vs cost
fig,ax=plt.subplots(1,2,figsize=(11,4.3))
for key,lab in (('nfe','RHS evaluations'),('nlu','sparse LU decompositions')):
    a=ax[0] if key=='nfe' else ax[1]
    cp=sorted([r for r in W if r['variant']=='coupled'],key=lambda r:fl(r,key))
    a.loglog([fl(r,key) for r in cp],[max(fl(r,'max_jco'),1e-12) for r in cp],'k-o',label='coupled (rtol 1e-4..1e-8)')
    for r in cp: a.annotate(f"{fl(r,'rtol'):.0e}",(fl(r,key),max(fl(r,'max_jco'),1e-12)),fontsize=8)
    for v,c in (('A','C0'),('B','C1')):
        for rt,m in ((1e-4,'^'),(1e-6,'v'),(1e-8,'o')):
            ss=sorted([r for r in W if r['variant']==v and r['order']=='CDC' and abs(fl(r,'rtol')/rt-1)<1e-6],key=lambda r:-fl(r,'dt_yr'))
            if ss: a.loglog([fl(r,key) for r in ss],[fl(r,'max_jco') for r in ss],color=c,marker=m,ls='-',lw=.8,label=f'split {v} CDC, rtol {rt:.0e}')
    a.set_xlabel(lab); 
ax[1].xaxis.set_major_locator(
    FixedLocator([2e3, 1e4, 2e4, 4e4])
)

ax[1].xaxis.set_major_formatter(
    FixedFormatter([
        r'$2\times10^3$',
        r'$10^4$',
        r'$2\times10^4$',
        r'$4\times10^4$'
    ])
)

ax[1].xaxis.set_minor_formatter(NullFormatter())
ax[0].set_ylabel('per-bin CO ice: max rel. error vs coupled (rtol 1e-10)')
ax[0].legend(fontsize=8); plt.tight_layout(); plt.savefig('fig_workprecision.pdf',dpi=130)

# fig.suptitle(r'$T_{\rm gas}=21.23$ K, $n_{\rm H}=2.2\times10^8$, $T_d$ 25$\to$13 K',fontsize=9)