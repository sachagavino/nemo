import sys, numpy as np
sys.argv=[sys.argv[0]]; from analyze import analyse
PI=3.1415926535898
def summ(d,t,dt=1e3):
    r=analyse(d,t); a=r['a']; n=r['n']; j=r['jco']
    theta=j/(n*4*PI*a**2*8e14)                     # CO monolayers per grain
    snow = theta/theta.max()
    kw = np.where(snow>=0.5)[0].min()               # warmest bin still >=half max coverage
    fdes = np.minimum(1, dt/r['tau_des']); ftr = np.minimum(1, dt/r['tau_ice_net'])
    P_gas = (ftr*j*fdes).sum()/max(r['gco'],1e-30)
    P_ice = (ftr*j*fdes).sum()/j.sum()
    band = slice(max(kw-1,0), min(kw+2,len(a)))
    print(f"{d:9s} t={r['t']:8.2e}  COgas/COtot={r['gco']/(r['gco']+j.sum()):.2e}  "
          f"theta(b1)/max={snow[0]:.1e} snow_bin={kw+1:2d}(Td={r['Td'][kw]:.1f})  "
          f"tau_ice_net@snow={r['tau_ice_net'][kw]:.1e} tau_des@snow={r['tau_des'][kw]:.1e} tau_fo={r['tau_fo']:.1e}  "
          f"proxy_gas={P_gas:.1e} proxy_ice={P_ice:.1e}")
runs=[('n2e8_21',1e4),('n2e8_23',1e4),('n2e8_25',1e4),('n2e8_27',1e4),
      ('n2e7_21',1e5),('n2e7_23',1e5),('n2e7_25',1e5),('n2e7_27',1e5),
      ('n2e6_21',1e6),('n2e6_23',1e6),('n2e6_25',1e6)]
for d,T in runs:
    for t in (T/10, T): summ(d,t)
