import csv,sys
keys=sys.argv[2].split(',') if len(sys.argv)>2 else ['run','dt_yr','max_chem','worst','max_grain','max_jco','jco_t1e3','co_gas','co_t1e3','jco_tot','L1','nfe','nje','nlu','cold','fail','cpu']
r=list(csv.DictReader(open(sys.argv[1]),delimiter='\t'))
def f(k,v):
    try: x=float(v); return f'{x:.0f}' if k in('nfe','nje','nlu','cold','fail','dt_yr','nst','macro') else f'{x:.2e}'
    except: return v
w={k:max(len(k),*(len(f(k,x[k])) for x in r)) for k in keys}
print('  '.join(k.rjust(w[k]) for k in keys))
for x in r: print('  '.join(f(k,x[k]).rjust(w[k]) for k in keys))
