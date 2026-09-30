import sys, numpy as np
src, out, tw, tc = sys.argv[1], sys.argv[2], float(sys.argv[3]), float(sys.argv[4])
rows=[l.split() for l in open(src) if l.strip() and not l.lstrip().startswith('!')]
a=np.array([float(r[0]) for r in rows])
td=tw+(tc-tw)*np.log(a/a[0])/np.log(a[-1]/a[0])
with open(out,'w') as f:
    f.write(f"! parametric Td: linear in log a, {tw} K at a_min -> {tc} K at a_max\n")
    for r,t in zip(rows,td):
        f.write(f"{r[0]}  {r[1]}  {t:.16e}  {r[3]}  {r[4]}\n")
