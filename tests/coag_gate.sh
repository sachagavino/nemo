#!/usr/bin/env bash
#==============================================================================
# Rung 1 coagulation analytic gate (constant kernel).
#
# Runs the dust-only fixture with a MONODISPERSE initial condition (all grains in
# bin 1) and a constant coagulation kernel K0, sized so the run spans the knee of
# the closed-form Smoluchowski total-number decay, then asserts:
#
#   Gate (a)  mass  Sum_k m_k n_k conserved to machine precision.
#   Gate (b1) total number  N(t) = N0 / (1 + K0 * nH * N0 * (t-t0) / 2)
#             reproduced to within 5% across the whole trajectory.
#
# The 1/2 is the constant-kernel self-collision factor. Overflow (top-bin) pairs
# are not generated in Rung 1; the grid is sized so they stay unpopulated, which
# is exactly what keeps the boundary-free analytic solution valid. If a future
# change breaks mass conservation, the redistribution weights, or the kernel
# normalisation, this gate fails.
#==============================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "[coag_gate] building (-O0)..."
make OPT=-O0 >/dev/null 2>&1

RUN="$(mktemp -d)"
trap 'rm -rf "$RUN"' EXIT
cp tests/fixture_dust_coag/*.in "$RUN"/
# Monodisperse IC + a window that spans the analytic knee (~2e3 yr for this grid/K0).
sed -i 's/^dust_ic = MRN/dust_ic = monodisperse/; s/^stop_time = .*/stop_time =  2.000E+03/; s/^nb_outputs = .*/nb_outputs = 25/' "$RUN"/parameters.in

echo "[coag_gate] running in $RUN ..."
( cd "$RUN" && "$ROOT/bin/nmgc" run >run.log 2>&1 )

python3 - "$RUN" <<'PYEOF'
import sys, collections
RUN = sys.argv[1]
rows = [l.split() for l in open(RUN+"/dust_distribution.out")
        if not l.startswith("!") and l.strip()]
byt = collections.OrderedDict()
for p in rows:
    if len(p) < 9: continue
    t = float(p[0]); d = byt.setdefault(t, {'num':0.0, 'mass':0.0})
    d['num']  += float(p[4])   # n_k [/H]
    d['mass'] += float(p[5])   # m_k*n_k [g/H]

P = {}
for l in open(RUN+"/parameters.in"):
    if '=' in l and not l.strip().startswith('!'):
        k, v = l.split('=', 1); P[k.strip()] = v.split('!')[0].strip()
K0 = float(P['constant_kernel_k0']); nH = float(P['initial_gas_density'])
YR = 3.1556952e7

times = list(byt)
t0 = times[0]; N0 = byt[t0]['num']; m0 = byt[t0]['mass']
mass_worst = max(abs(byt[t]['mass'] - m0)/m0 for t in times)
fit_worst = 0.0
for t in times:
    Nana = N0 / (1.0 + 0.5*K0*nH*N0*(t - t0)*YR)
    fit_worst = max(fit_worst, abs(byt[t]['num'] - Nana)/Nana)

print("=================================================================")
print(" Rung 1 coagulation analytic gate (constant kernel)")
print("=================================================================")
print(" outputs                = %d over %.4g..%.4g yr" % (len(times), times[0], times[-1]))
print(" N/N0 at end            = %.4f" % (byt[times[-1]]['num']/N0))
print(" gate (a) mass rel-err  = %.3e   (tol 1e-10)" % mass_worst)
print(" gate (b1) decay rel-err= %.3e   (tol 5e-2)"  % fit_worst)
ok = (mass_worst < 1e-10) and (fit_worst < 5e-2)
print(" RESULT: %s" % ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
PYEOF
rc=$?
echo "[coag_gate] $([ $rc -eq 0 ] && echo PASS || echo FAIL)"
exit $rc
