#!/usr/bin/env bash
# ===========================================================================
# tests/coag_brownian_gate.sh
#
# Integration gate for the free-molecular Brownian kernel. The Brownian kernel has
# no clean closed-form total-number decay (unlike the constant kernel), so instead
# of a full-trajectory analytic fit we assert:
#
#   Gate (a)  mass  Sum_k m_k n_k conserved to machine precision.
#   Early slope  for a MONODISPERSE start (all grains in bin 1), the only initial
#                channel is bin1+bin1, so the exact initial total-number rate is
#                   dN/dt|_0 = -1/2 * K11_phys * nH * N0^2,
#                   K11_phys = (2 a1)^2 * sqrt( 8 pi k_B T / (m1/2) ).
#                Measured from the first two outputs; must match within 5%.
#
# This checks the kernel is correctly WIRED into the RHS (magnitude + the 1/2
# self-pair factor + the nH density factor), complementing coag_kernel_check.sh,
# which checks the kernel VALUE. Run from repo root; exit 0 = PASS.
# ===========================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "[coag_brownian_gate] building (-O0)..."
make OPT=-O0 >/dev/null 2>&1

RUN="$(mktemp -d)"
trap 'rm -rf "$RUN"' EXIT
cp tests/fixture_dust_coag/*.in "$RUN"/
# Brownian kernel, monodisperse IC, early first outputs so the initial slope is clean.
sed -i 's/^coagulation_kernel = constant/coagulation_kernel = brownian/; s/^dust_ic = MRN/dust_ic = monodisperse/; s/^start_time = .*/start_time =  1.000E-01/; s/^stop_time = .*/stop_time =  5.000E+03/; s/^nb_outputs = .*/nb_outputs = 25/' "$RUN"/parameters.in

echo "[coag_brownian_gate] running in $RUN ..."
( cd "$RUN" && "$ROOT/bin/nmgc" run >run.log 2>&1 )

python3 - "$RUN" <<'PYEOF'
import sys, math, collections
RUN = sys.argv[1]
rows = [l.split() for l in open(RUN+"/dust_distribution.out")
        if not l.startswith("!") and l.strip()]
byt = collections.OrderedDict()
a1 = m1 = None
for p in rows:
    if len(p) < 9: continue
    t = float(p[0]); d = byt.setdefault(t, {'num':0.0, 'mass':0.0})
    d['num']  += float(p[4]); d['mass'] += float(p[5])
    if int(p[1]) == 1 and a1 is None:
        a1 = float(p[2]); m1 = float(p[3])

P = {}
for l in open(RUN+"/parameters.in"):
    if '=' in l and not l.strip().startswith('!'):
        k, v = l.split('=', 1); P[k.strip()] = v.split('!')[0].strip()
T  = float(P['initial_gas_temperature']); nH = float(P['initial_gas_density'])
kB = 1.3806488e-16; PI = math.pi; YR = 3.1556952e7

times = list(byt)
t0, t1 = times[0], times[1]
N0, N1 = byt[t0]['num'], byt[t1]['num']
m0 = byt[t0]['mass']
mass_worst = max(abs(byt[t]['mass'] - m0)/m0 for t in times)

# exact initial slope for a monodisperse (bin-1-only) start
K11phys = (2*a1)**2 * math.sqrt(8*PI*kB*T/(m1/2.0))
slope_analytic = -0.5 * K11phys * nH * N0**2
slope_measured = (N1 - N0)/((t1 - t0)*YR)
slope_relerr = abs(slope_measured - slope_analytic)/abs(slope_analytic)

print("=================================================================")
print(" Brownian coagulation gate (free molecular)")
print("=================================================================")
print(" gas T = %.2f K   nH = %.3e cm^-3   N/N0 end = %.4f"
      % (T, nH, byt[times[-1]]['num']/N0))
print(" gate (a) mass rel-err     = %.3e   (tol 1e-10)" % mass_worst)
print(" early dN/dt|0 analytic    = %.6e /H/s" % slope_analytic)
print(" early dN/dt|0 measured     = %.6e /H/s" % slope_measured)
print(" early-slope rel-err        = %.3e   (tol 5e-2)" % slope_relerr)
ok = (mass_worst < 1e-10) and (slope_relerr < 5e-2)
print(" RESULT: %s" % ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
PYEOF
rc=$?
echo "[coag_brownian_gate] $([ $rc -eq 0 ] && echo PASS || echo FAIL)"
exit $rc
