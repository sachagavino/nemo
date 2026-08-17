#!/usr/bin/env bash
# ===========================================================================
# tests/coag_ice_gate.sh
#
# Integration gate for ice transport. Runs the ice fixture with MONODISPERSE grains
# and ice in bin 1 (so the early-slope is a clean self-collision), constant kernel,
# and asserts:
#
#   Gate (a)  total ice  sum_k Y(J_k X) conserved to machine precision (ice is inert
#             here: coagulation only moves it between bins).
#   Dust-mass sum_k m_k n_k still conserved -- retained Rung-1 gate; the GRAIN_j
#             catalyst encoding must not perturb the grains (issue-1 failure mode).
#   Rate      ice early-slope for the monodisperse bin-1 start is the self-collision
#               dY(J01X)/dt|_0 = -K0 * nH * Y(GRAIN01)_0 * Y(J01X)_0
#             with the BARE K_11 (full, no 1/2). The 1/2 error would read ~2x here.
#
# Run from repo root; exit 0 = PASS.
# ===========================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "[coag_ice_gate] building (-O0)..."
make OPT=-O0 >/dev/null 2>&1

RUN="$(mktemp -d)"
trap 'rm -rf "$RUN"' EXIT
cp tests/fixture_ice_transport/*.in "$RUN"/
sed 's/^dust_ic = MRN/dust_ic = monodisperse/; s/^start_time = .*/start_time =  1.000E-02/; s/^stop_time = .*/stop_time =  1.000E+02/; s/^nb_outputs = .*/nb_outputs = 20/' \
    "$RUN"/parameters.in > "$RUN"/parameters.in.tmp && mv "$RUN"/parameters.in.tmp "$RUN"/parameters.in

echo "[coag_ice_gate] running in $RUN ..."
( cd "$RUN" && "$ROOT/bin/nmgc" run >run.log 2>&1 )

python3 - "$RUN" <<'PYEOF'
import sys, struct, re, collections
RUN = sys.argv[1]

# species index -> name
names = {}
for line in open(RUN+"/species.out"):
    for m in re.finditer(r"(\d+)\)\s+(\S+)", line):
        names[int(m.group(1))] = m.group(2)
nb_species = max(names)
ice_idx = sorted(i for i,nm in names.items() if re.match(r"^J\d\d", nm))       # all J_k X
j01     = sorted(i for i,nm in names.items() if re.match(r"^J01", nm))[0]       # a bin-1 ice
grain01 = [i for i,nm in names.items() if nm == "GRAIN01"][0]

P = {}
for l in open(RUN+"/parameters.in"):
    if '=' in l and not l.strip().startswith('!'):
        k,v = l.split('=',1); P[k.strip()] = v.split('!')[0].strip()
K0 = float(P['constant_kernel_k0']); nH = float(P['initial_gas_density'])
YR = 3.1556952e7

# abundances.out: per output -> record(time), record(temps+scalars), record(abundances)
data = open(RUN+"/abundances.out","rb").read(); pos = 0
def rec(d,p):
    (n,)=struct.unpack_from("<i",d,p); p+=4; pl=d[p:p+n]; p+=n
    (n2,)=struct.unpack_from("<i",d,p); p+=4
    assert n==n2; return pl,p
outs=[]
while pos < len(data):
    tp,pos=rec(data,pos); (t,)=struct.unpack("<d",tp)
    _,pos=rec(data,pos)                                   # skip temps/scalars via marker
    ap,pos=rec(data,pos); ab=struct.unpack("<%dd"%nb_species,ap)
    outs.append((t,ab))

# gate (a): total ice conservation
tot0 = sum(outs[0][1][i-1] for i in ice_idx)
ice_worst = max(abs(sum(ab[i-1] for i in ice_idx)-tot0)/tot0 for _,ab in outs)

# dust-mass conservation from dust_distribution.out
byt = collections.OrderedDict()
for p in (l.split() for l in open(RUN+"/dust_distribution.out") if not l.startswith("!") and l.strip()):
    if len(p) < 9: continue
    byt.setdefault(float(p[0]),0.0); byt[float(p[0])] += float(p[5])
m0 = byt[list(byt)[0]]
mass_worst = max(abs(v-m0)/m0 for v in byt.values())

# rate gate: ice early-slope (self-collision, bare K11)
t0,ab0 = outs[0]; t1,ab1 = outs[1]
Y1_0 = ab0[grain01-1]; Yice0 = ab0[j01-1]; Yice1 = ab1[j01-1]
slope_meas = (Yice1 - Yice0)/((t1 - t0))               # t already in s
slope_ana  = -K0 * nH * Y1_0 * Yice0
rate_relerr = abs(slope_meas - slope_ana)/abs(slope_ana)

print("=================================================================")
print(" Ice-transport gate (constant kernel, monodisperse bin-1 start)")
print("=================================================================")
print(" ice species  = %d (%s..)" % (len(ice_idx), names[ice_idx[0]]))
print(" gate (a) ice conservation rel-err = %.3e   (tol 1e-10)" % ice_worst)
print(" dust-mass conservation    rel-err = %.3e   (tol 1e-8)"  % mass_worst)
print(" ice early-slope           rel-err = %.3e   (tol 5e-2)"  % rate_relerr)
ok = (ice_worst < 1e-10) and (mass_worst < 1e-8) and (rate_relerr < 5e-2)
print(" RESULT: %s" % ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
PYEOF
rc=$?
echo "[coag_ice_gate] $([ $rc -eq 0 ] && echo PASS || echo FAIL)"
exit $rc
