#!/usr/bin/env bash
# ===========================================================================
# tests/coag_ice_scale.sh
#
# Scaling check for ice transport on the fiducial grid (16 bins, several ice
# species). Confirms:
#   * gate (a) total ice sum_k Y(J_k X) conserved to machine precision for EACH ice
#     species (transport is per-base and must not mix species);
#   * the term count stays tractable -- ice transport is O(N^2 N_ice); for the reduced
#     network (N_ice ~ tens) it must be ~1e4, not ~1e5. Printed by the injector.
#
# Run from repo root; exit 0 = PASS.
# ===========================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "[coag_ice_scale] building (-O0)..."
make OPT=-O0 >/dev/null 2>&1

RUN="$(mktemp -d)"; trap 'rm -rf "$RUN"' EXIT
cp tests/fixture_ice_fiducial/*.in "$RUN"/
echo "[coag_ice_scale] running 16-bin fiducial (3 ice species) in $RUN ..."
( cd "$RUN" && "$ROOT/bin/nmgc" run > run.log 2>&1 )
grep "reactions:" "$RUN"/run.log | sed 's/^/   /' || true

python3 - "$RUN" <<'PYEOF'
import sys, struct, re, collections
RUN = sys.argv[1]
names = {}
for line in open(RUN+"/species.out"):
    for m in re.finditer(r"(\d+)\)\s+(\S+)", line): names[int(m.group(1))] = m.group(2)
nb = max(names)
bases = collections.defaultdict(list)
for i,nm in names.items():
    mm = re.match(r"^J\d\d(.+)$", nm)
    if mm: bases[mm.group(1)].append(i)

data = open(RUN+"/abundances.out","rb").read(); pos = 0
def rec(d,p):
    (n,)=struct.unpack_from("<i",d,p); p+=4; pl=d[p:p+n]; p+=n
    (n2,)=struct.unpack_from("<i",d,p); p+=4; return pl,p
outs=[]
while pos < len(data):
    _,pos=rec(data,pos); _,pos=rec(data,pos); ap,pos=rec(data,pos)
    outs.append(struct.unpack("<%dd"%nb,ap))

# term count from the injector line
rc_line = [l for l in open(RUN+"/run.log") if "reactions:" in l]
n_total = None
if rc_line:
    m = re.search(r"=\s*(\d+)\s*total", rc_line[0]); n_total = int(m.group(1)) if m else None

print("=================================================================")
print(" Ice-transport scaling gate (16-bin fiducial)")
print("=================================================================")
worst = 0.0
for base in sorted(bases):
    tot0 = sum(outs[0][i-1] for i in bases[base])
    w = max(abs(sum(ab[i-1] for i in bases[base])-tot0)/tot0 for ab in outs)
    worst = max(worst, w)
    print(" J%-4s (%2d bins): total=%.4e  conservation rel-err=%.2e" % (base, len(bases[base]), tot0, w))
print(" ice species = %d   coag+ice reactions = %s" % (len(bases), n_total))
ok_cons  = worst < 1e-10
ok_count = (n_total is not None) and (n_total < 100000)   # ~1e4 regime, not 1e5
print(" gate (a) worst rel-err = %.2e   (tol 1e-10)  -> %s" % (worst, "OK" if ok_cons else "FAIL"))
print(" term count %s < 1e5 (tractable)               -> %s" % (n_total, "OK" if ok_count else "FAIL"))
print(" RESULT: %s" % ("PASS" if (ok_cons and ok_count) else "FAIL"))
sys.exit(0 if (ok_cons and ok_count) else 1)
PYEOF
rc=$?
echo "[coag_ice_scale] $([ $rc -eq 0 ] && echo PASS || echo FAIL)"
exit $rc
