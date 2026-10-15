#!/usr/bin/env bash
# ===========================================================================
# tests/coag_rung3_td_check.sh  --  P3: dust-temperature wiring guard (kept in harness)
#
# Rung 3 is the first load-bearing use of Td(a). Confirm the derived-grid path fills
# grain_temp(:) from the a^(-1/6) prescription, by checking dust_grid_active.out's Td
# column against  reference_dust_temperature * (a/reference_grain_radius)^(-1/6).
# A guard so a future grid/parameter change can't silently break dust temperature.
# Run from repo root; exit 0 = PASS.
# ===========================================================================
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
make OPT=-O0 >/tmp/td_b.log 2>&1 || { echo build FAIL; tail -10 /tmp/td_b.log; exit 3; }
R="$(mktemp -d)"; trap 'rm -rf "$R"' EXIT
cp networks/reduced_CHO/*.in "$R"/; cp inputs/element.in "$R"/
cp tests/fixture_rung3/parameters.in "$R"/parameters.in     # derived grid, size_scaled Td
( cd "$R" && "$ROOT/bin/nmgc" run >run.log 2>&1 ) || { echo "run FAIL"; tail -5 "$R"/run.log; exit 1; }
python3 - "$R" <<'PY'
import sys, re
R=sys.argv[1]
P={}
for l in open(R+"/parameters.in"):
    if '=' in l and not l.strip().startswith('!'):
        k,v=l.split('=',1); P[k.strip()]=v.split('!')[0].strip()
Tref=float(P.get('reference_dust_temperature','10'))
aref=float(P.get('reference_grain_radius','1.0e-5'))
worst=0.0
print(" bin   a[cm]        Td(dumped)     td_of_a        rel.err")
for l in open(R+"/dust_grid_active.out"):
    if l.startswith("!") or not l.strip(): continue
    p=l.split(); a=float(p[0]); Td=float(p[2])
    ana=Tref*(a/aref)**(-1.0/6.0); e=abs(Td-ana)/ana; worst=max(worst,e)
    print("  %s   %.4e   %.6f     %.6f     %.1e"%(p[4],a,Td,ana,e))
print(" P3 worst rel-err = %.2e"%worst)
import sys as _s; _s.exit(0 if worst<1e-6 else 1)
PY
rc=$?; echo "[td_check] $([ $rc -eq 0 ] && echo PASS || echo FAIL)"; exit $rc
