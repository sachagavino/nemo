#!/usr/bin/env bash
# ===========================================================================
# tests/coag_rung4_sparsity.sh  --  Rung 4 acceptance: symbolic pattern completeness
#                                   + FD-complete oracle correctness at the 20-bin grid
#
# Rung 4 is plumbing: complete the symbolic Jacobian pattern (the live-divisor block),
# wire MOSS=0, and keep the FD-complete path (mf=022) as the Rung 5 verification oracle.
# Production ships mf=121 (analytic); the omitted live-divisor entry is BOUNDED (floored
# SUMLAY), so it is harmless at Rung 4. This gate checks:
#
#   (a) structural completeness  -- the live-divisor entries are in the symbolic pattern
#                                   (the code's own static assert_live_divisor_in_pattern).
#   (b) correctness at scale     -- production (121) and the FD-complete oracle (022)
#                                   integrate the SAME RHS to the SAME solution at 20 bins:
#                                   same answer for meaningful species (<1% for ab>1e-8),
#                                   and neither run reports an elemental-conservation drift.
#   4-IV tractability            -- report both wall timings (022 is the slow upper bound).
#
# The FD-complete oracle is reached with NEMO_ORACLE_MF=22 (never a production default).
# Run from repo root; exit 0 = PASS. NOTE: the 022 oracle run is ~1-2 min at 20 bins.
#
# (Gate (c) "load-bearing convergence improvement" is retracted: the live divisor is
#  bounded by the SUMLAY floor, so it is not the stiff -1/Y^2 reciprocal -- that is the
#  dynamic surface-rate GTODN, which is Rung 5. See docs/PhaseII_rung4_sparsity_scale.md
#  section 8. 3b/3c physics gates run separately on production mf=121.)
# ===========================================================================
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
echo "[r4] building (-O0)..."; make OPT=-O0 >/tmp/r4_b.log 2>&1 || { echo build FAIL; tail -20 /tmp/r4_b.log; exit 3; }

WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT

# 20-bin fiducial grid (a_max = 0.5 um) on the Rung 3 fixture; coag on, symbolic.
mkcfg() { # $1=dir  $2..=extra sed exprs
  local d="$1"; shift; mkdir -p "$d"
  cp networks/reduced_CHO/*.in "$d"/; cp inputs/element.in "$d"/
  sed -e 's/^a_max = .*/a_max =  5.000E-05/' \
      -e 's/^coagulation = .*/coagulation = 1/' \
      -e 's/^sparsity = .*/sparsity = symbolic/' \
      tests/fixture_rung3/parameters.in > "$d"/parameters.in
}
runwall() { # $1=dir  $2=env  -> echoes wall seconds, leaves run.log
  local d="$1" e="$2" t0 t1
  t0=$(date +%s); ( cd "$d" && env $e "$ROOT/bin/nmgc" run >run.log 2>&1 ); local rc=$?
  t1=$(date +%s); echo $((t1-t0))
  return $rc
}

echo "[r4] production run (mf=121) at 20 bins..."
mkcfg "$WORK/prod"
w121=$(runwall "$WORK/prod" "DUMMY=1") || { echo "prod run FAIL"; tail -5 "$WORK/prod/run.log"; exit 1; }

echo "[r4] FD-complete oracle run (NEMO_ORACLE_MF=22) at 20 bins... (~1-2 min)"
mkcfg "$WORK/oracle"
w022=$(runwall "$WORK/oracle" "NEMO_ORACLE_MF=22") || { echo "oracle run FAIL"; tail -5 "$WORK/oracle/run.log"; exit 1; }

FAIL=0

echo "=================================================================="
# ---- (a) structural completeness: the code's own static assert ----
if grep -q "(gate a) live-divisor structural completeness: PASS" "$WORK/oracle/run.log"; then
  echo " (a) structural completeness : PASS  ($(grep -o '[0-9]\+ entries' "$WORK/oracle/run.log" | head -1))"
else
  echo " (a) structural completeness : FAIL  (assert_live_divisor_in_pattern did not pass)"; FAIL=1
fi

# ---- confirm the two runs used the intended methods ----
m121=$(grep -m1 "method flag mf" "$WORK/prod/run.log"   | grep -o '[0-9]\+$')
m022=$(grep -m1 "method flag mf" "$WORK/oracle/run.log" | grep -o '[0-9]\+$')
echo " methods                     : production mf=${m121:-?} , oracle mf=${m022:-?}"
[ "${m121:-0}" = "121" ] || { echo "   production is not 121!"; FAIL=1; }
[ "${m022:-0}" = "22"  ] || { echo "   oracle is not 22!"; FAIL=1; }

# ---- conservation: neither run may report an elemental drift ----
if grep -qiE "not conserved" "$WORK/prod/run.log" "$WORK/oracle/run.log"; then
  echo " conservation (3a)           : FAIL  (elemental drift reported)"; FAIL=1
else
  echo " conservation (3a)           : PASS  (no elemental drift in either run)"
fi

# ---- (b) same answer: production 121 vs oracle 022 ----
python3 - "$WORK/prod/abundances.out" "$WORK/oracle/abundances.out" <<'PY'
import sys, struct
def load(p):
    d=open(p,'rb').read();pos=0;o=[]
    def r(d,p):
        (n,)=struct.unpack_from("<i",d,p);p+=4;pl=d[p:p+n];p+=n;p+=4;return pl,p
    while pos<len(d):
        _,pos=r(d,pos);_,pos=r(d,pos);a,pos=r(d,pos);o.append(struct.unpack("<%dd"%(len(a)//8),a))
    return o
a=load(sys.argv[1]); b=load(sys.argv[2])
mx=0.0; nd=0
for ra,rb in zip(a,b):
    for x,y in zip(ra,rb):
        s=max(abs(x),abs(y))
        if s>1e-8:
            rr=abs(x-y)/s
            if rr>0.01: nd+=1
            if rr>mx: mx=rr
ok = (mx < 0.01)
print(" (b) same answer (ab>1e-8)   : %s  (max rel diff %.2e vs 121 ; species >1%% = %d)"
      % ("PASS" if ok else "FAIL", mx, nd))
sys.exit(0 if ok else 7)
PY
[ $? -eq 0 ] || FAIL=1

# ---- 4-IV tractability ----
echo " 4-IV timing (20 bins)       : production mf=121 ~${w121}s ; FD-complete oracle mf=022 ~${w022}s (upper bound)"
echo "=================================================================="

if [ "$FAIL" -eq 0 ]; then echo "[r4] PASS"; exit 0; else echo "[r4] FAIL"; exit 1; fi
