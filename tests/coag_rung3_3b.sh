#!/usr/bin/env bash
# ===========================================================================
# tests/coag_rung3_3b.sh  --  Rung 3 gate 3b: slow-coag reduction (primary acceptance)
#
# As coagulation is scaled into the slow regime (t_coag = 1/(K0 N) >> t_run), the coupled
# run must reduce to the coag-OFF reference. We run the reference (coag off) and TWO slow
# K0 values 10x apart, and require:
#   * the coupled run matches the reference to a small rtol, and
#   * rtol TRACKS ~linearly with K0 (10x smaller K0 -> ~10x tighter match) -- one passing
#     point is not the test; the linear collapse is,
#   * the divisor FLOOR is never triggered in the slow runs (a K0-independent floor hit is
#     the prime suspect for a 3b plateau).
# Identical stop_time / cadence across all three runs (from fixture_rung3/parameters.in).
# Run from repo root; exit 0 = PASS.
#
# GRID-AWARE K0 (why the slow points are computed, not hardcoded)
# --------------------------------------------------------------
# "Slow" is a statement about t_run/t_coag, not about a magic K0. Since
#   t_coag = 1/(K0 * n_grain),   n_grain = (sum of grain abundances/H) * n_gas  [cm^-3],
# the same physical slowness is  t_run/t_coag = t_run * K0 * n_grain = TARGET, so the K0
# that realises a given TARGET scales with the grid:  K0 = TARGET / (t_run * n_grain).
# n_grain depends on the bin count (at fixed dust mass, more/smaller bins => more grains),
# so a hardcoded K0 that is correctly "slow" on the 3-5 bin fixture silently falls BELOW
# the solver noise floor (~rtol) on the 20-bin production grid -- the linearity ratio then
# reads FAIL for a reason unrelated to the physics. We therefore MEASURE n_grain from the
# reference run and set K0 = TARGET/(t_run*n_grain). TARGET_A is chosen ~100x above the
# rtol=1e-4 noise floor and deep in the linear regime (t_run/t_coag << 1); TARGET_B =
# TARGET_A/10. This keeps 3b honest at any resolution (the convergence study runs
# non-fiducial grids), and keys the slow points off the grid's own grain number rather
# than two magic constants that only hold on one grid.
# ===========================================================================
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
echo "[3b] building (-O0)..."; make OPT=-O0 >/tmp/3b_b.log 2>&1 || { echo build FAIL; tail -20 /tmp/3b_b.log; exit 3; }

WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
run() {  # $1=label  $2=coagulation  $3=K0
  local d="$WORK/$1"; mkdir -p "$d"
  cp networks/reduced_CHO/*.in "$d"/; cp inputs/element.in "$d"/
  sed -e "s/^coagulation = .*/coagulation = $2/" -e "s/^constant_kernel_k0 = .*/constant_kernel_k0 = $3/" \
      tests/fixture_rung3/parameters.in > "$d"/parameters.in
  ( cd "$d" && "$ROOT/bin/nmgc" run >run.log 2>&1 ) || { echo "run $1 FAIL"; tail -5 "$d"/run.log; exit 1; }
}

# --- reference first (coag off; K0 irrelevant) ---
echo "[3b] reference (coag off)..."
run ref 0 "1.000E-14"

# --- derive grid-aware slow K0 from the reference's own grain number density ---
TARGET_A="1.0e-2"   # t_run/t_coag for the slow-A point: ~100x above rtol noise, linear regime
KVALS=$(python3 - "$WORK/ref" "tests/fixture_rung3/parameters.in" "$TARGET_A" <<'PY'
import sys, struct, re
refdir, params, target_a = sys.argv[1], sys.argv[2], float(sys.argv[3])
p={}
for line in open(params):
    m=re.match(r'\s*([A-Za-z0-9_]+)\s*=\s*([0-9.eE+-]+)', line)
    if m: p[m.group(1)]=m.group(2)
n_gas=float(p['initial_gas_density']); stop_yr=float(p['stop_time'])
YEAR=3.15576e7; t_run=stop_yr*YEAR
# grain number density from the reference at t0: sum of GRAIN abundances/H * n_gas
names={}
for line in open(refdir+"/species.out"):
    for m in re.finditer(r"(\d+)\)\s+(\S+)", line): names[int(m.group(1))]=m.group(2)
nb=max(names); gidx=[i-1 for i,nm in names.items() if nm.startswith("GRAIN")]
d=open(refdir+"/abundances.out","rb").read(); pos=0
def rec(x,pp):
    (n,)=struct.unpack_from("<i",x,pp);pp+=4;pl=x[pp:pp+n];pp+=n;pp+=4;return pl,pp
_,pos=rec(d,pos);_,pos=rec(d,pos);ap,pos=rec(d,pos)
row0=struct.unpack("<%dd"%nb, ap)
g_sum=sum(row0[i] for i in gidx)
n_grain=g_sum*n_gas
K0A=target_a/(t_run*n_grain); K0B=K0A/10.0
print("%.6e %.6e %.6e %.6e"%(K0A, K0B, n_grain, t_run))
PY
)
read K0A K0B NGRAIN TRUN <<< "$KVALS"
echo "[3b] grid-aware slow K0 from n_grain=${NGRAIN} cm^-3, t_run=${TRUN} s:"
echo "       K0_A=${K0A}  (t_run/t_coag=${TARGET_A})   K0_B=${K0B} (=K0_A/10)"
run slowA 1 "$K0A"
run slowB 1 "$K0B"

python3 - "$WORK" "$K0A" "$K0B" <<'PY'
import sys, struct, re
W=sys.argv[1]; K0A=float(sys.argv[2]); K0B=float(sys.argv[3])
def load(d):
    names={}
    for line in open(d+"/species.out"):
        for m in re.finditer(r"(\d+)\)\s+(\S+)", line): names[int(m.group(1))]=m.group(2)
    nb=max(names)
    data=open(d+"/abundances.out","rb").read();pos=0;outs=[]
    def rec(x,p):
        (n,)=struct.unpack_from("<i",x,p);p+=4;pl=x[p:p+n];p+=n;(n2,)=struct.unpack_from("<i",x,p);p+=4;return pl,p
    while pos<len(data):
        tp,pos=rec(data,pos);(t,)=struct.unpack("<d",tp);_,pos=rec(data,pos);ap,pos=rec(data,pos)
        outs.append((t,struct.unpack("<%dd"%nb,ap)))
    return names,outs
def worst_reldiff(a,b,floor=1e-10):
    # Reduction is measured on species above an abundance floor. Rationale: atol=1e-99
    # (negligible) and satol=max(atol,1e-16*init)~0 for ice, so trace species are
    # controlled PURELY at rtol=1e-4 relative -- their run-to-run rel-diff sits at the
    # ~rtol jitter floor with no K0 dependence. Species above 1e-10 carry a physical
    # reduction signal >> rtol jitter and collapse ~linearly in K0. The floor separates
    # signal from rtol jitter; it is NOT masking a real deviation (that would scale with
    # K0). The grid-aware K0 above keeps the >1e-10 signal above this floor at any grid.
    w=0.0
    for (ta, va),(tb,vb) in zip(a,b):
        for x,y in zip(va,vb):
            s=max(abs(x),abs(y))
            if s>floor: w=max(w, abs(x-y)/s)
    return w
def grain_min_total(d):  # min over run of per-bin (Y0+Y-), for floor check
    names,outs=load(d); idx={v:k for k,v in names.items()}
    import re as _re
    g0=[n for n in names.values() if _re.match(r"GRAIN\d+$",n)]; gm=[n for n in names.values() if _re.match(r"GRAIN\d+-$",n)]
    g0=sorted(g0,key=lambda s:int(_re.search(r"\d+",s).group())); gm=sorted(gm,key=lambda s:int(_re.search(r"\d+",s).group()))
    return min(min(ab[idx[a]-1]+ab[idx[b]-1] for a,b in zip(g0,gm)) for _,ab in outs)

_,ref = load(W+"/ref")
_,sA  = load(W+"/slowA")
_,sB  = load(W+"/slowB")
rA = worst_reldiff(ref, sA); rB = worst_reldiff(ref, sB)
floor_thr = 2.990e-39/2.0e8
gA = grain_min_total(W+"/slowA"); gB = grain_min_total(W+"/slowB")
print("=================================================================")
print(" Rung 3 gate 3b: slow-coag reduction to the coag-off reference")
print("=================================================================")
print(" K0_A = %.3e : worst rel-diff vs reference = %.3e"%(K0A,rA))
print(" K0_B = %.3e : worst rel-diff vs reference = %.3e"%(K0B,rB))
print(" linearity  : rA/rB = %.2f  (expect ~10 as K0 drops 10x)"%(rA/rB if rB>0 else float('inf')))
print(" floor      : min(Y0+Y-) slowA=%.2e slowB=%.2e  vs floor %.2e  (>> => never triggered)"%(gA,gB,floor_thr))
lin = 5.0 <= (rA/rB if rB>0 else 0) <= 20.0        # ~10x, allow 2x slack
small = rA < 5e-2
floor_ok = min(gA,gB) > floor_thr*1e3
ok = lin and small and floor_ok
print(" RESULT: %s%s"%("PASS" if ok else "FAIL",
      "" if ok else "  [%s%s%s]"%("small " if not small else "", "linearity " if not lin else "", "floor" if not floor_ok else "")))
sys.exit(0 if ok else 1)
PY
rc=$?; echo "[3b] $([ $rc -eq 0 ] && echo PASS || echo FAIL)"; exit $rc
