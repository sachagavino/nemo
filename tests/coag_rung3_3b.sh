#!/usr/bin/env bash
# ===========================================================================
# tests/coag_rung3_3b.sh  --  Rung 3 gate 3b: slow-coag reduction (primary acceptance)
#
# As coagulation is scaled into the slow regime (t_coag = 1/(K N) >> t_run), the coupled
# run must reduce to the coag-OFF reference. We run the reference (coag off) and TWO slow
# K0 values 10x apart, and require:
#   * the coupled run matches the reference to a small rtol, and
#   * rtol TRACKS ~linearly with K0 (10x smaller K0 -> ~10x tighter match) -- one passing
#     point is not the test; the linear collapse is,
#   * the divisor FLOOR is never triggered in the slow runs (a K0-independent floor hit is
#     the prime suspect for a 3b plateau).
# Identical stop_time / cadence across all three runs (from fixture_rung3/parameters.in).
# Run from repo root; exit 0 = PASS.
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
echo "[3b] reference (coag off) + two slow K0 (5e-18, 5e-19)..."
run ref 0 "1.000E-14"
run slowA 1 "5.000E-18"
run slowB 1 "5.000E-19"

python3 - "$WORK" <<'PY'
import sys, struct, re
W=sys.argv[1]
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
    # ~rtol jitter floor with no K0 dependence (verified: J04H's slow-vs-ref ratio is
    # non-monotonic 0.2-10.4, abundances identical to ~4 sig figs). Species above 1e-10
    # carry a physical reduction signal >> rtol jitter and collapse to exactly 10.00 at
    # every output (e.g. GRAIN01). The floor separates signal from rtol jitter; it is
    # NOT masking a real deviation (that would scale with K0, which J04H's does not).
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
print(" K0 = 5e-18 : worst rel-diff vs reference = %.3e"%rA)
print(" K0 = 5e-19 : worst rel-diff vs reference = %.3e"%rB)
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
