#!/usr/bin/env bash
# ===========================================================================
# tests/coag_rung3_3a.sh   --   Rung 3 gate 3a: conservation with chemistry on
#
# On the coupled reduced_CHO run (coag ON, chemistry + charging ON), assert:
#   * dust-core mass  Sum_k m_k*[Y(GRAIN_k0)+Y(GRAIN_k-)]  conserved to machine
#     precision (SUMMED OVER BOTH CHARGE STATES -- a neutral-only sum drifts on the
#     charging that moves grains between states at fixed m_k);
#   * ice-transport operator net-zero per species (structural, all ice-transport
#     reactions conserve their ice base with w1+w2=1) -- since total ice is NOT
#     conserved with chemistry on, this is the operator-level check.
# Run from repo root; exit 0 = PASS.
# ===========================================================================
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
echo "[3a] building (-O0)..."; make OPT=-O0 >/tmp/3a_b.log 2>&1 || { echo build FAIL; tail -20 /tmp/3a_b.log; exit 3; }
gfortran -O0 -ffree-line-length-none -Jbuild -Ibuild -c tests/coag_ice_operator_check.f90 -o build/coag_ice_operator_check.o 2>>/tmp/3a_b.log \
  || { echo "operator-check compile FAIL"; tail -20 /tmp/3a_b.log; exit 3; }
OBJS="$(ls build/*.o build/dust/*.o 2>/dev/null | grep -vE 'build/main\.o|_check\.o' | tr '\n' ' ')"
gfortran -O0 build/coag_ice_operator_check.o $OBJS -o bin/coag_ice_operator_check 2>>/tmp/3a_b.log || { echo "operator-check link FAIL"; exit 3; }

R="$(mktemp -d)"; trap 'rm -rf "$R"' EXIT
cp networks/reduced_CHO/*.in "$R"/; cp inputs/element.in "$R"/
sed 's/^coagulation = .*/coagulation = 1/' tests/fixture_rung3/parameters.in > "$R"/parameters.in
echo "[3a] coupled run (coag on, chemistry+charging on)..."
( cd "$R" && "$ROOT/bin/nmgc" run >run.log 2>&1 ) || { echo "run FAIL"; tail -5 "$R"/run.log; exit 1; }

echo "[3a] operator net-zero (structural)..."
( cd "$R" && "$ROOT/bin/coag_ice_operator_check" >op.log 2>&1 )
grep -q "RESULT: PASS" "$R"/op.log && OPOK=1 || OPOK=0
grep -E "reactions checked|worst|violations" "$R"/op.log | sed 's/^/   /'

python3 - "$R" "$OPOK" <<'PY'
import sys, struct, re, math
R=sys.argv[1]; opok=int(sys.argv[2])
names={}
for line in open(R+"/species.out"):
    for m in re.finditer(r"(\d+)\)\s+(\S+)", line): names[int(m.group(1))]=m.group(2)
nb=max(names); idx={v:k for k,v in names.items()}
g0=sorted((n for n in names.values() if re.match(r"GRAIN\d+$",n)),key=lambda s:int(re.search(r"\d+",s).group()))
gm=sorted((n for n in names.values() if re.match(r"GRAIN\d+-$",n)),key=lambda s:int(re.search(r"\d+",s).group()))
rad={}
for l in open(R+"/dust_grid_active.out"):
    if l.startswith("!") or not l.strip(): continue
    p=l.split(); rad[int(float(p[4]))]=float(p[0])
m=[(4/3)*math.pi*rad[k]**3*3.0 for k in range(1,len(g0)+1)]
data=open(R+"/abundances.out","rb").read();pos=0
def rec(d,p):
    (n,)=struct.unpack_from("<i",d,p);p+=4;pl=d[p:p+n];p+=n;(n2,)=struct.unpack_from("<i",d,p);p+=4;return pl,p
outs=[]
while pos<len(data):
    tp,pos=rec(data,pos);_,pos=rec(data,pos);ap,pos=rec(data,pos);outs.append(struct.unpack("<%dd"%nb,ap))
def dmass(ab,both=True): return sum(m[k]*(ab[idx[n0]-1]+(ab[idx[nm]-1] if both else 0)) for k,(n0,nm) in enumerate(zip(g0,gm)))
m0=dmass(outs[0]); wb=max(abs(dmass(ab)-m0)/m0 for ab in outs)
wn=max(abs(dmass(ab,False)-dmass(outs[0],False))/dmass(outs[0],False) for ab in outs)
print("=================================================================")
print(" Rung 3 gate 3a: conservation with chemistry on")
print("=================================================================")
print(" dust-core mass (both charges) worst rel-err = %.2e   (tol 1e-10)"%wb)
print(" (neutral-only would read %.2e -- why both-charge sum is required)"%wn)
print(" ice-transport operator net-zero            = %s"%("PASS" if opok else "FAIL"))
ok = wb<1e-10 and opok
print(" RESULT: %s"%("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
PY
rc=$?; echo "[3a] $([ $rc -eq 0 ] && echo PASS || echo FAIL)"; exit $rc
