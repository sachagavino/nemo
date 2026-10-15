#!/usr/bin/env bash
# ===========================================================================
# tests/coag_rung3_3c.sh  --  Rung 3 gate 3c: nominal-coag difference (signed)
#
# At the nominal kernel (t_coag ~ t_run), the coupled run must differ from the coag-off
# reference IN THE EXPECTED DIRECTION: dust-core mass migrates to larger bins (grains
# grow) and surface ice migrates with it (larger-bin J_k X populations rise). Measured
# as the mass-weighted and ice-weighted MEAN BIN INDEX at the final time -- both must be
# larger in the coupled run than in the reference. A directional check, not "differs".
# (Coverage-switch (silicate<->ASW) structure from the live SUMLAY is expected, not a bug.)
# Identical stop_time/cadence to the reference (fixture_rung3/parameters.in).
# Run from repo root; exit 0 = PASS.
# ===========================================================================
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
echo "[3c] building (-O0)..."; make OPT=-O0 >/tmp/3c_b.log 2>&1 || { echo build FAIL; tail -20 /tmp/3c_b.log; exit 3; }
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
run(){ d="$WORK/$1"; mkdir -p "$d"; cp networks/reduced_CHO/*.in "$d"/; cp inputs/element.in "$d"/;
  sed -e "s/^coagulation = .*/coagulation = $2/" -e "s/^constant_kernel_k0 = .*/constant_kernel_k0 = $3/" \
      tests/fixture_rung3/parameters.in > "$d"/parameters.in
  ( cd "$d" && "$ROOT/bin/nmgc" run >run.log 2>&1 ) || { echo "run $1 FAIL"; tail -5 "$d"/run.log; exit 1; }; }
echo "[3c] reference (coag off) + nominal coag (K0=1e-14)..."
run ref 0 "1.000E-14"
run nom 1 "1.000E-14"

python3 - "$WORK" <<'PY'
import sys, struct, re, math
W=sys.argv[1]
def load(d):
    names={}
    for line in open(d+"/species.out"):
        for m in re.finditer(r"(\d+)\)\s+(\S+)",line): names[int(m.group(1))]=m.group(2)
    nb=max(names); data=open(d+"/abundances.out","rb").read();pos=0;outs=[]
    def rec(x,p):
        (n,)=struct.unpack_from("<i",x,p);p+=4;pl=x[p:p+n];p+=n;struct.unpack_from("<i",x,p);p+=4;return pl,p
    while pos<len(data):
        tp,pos=rec(data,pos);_,pos=rec(data,pos);ap,pos=rec(data,pos);outs.append(struct.unpack("<%dd"%nb,ap))
    rad={}
    for l in open(d+"/dust_grid_active.out"):
        if l.startswith("!") or not l.strip(): continue
        p=l.split(); rad[int(float(p[4]))]=float(p[0])
    return names,outs,rad
def mean_bins(d):
    names,outs,rad=load(d); idx={v:k for k,v in names.items()}
    g0=sorted((n for n in names.values() if re.match(r"GRAIN\d+$",n)),key=lambda s:int(re.search(r"\d+",s).group()))
    gm=sorted((n for n in names.values() if re.match(r"GRAIN\d+-$",n)),key=lambda s:int(re.search(r"\d+",s).group()))
    ice=[n for n in names.values() if re.match(r"J\d\d",n)]
    m=[(4/3)*math.pi*rad[k]**3*3.0 for k in range(1,len(g0)+1)]
    ab=outs[-1]  # final time
    # mass-weighted mean bin index
    mass=[m[k]*(ab[idx[a]-1]+ab[idx[b]-1]) for k,(a,b) in enumerate(zip(g0,gm))]
    mbin=sum((k+1)*mass[k] for k in range(len(mass)))/sum(mass)
    # ice-weighted mean bin index (total ice, all species)
    icebin_num=0.0; icebin_den=0.0
    for nmk in ice:
        k=int(re.match(r"J(\d\d)",nmk).group(1)); icebin_num+=k*ab[idx[nmk]-1]; icebin_den+=ab[idx[nmk]-1]
    ibin=icebin_num/icebin_den if icebin_den>0 else 0.0
    return mbin, ibin
mR,iR=mean_bins(W+"/ref"); mN,iN=mean_bins(W+"/nom")
print("=================================================================")
print(" Rung 3 gate 3c: nominal-coag directional migration (final time)")
print("=================================================================")
print(" mass-weighted mean bin:  reference %.4f -> coupled %.4f  (%+.2e)"%(mR,mN,mN-mR))
print(" ice-weighted  mean bin:  reference %.4f -> coupled %.4f  (%+.2e)"%(iR,iN,iN-iR))
mass_up = mN > mR + 1e-6
ice_up  = iN > iR + 1e-6
print(" dust-core mass migrated to larger bins : %s"%("YES" if mass_up else "NO"))
print(" surface ice migrated to larger bins    : %s"%("YES" if ice_up else "NO"))
ok = mass_up and ice_up
print(" RESULT: %s"%("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
PY
rc=$?; echo "[3c] $([ $rc -eq 0 ] && echo PASS || echo FAIL)"; exit $rc
