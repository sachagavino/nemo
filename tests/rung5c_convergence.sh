#!/usr/bin/env bash
# Rung 5c: 021-vs-121 (and 022 for mechanism) convergence/step-stat comparison in a
# bin-depleting coagulation regime. Read-only wrt src. Usage: rung5c_convergence.sh K0 STOP NB
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
K0="${1:-1.000E-10}"; STOP="${2:-1.000E+07}"; NB="${3:-30}"
make OPT=-O0 >/tmp/5c_b.log 2>&1 || { echo build FAIL; tail -20 /tmp/5c_b.log; exit 3; }
R="$(mktemp -d)"; trap 'rm -rf "$R"' EXIT
cp networks/reduced_CHO/*.in "$R"/; cp inputs/element.in "$R"/
sed -e 's/^coagulation = .*/coagulation = 1/' \
    -e "s/^constant_kernel_k0 = .*/constant_kernel_k0 = ${K0}/" \
    -e "s/^stop_time = .*/stop_time = ${STOP}/" \
    -e "s/^nb_outputs = .*/nb_outputs = ${NB}/" \
    tests/fixture_rung3/parameters.in | sed 's/^sparsity = .*/sparsity = symbolic/' > "$R"/parameters.in
echo "regime: K0=${K0}  stop_time=${STOP} yr  nb_outputs=${NB}"
printf '%-6s %10s %10s %8s %8s %10s %14s\n' mf NST NFE NJE fails wall_s bin1_final
for MF in 121 21 22; do
  ( cd "$R" && t0=$(date +%s.%N); NEMO_ORACLE_MF=$MF "$ROOT/bin/nmgc" run > run_$MF.log 2>&1; t1=$(date +%s.%N); echo "$t1 $t0" > .t )
  nst=$(grep -E 'NST' "$R"/run_$MF.log | grep -oE '[0-9]+$'); nfe=$(grep -E 'NFE' "$R"/run_$MF.log | grep -oE '[0-9]+$')
  nje=$(grep -E 'NJE' "$R"/run_$MF.log | grep -oE '[0-9]+$'); fails=$(grep -E 'restarts' "$R"/run_$MF.log | grep -oE '[0-9]+$')
  wall=$(awk '{printf "%.3f", $1-$2}' "$R"/.t)
  b1=$(awk '$2==1{v=$5} END{printf "%.4e", v}' "$R"/dust_distribution.out)
  cp "$R"/abundances.out "$R"/ab_$MF.out
  printf '%-6s %10s %10s %8s %8s %10s %14s\n' "$MF" "$nst" "$nfe" "$nje" "$fails" "$wall" "$b1"
done
# cross-mf final-state agreement (max rel diff on final abundances row)
echo "--- final-state agreement (max rel diff, last abundances row) ---"
for PAIR in "21 121" "22 21"; do set -- $PAIR
  md=$(paste <(tail -1 "$R"/ab_$1.out | tr -s ' ' '\n') <(tail -1 "$R"/ab_$2.out | tr -s ' ' '\n') \
      | awk 'NF==2 && $1+0==$1 {a=$1;b=$2; d=(a-b); m=(a<0?-a:a); n=(b<0?-b:b); base=(m>n?m:n); if(base>1e-30){r=(d<0?-d:d)/base; if(r>mx)mx=r}} END{printf "%.3e", mx}')
  echo "  mf=$1 vs mf=$2 : max rel diff = $md"
done
