#!/usr/bin/env bash
# ===========================================================================
# tests/coag_ice_check.sh
#
# Hand-check (b): builds tests/coag_ice_check.f90, inits the ice-transport fixture,
# and verifies that the ice-transport reactions reuse the grain product weights
# (same eps) for a genuine eps=0.5 split pair (1,2) and the eps=0 clean merge (1,1).
# Run from repo root: ./tests/coag_ice_check.sh  (exit 0 = PASS).
# ===========================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "[coag_ice_check] building modules (-O0)..."
make OPT=-O0 >/tmp/cic_build.log 2>&1 || { echo "module build FAILED"; tail -30 /tmp/cic_build.log; exit 3; }

echo "[coag_ice_check] compiling + linking driver..."
gfortran -O0 -ffree-line-length-none -Jbuild -Ibuild -c tests/coag_ice_check.f90 -o build/coag_ice_check.o 2>>/tmp/cic_build.log \
  || { echo "driver compile FAILED"; tail -30 /tmp/cic_build.log; exit 3; }
OBJS="$(ls build/*.o build/dust/*.o 2>/dev/null | grep -vE 'build/main\.o|_check\.o' | tr '\n' ' ')"
gfortran -O0 build/coag_ice_check.o $OBJS -o bin/coag_ice_check 2>>/tmp/cic_build.log \
  || { echo "driver link FAILED"; tail -30 /tmp/cic_build.log; exit 3; }

RUNDIR="$(mktemp -d)"
trap 'rm -rf "$RUNDIR"' EXIT
cp tests/fixture_ice_transport/*.in "$RUNDIR"/

echo "[coag_ice_check] running in $RUNDIR ..."
( cd "$RUNDIR" && "$ROOT/bin/coag_ice_check" ) | tee /tmp/cic_run.log
grep -q "RESULT: PASS" /tmp/cic_run.log
echo "[coag_ice_check] PASS"
