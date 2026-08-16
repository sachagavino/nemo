#!/usr/bin/env bash
# ===========================================================================
# tests/coag_kernel_check.sh
#
# Builds and runs the Brownian-kernel absolute check (tests/coag_kernel_check.f90):
# links the driver against the compiled module objects (everything in build/ except
# main.o), initialises the dust-only fixture with coagulation_kernel=brownian, and
# checks coag_kernel(i,j) for two pairs against hand-computed values. Run from repo
# root:  ./tests/coag_kernel_check.sh   (exit 0 = PASS).
# ===========================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "[coag_kernel_check] building modules (-O0)..."
make OPT=-O0 >/tmp/ckc_build.log 2>&1 || { echo "module build FAILED"; tail -30 /tmp/ckc_build.log; exit 3; }

echo "[coag_kernel_check] compiling + linking driver..."
gfortran -O0 -ffree-line-length-none -Jbuild -Ibuild -c tests/coag_kernel_check.f90 -o build/coag_kernel_check.o 2>>/tmp/ckc_build.log \
  || { echo "driver compile FAILED"; tail -30 /tmp/ckc_build.log; exit 3; }
OBJS="$(ls build/*.o build/dust/*.o 2>/dev/null | grep -vE 'build/main\.o|_check\.o' | tr '\n' ' ')"
gfortran -O0 build/coag_kernel_check.o $OBJS -o bin/coag_kernel_check 2>>/tmp/ckc_build.log \
  || { echo "driver link FAILED"; tail -30 /tmp/ckc_build.log; exit 3; }

RUNDIR="$(mktemp -d)"
trap 'rm -rf "$RUNDIR"' EXIT
cp tests/fixture_dust_coag/*.in "$RUNDIR"/
sed -i 's/^coagulation_kernel = constant/coagulation_kernel = brownian/' "$RUNDIR"/parameters.in

echo "[coag_kernel_check] running in $RUNDIR ..."
( cd "$RUNDIR" && "$ROOT/bin/coag_kernel_check" ) | tee /tmp/ckc_run.log
grep -q "RESULT: PASS" /tmp/ckc_run.log
echo "[coag_kernel_check] PASS"
