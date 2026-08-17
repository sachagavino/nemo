#!/usr/bin/env bash
# ===========================================================================
# tests/jac_weight_check_debug.sh
#
# DIAGNOSTIC build of the weight test: rebuilds every module AND the driver with
# -g -fbacktrace -fcheck=all so a crash reports the exact array/index/file:line
# instead of an opaque address backtrace. Use this to pinpoint the macOS segfault.
# Same optional input-dir arg as jac_weight_check.sh (default: inputs/).
# Run from repo root:  ./tests/jac_weight_check_debug.sh   (or with a fixture dir).
# ===========================================================================
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
DBG="-O0 -g -fbacktrace -fcheck=all -ffree-line-length-none"

echo "[debug] rebuilding all modules with bounds checking (-fcheck=all)..."
make clean >/dev/null 2>&1
make debug >/tmp/jwcd_build.log 2>&1 || { echo "module debug build FAILED"; tail -30 /tmp/jwcd_build.log; exit 3; }

echo "[debug] compiling + linking driver with checks..."
gfortran $DBG -Jbuild -Ibuild -c tests/jac_weight_check.f90 -o build/jac_weight_check.o 2>>/tmp/jwcd_build.log \
  || { echo "driver compile FAILED"; tail -30 /tmp/jwcd_build.log; exit 3; }
OBJS="$(ls build/*.o build/dust/*.o 2>/dev/null | grep -vE 'build/main\.o|_check\.o' | tr '\n' ' ')"
gfortran -O0 -g -fbacktrace build/jac_weight_check.o $OBJS -o bin/jac_weight_check 2>>/tmp/jwcd_build.log \
  || { echo "driver link FAILED"; tail -30 /tmp/jwcd_build.log; exit 3; }

RUNDIR="$(mktemp -d)"
INPUTS_DIR="${1:-inputs}"
cp "$INPUTS_DIR"/*.in "$RUNDIR"/ 2>/dev/null || { echo "no $INPUTS_DIR/*.in"; exit 3; }
if [ -f "$RUNDIR/parameters.in" ]; then
  sed 's/^preliminary_test *=.*/preliminary_test = 0/' "$RUNDIR/parameters.in" > "$RUNDIR/parameters.in.tmp" \
    && mv "$RUNDIR/parameters.in.tmp" "$RUNDIR/parameters.in"
fi

echo "[debug] running (bounds-checked) in $RUNDIR ..."
echo "        -- the FIRST 'At line N of file src/...' below is the crash site --"
( cd "$RUNDIR" && "$ROOT/bin/jac_weight_check" ) 2>&1 | tail -25
echo "[debug] (rundir kept for inspection: $RUNDIR)"
