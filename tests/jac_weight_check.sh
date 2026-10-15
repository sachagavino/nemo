#!/usr/bin/env bash
# ===========================================================================
# tests/jac_weight_check.sh
#
# Builds and runs the Jacobian/RHS weight-plumbing test (tests/jac_weight_check.f90).
# It links the driver against the compiled module objects (everything in build/
# except main.o), initialises the full network from inputs/, and checks that the
# per-product weight reaches every get_jacobian deposit consistently with the RHS.
#
# Requires: gfortran, a populated inputs/ directory (the default network is fine;
# the test is network-agnostic -- it only needs a valid init). Run from repo root:
#     ./tests/jac_weight_check.sh
# Exit 0 = PASS.
# ===========================================================================
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "[jac_weight_check] building modules (-O0)..."
make OPT=-O0 >/tmp/jwc_build.log 2>&1 || { echo "module build FAILED"; tail -30 /tmp/jwc_build.log; exit 3; }

echo "[jac_weight_check] compiling + linking driver..."
gfortran -O0 -ffree-line-length-none -Jbuild -Ibuild -c tests/jac_weight_check.f90 -o build/jac_weight_check.o 2>>/tmp/jwc_build.log \
  || { echo "driver compile FAILED"; tail -30 /tmp/jwc_build.log; exit 3; }
OBJS="$(ls build/*.o build/dust/*.o 2>/dev/null | grep -vE 'build/main\.o|_check\.o' | tr '\n' ' ')"
gfortran -O0 build/jac_weight_check.o $OBJS -o bin/jac_weight_check 2>>/tmp/jwc_build.log \
  || { echo "driver link FAILED"; tail -30 /tmp/jwc_build.log; exit 3; }

RUNDIR="$(mktemp -d)"
trap 'rm -rf "$RUNDIR"' EXIT
# Optional first arg: input directory to stage (default: the full network in inputs/).
# Pass a coagulation fixture (e.g. tests/fixture_ice_transport) to exercise the
# grain- and ice-transport product-weight plumbing into get_jacobian.
INPUTS_DIR="${1:-inputs}"
cp "$INPUTS_DIR"/*.in "$RUNDIR"/ 2>/dev/null || { echo "no $INPUTS_DIR/*.in to seed the run"; exit 3; }
# The weight test exercises the RHS/Jacobian, not the network sanity checks, so skip
# preliminary_tests: it is unnecessary here, slow on the full network, and its
# full-network pass segfaults on some platforms (macOS). Portable in-place edit.
if [ -f "$RUNDIR/parameters.in" ]; then
  sed 's/^preliminary_test *=.*/preliminary_test = 0/' "$RUNDIR/parameters.in" > "$RUNDIR/parameters.in.tmp" \
    && mv "$RUNDIR/parameters.in.tmp" "$RUNDIR/parameters.in"
fi

echo "[jac_weight_check] running in $RUNDIR (inputs: $INPUTS_DIR) ..."
( cd "$RUNDIR" && "$ROOT/bin/jac_weight_check" ) | tee /tmp/jwc_run.log
grep -q "RESULT: PASS" /tmp/jwc_run.log
echo "[jac_weight_check] PASS"