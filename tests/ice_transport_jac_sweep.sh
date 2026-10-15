#!/usr/bin/env bash
# Rung 5b DIAGNOSIS (read-only): exhaustive per-reaction-analytic sweep of the type-50 Jacobian.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
make OPT=-O0 >/tmp/sw_b.log 2>&1 || { echo build FAIL; tail -20 /tmp/sw_b.log; exit 3; }
gfortran -O0 -ffree-line-length-none -Jbuild -Ibuild -c tests/ice_transport_jac_sweep.f90 -o build/ice_transport_jac_sweep.o 2>>/tmp/sw_b.log \
  || { echo "sweep compile FAIL"; tail -40 /tmp/sw_b.log; exit 3; }
OBJS="$(ls build/*.o build/dust/*.o 2>/dev/null | grep -vE 'build/main\.o|_check\.o|_diag\.o|_sweep\.o' | tr '\n' ' ')"
gfortran -O0 build/ice_transport_jac_sweep.o $OBJS -o bin/ice_transport_jac_sweep 2>>/tmp/sw_b.log \
  || { echo "sweep link FAIL"; tail -40 /tmp/sw_b.log; exit 3; }
R="$(mktemp -d)"; trap 'rm -rf "$R"' EXIT
cp networks/reduced_CHO/*.in "$R"/; cp inputs/element.in "$R"/
sed -e 's/^coagulation = .*/coagulation = 1/' -e 's/^sparsity = .*/sparsity = symbolic/' \
    tests/fixture_rung3/parameters.in > "$R"/parameters.in
( cd "$R" && "$ROOT/bin/ice_transport_jac_sweep" )
