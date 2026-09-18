#!/usr/bin/env bash
# Rung 5b DIAGNOSIS (read-only): term-by-term Jacobian of one ice-transport reaction.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
make OPT=-O0 >/tmp/it_b.log 2>&1 || { echo build FAIL; tail -20 /tmp/it_b.log; exit 3; }
gfortran -O0 -ffree-line-length-none -Jbuild -Ibuild -c tests/ice_transport_jac_diag.f90 -o build/ice_transport_jac_diag.o 2>>/tmp/it_b.log \
  || { echo "diag compile FAIL"; tail -40 /tmp/it_b.log; exit 3; }
OBJS="$(ls build/*.o build/dust/*.o 2>/dev/null | grep -vE 'build/main\.o|_check\.o|_diag\.o' | tr '\n' ' ')"
gfortran -O0 build/ice_transport_jac_diag.o $OBJS -o bin/ice_transport_jac_diag 2>>/tmp/it_b.log \
  || { echo "diag link FAIL"; tail -40 /tmp/it_b.log; exit 3; }
R="$(mktemp -d)"; trap 'rm -rf "$R"' EXIT
cp networks/reduced_CHO/*.in "$R"/; cp inputs/element.in "$R"/
sed -e 's/^coagulation = .*/coagulation = 1/' -e 's/^sparsity = .*/sparsity = symbolic/' \
    tests/fixture_rung3/parameters.in > "$R"/parameters.in
( cd "$R" && "$ROOT/bin/ice_transport_jac_diag" )
