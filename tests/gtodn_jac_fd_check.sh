#!/usr/bin/env bash
# Rung 5 entry-level FD check: analytic grain-column Jacobian vs finite-difference RHS.
# Run from repo root. Exit 0 = PASS.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"

echo "[fd] building library (-O0)..."; make OPT=-O0 >/tmp/fd_b.log 2>&1 || { echo build FAIL; tail -20 /tmp/fd_b.log; exit 3; }
gfortran -O0 -ffree-line-length-none -Jbuild -Ibuild -c tests/gtodn_jac_fd_check.f90 -o build/gtodn_jac_fd_check.o 2>>/tmp/fd_b.log \
  || { echo "check compile FAIL"; tail -30 /tmp/fd_b.log; exit 3; }
OBJS="$(ls build/*.o build/dust/*.o 2>/dev/null | grep -vE 'build/main\.o|_check\.o' | tr '\n' ' ')"
gfortran -O0 build/gtodn_jac_fd_check.o $OBJS -o bin/gtodn_jac_fd_check 2>>/tmp/fd_b.log \
  || { echo "check link FAIL"; tail -30 /tmp/fd_b.log; exit 3; }

R="$(mktemp -d)"; trap 'rm -rf "$R"' EXIT
cp networks/reduced_CHO/*.in "$R"/; cp inputs/element.in "$R"/
sed -e 's/^coagulation = .*/coagulation = 1/' \
    -e 's/^sparsity = .*/sparsity = symbolic/' \
    tests/fixture_rung3/parameters.in > "$R"/parameters.in
# Optional: shrink the coagulation kernel to isolate the dynamic-GTODN chemistry
# Jacobian from the (separately owned) type-50 ice-transport Jacobian.
if [ -n "${NEMO_TEST_KERNEL_K0:-}" ]; then
  sed -i -e "s/^constant_kernel_k0 = .*/constant_kernel_k0 = ${NEMO_TEST_KERNEL_K0}/" "$R"/parameters.in
  echo "[fd] kernel override: constant_kernel_k0 = ${NEMO_TEST_KERNEL_K0}"
fi

echo "[fd] running FD check (coagulation on)..."
( cd "$R" && "$ROOT/bin/gtodn_jac_fd_check" ); rc=$?
exit $rc
