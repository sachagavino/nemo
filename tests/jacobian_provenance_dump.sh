#!/usr/bin/env bash
# ===========================================================================
# tests/jacobian_provenance_dump.sh -- Phase III Jacobian provenance dump.
#
# Builds the read-only dump harness against the library objects, runs it on the
# committed fixture_jac_figure (init only, no integration), prints the gate
# report, and copies jac_pattern.tsv / jac_species.tsv back into the fixture dir
# for reproducibility. Run from anywhere. Exit 0 = gate PASS.
# ===========================================================================
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
FIX="tests/fixture_jac_figure"
BLOG=/tmp/jacdump_build.log

echo "[dump] building library (-O0)..."
make OPT=-O0 >"$BLOG" 2>&1 || { echo "library build FAIL"; tail -30 "$BLOG"; exit 3; }

echo "[dump] compiling harness..."
gfortran -O0 -ffree-line-length-none -Jbuild -Ibuild -c tests/jacobian_provenance_dump.f90 \
  -o build/jacobian_provenance_dump.o 2>>"$BLOG" \
  || { echo "harness compile FAIL"; tail -40 "$BLOG"; exit 3; }

# library objects, excluding the program mains and any test harness objects
OBJS="$(ls build/*.o build/dust/*.o 2>/dev/null | grep -vE 'build/main\.o|_check\.o|_dump\.o' | tr '\n' ' ')"
gfortran -O0 build/jacobian_provenance_dump.o $OBJS -o bin/jacobian_provenance_dump 2>>"$BLOG" \
  || { echo "harness link FAIL"; tail -40 "$BLOG"; exit 3; }

R="$(mktemp -d)"; trap 'rm -rf "$R"' EXIT
cp "$FIX"/*.in "$R"/ || { echo "fixture copy FAIL (is $FIX present?)"; exit 3; }

echo "[dump] running (init only, no integration) in $R ..."
( cd "$R" && "$ROOT/bin/jacobian_provenance_dump" ); rc=$?

if [ $rc -ne 0 ]; then
  echo "[dump] harness exited $rc (gate FAIL or build/setup error) -- see output above."
  exit $rc
fi

cp "$R"/jac_pattern.tsv "$R"/jac_species.tsv "$FIX"/ \
  && echo "[dump] committed TSVs to $FIX/ (jac_pattern.tsv, jac_species.tsv)"
echo "[dump] gate PASS."
