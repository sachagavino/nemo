#!/usr/bin/env bash
# ===========================================================================
# tests/reaction_census.sh -- per-type reaction/species census at N=4 and N=20.
#
# Same config (reduced_CHO network, coagulation on, 2-phase, surface chemistry on,
# mass_ratio=2); only a_max changes to set the bin count (N=4: a_max=1.2e-6; N=20:
# a_max=5.0e-5, the v1 fiducial). Prints the census from init_gasgrain for each N so
# the paper's dimensionality numbers can be quoted from one run.
# ===========================================================================
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
FIX="tests/fixture_jac_figure"          # reduced_CHO + coag + 2-phase + grain reactions
BLOG=/tmp/census_build.log

echo "[census] building (-O0)..."
make OPT=-O0 >"$BLOG" 2>&1 || { echo "library build FAIL"; tail -25 "$BLOG"; exit 3; }
gfortran -O0 -ffree-line-length-none -Jbuild -Ibuild -c tests/reaction_census.f90 \
  -o build/reaction_census.o 2>>"$BLOG" || { echo "harness compile FAIL"; tail -40 "$BLOG"; exit 3; }
OBJS="$(ls build/*.o build/dust/*.o 2>/dev/null | grep -vE 'build/main\.o|_check\.o|_dump\.o|_census\.o' | tr '\n' ' ')"
gfortran -O0 build/reaction_census.o $OBJS -o bin/reaction_census 2>>"$BLOG" \
  || { echo "harness link FAIL"; tail -40 "$BLOG"; exit 3; }

run_one () {   # $1 = N label, $2 = a_max
  local R; R="$(mktemp -d)"
  cp "$FIX"/*.in "$R"/
  sed -i -e "s/^a_max =.*/a_max =  $2/" "$R"/parameters.in
  echo
  echo "########## N=$1  (a_max=$2) ##########"
  ( cd "$R" && "$ROOT/bin/reaction_census" ) | grep -vaE 'writing|hdf5'
  rm -rf "$R"
}

run_one 4  "1.2000E-06"
run_one 20 "5.0000E-05"
echo
echo "[census] done."
