#!/usr/bin/env bash
# ===========================================================================
# tests/coag_ice_ic.sh
#
# Gate (d) for the ice initial-condition machinery, in two parts:
#
#   1. round-trip + shape (tests/coag_ice_ic_check.f90): area-weighted default
#      distributes each base-ice total X_total across bins as X_k ~ n_k a_k^2, so
#      sum_k X_k = X_total (round-trip) AND X_k/(n_k a_k^2) is constant (shape).
#
#   2. coverage-guard negative control: with a per-species OVERRIDE
#      (surface_ice_distribution.in) dumping all the ice onto a near-empty top bin,
#      the coverage guard must WARN; the area-weighted default on the same grid must
#      be SILENT. (theta is not capped; the guard only flags gross implausibility.)
#
# Run from repo root; exit 0 = PASS.
# ===========================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "[coag_ice_ic] building modules + driver (-O0)..."
make OPT=-O0 >/tmp/ciic_build.log 2>&1 || { echo "module build FAILED"; tail -30 /tmp/ciic_build.log; exit 3; }
gfortran -O0 -ffree-line-length-none -Jbuild -Ibuild -c tests/coag_ice_ic_check.f90 -o build/coag_ice_ic_check.o 2>>/tmp/ciic_build.log \
  || { echo "driver compile FAILED"; tail -30 /tmp/ciic_build.log; exit 3; }
OBJS="$(ls build/*.o build/dust/*.o 2>/dev/null | grep -vE 'build/main\.o|_check\.o' | tr '\n' ' ')"
gfortran -O0 build/coag_ice_ic_check.o $OBJS -o bin/coag_ice_ic_check 2>>/tmp/ciic_build.log \
  || { echo "driver link FAILED"; tail -30 /tmp/ciic_build.log; exit 3; }

fail=0

# ---- part 1: round-trip + shape (MRN grains so bins are all populated) --------
R1="$(mktemp -d)"; cp tests/fixture_ice_transport/*.in "$R1"/
( cd "$R1" && "$ROOT/bin/coag_ice_ic_check" ) | tee /tmp/ciic_run.log
grep -q "RESULT: PASS" /tmp/ciic_run.log || fail=1
rm -rf "$R1"

# ---- part 2: coverage-guard negative control (monodisperse: bins 2,3 near-empty) --
mono() { sed 's/^dust_ic = MRN/dust_ic = monodisperse/; s/^stop_time = .*/stop_time =  1.000E+00/; s/^nb_outputs = .*/nb_outputs = 1/' "$1"/parameters.in > "$1"/p && mv "$1"/p "$1"/parameters.in; }

echo "[coag_ice_ic] negative control A: override dumps CO on near-empty bin 3 -> expect WARN"
RA="$(mktemp -d)"; cp tests/fixture_ice_transport/*.in "$RA"/; mono "$RA"
printf '! base   f1 f2 f3\nCO   0.0 0.0 1.0\n' > "$RA"/surface_ice_distribution.in
( cd "$RA" && "$ROOT/bin/nmgc" run > run.log 2>&1 || true )
if grep -q "Warning (ice IC):" "$RA"/run.log; then echo "   WARN fired (correct)"; else echo "   FAIL: no warning on implausible override"; fail=1; fi
rm -rf "$RA"

echo "[coag_ice_ic] negative control B: area-weighted default (no override) -> expect SILENT"
RB="$(mktemp -d)"; cp tests/fixture_ice_transport/*.in "$RB"/; mono "$RB"
( cd "$RB" && "$ROOT/bin/nmgc" run > run.log 2>&1 || true )
if grep -q "Warning (ice IC):" "$RB"/run.log; then echo "   FAIL: default path warned"; fail=1; else echo "   silent (correct)"; fi
rm -rf "$RB"

if [ "$fail" -eq 0 ]; then echo "[coag_ice_ic] PASS"; else echo "[coag_ice_ic] FAIL"; fi
exit $fail
