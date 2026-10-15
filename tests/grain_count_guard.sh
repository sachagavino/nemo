#!/usr/bin/env bash
# ===========================================================================
# tests/grain_count_guard.sh -- v1 guard: reject grids with more than 99 bins.
#
# Grain/ice species-name indices are written with a 2-digit format (I2.2) into a
# character(2) buffer (input_output.f90, gasgrain.f90). An index >= 100 overflows
# the field to '**', so bins 100+ collapse onto duplicate names ('GRAIN**') and
# coagulation products are silently routed to the wrong bin -- mass is NOT
# conserved (the x2.37 gain seen on the 149-bin grid). Until the name builders AND
# parsers are widened to 3 digits (v2), init must REFUSE any grid with > 99 bins,
# on BOTH the derived path (dust_grid_count_bins) and the tabulated path.
#
# This test asserts:
#   A. derived    > 99 bins -> nmgc fails at init with the guard message.
#   B. tabulated  > 99 rows -> nmgc fails at init with the guard message.
#   C. derived   <= 99 bins -> the guard does NOT fire (grid is accepted).
# Exit 0 = all three behave. Run from anywhere.
# ===========================================================================
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
FIX="tests/fixture_dust_coag"
NMGC="$ROOT/bin/nmgc"
GUARD_RE='exceeding the hard limit'
rc_all=0

echo "[guard] building nmgc..."
make >/tmp/grainguard_build.log 2>&1 || { echo "build FAIL"; tail -20 /tmp/grainguard_build.log; exit 3; }

# setk <file> <key> <value> : replace the whole "key = ..." line in place.
setk() {
  local f="$1" key="$2" val="$3"
  awk -v k="$key" -v v="$val" '
    { if ($1==k && $2=="=") { print k" = "v } else { print } }' "$f" >"$f.tmp" && mv "$f.tmp" "$f"
}

newrun() {  # -> echoes a fresh temp run dir seeded from the fixture
  local r; r="$(mktemp -d)"; cp "$FIX"/*.in "$r"/ || { echo "fixture copy FAIL" >&2; exit 3; }
  echo "$r"
}

# --- Case A: derived, > 99 bins (mass_ratio 1.08 on the fixture range -> ~136) ---
A="$(newrun)"
setk "$A/parameters.in" dust_grid_source derived
setk "$A/parameters.in" dust_ic          MRN
setk "$A/parameters.in" mass_ratio        "1.080E+00"
outA="$( cd "$A" && "$NMGC" run 2>&1 )"; rcA=$?
rm -rf "$A"
if [ $rcA -ne 0 ] && grep -q "$GUARD_RE" <<<"$outA"; then
  echo "[guard] A derived >99  : PASS (nmgc exited $rcA, guard fired)"
else
  echo "[guard] A derived >99  : FAIL (rc=$rcA, guard message present: $(grep -qc "$GUARD_RE" <<<"$outA" && echo yes || echo no))"
  echo "$outA" | tail -5; rc_all=1
fi

# --- Case B: tabulated, 100 rows (> 99) ---
B="$(newrun)"
python3 - "$B/dust_grid_table.in" <<'PY'
import sys
with open(sys.argv[1], "w") as f:
    f.write("! test tabulated grid: 100 rows -> must fail the >99 guard\n")
    f.write("! columns: radius[cm]  1/n_k  T_d[K]  T_CR[K]\n")
    for i in range(100):
        a = 1.0e-6 * (1.1 ** i)
        f.write(f"{a:.15E}  1.000000E+10  1.000000E+01  1.500000E+01\n")
PY
setk "$B/parameters.in" dust_grid_source tabulated
setk "$B/parameters.in" dust_ic          tabulated
outB="$( cd "$B" && "$NMGC" run 2>&1 )"; rcB=$?
rm -rf "$B"
if [ $rcB -ne 0 ] && grep -q "$GUARD_RE" <<<"$outB"; then
  echo "[guard] B tabulated>99 : PASS (nmgc exited $rcB, guard fired)"
else
  echo "[guard] B tabulated>99 : FAIL (rc=$rcB, guard message present: $(grep -qc "$GUARD_RE" <<<"$outB" && echo yes || echo no))"
  echo "$outB" | tail -5; rc_all=1
fi

# --- Case C: derived, <= 99 bins (mass_ratio 1.15 -> ~75): guard must NOT fire ---
C="$(newrun)"
setk "$C/parameters.in" dust_grid_source derived
setk "$C/parameters.in" dust_ic          MRN
setk "$C/parameters.in" mass_ratio        "1.150E+00"
outC="$( cd "$C" && timeout 120 "$NMGC" run 2>&1 )"; rcC=$?
rm -rf "$C"
if ! grep -q "$GUARD_RE" <<<"$outC" && grep -q "NB OF GRAINS" <<<"$outC"; then
  echo "[guard] C derived<=99  : PASS (grid accepted, guard silent; $(grep 'NB OF GRAINS' <<<"$outC" | tail -1 | tr -s ' '))"
else
  echo "[guard] C derived<=99  : FAIL (guard fired on a <=99 grid, or init never reached; rc=$rcC)"
  echo "$outC" | tail -8; rc_all=1
fi

echo
if [ $rc_all -eq 0 ]; then echo "[guard] ALL PASS"; else echo "[guard] FAILURES ABOVE"; fi
exit $rc_all
