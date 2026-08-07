#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Stage 3 guard: the derived dust-grid path == the tabulated path, bit for bit.
#
#   ./tests/roundtrip_derived.sh            # default build (-O2)
#   OPT=-O0 ./tests/roundtrip_derived.sh    # match the equivalence build
#
# WHY THIS EXISTS
# ---------------
# equivalence.sh only exercises the *tabulated* grid path: its reference is
# nmgc-2.0's 2-bin, ~8-mass-decade grid, which no mass_ratio<=2 derivation can
# reproduce (see docs/STAGE3_GRID.md). The *derived* path is the new default and
# had no automated guard. This test provides one, and -- unlike equivalence.sh --
# it needs no external reference binary, so it runs anywhere (CI included).
#
# WHAT IT PROVES
# --------------
# get_grain_radii() fills the same per-bin arrays (radii, masses, temperatures,
# 1/abundance) whether they come from the derivation (dust_grid.f90) or from the
# tabulated reader, and the chemistry consumes them identically. Method:
#
#   1. run a DERIVED + MRN grid; it exports its active grid to
#      dust_grid_active.out (full es24.16e3 precision, round-trip-exact for a
#      double) at init;
#   2. feed that export back in as dust_grid_table.in on a second run with
#      dust_grid_source = tabulated, dust_ic = tabulated, every other input held
#      fixed;
#   3. require abundances.out (the whole trajectory: time, physical state,
#      abundances at every output) to be byte-identical between the two runs.
#
# A small grid is used on purpose: bin count changes only the species count and
# hence the runtime, never which code paths are taken, so a 4-bin grid guards
# exactly what a 20-bin science grid would, in a fraction of the time. The knobs
# below are overridable from the environment.
# ---------------------------------------------------------------------------
set -euo pipefail

# --- test knobs (kept small so the guard is cheap; override via env) --------
RT_AMIN=${RT_AMIN:-5.000E-07}   # a_min [cm]; with a_max below and mass_ratio=2
RT_AMAX=${RT_AMAX:-1.100E-06}   # a_max [cm] -> 4 mass-doubling bins
RT_MRATIO=${RT_MRATIO:-2.000E+00}
RT_STOP=${RT_STOP:-1.000E+03}   # stop_time [yr]
RT_NOUT=${RT_NOUT:-10}          # number of outputs

HERE=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
BIN="$HERE/bin/nmgc"

echo "==> building NEMO (${OPT:--O2})"
make -C "$HERE" clean >/dev/null
make -C "$HERE" ${OPT:+OPT=$OPT} >/dev/null

# portable in-place sed: BSD/macOS `sed -i` needs a backup-suffix argument and
# eats the following -e, so avoid -i entirely and go through a temp file. No
# arrays either, so this is safe on the bash 3.2 that macOS ships as /bin/bash.
sedi() { # sedi <file> <expr> [<expr> ...]
  local f=$1; shift
  local tmp="$f.sedi$$" e
  cp "$f" "$tmp"
  for e in "$@"; do
    sed "$e" "$tmp" > "$tmp.n" && mv "$tmp.n" "$tmp"
  done
  mv "$tmp" "$f"
}

# helper: rewrite the grid + run-length keys in a parameters.in to the test values
rewrite() { # file  (edits in place)
  local f=$1
  sedi "$f" \
    "s|^\(a_min[[:space:]]*=[[:space:]]*\)[^!]*|\1$RT_AMIN |" \
    "s|^\(a_max[[:space:]]*=[[:space:]]*\)[^!]*|\1$RT_AMAX |" \
    "s|^\(mass_ratio[[:space:]]*=[[:space:]]*\)[^!]*|\1$RT_MRATIO |" \
    "s|^\(stop_time[[:space:]]*=[[:space:]]*\)[^!]*|\1$RT_STOP |" \
    "s|^\(nb_outputs[[:space:]]*=[[:space:]]*\)[^!]*|\1$RT_NOUT |"
}

# --- side 1: the DERIVED + MRN grid (the new default path) ------------------
echo "==> running the derived grid (dust_grid_source=derived, dust_ic=MRN)"
mkdir -p "$WORK/derived"
cp "$HERE"/inputs/*.in "$WORK/derived/"
rewrite "$WORK/derived/parameters.in"
sedi "$WORK/derived/parameters.in" \
       "s|^\(dust_grid_source[[:space:]]*=[[:space:]]*\)[^!]*|\1derived |" \
       "s|^\(dust_ic[[:space:]]*=[[:space:]]*\)[^!]*|\1MRN |"
( cd "$WORK/derived" && "$BIN" run > run.log 2>&1 )
NBINS=$(grep -vc '^!' "$WORK/derived/dust_grid_active.out")
echo "    derived grid: $NBINS bins exported to dust_grid_active.out"

# --- side 2: feed that export back through the TABULATED path ---------------
echo "==> running the tabulated round-trip (fed the derived export back in)"
mkdir -p "$WORK/tabulated"
cp "$HERE"/inputs/*.in "$WORK/tabulated/"
rewrite "$WORK/tabulated/parameters.in"
sedi "$WORK/tabulated/parameters.in" \
       "s|^\(dust_grid_source[[:space:]]*=[[:space:]]*\)[^!]*|\1tabulated |" \
       "s|^\(dust_ic[[:space:]]*=[[:space:]]*\)[^!]*|\1tabulated |"
# the derived run's active grid IS the tabulated input
cp "$WORK/derived/dust_grid_active.out" "$WORK/tabulated/dust_grid_table.in"
( cd "$WORK/tabulated" && "$BIN" run > run.log 2>&1 )

# --- compare ----------------------------------------------------------------
echo "==> comparing"
fail=0
if cmp -s "$WORK/derived/abundances.out" "$WORK/tabulated/abundances.out"; then
  echo "    abundances.out          : byte-identical"
else
  echo "    abundances.out          : DIFFER"; fail=1
fi
if cmp -s "$WORK/derived/elemental_abundances.out" "$WORK/tabulated/elemental_abundances.out"; then
  echo "    elemental_abundances.out: byte-identical"
else
  echo "    elemental_abundances.out: DIFFER"; fail=1
fi

if [ "$fail" -ne 0 ]; then
  echo
  echo "==> FAIL: the derived and tabulated paths diverged."
  echo "    Converting to ASCII for a species-level diagnostic..."
  ( cd "$WORK/derived"   && "$BIN" outputs >/dev/null 2>&1 ) || true
  ( cd "$WORK/tabulated" && "$BIN" outputs >/dev/null 2>&1 ) || true
  python3 - "$WORK" <<'PY' || true
import glob, os, sys
work = sys.argv[1]
worst, where = 0.0, None
for f in sorted(glob.glob(os.path.join(work, 'derived/ab/*.ab'))):
    base = os.path.basename(f)
    if base == 'gas_phase.ab':   # species-name column, not a trajectory
        continue
    g = os.path.join(work, 'tabulated/ab', base)
    if not os.path.exists(g):
        continue
    A = [l.split() for l in open(f) if not l.startswith('!')]
    B = [l.split() for l in open(g) if not l.startswith('!')]
    for ra, rb in zip(A, B):
        for x, y in zip(ra, rb):
            x, y = float(x), float(y)
            if abs(x) > 1e-30:
                r = abs(x - y) / abs(x)
                if r > worst:
                    worst, where = r, (base, x, y)
print(f'    worst rel. diff : {worst:.3e}' + (f'   at {where}' if where else ''))
PY
  echo "    work kept for inspection: $WORK"
  exit 1
fi

echo "==> PASS: the derived path round-trips through the tabulated path bit-identical"
rm -rf "$WORK"
