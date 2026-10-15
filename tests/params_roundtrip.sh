#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Guard: the parameters.in rewrite is lossless and idempotent.
#
#   ./tests/params_roundtrip.sh
#
# WHY THIS EXISTS
# ---------------
# At start-up NEMO reads parameters.in and then REWRITES it in canonical form
# (write_parameters; the original is kept as parameters.in.bak). Re-running in
# the same directory therefore reads the REWRITTEN file. Two silent failure modes
# existed / can reappear:
#   (a) a key accepted by read_parameters_in but not emitted by write_parameters
#       is dropped on rewrite, so a re-run silently falls back to its default;
#   (b) reals were written with 4 significant digits (es10.3e2), so a re-run
#       silently used rounded values (e.g. fixture_rung3 a_max = 1.0001E-06 was
#       rewritten as 1.000E-06).
#
# WHAT IT CHECKS
# --------------
#   1. static: every key accepted by read_parameters_in (primary name of each
#      case(...) on an uncommented line) is written by write_parameters;
#   2. preservation: for every fixture, each numeric value in the ORIGINAL
#      parameters.in equals the value in the rewritten file (double compare);
#   3. idempotence: a second run in the same directory produces a byte-identical
#      abundances.out, and the second rewrite of parameters.in is byte-identical
#      to the first.
# No external reference binary is needed. Fixtures run in seconds.
# ---------------------------------------------------------------------------
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
SRC=$ROOT/src/input_output.f90
BIN=$ROOT/bin/nmgc
fail=0

# ---- 1. static key audit ---------------------------------------------------
reader_keys=$(awk '/^subroutine read_parameters_in\(\)/,/^end subroutine read_parameters_in/' "$SRC" \
  | grep -v '^[[:space:]]*!' | grep -o "case *( *'[A-Za-z0-9_]*'" | sed "s/case *( *'//; s/'//" | sort -u)
writer_keys=$(awk '/^subroutine write_parameters\(\)/,/^end subroutine write_parameters/' "$SRC" \
  | grep -v '^[[:space:]]*!' | grep -o "'[A-Za-z0-9_]* = '" | sed "s/'//g; s/ = //" | sort -u)
missing=$(comm -23 <(echo "$reader_keys") <(echo "$writer_keys") || true)
nr=$(echo "$reader_keys" | wc -l); nw=$(echo "$writer_keys" | wc -l)
if [ -n "$missing" ]; then
  echo "FAIL [static] keys read but never written (silently dropped on rewrite):"; echo "$missing" | sed 's/^/    /'
  fail=1
else
  echo "PASS [static] all $nr reader keys are written ($nw writer keys)"
fi

[ -x "$BIN" ] || { echo "build first: $BIN not found"; exit 2; }

# ---- 2 + 3. per-fixture preservation and idempotence ------------------------
numeric_pairs() {   # key value, for numeric values, comments stripped
  sed 's/!.*//' "$1" | awk -F'=' 'NF==2 {gsub(/[ \t]/,"",$1); gsub(/[ \t]/,"",$2);
       if ($2 ~ /^[-+]?[0-9]*\.?[0-9]+([eEdD][-+]?[0-9]+)?$/) { v=$2; gsub(/[dD]/,"e",v); print $1, v } }'
}
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
for fx in "$ROOT"/tests/fixture_*; do
  name=$(basename "$fx"); d=$TMP/$name; mkdir -p "$d"
  cp "$fx"/*.in "$d"/
  [ -f "$d"/gas_species.in ] || cp "$ROOT"/networks/reduced_CHO/*.in "$d"/   # fixtures that use reduced_CHO
  cp "$fx"/parameters.in "$d"/parameters.in                                   # fixture params win over the network's
  [ -f "$d"/element.in ] || cp "$ROOT"/inputs/element.in "$d"/
  cp "$d"/parameters.in "$d"/parameters.original
  ( cd "$d" && "$BIN" run > run1.log 2>&1 ) || { echo "FAIL [$name] first run crashed"; fail=1; continue; }
  cp "$d"/abundances.out "$d"/abundances.run1; cp "$d"/parameters.in "$d"/parameters.rewrite1

  # 2. every numeric value of the original survives the rewrite exactly
  bad=$(join <(numeric_pairs "$d"/parameters.original | sort) <(numeric_pairs "$d"/parameters.rewrite1 | sort) \
        | awk '{ if (($2+0) != ($3+0)) print "    " $1 ": " $2 " -> " $3 }')
  if [ -n "$bad" ]; then echo "FAIL [$name] values changed by the rewrite:"; echo "$bad"; fail=1
  else echo "PASS [$name] all numeric values preserved by the rewrite"; fi

  # 3. idempotence: re-run in place
  ( cd "$d" && "$BIN" run > run2.log 2>&1 ) || { echo "FAIL [$name] re-run crashed"; fail=1; continue; }
  if cmp -s "$d"/abundances.run1 "$d"/abundances.out && cmp -s "$d"/parameters.rewrite1 "$d"/parameters.in; then
    echo "PASS [$name] re-run in place: abundances.out and parameters.in byte-identical"
  else
    echo "FAIL [$name] re-run in place differs (abundances.out and/or second rewrite)"; fail=1
  fi
done

[ $fail -eq 0 ] && echo "ALL PASS" || { echo "SOME CHECKS FAILED"; exit 1; }
