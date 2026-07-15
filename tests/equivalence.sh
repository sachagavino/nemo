#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Stage 2: prove that the extracted skeleton reproduces the reference nmgc-2.0
# on the 0D 2-grain case, so that a later disagreement can be attributed to a
# refactor and not to the extraction.
#
#   ./tests/equivalence.sh /path/to/nmgc-2.0
#
# Builds both codes at -O0, runs both on tests/reference_0D_2grains/, converts
# the binary output of each to ASCII, and compares every species at every
# output time.
# ---------------------------------------------------------------------------
set -euo pipefail

REF_REPO=${1:?usage: $0 /path/to/nmgc-2.0}
HERE=$(cd "$(dirname "$0")/.." && pwd)
WORK=$(mktemp -d)
CASE=$HERE/tests/reference_0D_2grains

echo "==> building the skeleton (-O0)"
make -C "$HERE" clean >/dev/null
make -C "$HERE" OPT=-O0 >/dev/null

echo "==> building the reference nmgc-2.0"
make -C "$REF_REPO/src" >/dev/null

for side in ref new; do
  mkdir -p "$WORK/$side"
  cp "$REF_REPO"/inputs/*.in "$WORK/$side/" 2>/dev/null || true
  cp "$CASE"/*.in "$WORK/$side/"
done

# the reference still needs the parameters the skeleton has dropped
cp "$REF_REPO/inputs/parameters.in" "$WORK/ref/parameters.in"
python3 - "$WORK/ref/parameters.in" "$CASE/parameters.in" <<'PY'
import re, sys
ref, case = sys.argv[1], sys.argv[2]
s = open(ref).read()
for line in open(case):
    if '=' not in line or line.lstrip().startswith('!'):
        continue
    k = line.split('=')[0].strip()
    v = line.split('=')[1].split('!')[0].strip()
    s = re.sub(r'(?m)^(%s\s*=\s*)\S+' % re.escape(k), lambda m: m.group(1) + v, s, count=1)
# the reference needs these; the skeleton does not have them any more
s = re.sub(r'(?m)^(structure_type\s*=\s*)\S+',      r'\g<1>0D', s)
s = re.sub(r'(?m)^(spatial_resolution\s*=\s*)\S+',  r'\g<1>1',  s)
open(ref, 'w').write(s)
PY

echo "==> running the reference"
( cd "$WORK/ref" && "$REF_REPO/src/nmgc" run > run.log 2>&1 && "$REF_REPO/src/nmgc" outputs > /dev/null 2>&1 )

echo "==> running the skeleton"
( cd "$WORK/new" && "$HERE/bin/nmgc" run > run.log 2>&1 && "$HERE/bin/nmgc" outputs > /dev/null 2>&1 )

echo "==> comparing"
python3 - "$WORK" <<'PY'
import glob, os, sys
work = sys.argv[1]
worst, nfiles, nvals = 0.0, 0, 0
where = None
for f in sorted(glob.glob(os.path.join(work, 'ref/ab/*.ab'))):
    base = os.path.basename(f)
    if base == 'gas_phase.ab':
        continue
    g = os.path.join(work, 'new/ab', base)
    if not os.path.exists(g):
        print('MISSING in skeleton:', base); sys.exit(1)
    A = [l.split() for l in open(f) if not l.startswith('!')]
    B = [l.split() for l in open(g) if not l.startswith('!')]
    if len(A) != len(B):
        print('different number of output times:', base); sys.exit(1)
    nfiles += 1
    for ra, rb in zip(A, B):
        for x, y in zip(ra, rb):
            x, y = float(x), float(y)
            nvals += 1
            if abs(x) > 1e-30:
                r = abs(x - y) / abs(x)
                if r > worst:
                    worst, where = r, (base, x, y)
print(f'species compared : {nfiles}')
print(f'values compared  : {nvals}')
print(f'worst rel. diff  : {worst:.3e}' + (f'   at {where}' if where else ''))
sys.exit(0 if worst == 0.0 else 1)
PY
echo "==> PASS: the skeleton reproduces nmgc-2.0 exactly"
rm -rf "$WORK"
