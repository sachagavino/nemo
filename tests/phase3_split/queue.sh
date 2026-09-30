#!/usr/bin/env bash
# queue.sh <dir>... : run staged dirs sequentially; writes STATUS = ok / FAIL(reason).
# A split run whose log lacks the audit banner is marked FAIL (guards against a silently
# coupled run, e.g. split keys lost from parameters.in).
NEMO_ROOT=${NEMO_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}; BIN=$NEMO_ROOT/bin/nmgc
for d in "$@"; do
  ( cd $d
    split=$(grep -c "^split_mode = 1" parameters.in)
    $BIN run > run.log 2>&1; rc=$?
    if [ $rc -ne 0 ]; then echo "FAIL(rc=$rc)" > STATUS
    elif [ $split -eq 1 ] && ! grep -q "Phase III operator split" run.log; then echo "FAIL(split not active)" > STATUS
    else echo ok > STATUS; fi )
done
