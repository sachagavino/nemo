#!/usr/bin/env bash
# ===========================================================================
# tests/coag_diffusion_gate.sh -- Phase III Part B 1d.
#
# Numerical-diffusion convergence gate. Runs pure coagulation (chemistry off) on
# the tabulated grid+IC for the CONSTANT and ADDITIVE kernels over a mass_ratio
# sweep {2.0,1.5,1.3,1.15} with a_min,a_max fixed, and asserts that NEMO's error
# vs the analytic solution DECREASES MONOTONICALLY as the grid refines.
#
# Norm: CONTINUOUS L1 (L&L 2021 Eq.40, Gauss-integrated over each bin). TWO
# continuous-L1 norms are gated: the mass-weighted one on g=x f (the conservation-law
# norm) AND the surface-area-weighted one on the number density (weight a^2 ~ m^(2/3)),
# which is what the chemistry actually sees (accretion, monolayer site count ~ int a^2 n).
# Both must decrease monotonically under refinement. If the area-weighted norm also
# converges, the "over-diffused tail is chemically immaterial" claim is earned.
# We do NOT gate on the discrete L1 (Eq.41) or L2: both are POINT-evaluated and fragile
# in the sparse over-diffused large-mass tail -- non-monotone at evolved times, and
# their earlier apparent monotonicity was partly an artifact of the old [m_k,m_{k+1})
# IC offset (removed by the geometric-edge IC). They are reported as diagnostics only.
# NEMO is the k=0 Kovetz-Olund scheme (DustPy shows the same tail behaviour).
#
# GATE = monotone continuous-L1 AND discrete-L1 at the EVOLVED time (T=2 constant,
# tau=1 additive), both kernels. At that time ~50-63% of grains have coagulated,
# so the coagulation operator is genuinely exercised (unlike the tau=0.01 EOC
# point, which is IC-representation-dominated and kernel-independent).
#
# Also prints, as SUPPLEMENTARY context (not gated):
#   - spatial order at the L&L EOC time (dimless t=0.01);
#   - error-vs-time boundedness at fixed resolution.
#
# Exit 0 = gate PASS.
# ===========================================================================
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
BLOG=/tmp/coagdiff_build.log

echo "[coag_diffusion] building..."
make >"$BLOG" 2>&1 || { echo "BUILD FAIL"; tail -25 "$BLOG"; exit 3; }

STAGE="${1:-gate}"     # gate (default) | all
python3 scripts/coag_diffusion.py --stage "$STAGE"
rc=$?
[ $rc -eq 0 ] && echo "[coag_diffusion] gate PASS." || echo "[coag_diffusion] gate FAIL (rc=$rc)."
exit $rc
