#!/usr/bin/env bash
# ===========================================================================
# tests/coag_diffusion_gate.sh -- Phase III Part B 1d.
#
# Numerical-diffusion convergence gate. Runs pure coagulation (chemistry off) on
# the tabulated grid+IC for the CONSTANT and ADDITIVE kernels over a mass_ratio
# sweep {2.0,1.5,1.3,1.15} with a_min,a_max fixed, and asserts that NEMO's error
# vs the analytic solution DECREASES MONOTONICALLY as the grid refines.
#
# Norm: L1 only (continuous, L&L 2021 Eq.40; discrete, Eq.41) on the mass density
# g=x f. NEMO is the k=0 (piecewise-constant) Kovetz-Olund scheme in L&L's
# taxonomy. We DROPPED the L2 norm: it is not in L&L, it weights errors
# quadratically (dominated by the sparse large-mass tail rather than misplaced
# mass), and L1 is the natural norm for a conservation law -- the physically and
# chemically relevant grain quantities (surface area, site counts) are L1-type
# integrals. See the design-thread summary for the full rationale.
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
