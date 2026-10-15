# NEMO — Phase II Rung 5: amendment on type-50 verification and the FD noise floor

**Status:** correction note appended after the Rung 5 chemistry-Jacobian commit. Amends,
does not rewrite, the landed record. Companion to
`docs/PhaseII_rung5_dynamic_gtodn.md` and `docs/PhaseII_rung5_site_classification.md`.

## What the landed record says, and why it was provisional

The Rung 5 commit and its doc note that full-network `021≡022` is **"blocked on 5b"**,
on the working hypothesis that the type-50 ice-transport Jacobian had a missing gain-side
(non-reactant) entry analogous to the live divisor. That framing was recorded before the
diagnosis and is now **retracted**.

## The finding (diagnosis + exhaustive sweep)

There is **no type-50 Jacobian gap**. `get_jacobian` is entry-exact for the coagulation
subsystem — ice-transport and grain-grain, loss and gain, both charge columns.

- **Analytic completeness (no differencing).** `get_jacobian` (isolated to type-50) vs an
  independent per-reaction analytic hand-sum: 393 nonzero entries per state, checked in a
  light (ice 1e-10) and heavy (ice 1e-8) state → **786 comparisons, 0 mismatches**, max
  |A−B| = 2.0e-28, max relative = 2.3e-16 (machine epsilon). Covers grain columns
  (`∂/∂Y(GRAIN_k)`, the hypothesised "gap") and ice columns (`∂/∂Y(J_iX)`). Because the
  hand-sum loops each reaction independently, the exact match proves the loop drops no
  reaction and deposits every Podolak–Brauer weight and sign correctly.
- **Independent numerical anchor.** 10 subset-FD anchors on a background-free subset RHS
  (ice-transport and grain-grain, bins 2/3/4, both charge columns) each match the analytic
  (rel ≤ ~1e-10, truncation-limited). This catches any assumption the two analytic paths
  might share (e.g. a missing self-collision 2×) — none found.

**Why the hypothesis was wrong.** Ice transport carries the collision partner `GRAIN_j`
as an *explicit reactant*, so the reactant differentiator visits the reaction and deposits
`∂flux/∂Y(GRAIN_j)` to all compounds, including the redistributed product bins. The gain
entries are produced. This is genuinely unlike the live divisor, where the grain entered
only through a rate coefficient with no compound to trigger the column.

## The real cause of the `021≢022` failures: an FD-oracle resolution limit

The full-network `021≢022` failures on type-50 rows are **finite-difference cancellation
artifacts in the oracle, not Jacobian defects.** The coagulation coupling is ~1e-5–1e-6 of
a species' total (background-dominated) derivative — below the round-off floor of any RHS
differencing. Central differencing forms `hp − hm` in which the true signal is tens of
times below the subtraction's round-off floor; the library RHS and an independent hand RHS
produce byte-identical wrong differences because they accumulate in the same order
(deterministic cancellation, not noise). The same mechanism affects DLSODES's own internal
numerical Jacobian under `mf=121`.

**Consequence for the verification strategy.** The `021≡022` (analytic-vs-FD) gate built in
Rung 4 has a blind spot: **it cannot validate any Jacobian entry whose coupling is below
~1e-5 of the affected species' background derivative**, because the FD oracle is pure noise
there. This was invisible earlier only because the chemistry entries tested with it sit
above that floor. It is a property of finite differencing, not of NEMO.

## Decisions

1. **No Rung 5b.** Type-50 is entry-exact; there is nothing to fix. The "blocked on 5b"
   note is retracted.
2. **The gate is two-legged.** FD-oracle (`021≡022`) validates entries *above* the FD noise
   floor (the chemistry reciprocals — LH, accretion, photodesorption, sticking). A
   **per-reaction-analytic cross-check** validates entries *below* it (type-50 coagulation,
   and any future sub-background coupling) by constructing the derivative independently and
   comparing to `get_jacobian` with no differencing. Analytic-vs-independent-analytic is the
   stronger check where the signal is small; FD is the right check where it is large.
3. **Standing gate artifact.** `tests/ice_transport_jac_sweep.{f90,sh}` (the per-reaction
   sweep) is committed as the type-50 gate leg; `tests/ice_transport_jac_diag.*` is retained
   as the characterization tool.
4. **021 retarget is decided separately, on the 5c payoff, not on this finding.** This
   finding makes `021` *correct* on type-50 but not, by itself, *worth switching to* — those
   entries are dynamically negligible (which is why coagulation ran fine under `121` through
   Rung 4). The retarget is justified only if the *chemistry reciprocal* (the stiff LH
   `−1/Y²`) earns a convergence/robustness win at production scale (gate 5c). See the 5c
   run; hold the retarget until it is in.

## Scope of the proof

Verified for the fiducial `reduced_CHO` structure (4 bins, 11 ice bases, 2-phase). The
argument is general for that structure. A mantle-bearing or many-bin network would want its
own sweep run — the sweep tool now exists for it — but that is outside v1.
