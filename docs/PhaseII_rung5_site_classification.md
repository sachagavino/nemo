# NEMO — Phase II Rung 5: GTODN site enumeration & classification

**Status:** decisions LOCKED by the design thread (all six flags resolved). This is
the checked artifact for the Rung-5 site classification and the record of decisions.
**Branch/HEAD at analysis:** `phase2-coagulation` @ `8f546bc`.
**Scope:** v1 = **2-phase** (`is_3_phase = 0`), and — now enforced by a hard scope
guard (Flag ③④⑤) — `modify_rate_flag = 0`, `is_er_cir = 0` whenever coagulation is on.
**Method:** every read of `GTODN` / `GTODN_FIXED` / the derived `ab_lay`, `SUMLAY`
in the RHS builders was traced to its reaction type and reduced to the **net power
of the grain abundance** `Y_tot(k) = Y(GRAIN_k⁰)+Y(GRAIN_k⁻)` in the resulting
`reaction_rates(J)` coefficient, after algebraic cancellations.

---

## 0. The form (adopted) and the corrected floor rule

For any site whose coefficient obeys `coeff ∝ Y_tot(k)^ν` with `ν` locally
constant, the missing column-`GRAIN_k` Jacobian entry is

```
∂F/∂Y(GRAIN_k⁰) = ∂F/∂Y(GRAIN_k⁻) = ν · F / Y_tot(k)
```

`F` = reaction flux (the `tmp_value` `get_jacobian` already forms); both charge
columns identical because the divisor is the total. `ν = +1` linear, `ν = −1`
reciprocal. This is adopted as the canonical form (it subsumes the brief's
`−rate/Y²`).

**Floor rule (corrected — load-bearing).** The floor indicator is **not** keyed on
the sign of `ν`. It is keyed on **whether the coefficient routes through a
`max(Y_tot, floor)` divisor**:

> **Apply the floor (entry = 0 below floor) wherever `max(Y_tot, floor)` appears in
> the coefficient; otherwise no floor.**

- Accretion: bare `×Y_tot`, divides by nothing → **no floor** (`ν=+1`).
- LH: `×GTODN = ×1/max(Y_tot,floor)` → **floored** (`ν=−1`).
- Photodesorption: routes through the floored `SUMLAY` → **floored** even though
  `ν=+1`. (My earlier "floor only for `ν<0`" was correct for the A/B pair but
  breaks here; the divisor-routing rule is exact.)

Below the floor, `max(Y_tot,floor)` is constant, so any coefficient routed through
it has zero grain-derivative there.

---

## 1. Full enumeration

Lines are HEAD `8f546bc`, `src/ode_solver.f90` unless noted.

| # | line(s) | rxn type | how `Y_tot(k)` enters | ν | floor? | class | v1 | Rung-5 action |
|---|---------|----------|-----------------------|---|--------|-------|----|---------------|
| A | 1859–1860 | 99 accretion | explicit `/GTODN` ⇒ `×Y_tot` | +1 | no | **linear** | ✔ | live read + analytic `+F/Y_tot`, both cols |
| B | 2049 | 14 LH surface | explicit `×GTODN` ⇒ `×1/Y_tot` | −1 | yes | **reciprocal (stiff)** | ✔ | floored live read + analytic `−F/Y_tot`, 0 below floor, both cols; gated `modify_rate_flag==0` |
| C | 1884 | 66 photodesorp (ext UV) | via `MLAY/SUMLAY(k)` (floored), cap active when `SUMLAY≥MLAY` | +1 | yes (via SUMLAY) | **in-scope (Flag ②)** | ✔ | RHS read already live (Rung 3); **add** analytic `+F/Y_tot`, 0 below floor, 0 when cap inactive; both cols |
| D | 1903 | 67 photodesorp (CR UV) | same as C | +1 | yes (via SUMLAY) | **in-scope (Flag ②)** | ✔ | same as C |
| E | 1534 | (derived) `ab_lay = nb_sites/GTODN`, **frozen** | feeds A caps, ER | +1 | — | **frozen partial-coupling (Flag ①)** | ✔ | none; documented + cap-binding diagnostic |
| F | 1552–1556 | (derived) `SUMLAY/sumlaysurf/sumlaymant`, live+floored | Rung 3 | — | — | **Rung 3 RHS; entries via C/D** | ✔ | do **not** re-touch RHS read; entries supplied through C/D |
| G | 2099 | 30 ER, `ab_surf≤ab_lay` | `TSQ/GTODN/ab_lay`, `ab_lay=nb_sites/GTODN` ⇒ **cancels** | 0 | — | **cancels** | ✖ guard | none (scope guard) |
| H | 2102 | 30 ER, `ab_surf>ab_lay` | `TSQ/GTODN/ab_surf` ⇒ `×Y_tot` | +1 | no | **deferred v2** | ✖ guard | none (scope guard) |
| I | 2434,2439 | 14 (modify test) | `×GTODN` with `Y·GTODN<1` threshold | piecewise | — | **deferred v2** | ✖ guard | none (scope guard) |
| J | 2445,2451 | 14 (modify rewrite) | `/GTODN` into `DIFF` ⇒ cancels outer `×GTODN` | 0 (PICK) | — | **deferred v2** | ✖ guard | none (scope guard) |
| K | 2203–2204 | (derived) 3-phase `sumlay*` frozen | 3-phase | — | — | **deferred v2** | ✖ guard | none (scope guard) |
| L | 2256 | 40 surf→mant | `alpha ∝ ab_surf·GTODN/nb_sites` then `/ab_surf`, tangled | mixed | — | **deferred v2** | ✖ guard | none (scope guard) |
| M | 988,992 | 10/11 H2 single-grain | scalar `GTODN_FIXED`, `IS_GRAIN_REACTIONS==0` | n/a | — | **out of scope** | ✖ | none |

**In-scope set = three reaction classes, one form:** accretion (99, `ν=+1`, no
floor), photodesorption (66/67, `ν=+1`, floored-via-SUMLAY), LH (14, `ν=−1`,
floored). Sub-finding (LH mantle branch): when a reaction carries a mantle (`'K'`)
compound and `sumlaymant>1`, the LH `×GTODN` and the `/sumlaymant` divide (both
routed through the same floored divisor) **cancel to `ν=0`**. Mantle (`'K'`)
species exist only in the full/3-phase network — `reduced_CHO` and the fiducial
fixture have **zero** — so this is unreachable in v1 and is additionally blocked by
the `is_3_phase=0` guard. **Recommendation:** when coding B, mirror the RHS
`sumlaymant` branch (zero the entry when the divide is active) so 021 stays
FD-exact if a mantle-bearing network is ever run; cheap insurance, not a v1
requirement.

---

## 2. Decisions recorded against each flag

**Unifying form — ADOPTED.** `ν·F/Y_tot`, both charge columns identical,
`−F/Y_tot` (F = flux = `tmp_value`) for the reciprocal. Floor rule restated as
"floor wherever `max(Y_tot,floor)` appears" (§0).

**② SUMLAY Jacobian entries — 021 OWES THEM (forced, not optional).** Rung 4's
GRAIN_RANK superset already registered the 66/67 `(compound, GRAIN_k)` slots; the
FD oracle (022) fills them; production 021 is MITER=1 all-analytic with no
per-entry hybrid, so every registered, FD-filled entry must be analytic in 021 or
gate 5a fails. Photodesorption is a `ν=+1` site (same `+F/Y_tot` as accretion) but
routed through the floored `SUMLAY`, so it carries the floor indicator (and the
`SUMLAY≥MLAY` cap gate). Consequences recorded:
- **Landmine 7 corrected (design thread owns the erratum):** "monolayer count is
  already done" referred to the *RHS read* (Rung 3), which was conflated with the
  Jacobian entry. The entry never existed. Corrected rule: **do not re-touch the
  SUMLAY RHS read; DO supply its Jacobian entries.** Adding
  `d(photodesorption)/dY(GRAIN_k)` is not re-touching the count and is not a
  double-count (each reaction gets its own entry).
- **Gate 5a must include an ice-rich state**, not only the ice-free dust-only
  fiducial (where `SUMLAY≈0` leaves C/D unexercised). The ice-rich state is added
  to the 5a state set.

**③④⑤ confounders — ONE SCOPE GUARD.** Under coagulation, **hard-assert (stop, not
warn)** `modify_rate_flag==0 .and. is_er_cir==0 .and. is_3_phase==0`, with a
message naming the unsupported combination. This makes v1 scope explicit and
enforced (fails loudly rather than silently mis-computing a Jacobian), and matches
the paper's stated v1 scope. The LH reciprocal entry (B) is additionally gated on
`modify_rate_flag==0` (locally safe/redundant). Modified-LH / ER / 3-phase
Jacobians are **deferred to v2** explicitly (rows G–L).

**① ab_lay — FROZEN (confirmed), with two conditions.**
- (a) Recorded here as a **deliberate partial-coupling**: a frozen surface-site
  cap, the analogue of Rung 3's frozen surface area, and the v2 hook. Not silent.
- (b) The "second-order" claim is to be **checked, not asserted**: add a diagnostic
  that the accretion monolayer cap (`ab_lay`) rarely binds in the fiducial runs.
  Rationale: coagulation reduces total grain surface (`∝a²` falls faster than
  number rises), so a frozen `ab_lay` over-estimates available surface — harmless
  only if the cap seldom binds. If the diagnostic shows frequent binding,
  **escalate**; otherwise frozen is validated.

**⑥ grain_abundance_floor — PROMOTE TO A `parameters.in` KNOB (default 2.99e-39).**
Reasons: (a) enables the floor-insensitivity study (0.1× / 10×) without
recompiling — the actual proof the floor is chemically negligible; (b) one knob
enforces that the Rung-3 divisor floor and the Rung-5 rate-site floors share a
single value; (c) the paper implies tunability. Default unchanged, so existing
runs are unaffected. Paper wording to be corrected by the author: "set empirically
from the coupled models" → "set low enough to be chemically negligible — verified
by insensitivity of the results over [range] — and high enough to keep the
reciprocal terms finite" (it is a physically-motivated ~10 grains/AU³ threshold,
not a data fit). The insensitivity run is a validation-phase deliverable, noted now.

---

## 3. Net scope (as decided)

Three reaction classes (99 no-floor, 66/67 floored-via-SUMLAY, 14 floored), one
derivative form (`ν·F/Y_tot`, floor where a `max(Y_tot,floor)` divisor appears),
one confounder scope-guard (hard-assert under coag), `ab_lay` frozen + one
cap-binding diagnostic, and `grain_abundance_floor` promoted to a shared knob.
Gate 5a state set extended to include an ice-rich state.

---

## 4. Rung 5 implementation outcome (this commit)

**Status:** chemistry Jacobian **COMPLETE and FD-verified**. Type-50 ice-transport
Jacobian **VERIFIED entry-exact** by an exhaustive per-reaction-analytic sweep
(§4.3) — **Rung 5b is retracted; there is no gap to complete.** The residual
full-network `021≢022` on type-50 rows is a **finite-difference oracle artifact**
(§4.4), not a Jacobian defect. Production `mf` stays **121** pending the Rung 5c
convergence test (021-vs-121 in a bin-depleting regime); the retarget decision is
deferred to that evidence, not to any missing analytic work. **[UPDATED by §4.5:
Rung 5c is complete; the coagulation production path is retargeted to 021 on
symbolic sparsity.]**

### 4.1 Design refinement — recorded-δ (supersedes the static `ν·F/Y_tot` metadata)

The §0 form `∂F/∂Y(GRAIN_k) = ν·F/Y_tot` is correct for the clean sites (accretion
`ν=+1`, LH `ν=−1`, cap-active photodesorption `ν=+1`), but the **H/H₂ accretion
sticking coefficient is not a clean `ν=±1` site**: its rate is
`P · Y_tot · stick(SUMLAY(Y_tot))`, and `SUMLAY ∝ 1/Y_tot` (live, Rung 3), so

```
∂rate/∂Y_tot = rate·(1/Y_tot + slog·∂SUMLAY/∂Y_tot),  slog = ∂ln(stick)/∂SUMLAY
```

which carries a state-dependent correction on top of `ν=+1`. Rather than special-case
it, the implementation records, per in-scope reaction at its rate site (where SUMLAY,
`stick`, the floor, and the cap state are all in scope), the exact
`δ = ∂reaction_rates/∂Y_tot` into `gtodn_jac_dcoef(:)` each RHS evaluation;
`get_jacobian` then deposits `δ × (bilinear part)` into both grain charge columns,
`+w_p` to products and `−δ·bilinear` to reactants. This **subsumes** the ν-sign, the
`max(Y_tot,floor)` floor, the `SUMLAY≥MLAY` cap gate, and the sticking correction in
one place — no per-entry hybrid, no re-derivation of SUMLAY in `get_jacobian`. The
static `gtodn_jac_nu` / `gtodn_jac_floored` metadata is removed.

- `∂SUMLAY(k)/∂Y_tot = −SUMLAY/Y_tot` above the floor **and only where SUMLAY is live**
  (`coagulation` on); frozen SUMLAY (coag off / override) has zero derivative — this
  keeps the override path self-consistent.
- H/H₂ sticking: `slog = (stick_ice−stick_bare)/stick` for `SUMLAY≤1`, else 0
  (flat above 1 ML) — returned by a new out-arg of `sticking_special_cases`.
- LH mantle-multilayer branch (`HAS_MANTLE_COMPOUND .and. sumlaymant>1`): the two
  Y_tot dependences cancel, `δ=0` — coded per the §1 recommendation (v1-inert insurance).

### 4.2 Verification (tests/gtodn_jac_fd_check.*)

Analytic `get_jacobian` vs central-difference FD of `get_temporal_derivatives`, per
grain column, coagulation on. States: A zero ice (accretion), B trace ice on a real
LH reaction (reciprocal), C mild ice cap-off (**H/H₂ sticking active**), D heavy ice
(**photodesorption cap active**). Worst grain-column relative error:

| state | full run (kernel 1e-14) | chemistry-isolated (kernel→0) |
|-------|-------------------------|-------------------------------|
| A accretion         | 3.4e-10 ✔ | 3.4e-10 ✔ |
| B LH reciprocal     | 0.35 ✗ (type-50 only) | 3.4e-10 ✔ |
| C sticking, cap off | 1.03 ✗ (type-50 only) | 1.6e-8 ✔ |
| D cap on            | 3.7e-7 ✔ | 3.7e-7 ✔ |

Isolating the coagulation kernel (so type-50 ice-transport fluxes and their Jacobian
vanish while dynamic GTODN stays live) makes **all four states pass < 1e-5** — this
is the proof that the chemistry-coupling entries close gate 5a and that the full-run
failures are external to Rung 5.

### 4.3 Type-50 ice-transport Jacobian — VERIFIED entry-exact (Rung 5b retracted)

The working hypothesis for §4.3 was a gain-side non-reactant gap analogous to the
Rung-4 live divisor. A read-only diagnosis (tests/ice_transport_jac_diag.*) and an
exhaustive sweep (tests/ice_transport_jac_sweep.*) **refuted it**. `get_jacobian`
already produces the type-50 grain-column gain entries, because ice transport carries
the collision partner `GRAIN_j` as an **explicit reactant** — so the reactant
differentiator visits the reaction and deposits `∂flux/∂Y(GRAIN_j)` to every product
bin, including the redistributed non-reactant bins. (This is exactly what the live
divisor lacked: there the grain entered only through a rate *coefficient*, with no
compound to trigger the column.)

**Exhaustive proof (no finite differencing, so no cancellation).** For every column
and every row, `get_jacobian` isolated to type-50 was compared against an independent
hand-summed per-reaction analytic derivative of the RHS flux, in a light (ice 1e-10)
and heavy (ice 1e-8) state: **393 nonzero entries per state, 786 comparisons, 0
mismatches, agreement at machine epsilon** (max |Δ| 2.0e-28, max rel 2.3e-16). Because
the hand sum loops every type-50 reaction independently, the exact match proves the
loop drops no reaction, selects the correct partner reactant in each `∂flux/∂Y`, and
applies every Podolak/Brauer weight and sign correctly. Ten independent subset-FD
anchors (ice transport and grain-grain, bins 2/3/4, both charge columns) confirm the
analytic model itself against a background-free finite difference.

**Conclusion:** the type-50 (ice-transport + grain-grain) Jacobian is entry-exact.
There is no gap, no sign error, no weight error. **Rung 5b is not needed and is
retracted.**

### 4.4 Validation lesson — the FD noise floor and the two-legged gate

The full-network `021≢022` failures were finite-difference **cancellation artifacts**,
not defects. A coagulation coupling `∂(dY_i/dt)/∂Y(GRAIN_k)` can be a tiny fraction of
`dY_i/dt` when the species' derivative is dominated by other (grain-independent) terms:
J03CH2OH's ice-transport coupling to GRAIN01 is `~6e-6` of its chemistry-dominated
background (`dY/dt ≈ -9.5e-12`), GRAIN02's coupling `~2e-5` of its coagulation
background (`≈ -3.8e-8`). Central- or forward-differencing the RHS then subtracts two
nearly equal large numbers, and the signal is lost below `~ε·(background/signal)`. This
hits **any** FD Jacobian — the 022 oracle and the internal Jacobian under 121 alike.

**General rule (adopted): FD verification has a noise floor at roughly `1e-5` of the
background derivative.** Entries below that floor cannot be validated by differencing
the RHS. The equivalence gate is therefore **two-legged**:
- **above the floor** — FD (022 oracle / gtodn_jac_fd_check) remains the check; this
  covers the Rung-5 chemistry grain-column entries, which are the dominant coupling of
  the species they touch (verified `<1e-5`, §4.2).
- **below the floor** — the **per-reaction-analytic sweep** (ice_transport_jac_sweep)
  is the check: it compares `get_jacobian` against the closed-form derivative of the
  RHS flux and never differences the RHS, so it has no cancellation floor. This is the
  standing gate for the type-50 Jacobian.

The lesson generalizes beyond coagulation: whenever an analytic Jacobian entry is
expected to be much smaller than the row's dominant terms, prefer the analytic
cross-check over FD.


### 4.5 Rung 5c — solver-method retarget to 021 (coagulation path)

**Decision: production coagulation runs on 021 (MOSS=0, symbolic structure +
analytic get_jacobian) when `sparsity = symbolic` (the default); coag-on with
`sparsity = numerical` stays on 121; coag-off stays on 121.**

**Correction to the record (load-bearing).** In this DLSODES, `MF` decodes as
`MOSS = MF/100`, then `METH`, `MITER`. So **121 = MOSS=1, MITER=1** — it uses the
*analytic* `get_jacobian` for values (structure auto-probed), **not** finite
differences. FD is 022 (MITER=2). The earlier characterization of 121 as
"sparsity-only, values internally differenced" was wrong. Consequences:
- The analytic-Jacobian payoff (vs FD) was **always in production** under 121.
  Measured in a bin-depleting regime, 022 (FD) costs ~60x the RHS evaluations and
  ~50x the wall time of 121/021, with a convergence failure; 121 and 021 are the
  fast analytic path.
- Rung 5 was therefore an **analytic-completeness** milestone (the grain-column
  chemistry Jacobian is now entry-exact and gate-verifiable), **not** a
  convergence speedup: pre- vs post-Rung-5 production 121 in the strong regime is
  NST 9923 -> 10021, unchanged within noise, because the completed entries are
  sub-dominant (the same reason FD cannot resolve them, sec 4.4).

**Measured convergence (tests/rung5c_convergence.sh), all `sparsity=symbolic`:**

| regime | mf | Jacobian | NST | NFE | NJE | fails | wall |
|--------|----|----------|-----|-----|-----|-------|------|
| K0=1e-10, 10 Myr | 121 | analytic, auto-probe | 10021 | 13120 | 339 | 0 | 0.60s |
| K0=1e-10, 10 Myr | 021 | analytic, symbolic  | 9973  | 13161 | 341 | 0 | 0.61s |
| K0=1e-10, 10 Myr | 022 | **FD**, symbolic    | 32474 | 802127| 13022| 1 | 30.3s |
| K0=1e-11, 0.1 Myr| 121 | analytic, auto-probe | 6489  | 9014  | 247 | 0 | 0.41s |
| K0=1e-11, 0.1 Myr| 021 | analytic, symbolic  | 6505  | 9002  | 254 | 0 | 0.42s |

021 and 121 **tie** and give **bit-identical** results (max rel diff 0.000). So
the retarget is not chosen on speed.

**Why retarget anyway — the MOSS=1 single-state-probe hazard.** 121 infers the
sparsity structure by probing `get_jacobian`; the driver cold-starts each output
interval (`i_state=1`), so the probe is taken at each interval's *starting* state.
A coupling that is exactly zero at that state but nonzero later is dropped from the
structure, and Newton then runs without it. This project has three couplings that
are zero at plausible probe states: the live divisor before any ice forms and at a
floored bin, photodesorption at `SUMLAY~0`, and the LH reciprocal at `t=0`. The
first interval is ice-free, and any within-interval turn-on is missed for that
interval. The symbolic superset (021) carries every structural entry by
construction, closing the hazard. The `sparsity=symbolic` requirement is the
guarantee, not a burden; it is already the production default.

**Gates cleared before finalizing (both required):**
- (a) 021-on-symbolic completes a full fiducial science run (reduced_CHO, coag on,
  2 Myr) and is **bit-identical to 121** (max rel diff 0.000). Both mf hit 10
  recoverable `istate=-4` error-test-failure restarts on this config — identical
  across mf, so a pre-existing property of the run, not a retarget effect (flagged
  for awareness; worth a separate look at tolerances/first-interval span).
- (b) `equivalence.sh` **PASS** post-retarget (1414 species, 56560 values, worst
  rel diff 0.000): the coag-off / numerical path still reproduces nmgc-2.0 exactly,
  confirming the retarget does not disturb it. Chemistry FD gate and type-50 sweep
  gate both still green.