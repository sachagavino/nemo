# NEMO — Phase II Rung 5 execution brief: dynamic surface-rate GTODN (the stiff reciprocal)

**Git state.** Work on `phase2-coagulation`. Rungs 1–4 committed there. This is the **last
rung of Phase II**; a green Rung 5 makes v1 feature-complete. One bisectable commit. Do not
commit to `main`.

**Read first:**
- `docs/PhaseII_rung4_sparsity_scale.md` — this rung consumes what Rung 4 built: the
  `declare_jacobian_dependency` facility, the complete-pattern machinery, the MOSS=0 wiring,
  the `solver_method_flag` seam / `NEMO_ORACLE_MF` override, and — critically — the
  **FD-complete oracle (`mf=022`)** that verifies this rung's analytic entries.
- `docs/PhaseII_rung3_coupling.md` — the partial-coupling boundary this rung completes.
- The paper's `sec:gtodn` (the drafted "gas-to-dust ratio as a dynamic variable" text) —
  it is effectively this rung's design spec: the two-class distinction, the
  analytic-plus-FD-verification, the divisor floor, and the coag-off acceptance test.
- `src/ode_solver.f90` — `get_jacobian` (extended here for the reciprocal), the surface-rate
  sites (accretion ~1859, LH ~2049, diffusion ~2434, ER/CIR, CR desorption), and the frozen
  `GTODN` array those sites currently read.

---

## 0. What this rung is

Rung 3 made the *monolayer-counting* divisor live. Rung 5 makes **every surface-rate
GTODN(k) live** — the full dynamic substitution across the ~23 sites — and supplies the
**analytic reciprocal Jacobian entries** that this introduces, verified against the Rung 4
FD oracle, then switches production coag-on to the analytic path (`mf=021`).

### 0.1 The two classes (from `sec:gtodn`) — one is done, one is this rung

The grain abundance enters surface quantities two ways, and Rung 5 is only about the second:

- **Benign / bounded — already live since Rung 3.** The monolayer coverage
  `SUMLAY(k) ∝ 1/Y(GRAIN_k)` diverges, but the chemistry uses the *reactive fraction*
  `∝ MLAY/SUMLAY ∝ Y(GRAIN_k)` — the reciprocal inverts back to **linear**, so its Jacobian
  entry is bounded. Done. **Do not re-touch it.**
- **Demanding / stiff — this rung.** The Langmuir–Hinshelwood diffusive rates scale as
  `1/Y(GRAIN_k)` with **no compensating factor** (two adsorbed species meeting on one grain;
  spreading a fixed per-H abundance over more grains dilutes the encounter rate). So
  `∂(rate)/∂Y(GRAIN_k) ∝ −rate/Y(GRAIN_k)²` — the genuinely stiff reciprocal that stiffens
  as coagulation empties a bin. This, and any other ∝1/Y site, is Rung 5.

### 0.15 The site set is settled — three reaction classes, one form (see the classification doc)

The site enumeration is complete and locked in `docs/PhaseII_rung5_site_classification.md`
(the checked artifact). The brief's "~23 sites" was a raw textual-occurrence count; reduced
to rate coefficients that acquire a grain-abundance dependence **in the v1 2-phase config**,
the in-scope set is **three reaction classes under one derivative form**:

- **accretion (type 99): `ν=+1`, linear, no floor** — coefficient `∝ Y_tot`.
- **photodesorption (types 66/67): `ν=+1`, floored via SUMLAY** — routes through the floored
  monolayer count; entry is `+F/Y_tot` **only while the `SUMLAY≥MLAY` cap is active**, zero
  otherwise (piecewise — see landmine 8).
- **LH surface (type 14): `ν=−1`, reciprocal, floored** — the stiff one.

**The one form:** `∂F/∂Y(GRAIN_k⁰) = ∂F/∂Y(GRAIN_k⁻) = ν·F/Y_tot(k)`, where `F` is the flux
(the `tmp_value` `get_jacobian` already forms) and both charge columns are identical (the
divisor is the total). This subsumes the brief's `−rate/Y²` (that is this with `ν=−1`).

**The floor rule (corrected — keyed on the divisor, not on ν):** apply the floor (entry = 0
below floor) **wherever a `max(Y_tot,floor)` divisor appears in the coefficient**, otherwise
no floor. So accretion (bare `×Y_tot`) is unfloored despite being the "safe" sign, LH is
floored, and photodesorption is floored *even though `ν=+1`* because it routes through the
floored SUMLAY. Below the floor `max(Y_tot,floor)` is constant, so anything routed through it
has zero grain-derivative there.

Everything else (ER/CIR, modified-LH, 3-phase, single-grain `GTODN_FIXED`, the `ab_lay` cap)
is out of scope for v1 — cancelling, gated off, frozen by decision, or not per-bin. The
confounders are enforced off by a hard scope guard (decision 6). Full reasoning, line
numbers, and the per-flag decisions are in the classification doc; this brief does not
duplicate the table.

### 0.2 Why the reciprocal cannot come from the reactant differentiator

`get_jacobian` differentiates the bilinear `Y·Y` factor over a reaction's *reactants*,
assuming the rate is *proportional* to each reactant abundance. GRAIN_k enters the LH rate
as `1/Y(GRAIN_k)` — an *inverse*, not a product factor. Even if GRAIN_k appears as a
catalyst in the pseudo-reaction, the reactant differentiator would produce `+rate/Y`
(treating it as a product factor) instead of the correct `−rate/Y²`. So the reciprocal
entry **must be derived analytically and registered as a non-reactant dependency** via the
Rung 4 facility. This is the core of the rung, and the reason the FD oracle exists: a
hand-derived reciprocal is easy to get wrong in sign or form, and FD is the independent
check.

---

## 1. Decisions locked (do not relitigate)

1. **Dynamic GTODN is gated on the coagulation switch.** coag on → live GTODN at all sites +
   analytic reciprocal entries + `mf=021`. coag off → frozen GTODN + no reciprocal entries +
   `mf=121` = **exactly nmgc-2.0, bit-identical**. This is physically exact, not an
   approximation: coag off ⟹ `Y(GRAIN_k)` constant ⟹ dynamic ≡ frozen. Gating preserves the
   nmgc regression by construction (landmine 1). It mirrors the Rung 4 selector; Rung 5 flips
   the coag-on branch from `121` to `021`.
2. **Production coag-on switches to `mf=021`** (MOSS=0, analytic-complete `get_jacobian`).
   `mf=022` (FD-complete) is retained as the permanent verification oracle. `mf=121` remains
   the coag-off/nmgc path.
3. **Per-hydrogen abundances + divisor floor — not the θ formulation.** The occupation-per-
   grain `θ = Y_ice/Y(GRAIN)` reformulation removes the singularity but rewrites every
   surface rate; it is v2, deferred. Rung 5 keeps per-H and floors the divisor. Do not start
   a θ rewrite because a bin gets stiff — that is the temptation the floor exists to resist.
4. **The floor is on the divisor, not the state, and the floor indicator is keyed on the
   divisor, not on ν** (corrected — see §0.15). Apply the floor (entry = 0 below floor)
   wherever a `max(Y_tot,floor)` divisor appears in the coefficient. So LH and photodesorption
   are floored, accretion is not. Below the floor `max(·,floor)` is constant → zero
   grain-derivative (landmine 3).
5. **RHS conservation structure is untouched.** Dynamic GTODN changes rate *magnitudes*, not
   what is created/destroyed. Grain cores are untouched (GTODN is a rate factor), so 3a dust-
   mass conservation is unaffected. Verify, don't assume (gate 5d).
6. **One hard scope guard enforces the v1 confounder-off config.** When coagulation is on,
   **hard-assert (stop, not warn)** `modify_rate_flag==0 .and. is_er_cir==0 .and.
   is_3_phase==0`, with a message naming the unsupported combination. These are the three
   confounders whose GTODN power is branch-dependent or cancelling (classification doc rows
   G–L); v1 does not derive their Jacobians. The guard makes the code fail loudly rather than
   silently mis-compute an entry, and matches the paper's stated v1 scope. The LH entry is
   additionally gated `modify_rate_flag==0` (locally redundant with the guard, kept as
   defence in depth). Modified-LH / ER / 3-phase Jacobians are deferred to v2.
7. **`ab_lay` stays frozen** (a deliberate partial-coupling: a frozen surface-site cap, the
   analogue of Rung 3's frozen surface area, and the v2 hook). Its live variation is
   second-order vs the direct accretion/LH GTODN, and freezing keeps the A/B classification
   clean — but this is **checked, not asserted**: add a diagnostic that the accretion
   monolayer cap seldom binds in the fiducial runs (coagulation reduces total grain surface,
   so a frozen `ab_lay` over-estimates it — harmless only if the cap rarely binds). Escalate
   if it binds often.
8. **`grain_abundance_floor` is promoted to a `parameters.in` knob** (default `2.99e-39`,
   unchanged). One knob shared by the Rung-3 divisor floor and the Rung-5 rate-site floors, so
   they cannot drift apart; it also enables the floor-insensitivity validation run (0.1×/10×)
   without recompiling. It is a physically-motivated ~10 grains/AU³ threshold that results are
   *insensitive* to, not a fitted parameter.

---

## 2. The mechanism

**The classification is done and locked** in `docs/PhaseII_rung5_site_classification.md`
(§0.15). Do not re-derive it. The in-scope set is three reaction classes under the one form
`ν·F/Y_tot`: accretion (99, `ν=+1`, no floor), photodesorption (66/67, `ν=+1`, floored via
SUMLAY, piecewise on the `SUMLAY≥MLAY` cap — landmine 8), LH (14, `ν=−1`, floored). All other
GTODN reads are out of scope (cancelling, confounder-gated, frozen `ab_lay`, or single-grain
`GTODN_FIXED`) — enforced off by the decision-6 scope guard. The one open coding subtlety, the
LH `sumlaymant` cancellation (`ν=0`), is unreachable in v1 (mantle `'K'` species absent from
`reduced_CHO`/fixture, and 3-phase-guarded) but should be mirrored in the analytic entry as
cheap FD-insurance, with a one-line comment pointing at the v2 3-phase work.

**Make the in-scope sites live** (read instantaneous `Y_total(GRAIN_k) = Y(GRAIN_k⁰)+Y(GRAIN_k⁻)`,
floored where a `max(Y_tot,floor)` divisor appears), coag-gated per decision 1. **Centralize the gate into a
single predicate** — one `dynamic_gtodn_active` function that each site queries — rather than
re-testing the coagulation switch at each of the ~23 sites. This is what keeps the rung
cheap: the force-dynamic override (§4.4) and the coag gating then live in *one* place, a
one-line change, instead of being threaded through 23 sites where they would be error-prone
and where a single missed site would silently leave frozen GTODN in a live rate. If you find
yourself inlining the coag test per site, stop and centralize first.

**Register the analytic entries** through `declare_jacobian_dependency` (Rung 4's facility)
and compute their values in `get_jacobian`: `−rate/max(Y,floor)²·[Y>floor]` at reciprocal
sites, the bounded constant at linear sites. Both charge columns (`GRAIN_k⁰`, `GRAIN_k⁻`).

**Verify each against FD as you go, and keep the log.** For every registered entry, compare
analytic `get_jacobian` output against the FD oracle (`mf=022`) at a state set that
**straddles every piecewise boundary the entry has** — this is what catches a half-right
analytic form that passes on one side. The set must include: early/unfloored, mid,
late/floored (for the LH floor); an **ice-rich** state (so the photodesorption 66/67 entries
are exercised at all — they vanish at `SUMLAY≈0` in the ice-free dust-only fiducial); and,
for photodesorption, both a **cap-active** (`SUMLAY≥MLAY`, entry `+F/Y_tot`) and a
**cap-inactive** (entry 0) state. Catches sign errors, `+rate/Y` vs `−rate/Y²`, floor-
derivative mishandling, and the cap discontinuity. **Record each site's FD-agreement result in a
persisted per-site log** (site name, class, worst rel-diff vs FD, states tested), not just a
final pass/fail. Since this rung lands as one big commit (no per-site bisectability), that
log is how a later gate wobble is traced to the offending site without re-bisecting a commit
that no longer bisects. Build site-by-site (or class-by-class: accretion, LH, ER, CR),
verifying each against FD, then land as one bisectable rung commit when all pass.

---

## 3. Prerequisites — confirm inherited from Rung 4

- `declare_jacobian_dependency` present and general (Rung 4 shipped it reused-verbatim-ready).
- Complete symbolic pattern + MOSS=0 wiring consumable.
- FD oracle (`mf=022`) reachable via `NEMO_ORACLE_MF` and gate-verified at fiducial (Rung 4
  gate b). This rung *depends* on it — if the oracle run isn't in the permanent harness, add
  it first.
- `solver_method_flag` coag-on seam isolated (Rung 4 built it as the flip point).

---

## 4. What to build (one bisectable commit)

1. Site enumeration + classification table (reciprocal / linear / flagged), in the harness
   as a checked artifact — future readers must see which sites are which class.
2. **A single `dynamic_gtodn_active` predicate** gating live-vs-frozen GTODN, queried by
   every site — not a per-site coag test (§2). This is the seam the coag gating and the
   force-dynamic override both key on; centralizing it is what keeps both one-liners.
3. Live GTODN reads at all sites via that predicate; divisor floor at reciprocal sites.
4. Analytic non-reactant Jacobian entries via the facility (reciprocal + linear), both charge
   columns, floor-aware derivative.
5. `solver_method_flag`: coag-on → `021`; coag-off → `121`; `022` oracle via override. Plus a
   **force-dynamic override** (analogous to `NEMO_ORACLE_MF`) that flips
   `dynamic_gtodn_active` on with coag off, for the 5b acceptance test.
6. **Per-site FD-agreement log** (site, class, worst rel-diff vs FD, states tested),
   persisted as a checked artifact — the trace that substitutes for per-site bisectability
   under the one-commit choice.
7. Carry 3a diagnostics + step-stat/timing harness forward (dynamic GTODN on, fiducial).

---

## 5. The gate

**(5a) The Jacobian is verified two ways — FD above the noise floor, per-reaction-analytic
below it.** The `021≡022` (analytic-vs-FD) check has a hard resolution limit discovered
during this rung (see `docs/PhaseII_rung5_type50_amendment.md`): finite differencing cannot
resolve any Jacobian entry whose coupling is below ~1e-5 of the affected species' background
derivative — the signal falls under the round-off floor of the RHS subtraction, and the FD
oracle is pure noise there. So the gate has **two legs**:

- **FD leg (`021≡022`), for entries above the FD noise floor** — the Rung-5 chemistry
  reciprocals (accretion, LH, photodesorption, sticking). Analytic (`mf=021`) matches the FD
  oracle (`mf=022`) to tight tolerance at the piecewise-straddling state set of §2:
  early/mid/late including a **floored** state (LH floor derivative is zero, not
  `−rate/floor²`), an **ice-rich** state (or the photodesorption 66/67 entries are never
  exercised), and both a **cap-active and cap-inactive** photodesorption state (the entry is
  `+F/Y_tot` only while `SUMLAY≥MLAY`, zero otherwise). A mismatch here is a derivation
  error — fix the analytic term, do not loosen the tolerance.
- **Per-reaction-analytic leg, for entries below the FD noise floor** — the type-50
  coagulation Jacobian (ice-transport + grain-grain), whose coupling is ~1e-5–1e-6 of the
  background and therefore *unresolvable by FD by construction*. Validate these by
  constructing each entry from an independent per-reaction analytic hand-sum and comparing to
  `get_jacobian` with **no differencing** — machine-epsilon agreement is the bar
  (`tests/ice_transport_jac_sweep.{f90,sh}`, the standing gate leg). Type-50 is already
  verified entry-exact this way (786 comparisons, max rel 2.3e-16); keep the sweep as a
  permanent regression.

The per-site FD-agreement log (§4.6) is the checked artifact for the FD leg; the sweep output
is the artifact for the analytic leg. A later regression in either is traceable to a site
without re-bisecting. Note the general lesson for future entries: choose the leg by the
coupling's magnitude relative to the background derivative — FD only where the signal clears
its noise floor, independent analytic construction where it does not.

**(5b) Coag-off reduction to the frozen reference.** Using the force-dynamic override, run
dynamic GTODN with **coag off**; require agreement with the frozen-GTODN (uncoupled/nmgc)
reference to tight rtol. `Y(GRAIN_k)` is constant here, so the RHS is identical to machine
precision; the tight-rtol (not bit-identical) agreement is expected, arising from the extra
reciprocal Jacobian entries taking a different Newton path to the same solution. This is the
paper's stated acceptance test for the substitution.

**(5c) The load-bearing payoff — finally demonstrable (this is what Rung 4's gate c couldn't
show).** Same dynamic-GTODN physics, three Jacobian treatments at the physical fiducial in a
bin-depleting regime: (i) analytic-complete `mf=021`; (ii) FD-complete `mf=022`; (iii)
analytic with the **reciprocal entries omitted** (a test mode — register the rates live but
not the reciprocal Jacobian). Required outcome: (i) and (ii) converge and agree; (i) is
faster than (ii) (no FD probing); (iii) **degrades markedly or fails** as bins deplete. This
simultaneously proves the reciprocal is load-bearing (the paper's "matters at scale" claim),
the analytic form is correct (matches FD), and analytic is the right production choice
(fastest). If (iii) does **not** degrade, the reciprocal isn't actually stiff at the tested
config — deplete harder, or flag: it would mean the paper's central Jacobian argument needs
qualifying.

**(5d) Conservation at scale, dynamic GTODN on.** 3a holds: dust-core mass over both charges
to machine precision, ice-transport operator net-zero per species, charge conserved. Dynamic
GTODN changes rates, not conservation structure; confirm it.

**(5e) Regression + the production-retarget decision.** 3b/3c pass at fiducial. `equivalence.sh`
stays coag-off/`121` and **byte-identical to nmgc** (landmine 1 — verify explicitly). The
switch of production coag-on from `121` to `021` is **decided on the 5c payoff, not on the
type-50 finding.** Type-50 being entry-exact makes `021` *correct*, but those entries are
sub-noise and dynamically negligible (which is why coagulation ran fine under `121` through
Rung 4) — accuracy there is not a reason to switch. Retarget only if 5c shows the *chemistry
reciprocal* earning a convergence/robustness win at production scale. Run `021` vs `121` on
the fiducial in a bin-depleting regime, confirm trajectories agree (same RHS → same answer)
**and** compare step stats; retarget if `021` wins, otherwise keep `121` in production and
retain `021`/`022` as verification modes (a null result is honest — it means the reciprocal
isn't load-bearing at fiducial conditions, which is then stated rather than a payoff claimed).
Report `mf=021` timing vs the `mf=022` oracle (~118 s) either way.

---

## 6. Landmine checklist (Rung 5-specific)

1. **nmgc bit-identity is most at risk this rung.** Dynamic GTODN must be *fully* gated on
   coag — both the live rate reads **and** the reciprocal Jacobian entries. If the reciprocal
   entries leak into the coag-off `get_jacobian`, coag-off gains entries nmgc lacks and
   `equivalence.sh` breaks. Gate both; verify `equivalence.sh` green as an explicit gate.
2. **The reciprocal is not a reactant derivative** (§0.2). `+rate/Y` (reactant-style) is the
   wrong sign and form; the correct entry is `−rate/Y²`, analytic, non-reactant. The FD
   oracle (5a) catches this — do not skip the FD comparison to "save time."
3. **Floor derivative is zero below the floor.** The analytic entry is
   `−rate/max(Y,floor)²` only while `Y>floor`; at/below floor it is **0** (the floored
   divisor is constant). Getting this wrong gives a spurious large entry exactly as a bin
   empties — the worst place. FD at a floored state (5a) catches it.
4. **Use the class from the locked classification, don't re-derive it.** Accretion is `ν=+1`
   (bounded), LH is `ν=−1` (reciprocal), photodesorption is `ν=+1` but floored-via-SUMLAY.
   Applying `−rate/Y²` to accretion, or forgetting the floor on photodesorption, is wrong. The
   classification doc is the authority; the confounders are guarded off (decision 6).
5. **Both charge columns** per site (`GRAIN_k⁰`, `GRAIN_k⁻`) — the divisor is the total.
6. **Do not start the θ reformulation** (decision 3). Per-H + floor is v1. A stiffening bin
   is the floor's job, not a trigger to rewrite the surface rates.
7. **[CORRECTED] The monolayer-count RHS *read* is done (Rung 3); its Jacobian *entries* were
   never added — add them.** The original landmine said "monolayer count is already done,
   don't re-touch it," conflating the live RHS read (Rung 3) with the Jacobian entry (which
   `get_jacobian` has never produced for any GTODN site). Photodesorption (66/67) routes
   through SUMLAY, and `mf=021` owes those entries to match the FD oracle (flag ②). Correct
   rule: **do not re-touch the SUMLAY RHS read; DO add its Jacobian entries.** Each reaction
   gets its own entry — that is not a double-count.
8. **Photodesorption is piecewise on the cap.** Its entry is `+F/Y_tot` only while the
   `SUMLAY≥MLAY` cap is active, **zero** when inactive (no SUMLAY routing → no grain
   dependence). Test both sides in 5a (§2) — a one-sided test passes a half-wrong entry.
9. **The LH `sumlaymant` cancellation** (`ν=0`) is unreachable in v1 (mantle `'K'` species
   absent from `reduced_CHO`/fixture; 3-phase-guarded) but mirror the RHS branch in the entry
   as cheap FD-insurance, with a one-line comment pointing at the v2 3-phase work.

---

## 7. Standing conventions + hand-off

One bisectable commit. Site classification table is a checked artifact. Verify each analytic
entry against FD before committing. `equivalence.sh` coag-off/`121`/byte-identical. Escalate
(don't guess) on any site that won't classify, or if 5c shows no stiffness where the paper
claims it.

**Hand-off: this closes Phase II. v1 is feature-complete** — pure coagulation, ice
transport, chemistry coupling, complete symbolic Jacobian at fiducial scale, and dynamic
GTODN with analytic FD-verified reciprocal entries. What follows is **not another rung** but
the validation/science phase: the analytic-vs-FD agreement (5a) and the coag-off reduction
(5b) become permanent regressions; the resolution-convergence study and the diffusion-vs-DG
comparison (unlocked at Rung 4) get run; the simple worked example gets built; and the paper
sections that were drafted ahead of the code (`sec:gtodn`'s analytic claim, the `sec:charge`
machine-precision figure) become true and can be finalized. Coordinate the paper pass with
the design thread once 5a–5e are green.