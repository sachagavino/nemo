# NEMO — Phase II Rung 3 execution brief: couple chemistry with coagulation

**Git state.** Work on `phase2-coagulation`. Rungs 1 and 2 are committed there
(pure coagulation; ice transport). One commit for this rung, kept bisectable. Do
not commit to `main` (still pinned at `v0.4-skeleton`).

**Read first (committed artifacts this rung stands on):**
- `docs/PhaseII_execution_brief.md` — the Phase II map. **NB: its Rung 3 section is
  partially superseded by this brief** (see §0.1). Everything else in it still holds.
- `src/dust/dustevolution.f90` — the Rung 1/2 reaction-generation machinery you extend.
- `src/dust/dust_grid.f90`, `src/environment.f90` — grid + `td_of_a` (temperature).
- `networks/reduced_CHO/` — the network this rung runs on (see §3, P1).
- `docs/PhaseII_rung2_icetransport.md` — the ice-transport conservation gate you carry forward.

---

## 0. What this rung is

Turn the gas-grain chemistry back on and advance **coagulation + ice transport +
surface/gas chemistry in one implicit DLSODES solve** over a single state vector,
on the reduced C/H/O network, on a **small grid (3–5 bins)**. This is the first
genuinely coupled run. Scaling to the 20-bin fiducial and the symbolic-sparsity
plumbing are **Rung 4**, not here.

Rung 3 is a **deliberately partial** coupling. The grain population and the ice
loading are live; the *surface-reaction rates* still read a frozen grain surface
area. That inconsistency is intentional and is exactly what the gate exploits
(§5). Full dynamic GTODN across all surface-rate sites is **Rung 5**.

### 0.1 What this brief supersedes in the main Phase II brief

The main brief's Rung 3 bullet says "do NOT touch GTODN" and gates on "differs from
a chemistry-only run." Both are refined here:

- **GTODN is touched, surgically.** The monolayer-count *divisor* goes live and
  floored (§4c). The surface-*rate* GTODN stays frozen. This is a 1–2 site change,
  not the ~23-site Rung 5 sweep.
- **The gate is three sub-gates** (§5), not a single "differs" check. The "differs
  in the expected direction" check is one of the three (3c).

---

## 1. Decisions locked in the design thread (fold in; do not relitigate)

1. **Partial coupling.** Surface-reaction rates (accretion, Langmuir–Hinshelwood,
   Eley–Rideal/CIR) keep reading the **frozen** `GTODN(k)` array they already use in
   the uncoupled code — no change. Only the **monolayer/coverage count** switches to
   the **live** instantaneous grain abundance. Consequence: the chemistry sees
   correct ice loading and a stale surface area. Valid only when coagulation is slow
   over the run window — which is what gate 3b tests.

2. **Floor the DIVISOR, not the state.** The monolayer count `Y_ice / Y(GRAIN_k)`
   diverges as coagulation drains a bin. Floor the **denominator** at
   `grain_abundance_floor` (already in `parameters.in` = `2.990E-39` grains/cm³,
   ~10/AU³; converted to per-H in code). **Never floor the grain state variable** —
   that injects mass and breaks the dust-mass conservation gate.

3. **Grain charge — coagulate the TOTAL population; charging active.** Neutral-only
   coagulation is not physical and is not used. See §2 for the full charge scheme.

4. **Temperature: derived grid + parametrized `Td(a) ∝ a^(−1/6)`** via `td_of_a`
   (`environment.f90`). Not tabulated. The committed `parameters.in`
   (`dust_grid_source = derived`, `grain_temperature_type = fixed_to_dust_size`)
   already selects this path. **It has never been exercised where Td matters** —
   Rung 3 is its first load-bearing use, so verify it (P3).

5. **Reference curve = coag-OFF, charging-ON, chemistry-ON** reduced-network run.
   This is regime 1 (no coagulation), and it is the curve gate 3b must reduce to.
   Confirm it is green and saved labelled before comparing (P2).

6. **Grid:** mass-doubling (`mass_ratio = 2`), `a_min = 5 nm`, `a_max` set to yield
   **3–5 bins** (a very narrow size range is fine — this is a coupling proof, not a
   growth study). Do not run 20 bins here.

---

## 2. The grain-charge scheme (new v1 decision — implement exactly)

The code carries per-bin neutral and singly-negative grain states,
`INDGRAIN(k)` and `INDGRAIN_MINUS(k)` (init: all neutral, `GRAIN_k⁻ = 0`).

**Coagulate the total bin population**, both charge states. Charge is conserved on
every merge. With only `{0, −}` states available:

```
GRAIN_i⁰ + GRAIN_j⁰  →  w1 GRAIN_k1⁰ + w2 GRAIN_k2⁰
GRAIN_i⁰ + GRAIN_j⁻  →  w1 GRAIN_k1⁻ + w2 GRAIN_k2⁻            (product charge −1)
GRAIN_i⁻ + GRAIN_j⁻  →  w1 GRAIN_k1⁻ + w2 GRAIN_k2⁻  +  e⁻     (shed one electron)
```

- `(k1, k2, w1, w2)` come from the existing `coag_pair_target` map — **charge does
  not change the mass bracketing** (an electron's mass is negligible). Reuse it.
- The `(−)+(−)` product would be charge −2, which has no state; represent it as a
  singly-negative grain **plus a free electron**, so total charge is conserved
  exactly. The `e⁻` is a genuine product of that reaction.
- Rate coefficients: the **same** `coag_kernel` / `coag_kernel_bare` as Rung 1, with
  the diagonal-½ rule applied per charge-state self-pair exactly as today. The
  kernel is **charge-blind** (no Coulomb factor) — a v1 approximation, documented in
  the paper's charge caveat (τ = a k_B T/e² bounds it; error ≲10% only in the top
  bins, and Rung 3's grid never reaches them). Do not add Coulomb focusing.
- This multiplies the **grain** reaction count by ~4 (the charge-combo multiplicity).
  Still O(10⁴) at fiducial scale; trivial at 3–5 bins.

**Ice transport with charged partners.** Ice species `J_k X` are charge-agnostic
(ice sits on the core regardless of its charge). The transported-ice rate depends on
the **total** partner population, so generate the partner over both charges:

```
J_i X + GRAIN_j⁰  →  w1 J_k1 X + w2 J_k2 X + GRAIN_j⁰   (catalyst)
J_i X + GRAIN_j⁻  →  w1 J_k1 X + w2 J_k2 X + GRAIN_j⁻   (catalyst)
```

- Same `(w1, w2)` as the grain reaction (ice-follows-particle). Bare kernel on the
  diagonal, as in Rung 2. This doubles the Rung-2 ice-transport count (partner
  charge only — the ice-bearing side stays pooled in `J_i X`).

---

## 3. Prerequisites — verify/land BEFORE the coupling commit

**P0 — config-surface cleanup (recommended, separable).** Make `dust_grid_source`
the single grid switch; when `derived`, expose one temperature sub-switch
`{gas, size_scaled}` (`size_scaled` → `td_of_a`; `gas` → `get_grain_temperature_gas`).
When `tabulated`, Td comes from the table column, no sub-switch. Delete the free
top-level `grain_temperature_type` and its stale default `'fixed'` (which is not a
case in the `select` in `gasgrain.f90` and would hit the error branch). This can be
its own small prep commit; if you defer it, Rung 3 still runs as long as
`parameters.in` explicitly pins `derived` + `fixed_to_dust_size` (it currently does).

**P1 — confirm per-bin `GRAIN_k⁻` actually populates.** The grain charge-exchange
reactions live in `gas_reactions.in`, written against a *single* `GRAIN0`/`GRAIN⁻`
pair — this is the correct bin-agnostic network convention (the grid machinery
expands the network per bin, exactly as for surface reactions and accretion; the
network never names bins). The per-bin negative *species* exist (`INDGRAIN_MINUS(k)`,
init 0). **What the convention does not by itself guarantee, and what you must check,
is that the expansion writes a charging reaction into each per-bin `GRAIN_k⁻` and
scales its rate by grain cross-section.** If the single reaction fails to expand
per-bin, `GRAIN_k⁻` stays at zero forever, "coagulate the total population" (§2) is
silently identical to neutral-only, and the rung would "validate" a dormant charge
scheme.

This check is free — it rides on the P2 reference run. On that run, dump the per-bin
`GRAIN_k⁻` abundances and confirm:
- **non-zero** at some point (charging fires at all), and
- **varies across bins** in the right direction (larger grains more charged — the
  signature that the per-bin rate was scaled by cross-section).

Non-zero and bin-varying → P1 closed. All-zero → charging is inert; fix the
network/expansion **before** this rung, not during it (loop back to the design
thread). All-equal-across-bins despite differing sizes → charging populates but the
rate is not size-scaled — a physics bug worth surfacing even though the scheme would
still "run."

**P2 — generate and save the gate-3b reference (it does not exist yet).** This is
**not** `tests/equivalence.sh`. That script is the Stage-2 nmgc bit-identity
regression: 2-grain *tabulated* grid (`reference_0D_2grains`), whatever network that
case ships, coagulation off, `cmp` demanding worst-relative-diff exactly 0.0. It
proves the chemistry extraction is faithful (good — the chemistry engine under this
rung is trustworthy) but it is on the wrong grid, wrong network, and produces no
saved chemistry curve on the Rung 3 config.

The gate-3b reference is a distinct artifact: take the exact Rung 3 `parameters.in`
(derived 3–5 bin grid, `reduced_CHO`, charging on, a⁻¹ᐟ⁶ Td), set `coagulation = 0`,
run to the Rung 3 stop time, and **save the output labelled**. This gives 3b an
apples-to-apples comparison — coag-off vs coag-on-slow on identical grid/network/
temperature, differing only in the kernel. (This is the "chemistry-only reference"
from the design handoff, now explicitly charging-on so 3b compares like with like.)
Run P1's `GRAIN_k⁻` diagnostic on this same run.

**P3 — Td wiring (confirmed by design thread; keep the check in the loop).** The
executor has verified `td_of_a` fills `grain_temp(:)` correctly on the derived path —
implemented, wired, physically correct. Retain the one-line assertion (dump
`grain_temp(:)`, check against
`reference_dust_temperature · (a/reference_grain_radius)^(−1/6)`) in the Rung 3
harness so a future grid or parameter change can't silently break the first
load-bearing use of dust temperature.

---

## 4. What to build

### 4a. Charge-resolved grain coagulation
Extend the Rung-1 generation loops in `dustevolution.f90` to iterate over
`(bin, charge)` per reactant, emitting the three product-charge rules of §2,
including the `e⁻`-shedding product on the `(−)+(−)` channel. Reuse `coag_pair_target`
and the kernel functions unchanged.

### 4b. Ice transport over charged partners
Extend the Rung-2 ice-transport generation to emit the partner `GRAIN_j` in both
charge states (§2). Ice species stay charge-agnostic; do not create `J_k X⁻`.

### 4c. Live, floored monolayer divisor
Locate the monolayer/coverage computation (the site that divides total bin-k ice by
bin-k grain number to get monolayers per grain). Switch its divisor from the frozen
`GTODN(k)` to the **live** total grain abundance `Y(GRAIN_k⁰) + Y(GRAIN_k⁻)`, floored
at `grain_abundance_floor` (converted to per-H). **Leave every surface-*rate* site on
the frozen `GTODN(k)` array.** If you find yourself editing accretion / LH / ER-CIR
rate expressions, stop — that is Rung 5.

### 4d. Sparsity: run `sparsity = numerical` for this rung
The coagulation, ice-transport, and new live-divisor couplings are **not** in the
symbolic pattern yet (that append is Rung 4). Running `symbolic` here would give a
too-sparse Jacobian and silently degrade Newton. Run **numerical** sparsity for Rung
3 (it discovers the full pattern, including the new `d(rate)/dY(GRAIN_k)` from the
live divisor). The floor keeps that reciprocal entry finite. Symbolic-append + the
`subset` assert is Rung 4's job.

### 4e. Diagnostics — land these BEFORE turning the physics on
Each must read zero-drift on the reference before coupling adds complexity.
- **Dust-core mass:** `Σ_k m_k · [Y(GRAIN_k⁰) + Y(GRAIN_k⁻)]`. **Summed over both
  charge states** (see landmine 1). Chemistry does not touch refractory cores, so
  this is machine-precision conserved even with chemistry on.
- **Ice-transport operator net-zero:** instrument the coagulation-generated
  ice-transport terms only, and assert their net contribution to `Σ_bins Y(J_k X)`
  is zero per species X, every step. (Total ice is **not** conserved with chemistry
  on — surface reactions form/destroy it — so this is an operator-level check, not an
  end-to-end one.)
- **Dropped-flux (overflow) guard:** carry forward the Rung-1 diagnostic, now summing
  over the **total** bin population, both charges.

---

## 5. The gate — three sub-gates (all must pass)

Vocabulary (keep these distinct; do not conflate):
- **Regime 1 / no-coag:** `coagulation = 0`. The reference (P2). Already the Rung 1/2
  equivalence case.
- **Regime 2 / slow-coag:** coagulation ON but `t_coag = 1/(K·N) ≫ t_run`.
  *This is the new gate and does not exist yet.*
- **Regime 3 / nominal-coag:** coagulation acts significantly over the run.
  (The Rung-1 1000 yr / K0=1e-9 analytic test is regime 3 — it is **not** the slow
  limit; do not cite it as such.)

**3a — conservation with chemistry on.**
Dust-core mass conserved to machine precision (summed over both charges).
Ice-transport operator net-zero per species (§4d). Both hold with reactions active.

**3b — slow-coag reduction (the primary acceptance test).**
Scale the kernel into the slow regime (`t_coag ≫ t_run`) and require the coupled run
to match the regime-1 reference to **rtol ~ t_run / t_coag**. Make "slow" a number
you dial: scale `K0` down 10× and the match must tighten ~10×. Reach the limit by
**kernel scaling**, not by shortening the run (which couples it to output cadence).
Note this simultaneously tests the *whole apparatus* collapsing — live divisor,
frozen rate-GTODN, ice transport, charge-total coagulation all go quiescent as grains
stop moving. So a 3b failure is coarse (it does not localise which piece broke);
if it fires, bisect by disabling ice transport, then the live divisor, in turn.

**3c — nominal-coag difference (signed).**
At the nominal kernel, the coupled run differs from the reference **in the expected
direction**: dust-core mass moves to larger bins (grains grow) and surface ice moves
with it (population of larger-bin `J_k X` rises). A *directional* check, not just
"differs."

---

## 6. Landmine checklist (Rung 3-specific)

1. **The mass diagnostic must sum over charge.** Charging moves grains between
   `GRAIN_k⁰` and `GRAIN_k⁻` at fixed `m_k`. A neutral-only mass sum will read that
   as spurious coagulation drift and fail 3a for the wrong reason.
2. **Floor the divisor, not the state** (§1.2). The Rung-1 dust-mass diagnostic
   catches state-flooring as injected mass.
3. **Live divisor only at the monolayer site.** Surface-rate GTODN stays frozen. Do
   not drift into Rung 5.
4. **`e⁻` shedding conserves charge** — verify total charge balance holds across the
   `(−)+(−)` channel; a dropped electron there is a charge leak the conservation
   checks on dust mass will *not* see.
5. **`sparsity = numerical` for this rung** (§4d). Symbolic here = too-sparse Jacobian
   = silent Newton degradation.
6. **Td wiring** (P3): the assertion stays in the harness even though the parametrization
   is confirmed correct, so a later grid/parameter change can't break it silently.
7. **The gate-3b reference is not `equivalence.sh`** (P2). `equivalence.sh` is the
   nmgc bit-identity regression (2-grain tabulated, coag off, `cmp`). The 3b reference
   is a separate saved run: coag-OFF, charging-ON, chemistry-ON, on the Rung 3
   derived-grid `reduced_CHO` config. Do not conflate them, or 3b compares against the
   wrong curve on the wrong grid.
8. **The charge scheme can be silently inert** (P1). Per-bin `GRAIN_k⁻` species exist
   but a bin-agnostic charging reaction only populates them if the expansion writes it
   per bin. Confirm `GRAIN_k⁻` is non-zero and bin-varying on the P2 run before
   trusting any charge-related result.
9. **Charge-blind kernel is a documented approximation**, not a bug — but confirm the
   kernel functions are called identically for all charge combinations, so the
   approximation is uniform and not accidentally applied to only some channels.

---

## 7. Standing conventions (unchanged)

One commit for the rung, bisectable. Diagnostics exist and read zero-drift before the
physics they monitor. `equivalence.sh` stays pinned to `sparsity=numerical` for the
external nmgc regression; coagulation-off runs remain the reference. Keep coagulation
behind its enable flag so the regime-1 reference stays reproducible. Conceptual
questions go back to the design thread, not decided here.

**Hand-off out of this rung:** a green 3a/3b/3c on the small grid unlocks Rung 4
(symbolic-sparsity append + subset assert + scale to 20 bins), where the temperature
smear becomes quantitatively testable and the convergence figure gets its fiducial-grid
numbers.
