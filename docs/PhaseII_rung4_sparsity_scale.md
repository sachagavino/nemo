# NEMO — Phase II Rung 4 execution brief: symbolic sparsity + scale to fiducial grid

**Git state.** Work on `phase2-coagulation`. Rung 3 is committed there (charge-resolved
coagulation + ice transport + chemistry + charging, one implicit solve, gates 3a/3b/3c
green). One commit for this rung, bisectable. Do not commit to `main`.

**Read first:**
- `docs/PhaseII_rung3_coupling.md` — the coupled model this rung scales. Its physics is
  **frozen** for Rung 4; this rung changes only how the Jacobian pattern is built and how
  large the grid is.
- `docs/PhaseII_execution_brief.md` — the Phase II map (Rung 4 = "Jacobian/sparsity
  plumbing + scale to fiducial grid").
- The Stage-4 symbolic `(row,col)` builder (the "extensible (row,col) builder,
  `sparsity = numerical|symbolic` flag" from Stage 3–4) — the module this rung extends.
- `inputs/parameters.in` — already pins the 20-bin fiducial grid (`a_min = 5 nm`,
  `a_max = 0.5 µm`, `mass_ratio = 2`).

---

## 0. What this rung is

Two coupled jobs, both plumbing — **no physics changes**:

1. **Complete the symbolic Jacobian pattern.** Rung 3 ran on `sparsity = numerical`
   because the coagulation, ice-transport, and live-divisor couplings were not yet in the
   symbolic `(row,col)` list. Append them, so the symbolic pattern is provably complete,
   and make `symbolic` the production path with `numerical` retained as the verification
   reference.
2. **Scale from the 3–5 bin toy to the 20-bin fiducial grid.** This is where the term
   count reaches its real magnitude (O(10⁴) ice-transport terms), where the RHS hot-loop
   and Jacobian cost are first stressed, and where the per-bin temperature smear becomes
   physically active.

### 0.1 Why symbolic, not just "numerical works fine"

Numerical sparsity determines the pattern by probing the Jacobian at a *state* — so it
is state-dependent and can miss an entry that is near-zero at the probe point but
load-bearing later. The live-divisor entry is exactly this hazard: when a bin is at the
floor, `d SUMLAY(k)/dY(GRAIN_k) = 0`, so a probe taken while that bin is floored (or at
`t=0` before depletion) will not see the entry, yet it becomes stiff precisely as the bin
empties. A missing but load-bearing Jacobian entry degrades the Newton corrector where it
is stiffest. The symbolic pattern is state-independent and complete by construction; that
is the point of this rung.

---

## 1. Decisions locked (do not relitigate)

1. **Physics is frozen.** Same reactions, rates, kernel, charge scheme, live/frozen
   boundary as Rung 3. If a rate or reaction value changes, it is out of scope. The
   primary gate (§5, 4-I) is a bit-level equivalence between symbolic and numerical runs,
   which fails if physics drifts.
2. **Grid = committed fiducial** (20 bins, mass-doubling, 5 nm–0.5 µm). No new grid
   geometry; use `parameters.in` as committed.
3. **After this rung, `symbolic` is the default/production path; `numerical` is retained
   as the verification reference** and as the permanent subset-assert regression (§4c).
4. **`equivalence.sh` stays on `numerical`, coag off.** Do not wire the symbolic pattern
   into the nmgc bit-identity regression — see landmine 6.
5. **Dynamic surface-rate GTODN (the ~23-site reciprocal substitution) is Rung 5, not
   here.** Rung 4 adds only the live-divisor entries already present in the Rung 3 physics.

---

## 2. The one real design decision: how the live-divisor coupling enters the pattern

Most of the append is mechanical — coagulation and ice-transport pseudo-reactions are
bilinear, so their Jacobian entries derive from their reactants the same way two-body
chemistry does, and the existing builder already knows how to emit those. The
**non-mechanical** part is the live divisor.

In Rung 3, the monolayer count `SUMLAY(k)` reads the live total grain abundance
`Y(GRAIN_k⁰) + Y(GRAIN_k⁻)`. Any surface rate `R` that depends on `SUMLAY(k)` therefore
depends on `Y(GRAIN_k⁰)` and `Y(GRAIN_k⁻)` — **grain species that are not formal
reactants of `R`.** The standard builder derives Jacobian columns from a reaction's
reactant list; these entries are *not* in any reactant list. They are **declared
non-reactant dependencies**.

The decision: **the symbolic builder must accept explicit `(row, col)` additions that are
not derived from reactant lists.** The builder is already described as "extensible," so
this is expected to be a supported mode; if it currently only derives from reactants,
extend it to accept declared extra dependencies. For every surface rate `R` that reads
`SUMLAY(k)`, and every species `s` whose `dY_s/dt` includes `R`, emit
`(row = s, col = GRAIN_k⁰)` and `(row = s, col = GRAIN_k⁻)`.

**Build this mechanism generally, because Rung 5 reuses it.** The dynamic-GTODN
reciprocal entries of Rung 5 are the same class — a surface rate depending on a grain
abundance that is not its reactant. A clean "declare a non-reactant Jacobian dependency"
facility built here is exactly what Rung 5 needs at ~23 sites; a one-off hack for the
divisor is not. Design the facility, then use it for the divisor.

Note the columns are **both charge states** (the divisor is `Y⁰ + Y⁻`). Note also that
the structural entry exists whether or not the bin is currently floored — the sparsity
pattern is about which entries *can* be nonzero, not their value at any state — so it goes
in the pattern unconditionally.

---

## 3. Prerequisites — verify before building

**P1 — network expands cleanly to 20 bins.** The `reduced_CHO` network is bin-agnostic
(single `GRAIN0`/`GRAIN⁻`, single `J X` per species; the grid machinery expands per bin).
Confirm the expansion produces the expected 20-bin state vector with no hardcoded bin
count anywhere (no dimension pinned to the Rung 3 3–5 bin toy). Confirm the term counts
match the analytic expectation for 20 bins before trusting anything downstream.

**P2 — Rung 3 coupled model runs at 20 bins on `numerical` first.** Before touching the
symbolic append, confirm the *unchanged* Rung 3 physics integrates to completion on the
fiducial grid under `numerical` sparsity. This isolates "does it scale at all" from "is
the symbolic pattern complete." If it fails to integrate at 20 bins on numerical, that is
a scaling/tractability problem (→ §5 4-IV, possibly a design escalation), not a sparsity
problem, and must be understood before the append.

---

## 4. What to build

### 4a. Append coagulation + ice-transport entries to the symbolic pattern
Mechanical: these are bilinear, entries derive from reactants over both charge states.
The builder already emits this class for two-body chemistry; extend its coverage to the
generated coagulation and ice-transport pseudo-reactions.

### 4b. Append the live-divisor non-reactant entries
Per §2: the declared-dependency facility, then the `(s, GRAIN_k⁰)` / `(s, GRAIN_k⁻)`
entries for every rate reading `SUMLAY(k)`. This is the delicate part; build the facility
general for Rung 5 reuse.

### 4c. The subset-assert (permanent regression)
Compare the `numerical`-discovered pattern against the `symbolic` pattern and assert
**numerical ⊆ symbolic** — every entry numerical finds must be in the symbolic list.
Extra symbolic entries are safe (they are structural zeros; they cost a little
performance but never correctness). Missing symbolic entries are fatal.

**Probe numerical at multiple states, not one.** A single probe under-reports: probe
early (all bins above floor — the maximal pattern), and again late (some bins floored),
and take the **union** of numerical patterns for the assert. The early/unfloored probe is
the stringent one, because that is when `d SUMLAY/dY(GRAIN_k)` is nonzero and the
live-divisor entries actually appear numerically.

### 4d. Diagnostics carried forward at scale
The Rung 3 3a diagnostics (dust-core mass over both charges; ice-transport operator
net-zero per species; dropped-flux over both charges) must run at 20 bins. Add a
timing/characterization harness: wall-clock for RHS evaluation, Jacobian assembly, and
the sparse factorization, so §5 4-IV can be reported.

---

## 5. The gate

**4-I — symbolic ≡ numerical (the primary test).** Run the *same* coupled model on the
fiducial grid twice: once `sparsity = symbolic`, once `numerical`. Because the physics is
identical and a complete symbolic pattern contains every real Jacobian entry, the two
runs must agree. Ideal outcome is bit-identical; if the sparse-solver reordering differs
between the two pattern sources, agreement to a tight relative tolerance is acceptable
**only if** the discrepancy behaves like reordering roundoff (grows ~machine epsilon,
stays bounded) and **not** like a missing-entry drift (O(1) divergence in stiff regions,
worsening as bins deplete). A missing-entry signature fails the gate — go find the entry.
This test plus the 4c subset-assert together certify completeness two independent ways
(structural and dynamical).

**4-II — conservation at scale.** 3a holds at 20 bins: dust-core mass over both charges
to machine precision, ice-transport operator net-zero per species, charge conserved. No
new tolerance relaxation relative to Rung 3.

**4-III — coupling regression at scale.** Re-run 3b (slow-coag reduction, K0 scaling →
rtol ~ t_run/t_coag, ≥2 K0 values, floor-never-triggered) and 3c (nominal directional:
mass and ice migrate to larger bins) on the fiducial grid. The physics is unchanged, so
these should still pass; this is a regression that the scale-up did not break the
coupling. Expect 3c to show **more** structure than at 3–5 bins — the temperature smear
is now physically active across 20 bins — which is physics, not a failure.

**4-IV — tractability (characterize and report; escalate, don't paper over).** Report the
wall-clock breakdown (RHS, Jacobian, factorization) at 20 bins. This rung is the first
real test of the paper's claim that the reduced network keeps O(N²·N_ice) tractable as
generated reactions. If it is tractable, the numbers go in the paper's performance
discussion. **If it is not tractable, that is a design escalation back to this thread —
not a licence to start approximating the coupling.** The RHS hot-loop being slow is a
known, documented caveat; quantify it here.

---

## 6. Landmine checklist (Rung 4-specific)

1. **Live-divisor entries are non-reactant dependencies** (§2). The builder must accept
   declared `(row,col)` additions, not only reactant-derived ones. Build the facility
   general — Rung 5 reuses it.
2. **Both charge columns per live-divisor row** — the divisor is `Y⁰ + Y⁻`, so each such
   row couples to both `GRAIN_k⁰` and `GRAIN_k⁻`.
3. **Structural entries go in the pattern regardless of floor state.** Do not omit the
   live-divisor entry on the grounds that it is zero when floored — the pattern is about
   what *can* be nonzero. Omitting it reproduces exactly the numerical-probe hazard this
   rung exists to eliminate.
4. **Subset-assert must probe multiple states** (§4c). A single probe (especially at
   `t=0` or a floored state) under-reports and gives a false pass. Union the patterns;
   the early/unfloored probe is the stringent one.
5. **Superset, not exact match.** The requirement is numerical ⊆ symbolic. Do not trim
   symbolic to match numerical exactly — err toward extra (safe) entries over missing
   (fatal) ones.
6. **Do not wire symbolic into `equivalence.sh`.** The nmgc bit-identity regression needs
   the pattern that matches nmgc; the symbolic pattern lists coagulation entries (as
   structural zeros when coag is off) that can change the sparse LU ordering and break
   bit-identity with nmgc. Keep `equivalence.sh` on `numerical`, coag off. The symbolic
   path is a separate production/verification track.
7. **Physics must not move** (landmine's landmine). 4-I is bit-level precisely to catch
   an accidental physics change smuggled in with the plumbing. If 4-I shows a
   missing-entry-style drift, the first suspect is an incomplete pattern; the second is an
   unintended edit to a rate. Do not relax 4-I's tolerance to make it pass.
8. **Temperature smear is now real, not a bug.** At 20 bins adjacent bins differ in Td and
   desorption differs accordingly; 3c structure and any resolution-dependence of the
   chemistry are physics. Do not "fix" them.

---

## 7. Standing conventions + hand-off

One commit for the rung, bisectable. Diagnostics read zero-drift before the physics they
monitor. `equivalence.sh` stays `numerical`/coag-off (landmine 6). Conceptual questions
return to the design thread.

**What this rung unlocks (not Rung 4 deliverables, do not scope-creep into them):** with a
complete symbolic pattern and a working fiducial grid, the resolution-convergence study
(chemistry vs `mass_ratio`, quantifying the temperature-smear discretization error) and
the numerical-diffusion-vs-DG comparison become possible — these are the paper's
convergence figures and are science-phase analysis, not execution gates. A quick two-point
refinement smoke test (fiducial vs one finer grid, confirming refinement runs and results
move in the expected direction) is fine to include as a sanity check; the full study is
not this rung.

**Hand-off out of this rung:** green 4-I/4-II/4-III plus a tractability report unlocks
**Rung 5** (dynamic surface-rate GTODN across the ~23 sites), which reuses the
non-reactant-dependency facility built here for its reciprocal `−1/Y(GRAIN_k)²` entries.
