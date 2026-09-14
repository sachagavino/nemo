# NEMO — Phase II Rung 4 execution brief: complete the Jacobian + scale to fiducial grid

**Git state.** Work on `phase2-coagulation`. Rung 3 is committed there (charge-resolved
coagulation + ice transport + chemistry + charging, one implicit solve, gates 3a/3b/3c
green at the small-grid/slow-kernel config). One commit for this rung, bisectable. Do not
commit to `main`.

**Read first:**
- `docs/PhaseII_rung3_coupling.md` — the coupled model this rung scales. Its **RHS is
  frozen** for Rung 4. Note the `[Correction, Rung 4]` erratum in its §4d: the live-divisor
  Jacobian term was never produced by the analytic routine, which is the core of this rung.
- `docs/PhaseII_execution_brief.md` — the Phase II map.
- `src/ode_solver.f90` — `get_jacobian` (the analytic Jacobian, ~419–575) and
  `set_work_arrays` (builds IA/JA; recon confirmed it already writes the MOSS=0 slots).
- The Stage-4 symbolic `(row,col)` builder (`build_symbolic_sparsity`) — extended here.
- `inputs/parameters.in` — already pins the 20-bin fiducial grid.

---

## 0. What this rung is — and what recon changed about it

The original Rung 4 framing ("append couplings to the symbolic pattern; make symbolic the
production path") was **necessary but not sufficient**, for two reasons the Rung 4 recon
established (both confirmed by the design thread):

- **Finding 1 — the live-divisor Jacobian value does not exist.** With `mf=121`
  (MOSS=1/MITER=1), DLSODES takes both structure and values from the analytic
  `get_jacobian`, which differentiates only the bilinear `Y·Y` factor and treats the rate
  coefficient as constant. It never computes `d(rate)/dY(GRAIN_k)`. Adding a *pattern slot*
  for the live-divisor entry does nothing — a slot with no value is a structural zero.
- **Finding 2 — the symbolic pattern is not wired into the solve.** `mf=121` → MOSS=1, so
  ODEPACK ignores the IA/JA arrays `set_work_arrays` builds and rediscovers structure from
  `get_jacobian` on every restart. That is why numerical and symbolic runs are byte-identical
  today: the flag's only live effect was running the debug assert.

So Rung 4 has **three** jobs, not one:

1. **Wire the symbolic pattern into the solver** (`MOSS=0`), so the complete structure is
   actually used. Recon confirmed the IWORK slots are already written and the LENIW/LRW
   sizing is already generous — `mf` selection is the only change needed to activate it.
2. **Complete the pattern** — add the live-divisor structural entries via a general
   declared-dependency facility (§2), reused by Rung 5.
3. **Supply the missing values by finite difference** (`MITER=2`), which fills every pattern
   entry from the RHS — and the RHS *does* depend on `Y(GRAIN_k)` through SUMLAY, so the
   live-divisor values appear correctly with no analytic-reciprocal derivation.

**This rung does not touch the RHS and does not add analytic reciprocal Jacobian terms.**
Extending `get_jacobian` for the reciprocal `d(rate)/dY(GRAIN_k)` is **Rung 5** (see §1.2).

### 0.1 Why this is a hard blocker, not tidiness

At the physical fiducial config (`initial_dtg_mass_ratio = 1e-2`, 20 bins, coag on), the
incomplete Jacobian causes repeated `ISTATE=-4` error-test failures in **both** `mf=121`
and `mf=022` — because with the incomplete *pattern*, FD cannot fill entries the pattern
omits. So the physical fiducial config **cannot currently run**. Completing the Jacobian is
what unblocks it. This is the load-bearing "before" picture for the gate (§5, gate c).

---

## 1. Decisions locked (do not relitigate)

1. **The fork is resolved as a sequence: Rung 4 = finite-difference, Rung 5 = analytic.**
   Rung 4 completes the *structure* (symbolic, MOSS=0) and fills *values* by FD (MITER=2).
   FD is complete (no missing entries), so it does not violate the paper's objection, which
   is against *missing* entries, not FD ones. Rung 5 then extends `get_jacobian` analytically
   for all reciprocal terms at once (live divisor + the ~23 GTODN sites — same class),
   verifies each against Rung 4's FD path, and switches production to analytic for speed.
   The FD path is retained permanently as the FD-vs-analytic verification oracle.
2. **Analytic reciprocal Jacobian entries are Rung 5, not here.** Rung 4 does structure +
   MOSS=0 wiring + FD values only. No `get_jacobian` extension for the reciprocal term.
3. **`mf` is conditional on the coagulation switch, not global.** coag on → `022`
   (MOSS=0/METH=2/MITER=2 FD); coag off → `121` (MOSS=1/MITER=1 analytic, nmgc-compatible).
   This is semantically correct: coag off ⟹ grains constant ⟹ the live-divisor term is
   identically inactive ⟹ the complete Jacobian is not needed and `121` is both correct and
   sufficient. Implement as one small selector function so Rung 5 can retarget the coag-on
   branch to `021` and add a verification override in one place. **No explicit `mf` parameter
   yet** (it would let a user set coag-on + `121` and resurrect the ISTATE=-4 failures).
4. **`get_jacobian` is kept, dormant under coag-on** (FD fills values there), live under
   coag-off (the nmgc/`121` path). Rung 5 revives it for coag-on.
5. **`equivalence.sh` stays on coag-off → `121`** and remains bit-identical to nmgc-2.0
   (landmine 6). Do not put the coag-on FD path anywhere near the nmgc regression.
6. **Physics (RHS) is frozen.** Same reactions, rates, kernel, charge scheme, live/frozen
   boundary as Rung 3. The answer must not move (gate b).

### 1.2 Scope boundary with Rung 5 (read this before deciding anything looks in-scope)

The reciprocal `d(rate)/dY(GRAIN_k) ∝ −1/Y(GRAIN_k)²` is the mathematical class the paper
files under Rung 5. Rung 4 makes that term reach Newton **by finite difference over a
complete pattern**, which needs no analytic derivation. If you find yourself deriving or
coding an analytic reciprocal into `get_jacobian`, stop — that is Rung 5. Rung 4's job is
to make the *complete pattern real and used*, and let FD supply the values.

---

## 2. The pattern-completion facility (the one design piece)

Coagulation and ice-transport entries are bilinear — derived from reactants, the builder
already emits this class; extend its coverage to the generated pseudo-reactions over both
charge states.

The **live-divisor** entries are non-reactant dependencies: a surface rate `R` reading
`SUMLAY(k)` depends on `Y(GRAIN_k⁰)`+`Y(GRAIN_k⁻)`, which are not reactants of `R`. Build a
**general declared-dependency facility** in `build_symbolic_sparsity`: a reusable
`(row,col)` append merged into the key list before sort/dedup/CSC (grow the pass-1 upper
bound accordingly). Build it reusable — **Rung 5's ~23 GTODN sites are the same facility.**

Populate it for Rung 4 with the **GRAIN_RANK superset rule** (confirmed by the design
thread): for every reaction `r` with `GRAIN_RANK(r)=k>0`, for every compound `s` of `r`,
add `(s, INDGRAIN(k))` and `(s, INDGRAIN_MINUS(k))`. This is a clean superset of the true
SUMLAY-reading set (types 99-H/H2, 66, 67). Superset is the right call — extra structural
zeros are safe, and cheap under FD because those `GRAIN_k` columns are already probed for
coagulation, so no extra RHS evaluations. Never trim toward the exact set; missing entries
are fatal, extra ones are not.

---

## 3. Prerequisites — verify before building

- **P1 — 20-bin expansion + counts.** Confirmed by recon: 741 grain + 7942 ice-transport =
  8683 terms, 20 overflow pairs skipped, no hardcoded bin count. Keep the count assertion
  in the harness.
- **P2 — Rung 3 physics integrates at 20 bins.** Confirmed at toy dtg (~13 s, exit 0).
  Scaling itself is not the problem; the Jacobian is.
- **MOSS=0 sizing.** Confirmed a no-op — `set_work_arrays` already writes the MOSS=0 IWORK
  slots and the sizing formula is generous. Activating MOSS=0 is an `mf` change, not a
  sizing change.

---

## 4. What to build (one bisectable commit)

1. **`mf` selector** keyed on the coagulation switch: coag→`022`, else→`121` (§1, decision 3).
   Single function, structured for Rung 5's override.
2. **Declared-dependency facility** in `build_symbolic_sparsity` + the GRAIN_RANK superset
   population (§2). General and reusable.
3. **Diagnostics at scale:** carry the Rung 3 3a diagnostics (dust-core mass over both
   charges; ice-transport operator net-zero per species; dropped-flux over both charges) to
   20 bins. Add the DLSODES step-stat harness (`nst`/`nfe`/`nje` + error-test and
   convergence failure counts) and wall-clock for RHS / Jacobian / factorization.

---

## 5. The gate (the byte-identity 4-I is retracted — see §0)

**(a) Structural completeness — static.** Assert the live-divisor entries
(`(compound, INDGRAIN(k))` and `(compound, INDGRAIN_MINUS(k))` for the GRAIN_RANK rule) are
present in the symbolic pattern. Direct pattern check, not inferred from run agreement.
Plus the carried-forward subset-assert (numerical ⊆ symbolic) for the reactant-derived
entries, probed at an **early, unfloored** state (§ landmine 4). Together these replace the
retracted byte-identity gate: structural check certifies the live-divisor block, subset
certifies the reactant-derived block.

**(b) Complete-Jacobian correctness at scale.** With `mf=022` + the complete pattern, the
fiducial coupled model **integrates to completion** (the `ISTATE=-4` "before" is gone),
**conserves 3a** (dust-core mass both charges to machine precision; ice-transport operator
net-zero), and **gives the same answer as the prior run to tight rtol** (RHS unchanged, so
the solution must not move). If the answer moves past tight rtol, escalate — it means the
completed pattern changed the physics, which it must not.

**(c) Load-bearing demonstration.** Report DLSODES step stats at the **physical fiducial**
config, showing the complete pattern removes the `ISTATE=-4` failures / reduces step
rejections versus the incomplete pattern — **same answer, better convergence.** The "before"
is already measured (incomplete: `mf=121` 5×ISTATE=-4/27 s, `mf=022` 7×/109 s). The "after"
must converge cleanly. If it does **not** — if the complete FD pattern still fails to
converge at physical dtg — that is the escalation the design thread asked for: the FD-path
bet needs rethinking, stop and report.

**(d) Regression + tractability.** 3b (slow-coag reduction, ≥2 K0, floor-never-triggered)
and 3c (nominal directional; expect more temperature-smear structure at 20 bins — physics,
not failure) still pass at the fiducial grid. Report the FD-Jacobian wall-clock (the 4-IV
upper bound; recon: ~4× the RHS-only cost — analytic Rung 5 will beat it). If FD is
tractable, analytic certainly is; if FD is intractable, escalate before approximating.

---

## 6. Landmine checklist (Rung 4-specific)

1. **A pattern slot is not a Jacobian value.** The live-divisor entry needs its *value* to
   reach Newton; under Rung 4 that value comes from FD (MITER=2), not from `get_jacobian`.
   Do not "add the entry to the pattern" and assume it is live — verify via gate (c).
2. **MOSS=0 is mandatory, not optional.** MOSS=1 rediscovers structure from `get_jacobian`
   at the init state and drops entries that are value-zero at discovery (SUMLAY=0 before ice
   forms; floored bins). That state-dependence is the exact fragility the symbolic pattern
   eliminates. Wire MOSS=0.
3. **`mf` is conditional, never global** (§1, decision 3). A global `022` puts the coag-off nmgc
   regression on FD and breaks bit-identity (landmine 6).
4. **Subset-assert probes an early/unfloored state.** A probe at `t=0` or a floored state
   under-reports (live-divisor value is zero there) and gives a false pass. Union multiple
   probes; the unfloored one is the stringent one.
5. **Superset, not exact match** (§2). Numerical ⊆ symbolic; extra entries safe, missing
   entries fatal. Both charge columns per live-divisor row.
6. **Do not wire symbolic / FD into `equivalence.sh`.** The nmgc bit-identity regression
   needs nmgc's pattern and MITER=1; keep it coag-off → `121`.
7. **Numerical is NOT a valid coag-on reference.** It inherits `get_jacobian`'s
   incompleteness, so a numerical coag-on run is too sparse and FD would silently not fill
   the live-divisor entries. For coag-on, symbolic-complete is the only valid pattern.
   Numerical is retained only for the coag-off nmgc regression.
8. **RHS must not move.** Gate (b)'s tight-rtol match is the guard against smuggling a
   physics change in with the plumbing.
9. **Temperature smear is real at 20 bins**, not a bug — 3c structure and any
   resolution-dependence of the chemistry are physics.

---

## 7. Standing conventions + hand-off

One bisectable commit. Diagnostics read zero-drift before the physics they monitor.
`equivalence.sh` stays `numerical`/coag-off/`121`. Conceptual questions return to the
design thread. Escalate (don't approximate) if gate (c) fails with the complete pattern or
FD is intractable (§5).

**What this rung unlocks (not deliverables):** the resolution-convergence study (chemistry
vs `mass_ratio`, quantifying temperature-smear discretization error) and the
diffusion-vs-DG comparison — science-phase analysis, not gates. A two-point refinement
smoke test is fine as a sanity check; the full study is not this rung.

**Hand-off out of this rung:** green (a)–(d) plus a tractability report unlocks **Rung 5**
(analytic reciprocal Jacobian across the ~23 GTODN sites *and* the live divisor), which
reuses the declared-dependency facility built here and is verified against this rung's FD
path before production switches to `mf=021`.