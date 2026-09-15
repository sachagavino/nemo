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
3. **Build the FD-complete path (`MITER=2`) as the Rung 5 verification oracle.** FD fills
   every pattern entry from the RHS, and the RHS *does* depend on `Y(GRAIN_k)` through
   SUMLAY, so the live-divisor values appear correctly with no analytic derivation — an
   independent Jacobian to check Rung 5's analytic entries against. **This is built and
   gate-verified, but is not the production default** (§0.1): production coag-on stays on
   `mf=121`.

**This rung does not touch the RHS and does not add analytic reciprocal Jacobian terms.**
Extending `get_jacobian` for the reciprocal `d(rate)/dY(GRAIN_k)` is **Rung 5** (see §1.2).

### 0.1 What Rung 4 does and does NOT buy (corrected after the gate-c finding)

Rung 4 is **plumbing that pays off at Rung 5**, not a convergence win in itself. This was
established empirically after the plumbing was built (the gate-c experiment, below), and
it corrects an earlier misframing in this brief:

- The physical fiducial config is **not** blocked. `mf=121` (incomplete analytic) runs it
  to 10⁶ yr in ~28 s with ~4 recoverable restarts and the **correct answer**. The
  occasional `ISTATE=-4` are recovered restarts, not a wall.
- **Completing the live-divisor block gives no convergence improvement** at Rung 4, and FD
  is ~4× slower. Measured, nominal `k0=1e-9`, physical dtg, 10⁶ yr (smallest bin depletes
  ~10⁶×, genuinely load-bearing): FD-complete 11 restarts/118 s; FD-incomplete 7/116 s;
  `mf=121` 4/28 s — all same answer.
- **Why:** the Rung 4 live divisor is the *monolayer-counting* site, and it is **floored**
  (`1/max(Y⁰+Y⁻, floor)`). The floor caps it, so `d(SUMLAY)/dY(GRAIN_k) → 0` once a bin
  depletes — a **bounded** entry, not the stiff `−1/Y²` reciprocal. It feeds only bounded
  effects (photodesorption cap, saturating sticking). Every surface-*rate* GTODN is still
  frozen. The stiff reciprocal the paper says "matters at scale" is **Rung 5's dynamic
  GTODN**, not active here. The floor doing its job is exactly why there is no stiffness
  and no payoff at Rung 4 — a feature, not a defect.

So Rung 4's value is: the declared-dependency facility (reused by Rung 5's ~23 sites), the
complete symbolic pattern (needed by Rung 5's analytic path), the MOSS=0 wiring, and the
**FD-complete path built and gate-verified as the Rung 5 verification oracle.** The
convergence payoff arrives with Rung 5's analytic reciprocal, verified against this oracle.

---

## 1. Decisions locked (do not relitigate)

1. **The fork is resolved as a sequence: Rung 4 builds the FD path; Rung 5 makes it
   analytic. Production coag-on stays on `mf=121` until Rung 5.** Rung 4 completes the
   *structure* (symbolic, MOSS=0) and builds the FD-complete path (`mf=022`), gate-verified
   as the Rung 5 oracle — but does **not** make FD the production default, because it is a
   ~4× regression with no correctness or convergence benefit at Rung 4 (the live-divisor
   term is bounded; §0.1). An inexact analytic Jacobian (`mf=121`) still converges to the
   true RHS solution and does so faster and more robustly here. Rung 5 then extends
   `get_jacobian` analytically for all reciprocal terms at once (live divisor + the ~23
   GTODN sites — same class), verifies each against this rung's FD oracle, and switches
   production coag-on to `mf=021`. The FD path is retained permanently as the oracle.
2. **Analytic reciprocal Jacobian entries are Rung 5, not here.** Rung 4 does structure +
   MOSS=0 wiring + the FD-complete oracle only. No `get_jacobian` extension for the
   reciprocal term.
3. **`mf` selector keyed on the coagulation switch; production coag-on = `121` at Rung 4.**
   coag off → `121` (nmgc-compatible). coag on production → `121` for now (correct, fastest,
   most robust; the missing entry is bounded so harmless). The **`022` FD-complete path is
   invoked explicitly by the gate/oracle harness**, not the production default. Implement as
   one small selector so Rung 5 retargets production coag-on to `021` with `022` as the
   oracle, in one place. **No explicit production `mf` parameter** (it would let a user set
   coag-on + a wrong `mf`); the oracle mode is a test-harness selection, not a user knob.
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

1. **`mf` selector** keyed on the coagulation switch (§1, decision 3): coag-off → `121`;
   coag-on production → `121`; the `022` FD-complete path reachable as an explicit
   gate/oracle mode. Single function, structured for Rung 5's retarget to `021`.
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

**(b) Oracle correctness at scale.** With `mf=022` (FD-complete, the oracle mode), the
fiducial coupled model **integrates to completion**, **conserves 3a** (dust-core mass both
charges to machine precision; ice-transport operator net-zero), and **gives the same answer
as the production `mf=121` run to tight rtol** (same RHS, so the solution must not move).
This certifies the FD-complete path as a valid independent Jacobian — the property Rung 5
needs it for. If the `022` and `121` answers diverge past tight rtol, escalate: the
completed pattern would be perturbing the solution, which it must not.

**(c) [RETRACTED — do not gate on this.]** This gate originally required the complete
pattern to *improve convergence* versus the incomplete one. It does not, and should not, at
Rung 4: the live-divisor term is **floored and therefore bounded**, not the stiff `−1/Y²`
reciprocal, so completing it yields no convergence payoff (measured: FD-complete 11
restarts/118 s vs `mf=121` 4/28 s, same answer). That payoff belongs to **Rung 5's**
dynamic-GTODN reciprocal. The correct Rung 4 expectation is the opposite: **do not require
`022` to beat `121`** — it will not, and that is the floor working as designed. Instead,
**report** the step stats (`nst`/`nfe`/`nje` + failure counts) for `mf=121` (production) and
`mf=022` (oracle) side by side, as the documented baseline Rung 5 will improve upon.

**(d) Regression + tractability.** 3b (slow-coag reduction, ≥2 K0, floor-never-triggered)
and 3c (nominal directional; expect more temperature-smear structure at 20 bins — physics,
not failure) still pass at the fiducial grid, **on the production `mf=121` path**. Report
wall-clock for both `mf=121` (~28 s, production) and `mf=022` (~118 s, oracle upper bound);
analytic Rung 5 will beat the FD number. If `mf=121` is tractable (it is), production is
fine; the FD number bounds Rung 5 from above. Escalate only if the `022` oracle fails to
integrate or disagrees with `121` past tight rtol.

---

## 6. Landmine checklist (Rung 4-specific)

1. **A pattern slot is not a Jacobian value.** The live-divisor entry only reaches Newton
   when its *value* is supplied — by FD (`mf=022`, the oracle) or, at Rung 5, analytically.
   Under production `mf=121` the entry is absent, which is *acceptable at Rung 4* because the
   term is bounded (§0.1). Do not assume adding the pattern slot alone makes it live.
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