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
4. **The floor is on the divisor, not the state** (unchanged from Rung 3). At every ∝1/Y
   site the divisor is `max(Y(GRAIN_k⁰)+Y(GRAIN_k⁻), floor)`; the analytic derivative of
   `max(·,floor)` is **zero below the floor**, so the reciprocal Jacobian entry is zero when
   floored (landmine 3). Same `grain_abundance_floor`.
5. **RHS conservation structure is untouched.** Dynamic GTODN changes rate *magnitudes*, not
   what is created/destroyed. Grain cores are untouched (GTODN is a rate factor), so 3a dust-
   mass conservation is unaffected. Verify, don't assume (gate 5d).

---

## 2. The mechanism

**Enumerate and classify every GTODN site.** Find all sites reading the frozen `GTODN`
array (~23 per the paper: accretion, LH diffusive, ER/CIR, cosmic-ray desorption; the
monolayer count is already live and excluded). For each, classify how `Y(GRAIN_k)` enters:
- **reciprocal (∝ 1/Y)** — e.g. LH diffusive. Needs the floored divisor and the analytic
  `−rate/Y²` non-reactant entry. Stiff.
- **linear (∝ Y)** — e.g. accretion (more grains ⇒ more surface ⇒ faster accretion). Needs
  the live read and a **bounded** analytic non-reactant entry (constant in Y), but no floor
  and no stiffness.
- **anything else** — if a site's grain-abundance dependence is neither cleanly ∝Y nor ∝1/Y
  (an unusual ER or CR-desorption form), **flag it to the design thread** rather than
  guessing its derivative. Do not improvise a Jacobian for a form you can't cleanly classify.

**Make all sites live** (read instantaneous `Y_total(GRAIN_k) = Y(GRAIN_k⁰)+Y(GRAIN_k⁻)`,
floored at the reciprocal sites), coag-gated per decision 1. **Centralize the gate into a
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
analytic `get_jacobian` output against the FD oracle (`mf=022`) at several states
(early/unfloored, mid, late/floored) — this catches sign errors, the `+rate/Y` vs `−rate/Y²`
confusion, and floor-derivative mishandling. **Record each site's FD-agreement result in a
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

**(5a) Analytic ≡ FD — the correctness gate for every reciprocal.** For each registered
entry, analytic (`mf=021`) matches the FD oracle (`mf=022`) to tight tolerance at
early/mid/late states, including at least one **floored** state (confirms the floor
derivative is zero, not `−rate/floor²`). This is what the Rung 4 oracle was built for. The
per-site FD-agreement log (§4.6) is the checked artifact for this gate — it must show every
site passing, so a later regression is traceable to a site without re-bisecting. A mismatch
is a derivation error — fix the analytic term, do not loosen the tolerance.

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

**(5e) Regression + production switch.** 3b/3c pass at fiducial on the new production path
(`mf=021`). `equivalence.sh` stays coag-off/`121` and **byte-identical to nmgc** (landmine 1
— the gating is what preserves this; verify it explicitly, since this rung is the one most
likely to break it). Report `mf=021` production timing vs the `mf=022` oracle (021 should
beat the ~118 s FD number — the payoff over Rung 4).

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
4. **Classify before differentiating.** Linear sites (accretion) get a bounded entry, not a
   reciprocal one. Applying `−rate/Y²` to a ∝Y site is wrong. If a site won't classify
   cleanly, flag it, don't guess.
5. **Both charge columns** per site (`GRAIN_k⁰`, `GRAIN_k⁻`) — the divisor is the total.
6. **Do not start the θ reformulation** (decision 3). Per-H + floor is v1. A stiffening bin
   is the floor's job, not a trigger to rewrite the surface rates.
7. **Monolayer count is already done** (Rung 3, bounded). Don't re-touch it or double-count
   its entry.

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
