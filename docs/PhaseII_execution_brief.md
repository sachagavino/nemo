# NEMO — Phase II execution brief: coagulation coupling

**Git state:** all Stage 1–4 work is now on **`main`**, whose tip is tagged
**`v0.4-skeleton`** (commit `8f53b1a`, "Stage 4 Part 2: symbolic Jacobian
sparsity") — the last bit-identical commit before any coagulation physics. Phase II
work happens on the branch **`phase2-coagulation`**, already created from `main`.
Check it out with `git checkout phase2-coagulation` before starting. The
`v0.4-skeleton` tag remains the clean fallback point.

**How Phase II proceeds:** every rung commits on `phase2-coagulation` — one commit
per rung, so the ladder stays bisectable. Do not commit rung work to `main`;
`main` stays at the verified skeleton until a working coupled run is ready to merge
back.

**Stage 4 state:** symbolic sparsity + index maps complete; `sparsity=numerical|symbolic`
flag, default symbolic. **Read first:** `docs/EXTRACTION.md`, `src/dust/README.md`,
`docs/STAGE4_SPARSITY.md`.

This brief is the *map* for all of Phase II. You will be handed **one rung at a
time** as its own thread — do not implement rungs ahead of the one you are given.
Each rung ends with a commit and a green diagnostic; the next rung does not start
until the current diagnostic reads green. Conceptual questions go back to the
design thread, not decided here.

---

## 0. Notation and the physics→code translation (READ THIS FIRST)

The code works in abundances `Y_i = n_i / nH`, with `nH` constant (0D). The
existing two-body RHS builds, for a reaction `r1 + r2` with rate coefficient `R`:

```
contribution to dY(product)/dt  =  R * Y(r1) * Y(r2) * nH
```

Coagulation and ice transport are added as **generated pseudo-reactions** that
flow through this exact machinery — no new solver code, no dedicated dust
integrator. The whole task is: generate the right reactions with the right `R`,
and append their Jacobian couplings to the symbolic sparsity list.

**The discretized Smoluchowski equation, in code variables** (this is Eq. eq:smolu
of the paper; the target the RHS must reproduce):

```
dY_k/dt = (1/2) * sum_{i,j} M_ijk * K_ij * nH * Y_i * Y_j        [gain]
        -        Y_k * sum_j   K_kj * nH * Y_j                    [loss]
```

- `Y_k`      = abundance of GRAIN_k  ( = 1/GTODN(k) )
- `K_ij`     = coagulation kernel (cm^3 s^-1), see rung 1
- `M_ijk`    = Podolak/Brauer redistribution coefficient (dimensionless),
               frozen at init; nonzero only for the two bins bracketing m_i+m_j

### Reaction encoding: per-product weight field (DESIGN DECISION — supersedes earlier notation)

The current reaction arrays have **no stoichiometric-coefficient mechanism**: every
filled product slot receives the full `RATE`, with no per-product weight. The
redistribution needs fractional deposits (`epsilon` into bin k, `1-epsilon` into
k+1), which cannot be encoded as-is. Resolution, decided in the design thread:

**Extend the reaction record with a per-product weight field, default 1.0.** This is
NOT the "new solver code" the brief warns against — that warning is about the
*integrator* (no operator splitting, no separate dust ODE, DLSODES/Newton/Jacobian
machinery untouched). A per-product weight in the reaction *table*, read by the
existing RHS assembly loop, is a local data-structure extension. With weight = 1.0
for every existing chemical reaction, all current chemistry is reproduced
**bit-identically** (the equivalence tests must still pass). The weighted deposit is:

```
contribution to dY(product_p)/dt  =  w_p * R * Y(r1) * Y(r2) * nH
```

**One reaction per unordered pair. Loss counted once; gain split by weights.** Do
NOT split a pair into one reaction per product bin — that would double-count the
loss (the reactant depletion fires once per split reaction). Instead:

```
UNORDERED pair (i,j), product bins k and k+1:
    ONE generated reaction:
      reactants:  GRAIN_i , GRAIN_j          (loss: depletes each once, full rate)
      products:   GRAIN_k (weight epsilon) , GRAIN_{k+1} (weight 1-epsilon)   (gain)
    rate coefficient R:
      i != j :  R = K_ij           (loss counted once; NO 1/2)
      i == j :  R = (1/2) * K_ii    (explicit 1/2, COMBINATORIAL — see below)
```

The `M_ijk` of the Smoluchowski equation are realized as these product weights
(`epsilon`, `1-epsilon`), NOT folded into the rate. The rate carries only `K_ij`
(and the diagonal 1/2).

**The 1/2 lives ONLY on the diagonal (i == j)**, and it is the combinatorial factor
`n_i(n_i-1)/2 ~ n_i^2/2` (distinct pairs drawn from n_i identical grains), NOT a
double-counting correction. Off-diagonal pairs carry no 1/2 because each unordered
pair is written once. **Before writing any generated reaction, determine how the
existing RHS handles two identical reactants** (`GRAIN_i + GRAIN_i`, or the
chemistry's `H + H` surface reactions). If the code already applies a 1/2 (or an
`n(n-1)` vs `n^2` distinction) for identical reactants, the generated diagonal `R`
must NOT double it. Match the existing convention exactly. A factor-of-2 error here
is invisible until it shows up as dust-mass drift — exactly what the rung-1
diagnostic catches.

**The weight field MUST flow into the Jacobian, not just the RHS.** A weighted
product deposit `w_p * R * Y_i * Y_j` has partials `w_p * R * Y_j` and
`w_p * R * Y_i`. If the Jacobian assembly reads the same product records, it
inherits `w_p` for free; if there is a separate Jacobian code path, apply the weight
there too. Verify this explicitly — a weight that reaches the RHS but not the
Jacobian gives a correct derivative value with a wrong Jacobian, degrading Newton
silently.

**The same weight field serves ice transport** (ice-follows-particle uses the same
`epsilon`), so this is the general mechanism for all of Phase II, not a rung-1
special case.

### Mass conservation identity (the diagnostic's basis)

`M_ijk` conserves mass by construction: `sum_k M_ijk * m_k = m_i + m_j` for every
pair. Therefore `d/dt sum_k m_k Y_k = 0` exactly, if the gain/loss factors are
consistent. The dust-mass diagnostic monitors this sum; any drift = a factor bug.

---

## 1. Decisions locked in the design thread (do NOT re-litigate)

- **Grains are species; coagulation is generated bilinear pseudo-reactions.** One
  state vector, one implicit DLSODES solve. No operator splitting, no separate
  dust integrator.
- **Ice-follows-particle weighting for v1:** transported ice uses the SAME
  coefficient `M_ijk` (equivalently `epsilon`) as the grains. Conserves total ice
  by construction (ice and grains move in lockstep). Mass- and area-weighted
  variants are implemented behind the redistribution interface as DIAGNOSTICS
  only. Do not make mass/area the default.
- **Per-H abundance + floor**, NOT per-grain occupation theta. Theta removes the
  emptying-bin singularity but is a v2 rewrite. For v1: per-H with a floor on
  `Y(GRAIN_k)`.
- **Grid:** mass-doubling (`mass_ratio = 2`), fiducial 20 bins, 5 nm – 0.5 um.
  Grid already derived in `parameters.in` (Stage 3). Top-bin sink already hooked
  (`dust_mass_sink`).
- **Kernels:** constant (validation) + Brownian (physical v1). Behind a single
  interface returning `Delta v_ij`. Turbulent/drift kernels deferred.
- **Sticking p_stick = 1** (coagulation only). Fragmentation deferred (it is a
  DIFFERENT redistribution operator, products to smaller bins — not a v1 change).
- **2-phase only** (`is_3_phase = 0`). Do not enable 3-phase.
- **Reduced chemical network** for coupled runs (N_ice ~ tens, not ~950). Keeps
  the O(N^2 N_ice) ice-transport term count tractable (~10^4, not ~10^5).

---

## 2. The equations, in code variables (implement to match)

**Redistribution (frozen at init), Podolak/Brauer.** For the collision product
`m_+ = m_i + m_j`, find the bracketing bins `m_k <= m_+ < m_{k+1}` and set

```
epsilon = (m_{k+1} - m_+) / (m_{k+1} - m_k)
M_ijk   = epsilon          (into bin k)
M_ij,k+1 = 1 - epsilon      (into bin k+1)
```

Floating-point: compute the two mass contributions as `epsilon*m_k` and
`m_+ - epsilon*m_k` (NOT `(1-epsilon)*m_{k+1}` independently) so they sum to
exactly `m_+`. Self-collision `i=j` on the doubling grid gives `m_+ = 2 m_k =
m_{k+1}`, so `epsilon = 0`: entire product into bin k+1, no split. If
`m_+ >= m_N` (top bin), send the mass to `dust_mass_sink` (already hooked); it
does not re-enter the grid. Monitor it.

**Kernel.** `K_ij = sigma_ij * Delta v_ij * p_stick`, with `p_stick = 1`.
`sigma_ij = pi (a_i + a_j)^2`, frozen at init. Two `Delta v` models behind the
interface:

```
constant:   K_ij = K0                       (validation only; a fixed number)
Brownian:   Delta v_ij = sqrt( 8 kB T_gas (m_i + m_j) / (pi m_i m_j) )   (T_GAS, not T_dust)
```

Cross-check the Brownian prefactor against DustPy's source so the constant-kernel
validation and the Brownian implementation use consistent conventions.

**Ice transport (generated pseudo-reaction), ice-follows-particle.** For ice
species X on bin i, colliding with bin j, product bins k / k+1:

```
generated:  J_i X + GRAIN_j  ->  (epsilon) J_k X + (1-epsilon) J_{k+1} X
rate coeff: same K_ij as the grain reaction (and same 1/2-on-diagonal rule)
```

Total ice conserved because sum over product bins of the weight = 1 (epsilon +
(1-epsilon)).

**Dynamic GTODN (rung 5 only).** Replace frozen `GTODN(k)` with `1/Y(INDGRAIN(k))`
at every surface-rate site (~23 sites: accretion, Langmuir-Hinshelwood, ER/CIR,
monolayer counting). Add the analytic Jacobian entries
`d(rate)/dY(GRAIN_k)`; for diffusive rates `~ 1/Y_grain`, so the derivative
`~ -1/Y_grain^2`. Impose the floor on `Y(GRAIN_k)` to keep these finite as bins
empty. Value set empirically.

---

## 3. The rungs (each is a separate thread; each gates the next)

### Rung 1 — Pure grain coagulation
ADD: grains-as-species coagulation (redistribution + gain/loss + both kernels).
Chemistry OFF, ice OFF.

CONFIGURATION (how "chemistry off, ice off" is realized): **build a minimal,
self-contained input set** — a species file with only `GRAIN_k`, a reaction file
with only the generated coagulation reactions, no gas-phase or surface chemistry.
Do NOT rely on a runtime "disable chemistry" flag (one may not exist, and a flag
can leak). A minimal fixture is chemistry-free *by construction*, inspectable, and
becomes a permanent regression fixture for the coagulation sector. First confirm by
experiment the minimal species/reaction set the initializer will accept (strip to
grains, see what it demands — it may require a couple of dummy gas species to
initialize; if so, include them inert).

DIAGNOSTIC (both must pass):
  (a) **Total mass — machine precision.** `sum_k m_k Y_k` conserved throughout.
      This is the factor-of-2 trap: any drift = inconsistent gain/loss factor.
  (b) **Number decay vs analytic — discretization-limited, ~5%.** IC = monodisperse
      (all mass in bin 1, rest at floor). Kernel = constant, `K_ij = K0`. Compare
      the zeroth moment `M0(t) = sum_k n_k` against the analytic constant-kernel
      solution
          `M0(t) = M0(0) / (1 + (1/2) K0 M0(0) t)`
      at several times. Pass at rtol ~ 5% on the fiducial grid. Do NOT expect
      machine precision on (b) — the binned scheme approximates the continuous
      solution. THEN **demonstrate convergence**: refine the grid (smaller
      `mass_ratio`) and show the agreement tightening toward the analytic curve.
      The convergence demonstration is the real validation of the redistribution
      scheme and is a paper figure, not busywork.
      (Get the exact analytic form and its citation from Lombart & Laibe 2021,
      whose test cases these are — this also fixes the placeholder citation in the
      paper draft.)

This rung traps the self-collision factor-of-2 and validates redistribution.
Do NOT add ice or chemistry.

### Rung 2 — Ice transport on a 3-bin toy
ADD: generated ice-transport pseudo-reactions (ice-follows-particle). Still no gas
chemistry. Grid = 3 bins (cheap, hand-checkable).
DIAGNOSTIC: total ice per species conserved across bins to machine precision; the
clean-merge 3-grain example (1/3, 2/3 split) reproduces by hand.
Do NOT scale the grid or enable chemistry.

### Rung 3 — Couple with chemistry (reduced network, 3–5 bins)
ADD: turn gas-grain chemistry back on; coagulation + ice + chemistry in one
DLSODES call. Small grid.
DIAGNOSTIC: rung-1 and rung-2 conservation still hold WITH chemistry active; the
coupled run differs from a chemistry-only run in the expected direction (grains
grow, ice redistributes). First genuinely coupled run.
Do NOT scale to 20 bins yet; do NOT touch GTODN.

### Rung 4 — Jacobian/sparsity plumbing + scale to 20 bins
ADD: append coagulation and ice-transfer couplings to the Stage-4 symbolic
`(row,col)` list (grain columns go dense; cross-bin ice couplings appear). Scale to
the 20-bin fiducial.
DIAGNOSTIC: the debug `subset` assert (numerical ⊆ symbolic) STILL HOLDS with
coagulation live — if it fires, the append logic missed a coupling; solver
converges at 20 bins; conservation still machine-precision.

### Rung 5 — Dynamic GTODN
ADD: replace frozen GTODN at the ~23 sites; add reciprocal Jacobian entries;
impose the floor. LAST, because it modifies existing chemistry.
DIAGNOSTIC (in order):
  (a) coagulation OFF: dynamic GTODN reproduces the frozen-GTODN reference to a
      tight relative tolerance (Y(GRAIN) is constant, so it must reduce to frozen);
  (b) coagulation ON: runs stably through bin-emptying, no Inf/NaN, conservation
      holds.

---

## 4. Landmine checklist (the places design effort was spent)

0. **Product-weight field must reach the Jacobian, not just the RHS** (section 0).
   The weight defaults to 1.0 (existing chemistry bit-identical) and realizes the
   `epsilon`/`1-epsilon` redistribution. A weight applied in the RHS but missing
   from the Jacobian gives correct derivatives with a wrong Jacobian — silent Newton
   degradation. Verify both paths carry `w_p`. Foundational: this encoding underlies
   every generated reaction in Phase II (grains AND ice).
1. **Self-collision factor-of-2** (section 0). The `1/2` is diagonal-only and
   combinatorial. It rides on the rate `R`, NOT on the product weights. One reaction
   per unordered pair (loss once); do not split per product bin. Match the existing
   identical-reactant convention. Caught by rung-1 mass conservation.
2. **Ice-weight conservation.** Ice-follows-particle conserves by construction
   ONLY if the ice uses the same epsilon as the grains. Do not independently
   recompute an ice weight. Caught by rung-2 per-species ice conservation.
3. **Top-bin sink.** Mass with `m_+ >= m_N` goes to `dust_mass_sink`, monitored.
   Non-negligible accumulation => `a_max` too low (a finding, not a bug).
4. **Symbolic pattern must be superset-by-construction** (no rate-based pruning),
   as established in Stage 4. Coagulation couplings appended the same way. The
   `subset` assert is the guard — it must never fire.
5. **Reciprocal Jacobian stiffness (rung 5).** Missing `d/dY(GRAIN)` entries
   degrade Newton exactly where it is stiffest (emptying bins). Add them
   analytically; do not rely on solver robustness.
6. **GTODN floor (rung 5).** Keeps `1/Y_grain` finite. Chemically negligible,
   numerically safe. Empirical.
7. **Brownian prefactor** matches DustPy convention (section 2), else the
   constant-kernel validation and Brownian runs are inconsistent — and mass
   conservation will NOT catch it (it constrains the rate's structure, not its
   magnitude). Only the analytic constant-kernel timescale test catches it.

---

## 5. Standing conventions (unchanged from Stages 1–4)

- One commit per rung; separate commits are bisectable when the dust sector
  misbehaves.
- Diagnostics exist BEFORE the physics they monitor; each new diagnostic reads
  the known-correct (zero-drift) answer before the next rung adds complexity.
- `equivalence.sh` stays pinned to `sparsity=numerical` for the external nmgc
  regression; coagulation-off runs remain the reference.
- Keep the numerical/symbolic sparsity flag; keep coagulation behind an enable
  flag so chemistry-only reference runs stay reproducible.
```