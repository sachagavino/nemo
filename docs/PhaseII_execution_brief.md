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

### The rate coefficient of the generated reactions (THE landmine — get this exact)

**Convention: write each UNORDERED pair (i,j) once.** Under that convention:

```
OFF-DIAGONAL  (i != j):
    grain gain into k:   R = M_ijk * K_ij          (NO explicit 1/2)
    grain loss:          bin i and bin j each lose one grain per event
    => generated:  GRAIN_i + GRAIN_j  ->  (M_ijk) GRAIN_k + (M_ij,k+1) GRAIN_{k+1}

DIAGONAL      (i == j):
    grain gain into k:   R = (1/2) * M_iik * K_ii   (explicit 1/2, COMBINATORIAL)
    grain loss:          bin i loses TWO grains per event
    => generated:  GRAIN_i + GRAIN_i  ->  ...
```

The `1/2` lives **only on the diagonal**, and it is the combinatorial factor
`n_i(n_i-1)/2 ~ n_i^2/2` (distinct pairs drawn from n_i identical grains), NOT a
double-counting correction. Off-diagonal pairs carry no `1/2` because each
unordered pair is written once.

**Before writing any generated reaction, determine how the existing RHS handles a
reaction with two identical reactants** (`GRAIN_i + GRAIN_i`, or the chemistry's
`H + H` surface reactions). If the code already applies a `1/2` (or an `n(n-1)`
vs `n^2` distinction) for identical reactants, the generated diagonal `R` must NOT
double it. Match the existing convention exactly. A factor-of-2 error here is
invisible until it shows up as dust-mass drift — which is exactly what the rung-1
diagnostic is built to catch.

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
DIAGNOSTIC (both must pass):
  (a) dust mass `sum_k m_k Y_k` conserved to machine precision throughout;
  (b) constant-kernel run matches the analytic Smoluchowski solution.
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

1. **Self-collision factor-of-2** (section 0). The `1/2` is diagonal-only and
   combinatorial. Match the existing identical-reactant convention. Caught by
   rung-1 mass conservation.
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