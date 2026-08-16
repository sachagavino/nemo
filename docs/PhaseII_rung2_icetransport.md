# NEMO — Rung 2 execution brief: ice transport coupled to coagulation

**Repo:** github.com/sachagavino/nemo  **Branch:** `phase2-coagulation`
**Base:** the Rung 1b tip (Brownian kernel validated; patch `0005` applied and pushed).
**Read first:** the Phase II execution brief (`docs/PhaseII_execution_brief.md`),
`src/dust/README.md`, and the paper draft subsection `sec:icetransport` (the physics
this rung implements is written up there).

This is Rung 2 of the Phase II ladder. Rung 1 (pure grain coagulation) is green and
committed: constant + Brownian kernels validated against analytic solutions, mass
conserved to machine precision, coag-OFF bit-identical to nmgc-2.0 (verified directly
via `equivalence.sh`, worst rel. diff. 0). Rung 2 adds ice riding on the coagulating
grains. Do not start Rung 3 (chemistry coupling) or touch GTODN.

---

## 0. What this rung adds, and what it must NOT add

**ADD:** generated ice-transport pseudo-reactions, so that when grains coagulate the
ice mantles they carry move to the product bins. Ice-follows-particle weighting.
Plus: the surface-ice initial-condition machinery (area-weighted default + fractional
per-species override), because ice has to start *somewhere* on the grid before it can
be transported.

**Do NOT add:** gas-phase or surface *chemistry* (no reactions forming/destroying
ice — that is Rung 3). Ice in this rung is inert: it exists on grains, it rides along
when grains coagulate, and nothing else happens to it. This isolates the transport
mechanism so its conservation can be validated against a known answer before
chemistry complicates it.

**How "inert ice" is realized — extend the Rung-1 dust-only fixture, do NOT disable
chemistry in a full network.** Rung 1 established a chemistry-free minimal fixture
(grains-only species, empty reaction/surface files). Rung 2 extends *that* fixture
with a small set of inert ice species and the coagulation + ice-transport reactions —
nothing else. The ice is inert *by construction* (the fixture contains no adsorption,
desorption, or surface reactions), not by toggling flags. (In a full network,
`is_grain_reactions=0` would suppress adsorption and surface reactions, but there is
no flag to disable gas-phase reactions, and gas-phase reactions never write to the
`J_k X` ice species anyway. The fixture path avoids the question entirely and is
consistent with Rung 1.)

**Grid:** start on a 3-bin toy (cheap, hand-checkable), as the brief specifies. Scale
to the fiducial grid only after the toy conserves.

---

## 1. The physics, in code variables

When a bin-i grain and a bin-j grain coagulate, the ice species X on each is carried
to the product bins k and k+1. Ice-follows-particle: the ice uses the **same**
redistribution weight `epsilon` as the grains (Rung 1's product-weight field), so ice
and grains move in lockstep and total ice is conserved by construction.

For ice species X, the generated pseudo-reaction for the unordered pair (i,j) is:

```
J_i X + GRAIN_j  ->  (epsilon) J_k X  +  (1-epsilon) J_{k+1} X
```

- reactants: `J_i X` (the ice being transported) and `GRAIN_j` (the collision partner)
- products:  `J_k X`, `J_{k+1} X` with the SAME product weights epsilon, 1-epsilon as
  the grain reaction for that pair
- rate coefficient: the SAME `K_ij` as the grain reaction (constant or Brownian),
  including the 1/2 on the diagonal (i == j)

This is a bilinear term `~ K_ij * Y(J_i X) * Y(GRAIN_j) * nH`, evaluated by the same
RHS as everything else. It reuses the product-weight-field mechanism from Rung 1 — no
new solver machinery.

**Symmetry note (get this right):** a collision of bin i with bin j transports ice
from BOTH grains. `J_i X` rides with the bin-i grain into the product; `J_j X` rides
with the bin-j grain into the product. So each unordered pair (i,j) generates ice
transport for the ice on bin i AND the ice on bin j — two directions per pair per ice
species (except i==j, where they coincide). Do not transport only one side.

**Loss side:** the ice `J_i X` is depleted from bin i at the coagulation rate (it left
with the grain). This is the reactant-depletion the RHS does automatically — same as
the grain loss in Rung 1. Do not double-count it (one reaction per unordered pair per
direction; see the Phase II brief section 0 on loss-counted-once).

---

## 2. Surface-ice initial conditions (DESIGN DECISION — implement as specified)

Ice must be placed on the grid at init before it can be transported. The grid exists
first (from `parameters.in`, Stage 3); then ice ICs are distributed onto it.

**Default: area-weighted.** A user's `abundances.in` gives a single total surface
abundance per ice species (unchanged from current workflow, backward-compatible). At
init, distribute each species' total across bins in proportion to bin surface area:

```
X_k = X_total * (n_k a_k^2) / sum_j (n_j a_j^2)
```

This gives approximately uniform initial monolayer coverage across bins (ice coats
surface; more surface -> more ice). Near-empty bins (floored n_k) receive negligible
ice by construction — correct, no special-casing needed.

**Override: fractional, per-species.** An optional file `surface_ice_distribution.in`
lets the user specify, per species, a set of per-bin fractions (summing to 1) that
replace the area-weighting FOR THAT SPECIES ONLY:

```
CO   1.0 0.0 0.0 ... 0.0      # all CO ice on bin 1
```

- The TOTAL for each species still comes from `abundances.in` (one number). The
  override specifies only WHERE (dimensionless weights), never HOW MUCH.
- A species listed in the override uses its fractions; a species ABSENT stays
  area-weighted. The file is a set of per-species exceptions, NOT a global mode
  switch. (So the collisional-mixing experiment is a one-line file overriding only CO.)
- Fractions should sum to 1; validate and warn/normalize if they don't.

**Plausibility guard (override path only).** After applying override fractions and
multiplying by the total, check each bin's resulting coverage theta_k = X_k / sites_k.
If a bin's coverage is physically implausible (e.g. the user dumped all ice onto a
near-empty bin), WARN — do not hard-cap. In 2-phase (v1) there is no strict monolayer
ceiling (all ice is "surface", multilayers allowed), so the guard flags gross
implausibility, not theta > 1. The default (area-weighted) path cannot trigger this
and needs no guard.

**Fixture scaffolding — confirm empirically before building.** For the fixture to hold
ice at all it needs three things, and one of them is a potential snag:
  1. **Ice species defined** in the species list (`J_k X` for each ice species X and
     each bin k), with their properties (mass, and whatever surface parameters the
     placement needs). Defining a species is bookkeeping, NOT chemistry — a defined
     `J_k CO` with no reactions touching it is exactly the inert ice this rung wants.
  2. **Initial abundances**, which the executor writes into the fixture's
     `abundances.in` — one total per ice species. Choose SIMPLE, hand-traceable values
     for testability, not realism (e.g. CO total = 1e-4 per H, well above the floor).
     CO alone suffices for the conservation gate; add a second species only to confirm
     independent handling. For the hand-check (diagnostic b), a tabulated override
     placing all the ice on one bin makes the post-collision split a clean fraction to
     verify by arithmetic. The specific numbers are the executor's call (same as K0 and
     the monodisperse IC in Rung 1) — chosen to make each gate unambiguous.
  3. **Surface-structure parameters** — the coverage guard computes theta_k =
     X_k / sites_k, so it needs sites-per-bin (and area-per-bin, which the grid already
     provides). In the full code these ride along with the surface-chemistry setup that
     the fixture strips out. **Confirm empirically** (drive the initializer crash-by-
     crash, as in Rung 1) what surface scaffolding the ice-placement and coverage-guard
     paths require in a chemistry-free fixture. If sites-per-bin is not available
     without the chemistry setup, source it from the grid parameters directly. Flag
     back to the design thread if this turns out to need more than the grid provides.

---

## 3. Diagnostics (gate the rung)

**(a) Total ice per species conserved — machine precision.** For each ice species X,
`sum_k Y(J_k X)` must be constant throughout the run (ice is inert here — no chemistry
creates or destroys it, coagulation only moves it between bins). Any drift = a
transport bug (missing side, wrong weight, double-counted loss). This is the ice
analogue of Rung 1's dust-mass conservation and it is the primary gate.

**(b) Hand-checkable 3-grain example.** On the 3-bin toy with a simple IC, reproduce
the clean-merge case by hand: three grains in bin i each carrying theta molecules of X,
two coagulate to bin i+1 (self-collision, epsilon=0), giving
`N_i X = theta = (1/3) N_i X(t0)` and `N_{i+1} X = 2 theta = (2/3) N_i X(t0)`. The code
must reproduce this split. (This is the degenerate epsilon=0 case; also check a
non-degenerate split where ice divides between two product bins with epsilon != 0.)

**(c) Coag-OFF still bit-identical.** With coagulation disabled, the ice-transport
reactions are not generated, and the default build must still reproduce nmgc-2.0
bit-identically (`equivalence.sh`). The ice-IC machinery, when no override is present
and coagulation is off, must not perturb the reference case.

**(d) Ice-IC round-trip.** Area-weighted default and a fractional override both place
ice as specified; sum over bins recovers the `abundances.in` total for each species
(the distribution conserves the input total exactly).

---

## 4. Landmines specific to this rung

1. **Both directions per pair.** Each unordered pair (i,j) transports ice from bin i
   AND bin j. Missing one side silently loses half the ice transport — caught by
   diagnostic (a) as a slow drift, but easy to introduce. (i==j: one direction only.)
2. **Same epsilon as the grains.** Ice-follows-particle conserves BY CONSTRUCTION only
   if the ice reuses the grain reaction's product weight. Do not recompute an
   independent ice weight. (Mass-/area-weighted ice variants are diagnostics for a
   later comparison, NOT this rung — do not implement them now.)
3. **Weight into the Jacobian.** Same as Rung 1: the ice-transport product weights must
   reach `get_jacobian`, not just the RHS. The Rung-1 `jac_weight_check` fixture and
   the `freeze_dependent_rates` scaffold are the instruments — extend them to cover the
   ice-transport deposits. A weight in the RHS but not the Jacobian degrades Newton
   silently.
4. **Term count.** Ice transport is O(N^2 N_ice). On the 3-bin toy with a handful of
   ice species this is trivial; when scaling to the fiducial grid, confirm the reduced
   network keeps it tractable (~10^4, not ~10^5). Do NOT use a full network here.
5. **Total from abundances.in, not the override.** The override specifies fractions
   only. If an implementer reads absolute abundances from the override file, the total
   is double-specified and the two files can silently disagree. Fractions only; total
   always from `abundances.in`.
6. **Reaction injection ordering (from Rung 1).** Ice-transport reactions, like grain
   coagulation, must be in the reaction arrays before `relevant_reactions` is built, or
   their Jacobian columns will not exist. Same contiguous-type-block / tiling
   constraint. They can share the coagulation type block or form their own — either way
   the tiler must be satisfied and the ordering respected.

---

## 5. Suggested sequence within the rung

1. **Ice-IC machinery first** (area-weighted default + fractional override + guard +
   round-trip diagnostic (d)). Verify coag-OFF still bit-identical (c). This is a
   self-contained, testable piece with no transport yet.
2. **Ice-transport reaction generation** on the 3-bin toy, ice-follows-particle,
   reusing the Rung-1 product-weight field and kernel. Gate on (a) and (b).
3. **Extend `jac_weight_check`** to cover the ice-transport deposits (landmine 3).
4. **Scale to the fiducial grid** with the reduced network; confirm (a) still holds and
   the term count stays tractable (landmine 4).
5. Commit. One commit for the rung (or a clean split: IC machinery, then transport),
   bisectable.

Conceptual questions go back to the design thread. In particular: anything about
chemistry coupling, GTODN, or the mass-/area-weighted ice variants is OUT OF SCOPE for
Rung 2 — flag and defer rather than deciding here.