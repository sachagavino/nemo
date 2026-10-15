# Phase III (validation) — executor brief: Jacobian provenance dump + spy figures

**Owner of this brief:** design thread. **For:** executor thread.
**Status of decisions:** all design calls below are made and marked `DECIDED`. Where a
call is reversible cheaply, it says so. Do not silently re-decide; if a `DECIDED` item
looks wrong, stop and flag it back to the design thread (this is the §3 discipline —
executor-side inferences must not harden into facts).

**Goal.** Produce the two Jacobian figures for the paper (§5.3):
1. **Primary** — a block-structured *spy* plot of the coupled Jacobian, colored by
   provenance, showing gas-chemistry / surface-network / coagulation / ice-transport
   coexisting in one matrix (the "everything is one implicit Jacobian" thesis).
2. **Secondary** — a zoom on the non-reactant coupling entries (the grain-abundance
   couplings and the charge/electron-shedding entries), the visual companion to the
   Jacobian-verification subsection (§jacverif).

This is a **read-only** dump of the pattern the solver already builds, plus a classifier
and a plotting script. **No production code path is modified.** If you find yourself
editing `build_symbolic_sparsity`, `get_jacobian`, or `dust_coagulation_inject_static`,
stop — that is out of scope and a signal the approach has drifted.

---

## 0. Ground truth (already in the code — do not reinvent)

The pattern DLSODES receives is built once in `build_symbolic_sparsity()`
(`src/ode_solver.f90`) as CSC arrays:
- `IA_SYM(nb_species+1)`, `JA_SYM(NNZ_SYM)`, `NNZ_SYM` — the finished pattern.
  Column `c` occupies `JA_SYM(IA_SYM(c) : IA_SYM(c+1)-1)`, which lists its nonzero
  **rows**. (Key packing there is `col*(N+1)+row`; you do not need it — read IA/JA.)

Entries enter from exactly three sources (see the routine):
- **Reaction-derived** (pass 2): for each reaction `r`, for each reactant `col` and each
  compound (reactant *or* product) `row`, mark `(row,col)`. This covers gas chemistry,
  surface chemistry, coagulation grain redistribution, and ice transport — all of which
  are reactions (`REACTION_TYPE`, `REACTION_COMPOUNDS_ID`).
- **Declared non-reactant** (`build_live_divisor_dependencies()` → `DECLARED_JAC_ROW/COL`,
  `NB_DECLARED_JAC_DEPS`): the grain-abundance coupling block — every compound of every
  bin-`k` reaction coupled to `INDGRAIN(k)` and `INDGRAIN_MINUS(k)`.
- **Diagonal**: every `(i,i)`.

Value-side maps you will read for tagging (they add **no** pattern entries):
- `grain_col_bin(nb_species)` — species index → grain bin `k` for grain columns (else 0).
- `gtodn_jac_list(gtodn_jac_n)` — the reactions carrying the reciprocal `−rate/Y²` value
  (types 14 LH, 99 accretion, 66/67 photodesorption), built by
  `build_gtodn_jacobian_map()`.

Reaction typing you will use:
- `COAGULATION_TYPE` — coag + ice-transport pseudo-reactions (slots via
  `dust_coagulation_inject_static`). Grain-only vs ice-transport is distinguished by
  whether reactant slot 1 is a `J`-ice species (`ice_name` convention `Jkk X`) or a
  `GRAIN` species (`grain_name` convention).
- `e-` in product slot 6 of a coag reaction marks the `(−)+(−)→(−)+e⁻` electron-shed
  entries (`REACTION_COMPOUNDS_NAMES(6,slot)=='e-'`).
- Chemistry reactions are those with index `<= nb_chemistry_reactions`; generated coag
  reactions are `> nb_chemistry_reactions`. Surface vs gas within chemistry: surface
  types include 14/66/67 and accretion 99; treat everything else as gas. (Executor:
  confirm the surface-type set against `networks/reduced_CHO/grain_reactions.in` +
  `REACTION_TYPE`; do not assume beyond 14/66/67/99 without checking.)

**Co-location fact (do not fight it):** the live-divisor block and the reciprocal-GTODN
entries occupy the *same* `(row,col)` cells (GRAIN columns × bin-`k` reaction compounds).
They differ in value, not location. The spy plot therefore has **one** non-reactant-
coupling class; the GTODN-reciprocal subset is an optional *value-derived overlay*
(cells whose column bin has a reaction in `gtodn_jac_list`), never a separate color by
position.

---

## 1. Fixture

Create/confirm `tests/fixture_jac_figure/` = the `reduced_CHO` network with an
**ice-bearing** initial state and a **4-bin** derived grid:
- `coagulation = 1`, `sparsity = symbolic`, `coagulation_kernel = constant`.
- `mass_ratio = 2`, and `a_min`/`a_max` chosen so `dust_grid_count_bins` yields
  `nb_grains = 4` (from `nb_grains = int(log(m_max/m_min)/log(mass_ratio)) + 1`; pick
  `a_max` accordingly and print `nb_grains` at init to confirm — do not hard-code a bin
  count anywhere).
- Ice present at t0 so the ice-transport block is non-empty (reuse the ice IC path;
  `dust_ice_place` area-weighted default is fine — this figure is structural, values
  are irrelevant).

Rationale for 4 bins: legibility (§5.3). The figure carries a caption note that the
structure scales to the 20-bin fiducial; a 20-bin spy is too dense to read. **Do not**
also produce a 20-bin figure.

---

## 2. The dump (`tests/jacobian_provenance_dump.f90` + `.sh`)

Harness pattern: `call init_gasgrain()` on the fixture (this runs
`build_symbolic_sparsity`, so `IA_SYM/JA_SYM` are populated), then write two TSVs. No
integration, no solver step.

### 2a. `jac_pattern.tsv` — one row per nonzero entry
Columns: `row  col  mech_mask`

Walk every `(row,col)` in `IA_SYM/JA_SYM`. For each, compute a **bitmask** of ALL
contributing mechanisms (not just one — overlaps are the point):

| bit | mechanism        | test |
|-----|------------------|------|
| 0   | `diagonal`       | `row == col` |
| 1   | `gas_chem`       | ∃ chem reaction `r≤nb_chemistry_reactions`, non-surface type, with `col` a reactant and `row` a compound |
| 2   | `surface_net`    | as above but surface type (14/66/67/99, confirmed set) |
| 3   | `coag_grain`     | ∃ reaction `r>nb_chemistry_reactions`, `COAGULATION_TYPE`, reactant-1 a `GRAIN`, with `col` a reactant and `row` a compound |
| 4   | `ice_transport`  | as above but reactant-1 a `J`-ice species |
| 5   | `charge_electron`| the entry involves the `e-` row/col of a coag reaction whose slot-6 is `e-` |
| 6   | `grain_coupling` | `(row,col)` present in `DECLARED_JAC_ROW/COL` (the non-reactant block) |
| 7   | `gtodn_recip`    | bit 6 set AND `col`'s bin (`grain_col_bin(col)`) has a reaction in `gtodn_jac_list` whose compound == `row` (value-derived overlay flag) |

Notes:
- A cell may set several bits (e.g. `coag_grain | grain_coupling`). Emit the full mask;
  the plotter decides how to render overlaps.
- **Completeness gate (hard):** every entry in `IA_SYM/JA_SYM` must get a nonzero mask
  (diagonal counts). Any cell with `mask == 0` is a classifier bug — fail loudly with the
  offending `(row,col)` and species names. This is the analog of the "symbolic ⊇
  numerical" gate: an unclassifiable cell means the classifier does not understand the
  pattern.

### 2b. `jac_species.tsv` — one row per species (for the reorder + block boundaries)
Columns: `index  name  class  bin  charge`
- `class ∈ {gas, surface, grain, ice}` from the name (`GRAIN*`→grain, `J*`→ice, else via
  gas/surface species lists).
- `bin` = grain bin for grain/ice species (parse the 2-digit index in the name, or use
  `grain_col_bin` for grains), else 0.
- `charge ∈ {0,-}` for grains (from `grain_name` suffix), else blank.

### 2c. Cross-checks to print (cheap correctness gates)
- `NNZ_SYM` and `sum(bit set counts)` reconciled.
- count of `coag_grain` reactions and `ice_transport` reactions recovered from the tags
  vs. `nb_coag_grain_reactions` / `nb_ice_transport_reactions` (printed by
  `inject_static`) — these should be consistent.
- Commit both TSVs under `tests/fixture_jac_figure/` for reproducibility.

---

## 3. The plot (`scripts/plot_jacobian_spy.py`)

Consumes the two TSVs. Pure post-processing, no code units.

**Reorder (both figures):** build a permutation from `jac_species.tsv` grouping columns/
rows as `[ gas | surface | grain(bin 1..4, charge 0 then −) | ice(bin 1..4, per base) ]`.
Within grains, order by bin so the redistribution reads along the block diagonal. Draw
faint block-boundary lines and label the four families in the margins (an unlabeled spy
is abstract art — §5.3). Caption must state: "state vector reordered by species class for
legibility; production ordering differs."

**Primary figure (`figures/jacobian_spy_primary.pdf`):**
- Full 4-bin pattern, one marker per nonzero, colored by mechanism.
- `DECIDED` (reversible — a plotter flag): overlaps colored by **novelty-first
  precedence**: `charge_electron > gtodn_recip/grain_coupling > ice_transport >
  coag_grain > surface_net > gas_chem > diagonal`. Alternative renderings supported by
  the same dump: a distinct "mixed" color for multi-bit cells, or translucent overplot.
  Expose as `--overlap {precedence,mixed,overlay}`, default `precedence`.
- Legend names the physics families, not the bit numbers.

**Secondary figure (`figures/jacobian_coupling_zoom.pdf`):**
- Restrict to the columns/rows that carry the couplings: the GRAIN columns (the
  `grain_coupling` block) plus the coag/ice blocks and the `e-` row.
- Show `grain_coupling` cells, with the `gtodn_recip` subset marked (overlay: e.g. an
  outline or hatch on those cells), and the `charge_electron` entries highlighted.
- This is the panel that pairs with §jacverif: caption points to the two-legged
  verification (FD above the noise floor, per-reaction analytic below it).

**Do NOT:** draw a value heatmap (magnitudes span ~40 orders and would erase the
coupling cells — the opposite of the point); draw a fragmentation Jacobian (v2); mirror
DustPy's four-term (`C̃_ijk`) decomposition (that dissects the *borrowed* Brauer operator,
not NEMO's contribution).

---

## 4. Deliverables
- `tests/jacobian_provenance_dump.f90`, `tests/jacobian_provenance_dump.sh`
- `tests/fixture_jac_figure/` (parameters.in etc., + committed `jac_pattern.tsv`,
  `jac_species.tsv`)
- `scripts/plot_jacobian_spy.py`
- `figures/jacobian_spy_primary.pdf`, `figures/jacobian_coupling_zoom.pdf`

## 5. Acceptance
- Dump completeness gate (§2a) passes: no `mask==0` cell.
- Cross-checks (§2c) reconcile.
- Both figures render with legible, labeled blocks at 4 bins; the coupling block and the
  electron entries are visibly identifiable in the secondary.
- Report back to the design thread: `NNZ_SYM`, the per-mechanism cell counts, and a
  first-look PNG of each figure before finalizing PDF styling.

## 6. Open items to surface (not blockers)
- Confirm the exact surface-type set (beyond 14/66/67/99) against the reduced_CHO
  network before trusting the gas/surface split.
- If the primary figure reads as too sparse/too dense at 4 bins, report before restyling
  — bin count is a design-thread call, not an executor tweak.
