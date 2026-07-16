# Stage 3: derived dust grid and input restructuring

Stage 3 makes the dust size grid a **derived, first-class object** instead of a
hand-typed table. The grid geometry moves into `parameters.in`; NEMO computes
the representative masses, radii and per-bin temperatures at init.

## Two orthogonal concerns, two switches

The grid (geometry) and the initial distribution (population) are kept separate,
because coagulation later evolves the population *on* a fixed grid.

| switch | values | what it controls |
|---|---|---|
| `dust_grid_source` | `derived` \| `tabulated` | grid geometry: masses & radii of the bins |
| `dust_ic`          | `MRN` \| `tabulated`     | initial number per bin (the IC) |

* **Grid geometry** (`derived`): a grid geometric in mass,
  `m_k = m_min * mass_ratio^(k-1)`, from `a_min`, `a_max`, `mass_ratio` and
  `grain_density`. `nb_grains` is derived from `mass_ratio`, not typed in. The
  mass grid is the natural grid for coagulation: an `(i,j)` collision lands at
  `m_i + m_j`, trivially bracketed on a geometric mass grid.
* **Initial distribution** (`MRN`): `dn/da ~ a^dust_power_law_index`
  (default `-3.5`), integrated over each bin and normalised so that
  `sum_k n_k m_k = initial_dtg_mass_ratio * AMU` — the same convention as the
  single-grain GTODN of nmgc-2.0, to which it reduces for `nb_grains = 1`.
* **Temperatures**: on a derived grid, `Td(a)` and `T_CR,peak(a)` are named
  prescriptions in `environment.f90`, evaluated once on the grid. `Td(a)`
  defaults to the promoted `a^(-1/6)` scaling anchored at
  (`reference_grain_radius`, `reference_dust_temperature`). `fixed_to_dust_size`
  is now a real function on the derived path, not a table read.

The `mass_ratio <= 2` (mass-doubling) constraint is enforced in
`dust_grid_count_bins`: a coarser grid breaks the Podolak/Brauer redistribution
bracketing that coagulation will rely on.

## Why the regression runs through the *tabulated* path

The equivalence test reproduces nmgc-2.0 on the 0D 2-grain case bit-identical
(56560 values, worst relative difference 0). That case **cannot** be reproduced
by a physical derived grid, and this is by design, not a gap:

* its two representative radii, `2.3265e-6` and `1.0404e-3` cm, span a factor
  `~9e7` in mass — about 28 mass-doubling bins, not 2. Under `mass_ratio <= 2`
  no derivation yields a 2-bin grid there;
* its per-bin abundances are not MRN-normalised (their `sum_k n_k m_k` sits
  ~4.4x off the dtg-normalised value), so `dust_ic = MRN` would not reproduce
  them either.

So the reference grid is an *unphysical representative-size* setup, and the
regression's job is to prove the plumbing is unchanged. It does that through the
tabulated reader: `tests/reference_0D_2grains/` sets
`dust_grid_source = tabulated`, `dust_ic = tabulated`, reading the reference
radii, abundances and temperatures verbatim from `dust_grid_table.in` (the
renamed, byte-for-byte legacy `0D_grain_sizes.in`). `0D_grain_sizes.in` as a
name is retired; the tabulated reader stays, for restarts and for external
(DustPy / mcdust) hand-offs later.

## That the derived path feeds the chemistry identically

Proven directly. Every run exports its active grid to `dust_grid_active.out` at
full (round-tripping) precision. Feeding a **derived** grid's export back through
the **tabulated** path reproduces the run bit-identical (worst relative
difference 0), which shows the two entry points fill the same per-bin arrays and
the chemistry consumes them identically. (At a truncated 14-figure export the
round-trip drifts to ~1e-3 on trace species — the expected sensitivity of a
stiff network to a ~1e-13 input perturbation, not a code-path difference.)

The default derived science grid (`a_min = 5 nm`, `a_max = 0.5 um`,
`mass_ratio = 2`) yields 20 bins, radii stepping by `2^(1/3)`, with MRN mass
conservation exact to double precision.

## Top-bin sink (a grid property, decided now)

A top-bin self-collision scatters mass above `m_N`, which has no bin on a fixed
grid. v1 policy: it accumulates in `dust_mass_sink` (diagnostic, `[g per H]`),
which never re-collides. A non-negligible value at runtime means `a_max` is too
low. No dynamics are attached yet — the hook is baked into the grid so
coagulation does not have to retrofit it.

## New / changed `parameters.in` keys

Added: `dust_grid_source`, `a_min`, `a_max`, `mass_ratio`, `dust_ic`,
`dust_power_law_index`, `reference_dust_temperature`.
Renamed: `grain_radius` → `reference_grain_radius` (legacy names still accepted).
Unchanged and reused: `grain_density`, `initial_dtg_mass_ratio`.
