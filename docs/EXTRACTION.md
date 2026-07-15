# Extraction audit: nmgc-2.0 -> 0D skeleton

Reference commit: `sachagavino/nmgc-2.0`, `master`.

## Result

The skeleton reproduces the reference **exactly** on the 0D 2-grain 2-phase
case (see `tests/reference_0D_2grains/`):

| | reference nmgc-2.0 | skeleton |
|---|---|---|
| species | 1414 | 1414 |
| runtime to 1e5 yr, `-O0` | 158 s | 161 s |
| species x output times compared | — | 56560 values |
| worst relative difference | — | **0** |

`elemental_abundances.out` is byte-identical. The comparison is re-runnable:
`./tests/equivalence.sh /path/to/nmgc-2.0`.

## Kept, untouched

* DLSODES / ODEPACK (`src/odepack/`)
* `get_temporal_derivatives`, `get_jacobian`, `set_chemical_reactants`,
  `init_relevant_reactions`, `count_nonzeros`, `set_work_arrays`
* the multi-grain template-replication scheme in `read_species` / `read_reactions`
* `is_er_cir` and `is_crid`: dormant (default 0) but fully wired
* the 3-phase path (`is_3_phase = 1`): kept, not burned. v1 is 2-phase.
* `get_MRN_distribution` in `input_output.f90` — unused today, it is most of
  the analytic dust IC that stage 3 needs.

## Removed

| removed | why |
|---|---|
| `structure.f90` (1D grid, species diffusion, Crank-Nicholson, `init_1D_static`, `get_timestep_1D_diff`) | 1D |
| structure tables (`structure_evolution.dat`, `1D_static.dat`, `1D_grain_sizes.in`) | 1D / table-driven structure |
| `dust_temperature_module.f90` | Td will be a prescription, not a radiative-equilibrium solve |
| `photorates.f90` + `cross-sections/` | disk photorates |
| `major_reactions.f90`, `trace_major.f90`, `select_outputs.f90`, `rates.f90` | post-processing, belongs in `scripts/` |
| `is_h2_adhoc_form`, `photo_disk` | dead switches |
| the `x_i` spatial loop in `main.f90`, `reaction_rates_1D`, `NH_z`/`NH2_z`/... | 1D |

The 0D parts of `structure.f90` were carved out into **`environment.f90`**:
`get_structure_properties_fixed`, `get_grain_temperature_{fixed,fixed_to_dust_size,gas}`.
`get_timestep` and `structure_diffusion` disappeared entirely — in 0D the
integration step *is* the output step, and coagulation will live inside the ODE
system rather than beside it.

## Shape changes

`abundances(nb_species, spatial_resolution)` -> `abundances(nb_species)`.
`gas_temperature`, `visual_extinction`, `H_number_density` are now scalars;
`dust_temperature(nb_grains, spatial_resolution)` -> `dust_temperature(nb_grains)`.
This changes the layout of `abundances.out`; `outputs.f90` was rewritten to
match, and the python scripts in `scripts/` still need updating (they assume the
1D record layout).

## parameters.in migration

Dropped keys (silently ignored if left in the file):
`structure_type`, `spatial_resolution`, `is_structure_evolution`, `photo_disk`,
`is_h2_adhoc_form`, `height_h2formation`.

`grain_temperature_type` now accepts only `fixed`, `fixed_to_dust_size`, `gas`.

`height_h2formation` deserves a note: it was a *spatial* threshold above which
the Bron+2014 H2 formation rate was applied. In 0D with `x_i = 1` the test was
always true for any positive value, so the behaviour is preserved exactly by
gating on `is_h2_formation_rate` alone.

## Known issues, carried forward (not fixed here — stage 4)

1. **`set_work_arrays` infers the Jacobian sparsity numerically**, by evaluating
   the Jacobian and thresholding at `1e-99`. It must become a symbolic build
   from the reaction list. It is currently called on *every* output step.
2. **`set_dependant_rates` does ~90 string comparisons and 11
   `READ(c_i,'(I2)')` per RHS call**, to recover a grain rank that is already
   known at init. Precompute the index maps once.
3. **`FFLAGS` was empty in nmgc-2.0**, i.e. `-O0`. The new `Makefile` defaults to
   `-O2`; `make OPT=-O0` reproduces the reference build.

### One correction to the issue list

The reported `COND` bug — "`COND` is allocated `dim(nb_grains)` in
`gasgrain.f90` but shadowed by a local `cond(nb_reactions)` in
`set_dependant_rates`" — **is not a bug.** There is no global `COND`. The one in
`init_reaction_rates` is a *local* allocatable of that subroutine, and the
`cond(nb_reactions)` in `set_dependant_rates` is a separate local in a different
scope. They never alias. Both are correctly sized for their own use
(`nb_grains` for the per-bin geometric cross-section, `nb_reactions` for the
per-reaction one). The names are confusing and worth changing, but nothing is
being clobbered. Verified against the reference: results are bit-identical.

## Next

3. restructure the inputs (grid derived in `parameters.in`; dust IC in its own
   file; Td(a) and T_CR,peak(a) as named prescriptions; `0D_grain_sizes.in` dies)
4. clean the RHS hot loop, symbolic sparsity, `-O2`, re-verify against this same
   reference case

Coagulation only after all four are green.
