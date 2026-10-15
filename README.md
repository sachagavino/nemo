# NEMO

A 0D multi-grain gas–grain chemical code in which the **dust size distribution
evolves inside the ODE system**, by coagulation, alongside the chemistry — solved
together in one implicit step, not operator-split around it.

Grain size bins are carried as species and coagulation is encoded as generated
pseudo-reactions, so DLSODES advances the gas chemistry, the grain-surface
chemistry, the ice mantles and the dust size distribution as a single coupled
system with one analytic Jacobian. The ice on a grain rides with it when grains
coagulate.

The code is extracted from [NMGC-2.0](https://github.com/sachagavino/nmgc-2.0)
(Gavino et al.), itself descended from Nautilus (Hersant, Wakelam, Cossou,
Iqbal). The gas-phase network is `kida.uva.2014` (Wakelam et al. 2015) and the
surface network is that of Ruaud et al. (2015). Fortran core; python is for
post-processing only.

> **Status: v1, proof of concept.** Coagulation, collisional ice transport,
> charge-resolved grain merging, a symbolic-sparsity analytic Jacobian, and the
> coupled integration are implemented and validated (see **Validation** below).
> The coupling is demonstrated to be correct and cheaper than operator-splitting
> within one solver. Fragmentation, grain drift, a larger size range and
> time-dependent physical conditions are future work.

## Build and run

```sh
make                 # -O2, binary in bin/nmgc
make OPT=-O0         # matches the nmgc-2.0 reference build (bit-for-bit)
make debug           # bounds checking, warnings, backtrace

mkdir run && cd run
cp ../inputs/* .
../bin/nmgc run      # integrate
../bin/nmgc outputs  # binary -> ASCII in ab/, ml/, struct/
```

`inputs/` is the reference input set. To run the reduced C/H/O network used for
the coupled coagulation models, copy the species and reaction files from
`networks/reduced_CHO/` over the corresponding files in your run directory; the
full network lives in `networks/full/`.

## Coagulation

Coagulation is off by default. The relevant knobs in `parameters.in`:

```
coagulation         = 1          # master switch (0 = off)
coagulation_kernel  = constant   # constant | brownian
constant_kernel_k0  = 1.0E-09    # [cm^3/s], constant kernel only

dust_grid_source    = derived    # derived | tabulated
a_min               = 5.0E-07    # [cm] smallest representative radius (5 nm)
a_max               = 5.0E-05    # [cm] largest representative radius (0.5 um)
mass_ratio          = 2.0        # per-bin mass ratio (<= 2); nb_bins is derived
grain_density       = 3.0        # [g/cm^3] bulk grain material density

sparsity            = symbolic   # symbolic (default) | numerical
initial_dtg_mass_ratio = 1.0E-02 # dust-to-gas mass ratio; normalises the MRN dust IC
```

The size grid is **derived** from `(a_min, a_max, mass_ratio)` and the bulk
density — a geometric-in-mass grid, which is the natural grid for coagulation
(a collision product `m_i + m_j` is bracketed by adjacent bins). `mass_ratio` is
capped at 2 (mass-doubling): coarser grids break the Podolak/Brauer
redistribution. The tabulated path (`dust_grid_table.in`) is kept for the
bit-for-bit regression against nmgc-2.0 and for custom initial distributions.

Redistribution of a collision product between the two bracketing bins uses the
Kovetz–Olund/Brauer number weights `(eps, 1-eps)`, conserving mass to round-off.
Ice follows the grains with the *same* weights (ice-follows-particle), so each
ice species is conserved by construction. Grain merging is charge-resolved over
the `{0,-}` states, with a free electron shed in the `(-)+(-)` channel to
conserve charge. Collisions whose product exceeds the top bin are not generated;
the dropped-flux fraction is monitored and must stay negligible over a run
(~1.5% at the v1 fiducial).

The target v1 configuration is 5 nm – 0.5 µm, mass-doubling (~20 bins, top
representative radius ≈ 0.40 µm), 2-phase, reduced C/H/O network.

## Layout

```
src/
  global_variables.f90   state, parameters, procedure-pointer seams
  environment.f90        0D box properties; Td(a), T_CR,peak(a) prescriptions
  input_output.f90       readers/writers
  ode_solver.f90         RHS, Jacobian, rate coefficients, symbolic sparsity  <- the hot path
  gasgrain.f90           init, elemental conservation, solver method flag
  outputs.f90            binary -> ASCII
  main.f90               driver: output loop + DLSODES (no sub-stepping, on purpose)
  odepack/               DLSODES, untouched
  dust/                  dust_grid.f90 (size grid) + dustevolution.f90 (coagulation + ice transport)
networks/                full/ and reduced_CHO/ species + reaction sets
inputs/                  reference input set (parameters.in and the .in files)
tests/                   validation harnesses and fixtures (see below)
scripts/                 python post-processing
docs/                    design and validation notes (EXTRACTION.md, STAGE*, PhaseII/III)
```

## Validation

The coupled path is verified at several levels; the harnesses live in `tests/`.

- **nmgc-2.0 equivalence** (`equivalence.sh`): with coagulation off, numerical
  sparsity and `-O0`, NEMO reproduces nmgc-2.0 bit-for-bit (0 difference across
  56,560 values). This validates the extracted frozen path.
- **Coagulation vs analytic** (`coag_gate.sh`, `coag_diffusion_gate.sh`): the
  constant- and additive-kernel Smoluchowski solutions (Scott 1968; Golovin
  1963) are reproduced; the cell-averaged error isolating the coagulation
  operator converges monotonically under grid refinement.
- **Ice transport** (`coag_ice_*`): conservation of each ice species to ~1e-16;
  the transport operator is net-zero per species by construction.
- **Jacobian** (`ice_transport_jac_sweep.sh`, `gtodn_jac_fd_check.sh`): the
  analytic coagulation and chemistry-coupling entries agree with finite
  differences to the FD noise floor, verified on the fiducial 20-bin network.
- **Coupled vs operator-split** (`tests/phase3_split/`): the coupled solve
  converges to the operator-split result as the split step shrinks, and is both
  cheaper and more accurate than the split at any affordable cost — a
  method-fidelity-and-cost result.
- **Floor insensitivity, mass/charge conservation**: the grain-abundance divisor
  floor is never reached at the fiducial (32-decade margin; 0.1×/10× reruns
  byte-identical); dust-mass and charge drift are at round-off.

Numbers of record for the paper are in the project notes (`v1_validation_measurements.md`).

## Known limitations (v1)

- **Charge-blind coagulation kernel and charging rates.** The coagulation kernel
  ignores Coulomb/image-charge focusing, and the inherited grain-charging rates
  are purely geometric (which biases grains too negative). The two approximations
  act in opposite directions and largely cancel: a self-consistent focused-charging
  run (`charging_focusing` switch, default off) puts the net effect on the
  coagulation rate within a factor ~1.5 of unity. A charge-resolved kernel is
  future work.
- **Reduced network for the coupled runs.** Ice transport generates ~N²·N_ice
  terms, so the full network is intractable coupled; the C/H/O subset is used.
- **2-phase, constant conditions.** The 3-phase path is kept but not exercised
  under coagulation; fragmentation, drift, and time-dependent Td/T_gas/UV are not
  on this branch.
