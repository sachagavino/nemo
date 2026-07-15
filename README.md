# nmgc-ng

A 0D multi-grain gas–grain chemical code, being rebuilt so that the **dust size
distribution evolves inside the ODE system** alongside the chemistry —
coagulation first, other processes later.

The skeleton is extracted from [NMGC-2.0](https://github.com/sachagavino/nmgc-2.0)
(Gavino et al.), itself descended from Nautilus (Hersant, Wakelam, Cossou,
Iqbal). Fortran core, no python driver: python is for post-processing only.

> **Status: skeleton.** No new physics yet. The code currently does exactly what
> nmgc-2.0 does in 0D, and is verified to do so bit for bit
> (see [`docs/EXTRACTION.md`](docs/EXTRACTION.md)).

## Build and run

```sh
make                 # -O2, binary in bin/nmgc
make OPT=-O0         # matches the nmgc-2.0 reference build
make debug           # bounds checking, warnings, backtrace

mkdir run && cd run
cp ../inputs/* .
../bin/nmgc run      # integrate
../bin/nmgc outputs  # binary -> ASCII in ab/, ml/, struct/
```

## Layout

```
src/
  global_variables.f90   state, parameters, procedure-pointer seams
  environment.f90        0D box properties; Td(a), T_CR,peak(a) prescriptions
  input_output.f90       readers/writers
  ode_solver.f90         RHS, Jacobian, rate coefficients  <- the hot path
  gasgrain.f90           init, elemental conservation, preliminary tests
  outputs.f90            binary -> ASCII
  main.f90               driver: output loop + DLSODES
  odepack/               DLSODES, untouched
  dust/                  where coagulation will go -- read its README
inputs/                  reference input set
tests/                   equivalence test against nmgc-2.0
scripts/                 python post-processing
docs/EXTRACTION.md       what was kept, dropped, and changed, and why
```

## Roadmap

1. **strip to a 0D skeleton, no new physics** — done
2. **prove equivalence against reference nmgc-2.0** on a 0D 2-grain run — done,
   0 relative difference across 56,560 values
3. restructure the inputs: size grid *derived* from `(a_min, a_max, mass_ratio)`
   + bulk density (never a list); dust initial condition in its own file,
   parallel to `abundances.in`; `Td(a)` and `T_CR,peak(a)` as explicit switchable
   prescriptions evaluated once on the grid
4. clean the RHS hot loop, symbolic Jacobian sparsity, `-O2`, re-verify

Coagulation lands only after all four are green. Target v1 grid: 5 nm – 0.5 µm,
mass-doubling, ~20 bins, 2-phase.
