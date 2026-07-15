# `src/dust/` — where the dust size distribution will live

Nothing here is implemented yet, and nothing here is compiled. This directory
fixes the *shape* of the coagulation code so that adding it later does not
require touching the driver, the solver, or the chemistry.

## The state vector

The dust bins become part of the ODE state, integrated by DLSODES together with
the chemistry — **not** operator-split around it. `main.f90` therefore has no
sub-stepping loop, on purpose.

## Two separate interfaces, never one

Coagulation sits behind two independent objects. Keeping them apart is what
makes fragmentation a drop-in later rather than a rewrite.

1. **Kernel** — returns `K_ij`, the rate at which bins *i* and *j* collide.
   Depends on the collision physics (relative velocities, cross sections,
   sticking). One implementation = one physical prescription.

2. **Redistribution operator** — returns `C_ijk`, the mass fraction of an
   (i,j) collision that lands in bin *k*. Podolak/Brauer-style mass-conserving
   redistribution first; fragmentation later, as *another implementation of the
   same interface*, not as a special case bolted onto the first.

## Ice transfer is derived, never specified

When two grains coagulate, their ices go with them. The ice-transfer weights are

    f_ijk = C_ijk * m_k / (m_i + m_j)

They are **derived from the redistribution operator**. They must never be given
their own independent prescription: if `C_ijk` and `f_ijk` can disagree, dust
mass and ice mass will silently stop being conserved together, and the
diagnostics will not tell you which one is wrong.

## Diagnostics that must exist before the first coagulation term

* dust mass, `sum_k m_k Y_k`
* grain number budget
* total ice per species, summed over bins
* elemental conservation (already exists, see `check_conservation`)
