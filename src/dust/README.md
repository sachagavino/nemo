# `src/dust/` — where the dust size distribution lives

As of stage 3 this directory holds its first live, compiled module,
`dust_grid.f90`: the size grid is now a derived, first-class object (see
`docs/STAGE3_GRID.md`). The coagulation code — `dustevolution.f90` and the two
interfaces below — is still to come and still shapes the layout: adding it must
not require touching the driver, the solver, or the chemistry.

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
**derived from the redistribution operator** — they are never given their own
independent prescription. That decision is what keeps dust mass and ice mass
conserved together: if the ice weights could disagree with `C_ijk`, the two
budgets would silently drift apart and the diagnostics would not tell you which
one was wrong.

The default weighting is `C_ijk` itself: ice follows particle mass. Since
`sum_k C_ijk = 1` by construction, `C_ijk` already conserves ice — there is no
extra `m_k/(m_i+m_j)` factor to apply. Two alternative weightings are selectable
behind the same interface:

* **mass-weighted** — ice split in proportion to the mass each product bin
  receives (this *is* `C_ijk` when redistribution is mass-conserving);
* **area-weighted** — ice split in proportion to product surface area, which is
  the natural choice when the ice is re-condensing onto fresh surface rather
  than riding along with the mass (Houge & Krijt 2023 chose area-weighting for
  re-condensation for exactly this reason — worth a code comment when we get
  there).

All three conserve total ice. They differ only in *where* the ice lands across
the product bins, so the spread between them is a grid-resolution diagnostic: if
mass-weighted and area-weighted disagree noticeably, the grid is too coarse to
resolve the ice transfer, not a sign that one of them is wrong.

## Diagnostics that must exist before the first coagulation term

* dust mass, `sum_k m_k Y_k`
* grain number budget
* total ice per species, summed over bins
* elemental conservation (already exists, see `check_conservation`)
