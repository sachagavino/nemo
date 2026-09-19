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

2. **Redistribution operator** — returns `C_ijk`, the fraction of the (i,j)
   collision *product* assigned to bin *k*. Podolak/Brauer: the weights are
   `(eps, 1-eps)` on grain **number**, and they conserve mass by construction
   (`eps*m_k + (1-eps)*m_{k+1} = m_i + m_j`), not because they are themselves
   mass fractions. Fragmentation later, as *another implementation of the same
   interface*, not as a special case bolted onto the first.

## Ice transfer is derived, never specified

When two grains coagulate, their ices go with them. The ice-transfer weights are
**derived from the redistribution operator** — they are never given their own
independent prescription. That decision is what keeps dust mass and ice mass
conserved together: if the ice weights could disagree with `C_ijk`, the two
budgets would silently drift apart and the diagnostics would not tell you which
one was wrong.

The default weighting is `C_ijk` itself: **ice-follows-particle**. The ice is
split with the *same* `(eps, 1-eps)` as the grain number, so ice and grains move
in lockstep: total ice is conserved automatically (`eps + (1-eps) = 1`), no bin
can receive ice without the grains to hold it, and no separate conservation guard
is needed. This is the choice adopted for v1. Two alternative weightings are
selectable behind the same interface, kept as diagnostics:

* **mass-weighted** — ice split in proportion to the *mass* each product bin
  receives, `∝ eps*m_k`. This is NOT `C_ijk`; it decouples ice from grain number
  and needs an explicit conservation guard to avoid spurious drift;
* **area-weighted** — ice split in proportion to product surface area,
  `∝ eps*a_k^2` — the weighting that keeps monolayer coverage
  `theta ∝ N_ice/a^2` continuous across the split. Houge & Krijt (2023) use it
  for *re-condensing* ice, where it reflects a real gas-to-grain flux ∝ surface;
  here the ice is merely reassigned across the two product bins, so it is a pure
  bookkeeping choice, not a flux, and (like mass-weighting) needs the guard.

All three conserve total ice. They differ only in *where* the ice lands across
the product bins, so the spread between them is a grid-resolution diagnostic: if
particle- and area-weighted disagree noticeably, the grid is too coarse to
resolve the ice transfer, not a sign that one of them is wrong.

## Diagnostics that must exist before the first coagulation term

* dust mass, `sum_k m_k Y_k`
* grain number budget
* total ice per species, summed over bins
* elemental conservation (already exists, see `check_conservation`)
