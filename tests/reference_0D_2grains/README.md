# Reference case: 0D, 2 grain bins, 2-phase

The case used to prove that the extraction did not change the physics.

* 2 grain bins (2.33e-6 cm and 1.04e-3 cm), `multi_grain = 1`
* 2-phase (`is_3_phase = 0`)
* 1414 species, 20 log-spaced outputs, 1 yr -> 1e5 yr
* Td comes from the `dust_grid_table.in` column (tabulated grid; no temperature sub-switch)

Everything else (network, abundances, surface parameters) is the stock
nmgc-2.0 input set, taken unmodified.
