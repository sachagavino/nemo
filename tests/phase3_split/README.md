# Phase III — coupled vs operator-split (money-plot) tooling

Code: `split_mode` in `parameters.in` (default 0 = coupled production path, bit-identical to
`phase2-coagulation` @352b51c). Split keys: `split_mode = 1`, `split_dt` [yr],
`split_variant` (A: {coag+ice}|{chem}; B: {coag+chem}|{ice transport}),
`split_order` (CDC: chem(dt/2)-dust(dt)-chem(dt/2); DCD). Two-phase only.
Debug: `NEMO_SPLIT_SELFCHECK=1` prints RHS/Jacobian additivity (f_C+f_D=f_full,
J_C+J_D=J_full) and dust-block Jacobian vs central FD of the masked RHS, once per run.

## Files
- `mk.sh <dir> <rtol> [dt A|B CDC|DCD]` — stage the target regime (see header).
- `queue.sh <dir>...` — run sequentially; `STATUS` = ok / FAIL (fails a split run whose log
  lacks the split audit banner).
- `compare.py <ref> <run>` — per-time gate table + conservation drift.
  Gate set: gas+ice species with Y_ref > 1e-12, plus ALL grain bins (neutral+negative totals),
  plus per-bin JCO (JkCO), gas CO, total CO ice. L1/L2 over the above-cutoff gas+ice set.
- `collect.py <ref> <out.tsv> <dir>...` — error metrics (max over comparison times) + cost
  (NST, NFE, NJE, NLU, cold starts, macro-steps, restarts, CPU) per run.
- `show.py <tsv> [cols]`, `figs.py` — table / draft figures.
- `dust_grid_table_Td25-13.in` — 20-bin derived+MRN grid, Td linear in log a, 25 -> 13 K.
- `step0/` — Step-0 regime-scan analysis (`analyze.py`, `summary.py`, `mktable.py`) and the
  exported derived grid used to build the Td table.
- `data_Tgas21/` — **CURRENT** results at T_gas = 21.233896 K (area-weighted mean Td of the
  initial distribution): `sweep_all.tsv`, `wp_all.tsv`, figures, and `gate_order_c1_topbin.txt`
  (decided gate, convergence orders incl. t<=8 kyr, c1 smoothed-kink test, top-bin ice).
- `gate.py <ref> <run>` — the decided correctness gate (design note sec. 1): (i) per-bin CO ice,
  gas CO, total CO ice above Y>1e-12; (ii) all grain bins; (iii) conservation (< 1e-10 relative,
  executor-proposed bound); (iv) mask audit. All-species max reported as a diagnostic only.
- `order.py <ref> <prefix>` — measured local convergence orders, full window and t <= 8 kyr.
- `topbin.py <ref> <run>` — top-bin CO-ice share, monolayers per bin, per-bin split error profile.
- `data/` — SUPERSEDED first sweep at the fixture T_gas = 27.41 K (kept for the record).
- `data/sweep_all.tsv` — errors vs coupled rtol-1e-8 reference (A/B x CDC/DCD, dt 25..1000 yr).
- `data/wp_all.tsv` — work-precision: errors vs coupled rtol-1e-10 reference; coupled rtol
  1e-4..1e-8; split A/B CDC at rtol 1e-4, 1e-6, 1e-8.
- `data/fig_accuracy.png`, `data/fig_workprecision.png` — draft figures (unstyled).

Reproduction check: `mk.sh` + `queue.sh` reproduce the committed sweep runs bit-for-bit.
Runtime (1 core): coupled rtol 1e-8 ~21 s; split dt=1000 ~24 s, dt=25 ~270 s.

## T_gas = 21.233896 K rerun (design-thread decision note)
`data_Tgas21/` supersedes `data/` (which is at the fixture T_gas = 27.41 K; kept for the record).
- `sweep_all.tsv` (vs coupled rtol 1e-8), `wp_all.tsv` (vs coupled rtol 1e-10), figures.
- `gate_order_c1_topbin.txt`: decided gate (gate.py) for A CDC dt=25 yr, measured convergence
  orders incl. t <= 8 kyr (order.py), the c1 smoothed-kink test, top-bin ice analysis (topbin.py).
- `gate.py <ref> <run>`: (i) per-bin CO ice / gas CO / total CO ice above Y>1e-12, (ii) grain-bin
  totals, both < 1e-4; (iii) conservation drift < 1e-10 relative (executor-proposed bound);
  (iv) mask audit from the run log. The all-species max is printed as a diagnostic, not gated.
- The c1 test used a THROWAWAY build (smoothed H/H2 1-ML sticking switch, w = 0.25 ML, and
  smoothed photodesorption cap, w = 0.5 ML); that build is intentionally NOT committed.
