#!/usr/bin/env bash
# mk.sh <dir> <rtol> [<split_dt_yr> <A|B> <CDC|DCD>]
# Stages the Phase III target regime: reduced_CHO, chemistry on, Brownian kernel, 20 bins
# (5 nm - 0.40 um, mass ratio 2, MRN IC), parametric Td 25 -> 13 K (tabulated path; table =
# derived+MRN export with the Td column replaced), n_H = 2.212e8, T_gas = 27.41 K (fixture),
# dtg = 1e-2, T_run = 1e4 yr, linear outputs every 1e3 yr, symbolic sparsity (mf=21).
set -eu
NEMO_ROOT=${NEMO_ROOT:-$(cd "$(dirname "$0")/../.." && pwd)}; HERE=$(cd "$(dirname "$0")" && pwd)
d=$1; mkdir -p "$d"
cp $NEMO_ROOT/networks/reduced_CHO/*.in "$d"/; cp $NEMO_ROOT/inputs/element.in "$d"/
cp $HERE/dust_grid_table_Td25-13.in "$d"/dust_grid_table.in
sed -e "s/^coagulation = .*/coagulation = 1/" -e "s/^coagulation_kernel = .*/coagulation_kernel = brownian/" \
    -e "s/^a_max = .*/a_max =  5.000E-05/" -e "s/^sparsity = .*/sparsity = symbolic/" \
    -e "s/^sparsity_check_steps = .*/sparsity_check_steps = 0/" \
    -e "s/^initial_dtg_mass_ratio = .*/initial_dtg_mass_ratio =  1.000E-02/" \
    -e "s/^initial_gas_density = .*/initial_gas_density = 2.212E+08/" \
    -e "s/^dust_grid_source = .*/dust_grid_source = tabulated/" -e "s/^dust_ic = .*/dust_ic = tabulated/" \
    -e "s/^relative_tolerance = .*/relative_tolerance = $2/" -e "s/^output_type = .*/output_type = linear/" \
    -e "s/^start_time = .*/start_time = 0.0/" -e "s/^nb_outputs = .*/nb_outputs = 11/" \
    -e "s/^stop_time = .*/stop_time = 1.0E+04/" \
    $NEMO_ROOT/tests/fixture_rung3/parameters.in > "$d"/parameters.in
if [ $# -ge 5 ]; then printf "split_mode = 1\nsplit_dt = %s\nsplit_variant = %s\nsplit_order = %s\n" $3 $4 $5 >> "$d"/parameters.in; fi
