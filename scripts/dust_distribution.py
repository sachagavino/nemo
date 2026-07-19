#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Reduction of dust_distribution.out into the two views a human wants to look at:
#   * sigma(a) = m * dn/dln a   -- the mass-weighted size distribution
#   * n(a, t)                   -- raw grain number per bin over time
# plus the dust-mass / grain-number conservation table.
#
# NEMO writes dust_distribution.out as RAW STATE (number per bin per H), on
# purpose: the mass weighting and the per-log-interval density are a reduction,
# and reductions live here in scripts/, not in the Fortran core. This script is
# that reduction.
#
# Columns of dust_distribution.out (one row per (time, bin)):
#   time[yr]  bin  radius[cm]  mass[g]  n_k[/H]  m_k*n_k[g/H]
#   dust_mass[g/H]  grain_number[/H]  dust_mass_sink[g/H]
#
# Usage (key=value args, no spaces; matches the other scripts here):
#   scripts/dust_distribution.py                      # reads ./dust_distribution.out
#   scripts/dust_distribution.py file=path/to/out ext=png
#   scripts/dust_distribution.py times=1,1e3,1e5      # only these output times [yr]
#
# Outputs: dust_sigma.<ext>, dust_number.<ext>, and a conservation table to stdout.

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')            # headless: write files, never require a display
import matplotlib.pyplot as plt

INFILE = 'dust_distribution.out'
EXT = 'pdf'
WANT_TIMES = None                # None = all output times

HELP = """AIM: reduce dust_distribution.out to sigma(a)=m*dn/dln a and n(a,t),
and print the dust-mass / grain-number conservation table.

Run in a folder that contains a finished NEMO run (or pass file=...).
Arguments (key=value, no spaces):
  file=dust_distribution.out : input file
  ext=pdf                    : output image extension (pdf, png, ...)
  times=1,1e3,1e5            : restrict to these output times [yr] (nearest match)
  help                       : this message
"""

for arg in sys.argv[1:]:
    if arg == 'help':
        print(HELP); sys.exit(0)
    if '=' not in arg:
        print("unknown argument '%s'\n%s" % (arg, HELP)); sys.exit(1)
    key, val = arg.split('=', 1)
    if key == 'file':
        INFILE = val
    elif key == 'ext':
        EXT = val
    elif key == 'times':
        WANT_TIMES = [float(x) for x in val.split(',')]
    else:
        print("unknown key '%s'\n%s" % (key, HELP)); sys.exit(1)

if not os.path.exists(INFILE):
    print("No '%s' here. Run NEMO first, or pass file=... (see 'help')." % INFILE)
    sys.exit(1)

# --- load -------------------------------------------------------------------
# columns: 0 time, 1 bin, 2 radius, 3 mass, 4 n_k, 5 m_k*n_k, 6 dust_mass,
#          7 grain_number, 8 dust_mass_sink
data = np.loadtxt(INFILE, comments='!')
if data.ndim == 1:
    data = data[None, :]

times_all = np.unique(data[:, 0])
nbin = int(data[:, 1].max())

# reshape into [n_time, n_bin] using the (time, bin) grid
def grid(col):
    out = np.empty((times_all.size, nbin))
    for it, t in enumerate(times_all):
        rows = data[data[:, 0] == t]
        rows = rows[np.argsort(rows[:, 1])]     # by bin index
        out[it, :] = rows[:, col]
    return out

radius = grid(2)[0, :]          # [cm], time-independent grid
mass = grid(3)[0, :]            # [g]
n_kt = grid(4)                  # n_k(t) [/H]           -> the n(a,t) view
dust_mass = grid(6)[:, 0]       # [g/H] per time (repeated across bins)
grain_number = grid(7)[:, 0]    # [/H]  per time
sink = grid(8)[:, 0]            # [g/H] per time

# --- sigma(a) = m * dn/dln a ------------------------------------------------
# n_k is the NUMBER integrated over bin k; dn/dln a = n_k / (dln a)_k. On the
# geometric NEMO grid (dln a) is uniform; compute it per bin from neighbours so
# this also works for a tabulated / non-uniform grid. Base choice (ln here) only
# rescales sigma by a constant (ln 10 for log10); the shape is what matters.
lna = np.log(radius)
dlna = np.empty_like(lna)
dlna[1:-1] = 0.5 * (lna[2:] - lna[:-2])
dlna[0] = lna[1] - lna[0]
dlna[-1] = lna[-1] - lna[-2]
sigma_kt = n_kt * mass[None, :] / dlna[None, :]     # [g/H per ln a]

# --- which times to draw ----------------------------------------------------
if WANT_TIMES is None:
    sel = np.arange(times_all.size)
else:
    sel = [int(np.argmin(np.abs(times_all - t))) for t in WANT_TIMES]
    sel = sorted(set(sel))

# --- plot: sigma(a) ---------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))
for it in sel:
    ax.plot(radius, sigma_kt[it, :], marker='o', ms=3,
            label='%.3g yr' % times_all[it])
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel(r'grain radius $a$ [cm]')
ax.set_ylabel(r'$\sigma(a) = m\,\mathrm{d}n/\mathrm{d}\ln a$  [g per H per $\ln a$]')
ax.set_title(r'Mass-weighted size distribution')
if len(sel) <= 12:
    ax.legend(fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig('dust_sigma.%s' % EXT)
print("wrote dust_sigma.%s" % EXT)

# --- plot: n(a, t) ----------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))
for it in sel:
    ax.plot(radius, n_kt[it, :], marker='o', ms=3,
            label='%.3g yr' % times_all[it])
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel(r'grain radius $a$ [cm]')
ax.set_ylabel(r'$n_k$ (number per bin, per H)')
ax.set_title(r'Grain number distribution $n(a,t)$ (raw state)')
if len(sel) <= 12:
    ax.legend(fontsize=8, ncol=2)
fig.tight_layout()
fig.savefig('dust_number.%s' % EXT)
print("wrote dust_number.%s" % EXT)

# --- conservation table -----------------------------------------------------
# This file IS the dust-mass / grain-number conservation diagnostic. Print the
# drift relative to the first output; with a static distribution it should sit
# at the floating-point floor.
dm0, gn0 = dust_mass[0], grain_number[0]
print("\n%-14s  %-22s  %-22s  %-12s" %
      ('time[yr]', 'dust_mass[g/H]', 'grain_number[/H]', 'sink[g/H]'))
for it in range(times_all.size):
    print("%-14.4e  %-22.15e  %-22.15e  %-12.4e" %
          (times_all[it], dust_mass[it], grain_number[it], sink[it]))
print("\nmax |Δ dust_mass|/dust_mass[0]    : %.3e" %
      (np.max(np.abs(dust_mass - dm0)) / dm0))
print("max |Δ grain_number|/grain_number[0]: %.3e" %
      (np.max(np.abs(grain_number - gn0)) / gn0))
if sink.max() > 1e-3 * dm0:
    print("WARNING: dust_mass_sink is non-negligible -> a_max is too low.")
