#!/usr/bin/env python3
"""
scripts/plot_coag_diffusion.py -- Phase III Part B 1e figure (post-processing).

Pure plotting: reads the TSVs written by
  python3 scripts/coag_diffusion.py --stage figure-data
and renders the two-panel numerical-diffusion figure. No NEMO, no code units --
edit this freely and re-run to restyle without re-running the sims.

  Left  : additive kernel, mass density g=x f vs grain mass at tau=0.5,1,2, for a
          coarse (mass_ratio 2) and a fine (mass_ratio 1.3) grid vs the analytic
          solution -- shows k=0 over-diffusion into the large-mass tail shrinking
          with resolution.
  Right : relative L1 error (continuous + discrete, L&L 2021 Eq.40/41) vs bin
          count at the evolved time, both kernels, monotone -> converges.

Usage:  python3 scripts/plot_coag_diffusion.py [--datadir DIR] [--ext png|pdf]
"""
import argparse
import csv
import os
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

TAU_COLORS = {0.5: "#1b9e77", 1.0: "#d95f02", 2.0: "#7570b3"}


def _read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _bin_edges(m):
    e = np.empty(m.size + 1)
    e[1:-1] = np.sqrt(m[:-1] * m[1:])
    e[0] = m[0] / np.sqrt(m[1] / m[0])
    e[-1] = m[-1] * np.sqrt(m[-1] / m[-2])
    return e


def plot(datadir, ext):
    nemo = _read_tsv(os.path.join(datadir, "panelA_nemo.tsv"))
    ana = _read_tsv(os.path.join(datadir, "panelA_analytic.tsv"))
    pb = _read_tsv(os.path.join(datadir, "panelB.tsv"))

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 4.4))

    # ---------- Panel A ----------
    # analytic curves
    ana_by_tau = defaultdict(list)
    for r in ana:
        ana_by_tau[float(r["tau"])].append((float(r["mass_g"]), float(r["g_massdensity"])))
    gmax = 0.0
    for tau, pts in sorted(ana_by_tau.items()):
        pts = np.array(sorted(pts))
        axA.loglog(pts[:, 0], pts[:, 1], "-", color=TAU_COLORS[tau], lw=1.7,
                   label=f"analytic τ={tau:g}")
        gmax = max(gmax, pts[:, 1].max())
    # NEMO points: coarse = open circles, fine = filled small
    nemo_by = defaultdict(lambda: defaultdict(list))   # ratio -> tau -> [(m,g)]
    ratios_seen = set()
    for r in nemo:
        ratio = float(r["ratio"]); tau = float(r["tau"])
        nemo_by[ratio][tau].append((float(r["mass_g"]), float(r["g_massdensity"])))
        ratios_seen.add((ratio, int(r["nbins"])))
    styles = {}   # coarse vs fine markers
    ratios_sorted = sorted(ratios_seen, reverse=True)         # coarse (big ratio) first
    marker_for = {ratios_sorted[0][0]: ("o", "none"), ratios_sorted[-1][0]: (".", None)}
    for (ratio, nb) in ratios_sorted:
        mk, mfc = marker_for.get(ratio, ("x", None))
        for tau, pts in nemo_by[ratio].items():
            pts = np.array(sorted(pts))
            axA.loglog(pts[:, 0], pts[:, 1], mk, color=TAU_COLORS[tau], ms=5,
                       mfc=(TAU_COLORS[tau] if mfc is None else mfc), lw=0)
    # faint bin-edge gridlines for the COARSE grid (L&L style)
    coarse_ratio = ratios_sorted[0][0]
    coarse_m = np.array(sorted({float(r["mass_g"]) for r in nemo
                                if float(r["ratio"]) == coarse_ratio}))
    for e in _bin_edges(coarse_m):
        axA.axvline(e, color="0.9", lw=0.5, zorder=0)
    # marker legend (resolutions)
    from matplotlib.lines import Line2D
    res_handles = [
        Line2D([], [], color="0.3", marker="o", mfc="none", lw=0,
               label=f"NEMO coarse ({ratios_sorted[0][1]} bins)"),
        Line2D([], [], color="0.3", marker=".", lw=0,
               label=f"NEMO fine ({ratios_sorted[-1][1]} bins)"),
    ]
    lg1 = axA.legend(loc="lower left", fontsize=8, frameon=False)
    axA.add_artist(lg1)
    axA.legend(handles=res_handles, loc="upper right", fontsize=8, frameon=False)
    axA.set_xlabel("grain mass m [g]")
    axA.set_ylabel(r"mass density $g = x\,f$  [g / H]")
    axA.set_title("Additive kernel: over-diffusion shrinks with resolution")
    axA.set_ylim(gmax * 1e-6, gmax * 3)                       # clip the noise tail

    # ---------- Panel B ----------
    by_kernel = defaultdict(list)
    for r in pb:
        by_kernel[r["kernel"]].append((int(r["nbins"]), float(r["cont_L1"]), float(r["disc_L1"])))
    kstyle = {"constant": ("s", "#1f77b4"), "additive": ("^", "#2ca02c")}
    nb_all = []
    for kernel, rows in by_kernel.items():
        rows = np.array(sorted(rows))
        mk, col = kstyle[kernel]
        axB.loglog(rows[:, 0], rows[:, 1], mk + "-", color=col, label=f"{kernel} cont-L1")
        axB.loglog(rows[:, 0], rows[:, 2], mk + "--", color=col, mfc="none",
                   label=f"{kernel} disc-L1")
        nb_all += list(rows[:, 0])
        e0, nb0 = rows[0, 1], rows[0, 0]
    nb_ref = np.array(sorted(set(nb_all)), float)
    axB.loglog(nb_ref, e0 * (nb_ref / nb0)**-1.0, ":", color="0.5",
               label="slope −1 (first-order k=0)")
    axB.set_xlabel("number of bins")
    axB.set_ylabel("relative L1 error on g")
    axB.set_title("Convergence at evolved time (T=2 const, τ=1 add)")
    axB.legend(fontsize=8, frameon=False)

    fig.tight_layout()
    out = os.path.join(os.path.dirname(datadir) or ".", f"coag_diffusion.{ext}")
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", default="figures/coag_diffusion_data")
    ap.add_argument("--ext", default="png", choices=["png", "pdf"])
    args = ap.parse_args()
    plot(args.datadir, args.ext)


if __name__ == "__main__":
    main()
