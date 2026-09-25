#!/usr/bin/env python3
"""
scripts/plot_compare_dustpy.py -- Phase III Part C figures (post-processing).

Pure plotting from the TSVs written by
  python3 scripts/compare_dustpy.py
(figures/coag_dustpy_data/compare_{constant,additive}.tsv). No NEMO, no DustPy --
edit freely and re-run to restyle without re-running the codes.

Produces two figures:
  (1) coag_dustpy_l1.{png,pdf}        -- the L&L-style L1 comparison: each code's
      relative L1 error vs the exact analytic as a function of dimensionless time,
      both kernels. Quantitative "NEMO and DustPy reach the analytic to comparable
      error" statement.
  (2) coag_nemo_vs_dustpy.{png,pdf}   -- the broad-community figure: NEMO and DustPy
      grain-size distributions overlaid with the exact analytic, both kernels.
      "Two independent codes, same answer" at a glance.

L1 is computed here from the normalised distributions (sum_k |q_code - q_exact|),
the continuous-family metric used for the gate; it is convention-robust (unlike the
point-evaluated disc-L1/L2, which are diagnostics only).
"""
import argparse
import csv
import os
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
KERNELS = ("constant", "additive")
TAU_COLORS = {0.5: "#1b9e77", 1.0: "#d95f02", 2.0: "#7570b3", 0.01: "#999999"}
KCOL = {"constant": "#1f77b4", "additive": "#2ca02c"}


def load(datadir, kernel):
    """Return {t: (mass[], nemo[], dustpy[], exact[])} from compare_<kernel>.tsv."""
    by_t = defaultdict(lambda: ([], [], [], []))
    with open(os.path.join(datadir, f"compare_{kernel}.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            t = float(r["dimless_t"])
            m, n, d, e = by_t[t]
            m.append(float(r["mass_g"])); n.append(float(r["nemo_norm"]))
            d.append(float(r["dustpy_norm"])); e.append(float(r["exact_norm"]))
    return {t: tuple(np.array(a) for a in v) for t, v in by_t.items()}


def figure_l1(datadir, ext):
    """L1(code vs exact) as a function of dimensionless time, both kernels."""
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for kernel in KERNELS:
        data = load(datadir, kernel)
        ts = sorted(data)
        l1N, l1D = [], []
        for t in ts:
            m, n, d, e = data[t]
            l1N.append(np.sum(np.abs(n - e)))
            l1D.append(np.sum(np.abs(d - e)))
        ax.plot(ts, l1N, "o-", color=KCOL[kernel], label=f"{kernel}: NEMO")
        ax.plot(ts, l1D, "x--", color=KCOL[kernel], label=f"{kernel}: DustPy")
    ax.set_xlabel("dimensionless time  T (constant) / \u03c4 (additive)")
    ax.set_ylabel("relative L1 error vs exact analytic")
    ax.set_title("NEMO and DustPy reach the analytic to comparable error\n"
                 "(20-bin fiducial grid, k=0 scheme both codes)")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(True, alpha=0.3)
    out = os.path.join(ROOT, "figures", f"coag_dustpy_l1.{ext}")
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    print(f"wrote {out}")


def figure_overlay(datadir, ext, times=(0.5, 1.0, 2.0)):
    """NEMO vs DustPy grain-size distributions overlaid with exact analytic."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for ax, kernel in zip(axes, KERNELS):
        data = load(datadir, kernel)
        ymax = 0.0
        for t in times:
            if t not in data:
                continue
            m, n, d, e = data[t]
            c = TAU_COLORS.get(t, "#333333")
            ax.loglog(m, e, "-", color=c, lw=1.5, label=f"exact t*={t:g}")
            ax.loglog(m, n, "o", color=c, ms=5, mfc="none")
            ax.loglog(m, d, "x", color=c, ms=5)
            ymax = max(ymax, e.max())
        ax.set_xlabel("grain mass m [g]")
        ax.set_ylabel(r"normalised $N\,m^2$")
        ax.set_title(f"{kernel} kernel")
        ax.set_ylim(ymax * 1e-5, ymax * 3)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color="0.3", lw=1.5, label="exact analytic"),
               Line2D([], [], color="0.3", marker="o", mfc="none", lw=0, label="NEMO"),
               Line2D([], [], color="0.3", marker="x", lw=0, label="DustPy")]
    axes[0].legend(fontsize=8, frameon=False, loc="lower center")
    axes[1].legend(handles=handles, fontsize=8, frameon=False, loc="upper right")
    fig.suptitle("NEMO vs DustPy: same 0D coagulation problem, two independent codes")
    out = os.path.join(ROOT, "figures", f"coag_nemo_vs_dustpy.{ext}")
    fig.tight_layout(); fig.savefig(out, dpi=150); plt.close(fig)
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", default=os.path.join(ROOT, "figures", "coag_dustpy_data"))
    ap.add_argument("--ext", default="png", choices=["png", "pdf"])
    ap.add_argument("--overlay-times", default="0.5,1.0,2.0")
    args = ap.parse_args()
    figure_l1(args.datadir, args.ext)
    figure_overlay(args.datadir, args.ext,
                   times=tuple(float(s) for s in args.overlay_times.split(",")))


if __name__ == "__main__":
    main()
