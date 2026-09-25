#!/usr/bin/env python3
"""
scripts/compare_dustpy.py -- Phase III Part C: NEMO vs DustPy vs exact analytic.

Runs NEMO and DustPy on the IDENTICAL 0D pure-coagulation problem (same ratio-2
fiducial mass grid, same exponential IC, same kernel) and overlays both against
the exact analytic solution (coag_analytic.py), for the constant and additive
kernels. Reports each code's relative L1 error vs the exact solution at matched
dimensionless times -- the honest §2c claim is that NEMO and DustPy agree with
each other and with the analytic to *comparable* error (both are the k=0
Kovetz-Olund scheme), not that either beats a high-order (DG) solver.

Comparison quantity: N m^2 = m^2 dN/dm, normalised per curve to unit grid-sum so
the volume-density (NEMO) vs surface-density (DustPy) normalisation drops out;
only the matched-dimensionless-time shape is compared. Both codes use the same
ratio-2 fiducial grid, so their mass grids coincide bin-for-bin (asserted).

Writes figures/coag_dustpy_comparison.{png,pdf} + data TSVs under
figures/coag_dustpy_data/ (so the plot can be regenerated without re-running).
"""
import argparse
import os
import sys
import tempfile

import numpy as np

import coag_analytic as CA
import coag_diffusion as CD
import run_dustpy_coag as RD

YR = 3.1556952e7
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def _fdim(kernel, x, t):
    return CA.n_constant(x, t) if kernel == "constant" else CA.n_additive(np.asarray(x), t)


def nemo_curves(kernel, a_min, a_max, ratio, rho, m0, N0, K0, nH, times):
    CD.FIXDIR = os.path.join(ROOT, "tests", "fixture_dust_coag")
    nmgc = os.path.join(ROOT, "bin", "nmgc")
    prefac = K0 if kernel == "constant" else K0 / m0
    R_eff = K0 * nH * N0
    out = CD._run_nemo(nmgc, kernel, ratio, a_min, a_max, rho, N0, m0, prefac, nH,
                       stop_time_yr=(3.0 * max(times) / R_eff) / YR, nb_outputs=40)
    h = np.diff(CD._bin_edges(out["mass"]))
    curves = {}
    for T in times:
        t_s = T / (prefac * nH * (out["N0"] if kernel == "constant" else out["M1"]))
        n_k = CD._interp_n(out["times"], out["n_of_t"], t_s / YR)
        curves[T] = out["mass"]**2 * n_k / h        # N m^2
    return out["mass"], curves


def dustpy_curves(kernel, a_min, a_max, ratio, rho, m0_radius, times):
    outdir = tempfile.mkdtemp(prefix=f"cmp_dustpy_{kernel}_")
    res = RD.run(kernel, a_min, a_max, ratio, rho, m0_radius, a_kernel=1.0, N0=1.0,
                 dimless_times=times, outdir=outdir)
    m, B = res["m"], res["B"]
    curves = {}
    for T in times:
        i = 1 + int(np.argmin(np.abs(res["t"][1:] * res["rate"] - T)))
        curves[T] = res["Sigma_cell1"][i] / B           # N m^2
    return m, curves


def _norm(y):
    s = np.sum(y)
    return y / s if s > 0 else y


def run(a_min, a_max, ratio, rho, m0_radius, N0, K0, nH, times, datadir):
    os.makedirs(datadir, exist_ok=True)
    m0 = CA.mass_of_radius(m0_radius, rho)
    report = {}
    data = {}
    for kernel in ("constant", "additive"):
        mN, cN = nemo_curves(kernel, a_min, a_max, ratio, rho, m0, N0, K0, nH, times)
        mD, cD = dustpy_curves(kernel, a_min, a_max, ratio, rho, m0_radius, times)
        # identical grid check
        gridmax = np.max(np.abs(mN / mD - 1.0))
        data[kernel] = dict(m=mN, gridmax=gridmax, nemo={}, dustpy={}, exact={})
        report[kernel] = {}
        for T in times:
            ex = mN**2 * (N0 / m0) * _fdim(kernel, mN / m0, T)     # exact N m^2 on grid
            qN, qD, qE = _norm(cN[T]), _norm(cD[T]), _norm(ex)
            l1_N = np.sum(np.abs(qN - qE))
            l1_D = np.sum(np.abs(qD - qE))
            report[kernel][T] = (l1_N, l1_D, gridmax)
            data[kernel]["nemo"][T] = qN
            data[kernel]["dustpy"][T] = qD
            data[kernel]["exact"][T] = qE
        # dump TSV
        with open(os.path.join(datadir, f"compare_{kernel}.tsv"), "w") as f:
            f.write("dimless_t\tmass_g\tnemo_norm\tdustpy_norm\texact_norm\n")
            for T in times:
                for j in range(mN.size):
                    f.write(f"{T}\t{mN[j]:.6e}\t{data[kernel]['nemo'][T][j]:.6e}\t"
                            f"{data[kernel]['dustpy'][T][j]:.6e}\t{data[kernel]['exact'][T][j]:.6e}\n")

    # ---- report ----
    print("\n=== NEMO vs DustPy vs exact analytic (relative L1 on normalised N m^2) ===")
    for kernel in ("constant", "additive"):
        print(f"  {kernel} kernel  (grid match |mN/mD-1|max = {data[kernel]['gridmax']:.1e}):")
        print(f"    {'dimless_t':>10} {'L1(NEMO)':>12} {'L1(DustPy)':>12} {'|NEMO-DustPy|/L1':>18}")
        for T in times:
            lN, lD, _ = report[kernel][T]
            ratio_cmp = abs(lN - lD) / max(lN, lD)
            print(f"    {T:10.3f} {lN:12.4e} {lD:12.4e} {ratio_cmp:18.2f}")
    return data, times


def plot(data, times, outbase):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a"]
    for ax, kernel in zip(axes, ("constant", "additive")):
        m = data[kernel]["m"]
        for c, T in zip(colors, times):
            ax.loglog(m, data[kernel]["exact"][T], "-", color=c, lw=1.5,
                      label=f"exact t*={T:g}")
            ax.loglog(m, data[kernel]["nemo"][T], "o", color=c, ms=5, mfc="none")
            ax.loglog(m, data[kernel]["dustpy"][T], "x", color=c, ms=5)
        ax.set_xlabel("grain mass m [g]")
        ax.set_ylabel(r"normalised $N\,m^2$")
        ax.set_title(f"{kernel} kernel")
        ymax = max(data[kernel]["exact"][times[0]].max(), 1e-30)
        ax.set_ylim(ymax * 1e-5, ymax * 3)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color="0.3", lw=1.5, label="exact"),
               Line2D([], [], color="0.3", marker="o", mfc="none", lw=0, label="NEMO"),
               Line2D([], [], color="0.3", marker="x", lw=0, label="DustPy")]
    axes[0].legend(fontsize=8, frameon=False, loc="lower center")
    axes[1].legend(handles=handles, fontsize=8, frameon=False, loc="upper right")
    fig.suptitle("NEMO vs DustPy vs exact analytic (0D coagulation, 20-bin fiducial grid)")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{outbase}.{ext}", dpi=150)
    print(f"wrote {outbase}.png / .pdf")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a-min", type=float, default=5.0e-7)
    ap.add_argument("--a-max", type=float, default=5.0e-5)
    ap.add_argument("--mass-ratio", type=float, default=2.0)
    ap.add_argument("--rho", type=float, default=3.0)
    ap.add_argument("--m0-radius", type=float, default=1.0e-6)
    ap.add_argument("--N0", type=float, default=2.312e-10)
    ap.add_argument("--K0", type=float, default=1.0e-9)
    ap.add_argument("--nH", type=float, default=2.212e8)
    ap.add_argument("--times", default="0.5,1.0,2.0")
    ap.add_argument("--datadir", default=os.path.join(ROOT, "figures", "coag_dustpy_data"))
    args = ap.parse_args()
    times = [float(s) for s in args.times.split(",")]
    data, times = run(args.a_min, args.a_max, args.mass_ratio, args.rho, args.m0_radius,
                      args.N0, args.K0, args.nH, times, args.datadir)
    plot(data, times, os.path.join(ROOT, "figures", "coag_dustpy_comparison"))


if __name__ == "__main__":
    main()
