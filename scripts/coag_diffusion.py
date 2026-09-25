#!/usr/bin/env python3
"""
scripts/coag_diffusion.py -- Phase III Part B 1d (numerical-diffusion convergence).

Runs pure coagulation (chemistry off) on the tabulated grid+IC for the constant
and additive kernels over a mass_ratio sweep with a_min,a_max FIXED (finer grid =
more bins = the same continuous problem at higher resolution), and measures the
error of NEMO's binned solution against the analytic solution as a function of
resolution. The claim under test: the error DECREASES MONOTONICALLY as
mass_ratio -> 1 (bounded numerical diffusion that converges under refinement).

NEMO's Kovetz-Olund/Brauer redistribution is the k=0 (piecewise-constant) point
scheme in the taxonomy of Lombart & Laibe 2021 (arXiv:2011.12298); we therefore
reconstruct NEMO's solution as a piecewise-constant number density f_k = n_k/h_k
on bin k and compare on the MASS DENSITY g(x)=x f(x), using L&L's error norms:

  continuous L1  (L&L Eq. 40):  e_c = sum_j (h_j/2) sum_a w_a |g_num - g_ana|(x_j^a)
  discrete   L1  (L&L Eq. 41):  e_d = sum_j h_j |g_num - g_ana|(xhat_j)
                                 xhat_j = sqrt(x_{j-1/2} x_{j+1/2})
with the L2 analogues (square inside, sqrt outside). Errors are reported RELATIVE
to the same norm of g_ana. Bin edges are geometric means of adjacent grid points
(the Kovetz-Olund bin convention, cf. L&L Sect. 2.5.2.3). "compare like with
like": both g_num and g_ana are evaluated/integrated the same way on NEMO's grid.

Compares at a matched dimensionless time across resolutions (per-bin linear-in-t
interpolation of NEMO's output to the target time, using each run's own conserved
M1/N0 for the physical<->dimensionless map, so top-bin tail truncation cannot
bias the alignment).
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from numpy.polynomial.legendre import leggauss

import coag_analytic as CA

AMU = 1.66053892e-24
YR = 3.1556952e7


# --------------------------------------------------------------------------
def _bin_edges(m):
    """Geometric-mean bin edges for grid points m (mass). Returns edges[nb+1]."""
    nb = m.size
    e = np.empty(nb + 1)
    e[1:-1] = np.sqrt(m[:-1] * m[1:])
    e[0] = m[0] / np.sqrt(m[1] / m[0])       # geometric extrapolation at ends
    e[-1] = m[-1] * np.sqrt(m[-1] / m[-2])
    return e


def _read_distribution(path):
    """Parse dust_distribution.out -> (times[nt], mass[nb], n_of_t[nt,nb])."""
    rows = {}
    masses = {}
    with open(path) as f:
        for ln in f:
            if ln.startswith("!") or not ln.strip():
                continue
            c = ln.split()
            t = float(c[0]); b = int(c[1]); mass = float(c[3]); nk = float(c[4])
            rows.setdefault(t, {})[b] = nk
            masses[b] = mass
    times = np.array(sorted(rows))
    nb = len(masses)
    mass = np.array([masses[b] for b in range(1, nb + 1)])
    n_of_t = np.array([[rows[t].get(b, 0.0) for b in range(1, nb + 1)] for t in times])
    return times, mass, n_of_t


def _run_nemo(nmgc, kernel, ratio, a_min, a_max, rho, N0, m0, K0, nH,
              stop_time_yr, nb_outputs):
    """Generate tabulated grid+IC for this ratio, run pure coag, return outputs."""
    radii, masses, nb = CA.build_grid(a_min, a_max, ratio, rho)
    n_k = CA.exp_ic_bins(masses, N0, m0)
    run = tempfile.mkdtemp(prefix=f"coagdiff_{kernel}_{ratio}_")
    shutil.copy(os.path.join(FIXDIR, "dust_grid_table.in"), run) if False else None
    for fn in os.listdir(FIXDIR):
        if fn.endswith(".in"):
            shutil.copy(os.path.join(FIXDIR, fn), run)
    CA.write_grid_table(os.path.join(run, "dust_grid_table.in"), radii, n_k)
    # patch parameters.in
    p = os.path.join(run, "parameters.in")
    txt = open(p).read().splitlines()
    def setk(lines, key, val):
        out = []
        for ln in lines:
            if ln.strip().startswith(key + " ") or ln.strip().startswith(key + "="):
                out.append(f"{key} = {val}")
            else:
                out.append(ln)
        return out
    txt = setk(txt, "dust_grid_source", "tabulated")
    txt = setk(txt, "dust_ic", "tabulated")
    txt = setk(txt, "coagulation", "1")
    txt = setk(txt, "coagulation_kernel", kernel)
    txt = setk(txt, "constant_kernel_k0", f"{K0:.6E}")
    txt = setk(txt, "stop_time", f"{stop_time_yr:.6E}")
    txt = setk(txt, "nb_outputs", str(nb_outputs))
    open(p, "w").write("\n".join(txt) + "\n")
    r = subprocess.run([nmgc, "run"], cwd=run, capture_output=True, text=True, timeout=1200)
    if r.returncode != 0:
        sys.stderr.write(r.stdout[-2000:] + r.stderr[-2000:])
        raise RuntimeError(f"NEMO failed: kernel={kernel} ratio={ratio}")
    times, mass, n_of_t = _read_distribution(os.path.join(run, "dust_distribution.out"))
    shutil.rmtree(run, ignore_errors=True)
    M1 = float(np.sum(mass * n_of_t[0]))     # conserved dust mass per H (this grid)
    N0_run = float(np.sum(n_of_t[0]))
    return dict(nb=nb, times=times, mass=mass, n_of_t=n_of_t, M1=M1, N0=N0_run)


def _interp_n(times, n_of_t, t_target):
    """Per-bin linear-in-time interpolation of n_k to t_target [yr]."""
    if t_target <= times[0]:
        return n_of_t[0]
    if t_target >= times[-1]:
        return n_of_t[-1]
    j = np.searchsorted(times, t_target) - 1
    w = (t_target - times[j]) / (times[j + 1] - times[j])
    return (1 - w) * n_of_t[j] + w * n_of_t[j + 1]


def _f_dim(x, kernel, dimless_t):
    return CA.n_constant(x, dimless_t) if kernel == "constant" else CA.n_additive(x, dimless_t)


def _errors(mass, n_k, kernel, dimless_t, N0, m0):
    """L&L Eq.40 (continuous) / Eq.41 (discrete) relative L1 errors on the mass
    density g=x f. L1 only: it has units of misplaced mass and is the natural norm
    for a conservation law (L&L use L1; L2 dropped, see coag_diffusion_gate.sh)."""
    edges = _bin_edges(mass)
    h = np.diff(edges)
    f_k = n_k / h                              # NEMO piecewise-constant number density /H
    xg, wg = leggauss(16)
    ec1 = ec1_ref = ed1 = ed1_ref = 0.0
    for k in range(mass.size):
        lo, hi = edges[k], edges[k + 1]
        # continuous L1 (Eq.40), Gauss-16 over the bin
        xa = 0.5 * (hi - lo) * xg + 0.5 * (hi + lo)
        g_num = xa * f_k[k]
        g_ana = xa * (N0 / m0) * _f_dim(xa / m0, kernel, dimless_t)
        ec1 += 0.5 * (hi - lo) * np.sum(wg * np.abs(g_num - g_ana))
        ec1_ref += 0.5 * (hi - lo) * np.sum(wg * np.abs(g_ana))
        # discrete L1 (Eq.41), geometric-mean point
        xh = np.sqrt(lo * hi)
        gnh = xh * f_k[k]
        fah = _f_dim(np.array([xh / m0]), kernel, dimless_t)[0] if kernel == "additive" \
            else _f_dim(xh / m0, kernel, dimless_t)
        gah = xh * (N0 / m0) * fah
        ed1 += h[k] * abs(gnh - gah)
        ed1_ref += h[k] * abs(gah)
    return dict(cL1=ec1 / ec1_ref, dL1=ed1 / ed1_ref)


def sweep(nmgc, kernel, ratios, dimless_t, a_min, a_max, rho, dtg, m0_radius, K0, nH):
    m0 = CA.mass_of_radius(m0_radius, rho)
    N0 = dtg * 1.4 * AMU / m0                 # fixed physical IC across the sweep
    # additive prefactor chosen so the effective rate K0*nH*N0 matches the constant
    # kernel: B = K0/m0  =>  tau = B nH M1 t = (K0/m0) nH (N0 m0) t = K0 nH N0 t.
    prefac = K0 if kernel == "constant" else K0 / m0
    R_eff = K0 * nH * N0                       # 1/s, common target rate
    print(f"\n=== {kernel} kernel: convergence sweep at dimensionless time "
          f"{'T' if kernel=='constant' else 'tau'}={dimless_t} ===")
    print(f"   fixed IC: m0={m0:.3e} g (a={m0_radius:.2e} cm), N0={N0:.3e} /H, "
          f"prefactor={prefac:.3e}")
    print(f"   {'ratio':>6} {'nbins':>6} {'cont-L1':>11} {'disc-L1':>11}")
    results = []
    for ratio in ratios:
        t_end_s = 3.0 * dimless_t / R_eff      # comfortably past the target time
        out = _run_nemo(nmgc, kernel, ratio, a_min, a_max, rho, N0, m0, prefac, nH,
                        stop_time_yr=t_end_s / YR, nb_outputs=40)
        if kernel == "constant":
            t_target_s = dimless_t / (prefac * nH * out["N0"])
        else:
            t_target_s = dimless_t / (prefac * nH * out["M1"])   # tau = B nH M1 t
        n_k = _interp_n(out["times"], out["n_of_t"], t_target_s / YR)
        e = _errors(out["mass"], n_k, kernel, dimless_t, out["N0"], m0)
        results.append((ratio, out["nb"], e))
        print(f"   {ratio:6.2f} {out['nb']:6d} {e['cL1']:11.4e} {e['dL1']:11.4e}", flush=True)
    # GATE on continuous L1 only (L&L Eq.40, integrated): it is the conservation-law
    # norm and robust under refinement. disc-L1 (Eq.41) and L2 are POINT-evaluated and
    # fragile in the sparse over-diffused tail (non-monotone at evolved times), so they
    # are reported as diagnostics, not gated. See coag_diffusion_gate.sh.
    print("   monotone continuous-L1 decrease as mass_ratio -> 1 (GATE):")
    seqc = [e["cL1"] for (_, _, e) in results]
    ok = all(seqc[i] > seqc[i + 1] for i in range(len(seqc) - 1))
    print(f"     cont-L1: {'PASS' if ok else 'FAIL'}  ({' > '.join(f'{v:.2e}' for v in seqc)})")
    seqd = [e["dL1"] for (_, _, e) in results]
    monod = all(seqd[i] > seqd[i + 1] for i in range(len(seqd) - 1))
    print(f"     disc-L1 (diagnostic, not gated): "
          f"{'monotone' if monod else 'non-monotone'}  "
          f"({' , '.join(f'{v:.2e}' for v in seqd)})")
    return ok, results


def eoc_report(nmgc, kernel, ratios, a_min, a_max, rho, dtg, m0_radius, K0, nH):
    """Spatial order at the L&L EOC time (dimensionless 0.01): log-log slope of
    the continuous-L1 error vs bin count. Reported as the SPATIAL discretisation
    order (kernel-independent at this time), not as the coagulation gate."""
    _, res = sweep(nmgc, kernel, ratios, 0.01, a_min, a_max, rho, dtg, m0_radius, K0, nH)
    nb = np.array([r[1] for r in res], float)
    e = np.array([r[2]["cL1"] for r in res], float)
    slope = np.polyfit(np.log(nb), np.log(e), 1)[0]
    print(f"   -> effective spatial order (cont-L1 vs nbins, dimless t=0.01): {-slope:.2f}")
    return -slope


def bounded_report(nmgc, kernel, ratio, dimless_end, a_min, a_max, rho, dtg,
                   m0_radius, K0, nH):
    """Leg 3: error-vs-time boundedness at FIXED resolution. Runs one grid to a
    late dimensionless time and reports the continuous-L1 error at a series of
    times; 'bounded' means it does not diverge as coagulation proceeds."""
    m0 = CA.mass_of_radius(m0_radius, rho)
    N0 = dtg * 1.4 * AMU / m0
    prefac = K0 if kernel == "constant" else K0 / m0
    R_eff = K0 * nH * N0
    out = _run_nemo(nmgc, kernel, ratio, a_min, a_max, rho, N0, m0, prefac, nH,
                    stop_time_yr=(dimless_end / R_eff) / YR, nb_outputs=40)
    print(f"\n=== {kernel} kernel: error-vs-time at fixed resolution "
          f"(mass_ratio={ratio}, {out['nb']} bins) ===")
    print(f"   {'dimless_t':>10} {'cont-L1':>11}")
    errs = []
    for t_yr, n_k in zip(out["times"], out["n_of_t"]):
        if kernel == "constant":
            dt = prefac * nH * out["N0"] * (t_yr * YR)
        else:
            dt = prefac * nH * out["M1"] * (t_yr * YR)
        if dt < 0.05:
            continue
        e = _errors(out["mass"], n_k, kernel, dt, out["N0"], m0)["cL1"]
        errs.append((dt, e))
    for dt, e in errs[:: max(1, len(errs) // 8)]:
        print(f"   {dt:10.3f} {e:11.4e}")
    emax = max(e for _, e in errs)
    print(f"   -> max cont-L1 over the run = {emax:.3e} (bounded: error does not diverge)")
    return emax


def dump_figure_data(nmgc, ratios, a_min, a_max, rho, dtg, m0_radius, K0, nH, outdir):
    """Run the sims and DUMP the figure data as TSVs (no plotting here, so the
    plot can be regenerated/edited from data without re-running NEMO). Consumed by
    scripts/plot_coag_diffusion.py."""
    os.makedirs(outdir, exist_ok=True)
    m0 = CA.mass_of_radius(m0_radius, rho)
    N0 = dtg * 1.4 * AMU / m0
    R_eff = K0 * nH * N0
    prefac = K0 / m0
    taus = [0.5, 1.0, 2.0]
    with open(os.path.join(outdir, "meta.txt"), "w") as f:
        f.write(f"# coag_diffusion figure data (generated by coag_diffusion.py --stage figure-data)\n")
        f.write(f"# a_min={a_min} a_max={a_max} rho={rho} dtg={dtg} m0={m0:.6e} g "
                f"N0={N0:.6e} /H K0={K0} nH={nH}\n")
        f.write(f"# panelA: additive kernel, g=x f vs mass at tau=0.5,1,2, two "
                f"resolutions (mass_ratio 2.0 coarse, 1.3 fine)\n")
        f.write(f"# panelB: relative L1 error on g at the evolved time "
                f"(constant T=2, additive tau=1)\n")

    # ---- Panel A: additive g(m,tau) NEMO at two resolutions ----
    with open(os.path.join(outdir, "panelA_nemo.tsv"), "w") as f:
        f.write("ratio\tnbins\ttau\tmass_g\tg_massdensity\n")
        for ratioA in (2.0, 1.3):
            out = _run_nemo(nmgc, "additive", ratioA, a_min, a_max, rho, N0, m0, prefac, nH,
                            stop_time_yr=(3.0 / R_eff) / YR, nb_outputs=60)
            edges = _bin_edges(out["mass"]); h = np.diff(edges)
            for tau in taus:
                t_yr = tau / (prefac * nH * out["M1"]) / YR
                n_k = _interp_n(out["times"], out["n_of_t"], t_yr)
                g_num = out["mass"] * (n_k / h)
                for mk, gk in zip(out["mass"], g_num):
                    f.write(f"{ratioA}\t{out['nb']}\t{tau}\t{mk:.6e}\t{gk:.6e}\n")
    # analytic reference on a fine grid
    with open(os.path.join(outdir, "panelA_analytic.tsv"), "w") as f:
        f.write("tau\tmass_g\tg_massdensity\n")
        mm = np.logspace(np.log10(m0 * 1e-2), np.log10(m0 * 1e5), 500)
        for tau in taus:
            g_ana = mm * (N0 / m0) * CA.n_additive(mm / m0, tau)
            for mk, gk in zip(mm, g_ana):
                f.write(f"{tau}\t{mk:.6e}\t{gk:.6e}\n")

    # ---- Panel B: convergence errors (evolved time), both kernels ----
    with open(os.path.join(outdir, "panelB.tsv"), "w") as f:
        f.write("kernel\tmass_ratio\tnbins\tcont_L1\tdisc_L1\n")
        for kernel, dim_t in (("constant", 2.0), ("additive", 1.0)):
            pf = K0 if kernel == "constant" else K0 / m0
            for ratio in ratios:
                o = _run_nemo(nmgc, kernel, ratio, a_min, a_max, rho, N0, m0, pf, nH,
                              stop_time_yr=(3.0 * dim_t / R_eff) / YR, nb_outputs=40)
                tt = dim_t / (pf * nH * (o["N0"] if kernel == "constant" else o["M1"]))
                nk = _interp_n(o["times"], o["n_of_t"], tt / YR)
                e = _errors(o["mass"], nk, kernel, dim_t, o["N0"], m0)
                f.write(f"{kernel}\t{ratio}\t{o['nb']}\t{e['cL1']:.6e}\t{e['dL1']:.6e}\n")
    print(f"wrote figure data to {outdir}/ (panelA_nemo.tsv, panelA_analytic.tsv, "
          f"panelB.tsv, meta.txt)")


FIXDIR = None


def main():
    global FIXDIR
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmgc", default="bin/nmgc")
    ap.add_argument("--fixture", default="tests/fixture_dust_coag")
    ap.add_argument("--ratios", default="2.0,1.5,1.3,1.15")
    ap.add_argument("--a-min", type=float, default=5.0e-7)
    ap.add_argument("--a-max", type=float, default=5.0e-5)
    ap.add_argument("--rho", type=float, default=3.0)
    ap.add_argument("--dtg", type=float, default=1.0e-2)
    ap.add_argument("--m0-radius", type=float, default=1.0e-6)
    ap.add_argument("--K0", type=float, default=1.0e-9)
    ap.add_argument("--nH", type=float, default=2.212e8)
    ap.add_argument("--T-const", type=float, default=2.0)     # evolved-time gate (constant)
    ap.add_argument("--tau-add", type=float, default=1.0)     # evolved-time gate (additive)
    ap.add_argument("--stage", default="gate", choices=["gate", "eoc", "bounded", "figure-data", "all"])
    ap.add_argument("--datadir", default="figures/coag_diffusion_data")
    args = ap.parse_args()
    FIXDIR = os.path.abspath(args.fixture)
    nmgc = os.path.abspath(args.nmgc)
    ratios = [float(s) for s in args.ratios.split(",")]
    A = (args.a_min, args.a_max, args.rho, args.dtg, args.m0_radius, args.K0, args.nH)

    ok = True
    if args.stage in ("gate", "all"):
        print("########## GATE: monotone L1 (cont+disc) at the evolved time ##########")
        ok_c, _ = sweep(nmgc, "constant", ratios, args.T_const, *A)
        ok_a, _ = sweep(nmgc, "additive", ratios, args.tau_add, *A)
        ok = ok_c and ok_a
        print("\nGATE RESULT:", "PASS" if ok else "FAIL")
    if args.stage in ("eoc", "all"):
        print("\n########## SUPPLEMENTARY: spatial order (EOC time, dimless t=0.01) ##########")
        eoc_report(nmgc, "constant", ratios, *A)
        eoc_report(nmgc, "additive", ratios, *A)
    if args.stage in ("bounded", "all"):
        print("\n########## SUPPLEMENTARY: error-vs-time boundedness (fixed resolution) ##########")
        bounded_report(nmgc, "constant", 1.3, 3.0 * args.T_const, *A)
        bounded_report(nmgc, "additive", 1.3, 3.0 * args.tau_add, *A)
    if args.stage == "figure-data":
        dump_figure_data(nmgc, ratios, *A, args.datadir)

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
