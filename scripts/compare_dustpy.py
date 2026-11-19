#!/usr/bin/env python3
"""
scripts/compare_dustpy.py -- Phase III Part C: NEMO vs DustPy vs exact analytic.

Like-for-like rerun (claude/dustpy_rerun_brief.md, 2026-11):
  * both codes start from the SAME discrete bin contents (run_dustpy_coag.point_ic),
    asserted at t = 0;
  * DustPy keeps its (1 - delta_ij/2) same-bin factor (run_dustpy_coag.run);
  * the DustPy time step is swept (sim.t.cfl);
  * NEMO writes its outputs EXACTLY at the comparison times (output_type = linear,
    start_time = 0, step = smallest comparison time), at a tight relative tolerance,
    so NEMO is a time-converged reference and no output interpolation enters;
  * a "legacy" DustPy baseline (uniform kernel override + B m^2 dN/dm IC, the
    pre-rerun harness) can be scored with the same metrics for comparability.

Metrics, all on m^2 dN/dm at the grid points normalised to unit sum
(m^2 dN/dm = m^2 n_k / h_k, h_k the geometric-mean-edge width, the same for both
codes; the per-code normalisation removes the surface- vs volume-density units):
  L1(A, B)  = sum_k |q_A - q_B|
  Tfit/T    = argmin_{T'} L1(code, exact(T')) / T, bounded on [0.3 T, 3 T], xatol 1e-6
Time mapping: physical N0 for both codes (T = K0 nH N0 t; tau = B nH N0 m0 t).

Writes figures/coag_dustpy_data/compare_<kernel>_cfl<cfl>.tsv for every cfl (and
compare_<kernel>_legacy_cfl<cfl>.tsv for the baseline), copies the smallest-cfl
pair to compare_<kernel>.tsv (read by plot_compare_dustpy.py), and writes
metrics.tsv with every number of the report.

Usage:
  python3 scripts/compare_dustpy.py                                  # ratio 2, full sweep
  python3 scripts/compare_dustpy.py --mass-ratio 1.5 --cfls 0.1 \
      --no-legacy --datadir figures/coag_dustpy_data/ratio1.5
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from scipy.optimize import minimize_scalar

import coag_analytic as CA
import run_dustpy_coag as RD

YR = 3.1556952e7
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FIXDIR = os.path.join(ROOT, "tests", "fixture_dust_coag")
NMGC = os.path.join(ROOT, "bin", "nmgc")


# --------------------------------------------------------------------------
# metrics
# --------------------------------------------------------------------------
def _norm(y):
    s = np.sum(y)
    return y / s if s > 0 else y


def exact_shape(kernel, m, m0, T):
    """Exact m^2 dN/dm at the grid points, normalised to unit sum."""
    x = m / m0
    f = CA.n_constant(x, T) if kernel == "constant" else CA.n_additive(np.asarray(x), T)
    return _norm(f * m**2)


def tfit(kernel, m, m0, q, T):
    """Best-fit T'/T of a normalised code curve q against the exact shape."""
    r = minimize_scalar(lambda s: np.abs(q - exact_shape(kernel, m, m0, s)).sum(),
                        bounds=(0.3 * T, 3.0 * T), method="bounded",
                        options=dict(xatol=1e-6))
    return r.x / T


# --------------------------------------------------------------------------
# NEMO (own runner: exact output times, configurable tolerance; Part B's
# coag_diffusion._run_nemo is left untouched)
# --------------------------------------------------------------------------
def _setk(lines, key, val):
    out, hit = [], False
    for ln in lines:
        s = ln.strip()
        if s.startswith(key + " ") or s.startswith(key + "="):
            out.append(f"{key} = {val}"); hit = True
        else:
            out.append(ln)
    if not hit:
        raise KeyError(f"{key} not in fixture parameters.in")
    return out


def _read_distribution(path):
    rows, masses, dmass = {}, {}, {}
    with open(path) as f:
        for ln in f:
            if ln.startswith("!") or not ln.strip():
                continue
            c = ln.split()
            t = float(c[0]); b = int(c[1])
            rows.setdefault(t, {})[b] = float(c[4])
            masses[b] = float(c[3]); dmass[t] = float(c[6])
    times = np.array(sorted(rows))
    nb = len(masses)
    mass = np.array([masses[b] for b in range(1, nb + 1)])
    n = np.array([[rows[t][b] for b in range(1, nb + 1)] for t in times])
    return times, mass, n, np.array([dmass[t] for t in times])


def run_nemo(kernel, ratio, a_min, a_max, rho, n_k, prefac, nH, stop_yr, nb_outputs, rtol):
    radii, masses, nb = CA.build_grid(a_min, a_max, ratio, rho)
    run = tempfile.mkdtemp(prefix=f"nemo_{kernel}_{ratio}_")
    for fn in os.listdir(FIXDIR):
        if fn.endswith(".in"):
            shutil.copy(os.path.join(FIXDIR, fn), run)
    CA.write_grid_table(os.path.join(run, "dust_grid_table.in"), radii, n_k)
    p = os.path.join(run, "parameters.in")
    txt = open(p).read().splitlines()
    for key, val in [("dust_grid_source", "tabulated"), ("dust_ic", "tabulated"),
                     ("coagulation", "1"), ("coagulation_kernel", kernel),
                     ("constant_kernel_k0", f"{prefac:.15E}"),
                     ("output_type", "linear"), ("start_time", "0.0"),
                     ("stop_time", f"{stop_yr:.15E}"), ("nb_outputs", str(nb_outputs)),
                     ("relative_tolerance", f"{rtol:.3E}")]:
        txt = _setk(txt, key, val)
    open(p, "w").write("\n".join(txt) + "\n")
    r = subprocess.run([NMGC, "run"], cwd=run, capture_output=True, text=True, timeout=3600)
    if r.returncode != 0:
        sys.stderr.write(r.stdout[-2000:] + r.stderr[-2000:])
        raise RuntimeError(f"NEMO failed: kernel={kernel} ratio={ratio}")
    out = _read_distribution(os.path.join(run, "dust_distribution.out"))
    shutil.rmtree(run, ignore_errors=True)
    return out


def nemo_curves(kernel, ratio, a_min, a_max, rho, m0, N0, K0, nH, times, rtol):
    """NEMO N m^2 at exactly the comparison times, plus its t~0 state and mass drift."""
    _, masses, _ = CA.build_grid(a_min, a_max, ratio, rho)
    n_k = RD.point_ic(masses, N0, m0)
    prefac = K0 if kernel == "constant" else K0 / m0
    amp = N0 if kernel == "constant" else N0 * m0        # physical-N0 time mapping
    t_unit_yr = 1.0 / (prefac * nH * amp) / YR           # 1 dimensionless unit, in yr
    step = min(times)
    nsteps = [T / step for T in times]
    if any(abs(s - round(s)) > 1e-9 for s in nsteps):
        raise ValueError("comparison times must be integer multiples of the smallest")
    nmax = int(round(max(nsteps)))
    t_out, mass, n_of_t, dmass = run_nemo(kernel, ratio, a_min, a_max, rho, n_k, prefac,
                                          nH, max(times) * t_unit_yr, nmax + 1, rtol)
    h = np.diff(RD.bin_edges(mass))
    curves = {}
    for T in times:
        j = int(np.argmin(np.abs(t_out - T * t_unit_yr)))
        if abs(t_out[j] / (T * t_unit_yr) - 1.0) > 1e-8:
            raise RuntimeError(f"NEMO has no output at T={T}")
        curves[T] = mass**2 * n_of_t[j] / h
    # t ~ 0 state. NEMO writes only after an integration interval, so integrate to
    # T = 1e-30. It must be that small: an exponential-tail bin k is fed by collisions of
    # far more populated bins, so its RELATIVE change is ~ T * n_1 n_{k-1} / n_k, up to
    # ~1e10 * T for the bins kept in the check (> 1e-20 of the total).
    t0_out, _, n0, _ = run_nemo(kernel, ratio, a_min, a_max, rho, n_k, prefac, nH,
                                1e-30 * t_unit_yr, 2, rtol)
    return dict(m=mass, curves=curves, n_init=n0[0] / N0, n_ic=n_k / N0,
                mass_drift=float(np.max(np.abs(dmass / dmass[0] - 1.0))),
                t0_dimless=float(t0_out[0] / t_unit_yr))


def dustpy_curves(kernel, ratio, a_min, a_max, rho, m0_radius, times, cfl, legacy):
    outdir = tempfile.mkdtemp(prefix=f"dp_{kernel}_")
    res = RD.run(kernel, a_min, a_max, ratio, rho, m0_radius, a_kernel=1.0, N0=1.0,
                 dimless_times=times, outdir=outdir, cfl=cfl, legacy=legacy)
    curves = {T: RD.nm2(res, RD.snapshot_index(res, T)) for T in times}
    M = res["Sigma_cell1"].sum(axis=1)
    out = dict(m=res["m"], curves=curves, mass_ratio_t={T: float(M[RD.snapshot_index(res, T)] / M[0])
                                                         for T in times},
               n_init=res["Sigma_cell1"][0] / res["m"], floor_n=res["floor"] / res["m"],
               gridmax=res["gridmax"])
    shutil.rmtree(outdir, ignore_errors=True)
    return out


def check_t0(nemo, dp, kernel, tol=1e-12, sig=1e-20):
    """Assert the two codes' initial n_k (in units of N0) agree. Bins below `sig` of the
    total (exponential tail; DustPy holds them at its one-particle floor, ~1e-27 here)
    are compared only in absolute terms, relative to the total."""
    a, b = nemo["n_init"], dp["n_init"]
    tot = a.sum()
    mask = a > sig * tot
    rel = float(np.max(np.abs(b[mask] / a[mask] - 1.0)))
    absd = float(np.max(np.abs(b - a)) / tot)
    vs_ic = float(np.max(np.abs(a[mask] / nemo["n_ic"][mask] - 1.0)))
    ok = rel <= tol and absd <= tol
    print(f"  t=0 [{kernel}]: per-bin max rel diff {rel:.2e} over {mask.sum()}/{a.size} bins "
          f"(> {sig:g} of total); max abs diff / total {absd:.2e}; "
          f"NEMO t~0 vs IC {vs_ic:.2e} (at T={nemo['t0_dimless']:.1e}) -> {'OK' if ok else 'FAIL'}")
    return ok, rel, absd


# --------------------------------------------------------------------------
def write_tsv(path, m, times, qN, qD, qE):
    with open(path, "w") as f:
        f.write("dimless_t\tmass_g\tnemo_norm\tdustpy_norm\texact_norm\n")
        for T in times:
            for j in range(m.size):
                f.write(f"{T}\t{m[j]:.6e}\t{qN[T][j]:.6e}\t{qD[T][j]:.6e}\t{qE[T][j]:.6e}\n")


def run(a_min, a_max, ratio, rho, m0_radius, N0, K0, nH, times, cfls, legacy_cfl,
        nemo_rtol, datadir, nemo_rtol_check=None):
    os.makedirs(datadir, exist_ok=True)
    m0 = CA.mass_of_radius(m0_radius, rho)
    rows = []          # (kernel, case, cfl, T, L1_nemo_ex, L1_dp_ex, L1_nemo_dp, tfit_N, tfit_D, M/M0)
    t0_ok = True
    for kernel in ("constant", "additive"):
        nemo = nemo_curves(kernel, ratio, a_min, a_max, rho, m0, N0, K0, nH, times, nemo_rtol)
        m = nemo["m"]
        qE = {T: exact_shape(kernel, m, m0, T) for T in times}
        qN = {T: _norm(nemo["curves"][T]) for T in times}
        tN = {T: tfit(kernel, m, m0, qN[T], T) for T in times}
        print(f"\n[{kernel}] ratio {ratio}, nb={m.size}; NEMO rtol {nemo_rtol:g}, "
              f"dust-mass drift {nemo['mass_drift']:.1e}")
        if nemo_rtol_check is not None:
            alt = nemo_curves(kernel, ratio, a_min, a_max, rho, m0, N0, K0, nH, times,
                              nemo_rtol_check)
            for T in times:
                d = np.abs(_norm(alt["curves"][T]) - qN[T]).sum()
                print(f"  NEMO rtol {nemo_rtol_check:g} vs {nemo_rtol:g}: L1 at T={T:g} = {d:.2e}")
        cases = [("fixed", c, False) for c in cfls]
        if legacy_cfl is not None:
            cases.append(("legacy", legacy_cfl, True))
        for case, cfl, legacy in cases:
            dp = dustpy_curves(kernel, ratio, a_min, a_max, rho, m0_radius, times, cfl, legacy)
            if np.max(np.abs(dp["m"] / m - 1.0)) > 1e-12:
                raise RuntimeError("grids differ")
            if not legacy:
                ok, _, _ = check_t0(nemo, dp, kernel)
                t0_ok &= ok
            qD = {T: _norm(dp["curves"][T]) for T in times}
            for T in times:
                rows.append((kernel, case, cfl, T,
                             np.abs(qN[T] - qE[T]).sum(), np.abs(qD[T] - qE[T]).sum(),
                             np.abs(qN[T] - qD[T]).sum(), tN[T],
                             tfit(kernel, m, m0, qD[T], T), dp["mass_ratio_t"][T]))
            tag = f"legacy_cfl{cfl:g}" if legacy else f"cfl{cfl:g}"
            write_tsv(os.path.join(datadir, f"compare_{kernel}_{tag}.tsv"), m, times, qN, qD, qE)
        # canonical TSV (read by plot_compare_dustpy.py): smallest cfl, fixed setup
        shutil.copy(os.path.join(datadir, f"compare_{kernel}_cfl{min(cfls):g}.tsv"),
                    os.path.join(datadir, f"compare_{kernel}.tsv"))

    hdr = ("kernel\tcase\tcfl\tT\tL1_nemo_exact\tL1_dustpy_exact\tL1_nemo_dustpy"
           "\tTfit_nemo\tTfit_dustpy\tdustpy_M_over_M0")
    with open(os.path.join(datadir, "metrics.tsv"), "w") as f:
        f.write(f"# ratio {ratio}; NEMO rtol {nemo_rtol:g}; t0 check {'OK' if t0_ok else 'FAIL'}\n")
        f.write(hdr + "\n")
        for r in rows:
            f.write("\t".join([r[0], r[1], f"{r[2]:g}", f"{r[3]:g}"]
                              + [f"{v:.6e}" for v in r[4:]]) + "\n")
    print("\n=== relative L1 on normalised m^2 dN/dm at grid points; Tfit/T; DustPy M/M0 ===")
    print(f"{'kernel':9}{'case':8}{'cfl':>6}{'T':>6}{'N-ex':>9}{'D-ex':>9}{'N-D':>10}"
          f"{'TfN':>7}{'TfD':>7}{'M/M0-1':>11}")
    for r in rows:
        print(f"{r[0]:9}{r[1]:8}{r[2]:6g}{r[3]:6g}{r[4]:9.4f}{r[5]:9.4f}{r[6]:10.2e}"
              f"{r[7]:7.3f}{r[8]:7.3f}{r[9]-1:11.1e}")
    print(f"t=0 check overall: {'OK' if t0_ok else 'FAIL'}")
    return rows, t0_ok


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
    ap.add_argument("--times", default="0.01,0.5,1.0,2.0")
    ap.add_argument("--cfls", default="0.1,0.03,0.01")
    ap.add_argument("--legacy-cfl", type=float, default=0.1)
    ap.add_argument("--no-legacy", action="store_true")
    ap.add_argument("--nemo-rtol", type=float, default=1e-8)
    ap.add_argument("--nemo-rtol-check", type=float, default=None,
                    help="also run NEMO at this rtol and report its distance to the reference")
    ap.add_argument("--datadir", default=os.path.join(ROOT, "figures", "coag_dustpy_data"))
    args = ap.parse_args()
    times = [float(s) for s in args.times.split(",")]
    cfls = [float(s) for s in args.cfls.split(",")]
    _, ok = run(args.a_min, args.a_max, args.mass_ratio, args.rho, args.m0_radius, args.N0,
                args.K0, args.nH, times, cfls, None if args.no_legacy else args.legacy_cfl,
                args.nemo_rtol, args.datadir, args.nemo_rtol_check)
    print("data dumped; plot with: python3 scripts/plot_compare_dustpy.py")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
