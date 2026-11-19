#!/usr/bin/env python3
"""
scripts/run_dustpy_coag.py -- Phase III Part C: DustPy 0D coagulation matched to NEMO.

Follows DustPy's "Analytical Coagulation Kernels" recipe
(https://stammler.github.io/dustpy/test_analytical_coagulation_kernels.html):
Nr=3 radial cells (compare the middle cell [1]; the other two are boundaries),
gas evolution removed, dust radial advection off, pure sticking, and the
collision kernel overridden to the prescribed constant / additive form.

Like-for-like setup (2026-11 rerun, claude/dustpy_rerun_brief.md):
  * KERNEL keeps DustPy's combinatorial factor (1 - delta_ij/2), as DustPy's own
    production kernel does (dustpy/std/dust.f90, subroutine kernel:
    K = (1 - 0.5*kdelta(j,i)) * ...). DustPy's sticking matrix removes two grains per
    same-bin event and assumes that factor. The uniform override used before
    (and in DustPy's analytic-kernel notebook) doubled same-bin collisions, which made
    DustPy run ~20% fast on the ratio-2 grid.
  * IC: the SAME discrete bin contents as NEMO. n_k = point_ic(m) = dN/dm|_{m_k} times
    the geometric-mean-edge width (NEMO's ic="point"); DustPy gets Sigma_k = m_k * n_k.
  * DustPy time step: sim.t.cfl (default 0.1; dt = cfl * min Sigma/|S|, first-order
    implicit Euler).

legacy=True reproduces the pre-rerun harness (uniform kernel override,
Sigma = B m^2 dN/dm IC) for the comparability baseline.

Readout (nm2): N m^2 = m^2 n_k / h_k with n_k = Sigma_k / m_k and h_k the
geometric-mean-edge width, exactly as for NEMO (legacy mode keeps the old Sigma/B).
The caller normalises to unit sum, so the surface- vs volume-density units drop out.

Usage: python3 scripts/run_dustpy_coag.py --kernel constant --cfl 0.01 --outdir /tmp/dp
"""
import argparse
import os
import tempfile

import numpy as np

import coag_analytic as CA


def bin_edges(m):
    """Geometric-mean bin edges (identical to coag_diffusion._bin_edges)."""
    nb = m.size
    e = np.empty(nb + 1)
    e[1:-1] = np.sqrt(m[:-1] * m[1:])
    e[0] = m[0] / np.sqrt(m[1] / m[0])
    e[-1] = m[-1] * np.sqrt(m[-1] / m[-2])
    return e


def point_ic(m, N0, m0):
    """Shared discrete IC: exponential number density sampled at m_k, times the
    geometric-mean-edge width (NEMO's ic='point'). Number per bin."""
    return (N0 / m0) * np.exp(-m / m0) * np.diff(bin_edges(m))


def matched_grid_params(a_min, a_max, mass_ratio, rho):
    """NEMO's exact geometric mass grid, and DustPy ini params that reproduce it.
    DustPy: Nm = int(ceil(decades)*Nmbpd)+1, m = logspace(mmin, mmax, Nm). Choosing
    mmax = m_min*ratio^(nb-1) (NEMO's top bin) and Nmbpd in the integer window
    gives logspace ratio exactly = mass_ratio and Nm = nb."""
    radii, masses, nb = CA.build_grid(a_min, a_max, mass_ratio, rho)
    mmin, mmax = masses[0], masses[-1]
    decades = int(np.ceil(np.log10(mmax / mmin)))
    Nmbpd = (nb - 0.5) / decades
    return dict(masses=masses, nb=nb, mmin=mmin, mmax=mmax, Nmbpd=Nmbpd, decades=decades)


def run(kernel, a_min, a_max, mass_ratio, rho, m0_radius, a_kernel, N0,
        dimless_times, outdir, cfl=0.1, legacy=False):
    from dustpy import Simulation

    gp = matched_grid_params(a_min, a_max, mass_ratio, rho)
    m0 = CA.mass_of_radius(m0_radius, rho)

    sim = Simulation()
    sim.ini.grid.Nr = 3
    sim.ini.grid.rmin = 1.0 * 1.495978707e13               # 1 au (value irrelevant, 0D)
    sim.ini.grid.rmax = 10.0 * 1.495978707e13
    sim.ini.grid.mmin = gp["mmin"]
    sim.ini.grid.mmax = gp["mmax"]
    sim.ini.grid.Nmbpd = gp["Nmbpd"]
    sim.ini.dust.allowDriftingParticles = True
    sim.initialize()

    m = np.asarray(sim.grid.m)
    if m.size != gp["nb"]:
        raise RuntimeError(f"DustPy Nm={m.size} != NEMO nb={gp['nb']}")
    gridmax = float(np.max(np.abs(m / gp["masses"] - 1.0)))
    if gridmax > 1e-12:
        raise RuntimeError(f"DustPy/NEMO mass grids differ: max|mD/mN-1|={gridmax:.2e}")
    Nm = m.size
    A = np.mean(m[1:] / m[:-1])
    B = 2.0 * (A - 1.0) / (A + 1.0)                         # DustPy bin-width constant

    # ---- turn everything off except sticking coagulation ----
    del sim.integrator.instructions[1]                     # drop gas integration
    sim.gas.S.tot[...] = 0.0; sim.gas.S.tot.updater = None
    sim.dust.v.rad[...] = 0.0; sim.dust.v.rad.updater = None
    sim.dust.p.frag[...] = 0.0; sim.dust.p.frag.updater = None
    sim.dust.p.stick[...] = 1.0; sim.dust.p.stick.updater = None
    diag = np.ones((Nm, Nm)) if legacy else (1.0 - 0.5 * np.eye(Nm))
    if kernel == "constant":
        sim.dust.kernel[...] = a_kernel * diag[None, ...]
    elif kernel == "additive":
        sim.dust.kernel[...] = a_kernel * ((m[:, None] + m[None, :]) * diag)[None, ...]
    else:
        raise ValueError(kernel)
    sim.dust.kernel.updater = None

    # ---- IC in cell [1] ----
    if legacy:
        dNdm = (N0 / m0) * np.exp(-m / m0)
        Sig0 = B * m**2 * dNdm
    else:
        Sig0 = m * point_ic(gp["masses"], N0, m0)
    sim.dust.Sigma[...] = sim.dust.SigmaFloor[...]
    sim.dust.Sigma[1, :] = np.maximum(Sig0, sim.dust.SigmaFloor[1, :])
    sim.t = 1.0e-12                                         # tiny t0 (dimensionless ~0)
    sim.t.cfl = cfl
    sim.update()

    rate = a_kernel * N0 * (1.0 if kernel == "constant" else m0)
    snaps = np.array([dt / rate for dt in dimless_times])
    sim.t.snapshots = np.hstack([sim.t, snaps])

    os.makedirs(outdir, exist_ok=True)
    sim.writer.datadir = outdir
    sim.writer.overwrite = True
    sim.verbosity = 0
    sim.run()

    Sig = np.asarray(sim.writer.read.sequence("dust.Sigma"))[:, 1, :]   # (nt, Nm)
    tt = np.asarray(sim.writer.read.sequence("t"))
    floor = np.asarray(sim.dust.SigmaFloor)[1, :]
    return dict(kernel=kernel, m=m, B=B, N0=N0, m0=m0, a=a_kernel, cfl=cfl,
                legacy=legacy, t=tt, Sigma_cell1=Sig, floor=floor, Sig0=Sig0,
                dimless_times=list(dimless_times), rate=rate, nb=Nm, gridmax=gridmax)


def nm2(res, i):
    """DustPy N m^2 (= m^2 dN/dm, un-normalised) at snapshot i, cell [1]."""
    m = res["m"]
    if res["legacy"]:
        return res["Sigma_cell1"][i] / res["B"]             # old convention
    return m**2 * (res["Sigma_cell1"][i] / m) / np.diff(bin_edges(m))


def snapshot_index(res, T):
    """Snapshot at dimensionless time T; asserts the snapshot landed on it."""
    i = 1 + int(np.argmin(np.abs(res["t"][1:] * res["rate"] - T)))
    rel = abs(res["t"][i] * res["rate"] - T) / T
    if rel > 1e-9:
        raise RuntimeError(f"DustPy snapshot at {res['t'][i]*res['rate']} != T={T}")
    return i


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kernel", default="constant", choices=["constant", "additive"])
    ap.add_argument("--a-min", type=float, default=5.0e-7)
    ap.add_argument("--a-max", type=float, default=5.0e-5)
    ap.add_argument("--mass-ratio", type=float, default=2.0)
    ap.add_argument("--rho", type=float, default=3.0)
    ap.add_argument("--m0-radius", type=float, default=1.0e-6)
    ap.add_argument("--a-kernel", type=float, default=1.0)
    ap.add_argument("--N0", type=float, default=1.0)
    ap.add_argument("--times", default="0.5,1.0,2.0")
    ap.add_argument("--cfl", type=float, default=0.1)
    ap.add_argument("--legacy", action="store_true",
                    help="pre-rerun harness: uniform kernel override + B m^2 dN/dm IC")
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    dimless = [float(s) for s in args.times.split(",")]
    outdir = args.outdir or tempfile.mkdtemp(prefix=f"dustpy_{args.kernel}_")
    res = run(args.kernel, args.a_min, args.a_max, args.mass_ratio, args.rho,
              args.m0_radius, args.a_kernel, args.N0, dimless, outdir,
              cfl=args.cfl, legacy=args.legacy)
    M = res["Sigma_cell1"].sum(axis=1)
    print(f"{args.kernel}: nb={res['nb']} cfl={args.cfl} legacy={args.legacy} "
          f"M/M0 at snapshots: {np.array2string(M / M[0], precision=6)}")


if __name__ == "__main__":
    main()
