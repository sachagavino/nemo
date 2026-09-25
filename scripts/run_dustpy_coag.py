#!/usr/bin/env python3
"""
scripts/run_dustpy_coag.py -- Phase III Part C: DustPy 0D coagulation matched to NEMO.

Follows DustPy's canonical "Analytical Coagulation Kernels" recipe
(https://stammler.github.io/dustpy/test_analytical_coagulation_kernels.html):
Nr=3 radial cells (compare the middle cell [1]; the other two are boundaries),
gas evolution removed, dust radial advection off, pure sticking, and the
collision kernel overridden to the prescribed constant / additive form.

Two deliberate departures from DustPy's own test, both to match NEMO / raise rigor:
  * IC = the SAME exponential n(m,0) NEMO uses (coag_analytic.exp_ic_bins), placed
    in cell [1] -- not DustPy's monodisperse start.
  * validated against the EXACT Golovin / constant analytic (coag_analytic.py),
    not DustPy's Bessel-free large-mass asymptotic.

Convention mapping (NEMO  <->  DustPy):
  grid      m_k = m_min*ratio^(k-1)         logspace(mmin,mmax,Nm) with edges matched
  kernel    K_ij [cm3/s | cm3/s/g]          a | a*(m_i+m_j)   (a dimensionless here)
  state     n_k = nH*Y_k  [number/H]        Sigma_k = N_k*m_k^2*B  (B=2(A-1)/(A+1))
  compare   dimensionless T=a N0 t / tau=a N0 m0 t   (same on both sides)
Comparison quantity is N*m^2 = m^2 dN/dm (DustPy's convert(Sigma)=Sigma/B), so the
volume-density (NEMO) vs surface-density (DustPy) normalization drops out and only
the matched dimensionless-time shape is compared. Both codes are the k=0
Kovetz-Olund scheme, so both should show the same over-diffusion vs the exact curve.

Usage: python3 scripts/run_dustpy_coag.py --kernel constant --outdir /tmp/dp_const
"""
import argparse
import os
import sys
import tempfile

import numpy as np

import coag_analytic as CA


def matched_grid_params(a_min, a_max, mass_ratio, rho):
    """NEMO's exact geometric mass grid, and DustPy ini params that reproduce it.
    DustPy: Nm = int(decades*Nmbpd)+1, m = logspace(mmin, mmax, Nm). Choosing
    mmax = m_min*ratio^(nb-1) (NEMO's top bin) and Nmbpd in the integer window
    gives logspace ratio exactly = mass_ratio and Nm = nb."""
    radii, masses, nb = CA.build_grid(a_min, a_max, mass_ratio, rho)
    mmin, mmax = masses[0], masses[-1]                      # NEMO's first/last bin masses
    decades = int(np.ceil(np.log10(mmax / mmin)))
    Nmbpd = (nb - 0.5) / decades                            # midpoint of [(nb-1)/dec, nb/dec)
    return dict(masses=masses, nb=nb, mmin=mmin, mmax=mmax, Nmbpd=Nmbpd, decades=decades)


def run(kernel, a_min, a_max, mass_ratio, rho, m0_radius, a_kernel, N0,
        dimless_times, outdir):
    from dustpy import Simulation

    gp = matched_grid_params(a_min, a_max, mass_ratio, rho)
    m0 = CA.mass_of_radius(m0_radius, rho)

    sim = Simulation()
    sim.ini.grid.Nr = 3
    sim.ini.grid.rmin = 1.0 * 1.495978707e13               # 1 au in cm (value irrelevant, 0D)
    sim.ini.grid.rmax = 10.0 * 1.495978707e13
    sim.ini.grid.mmin = gp["mmin"]
    sim.ini.grid.mmax = gp["mmax"]
    sim.ini.grid.Nmbpd = gp["Nmbpd"]
    sim.ini.dust.allowDriftingParticles = True
    sim.initialize()

    m = np.asarray(sim.grid.m)
    if m.size != gp["nb"]:
        print(f"WARNING: DustPy Nm={m.size} != NEMO nb={gp['nb']} "
              f"(Nmbpd={gp['Nmbpd']:.4f}); grid edges may not match exactly.")
    A = np.mean(m[1:] / m[:-1])
    B = 2.0 * (A - 1.0) / (A + 1.0)                         # DustPy bin-width constant

    # ---- DustPy canonical recipe: turn everything off except sticking coagulation ----
    del sim.integrator.instructions[1]                     # drop gas integration
    sim.gas.S.tot[...] = 0.0; sim.gas.S.tot.updater = None
    sim.dust.v.rad[...] = 0.0; sim.dust.v.rad.updater = None
    sim.dust.p.frag[...] = 0.0; sim.dust.p.frag.updater = None
    sim.dust.p.stick[...] = 1.0; sim.dust.p.stick.updater = None
    if kernel == "constant":
        sim.dust.kernel[...] = a_kernel
    elif kernel == "additive":
        sim.dust.kernel[...] = a_kernel * (m[:, None] + m[None, :])[None, ...]
    else:
        raise ValueError(kernel)
    sim.dust.kernel.updater = None

    # ---- IC: same continuous exponential as NEMO, in cell [1] ----
    # DustPy convention: Sigma = B * m^2 * (dN/dm)  (so convert(Sigma)=Sigma/B = N m^2).
    dNdm = (N0 / m0) * np.exp(-m / m0)                     # exponential number density dN/dm
    sim.dust.Sigma[...] = sim.dust.SigmaFloor[...]
    sim.dust.Sigma[1, :] = np.maximum(B * m**2 * dNdm, sim.dust.SigmaFloor[1, :])
    sim.t = 1.0e-12                                         # tiny t0 (IC ~ dimensionless 0)
    sim.update()

    # dimensionless time -> physical: constant T=a N0 t ; additive tau=a N0 m0 t
    rate = a_kernel * N0 * (1.0 if kernel == "constant" else m0)
    snaps = np.array([dt / rate for dt in dimless_times])
    sim.t.snapshots = np.hstack([sim.t, snaps])

    os.makedirs(outdir, exist_ok=True)
    sim.writer.datadir = outdir
    sim.writer.overwrite = True
    sim.run()

    Sig = sim.writer.read.sequence("dust.Sigma")           # (nt, Nr, Nm)
    tt = np.asarray(sim.writer.read.sequence("t"))
    mm = np.asarray(sim.writer.read.sequence("grid.m"))[0]
    return dict(kernel=kernel, m=mm, B=B, N0=N0, m0=m0, a=a_kernel,
                t=tt, Sigma_cell1=np.asarray(Sig)[:, 1, :],
                dimless_times=dimless_times, rate=rate, nb=m.size)


def validate(res):
    """Overlay DustPy (N m^2) vs exact analytic at each snapshot; report mass
    conservation and a relative L1 error on the mass density."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    m, B, N0, m0 = res["m"], res["B"], res["N0"], res["m0"]
    kernel = res["kernel"]
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"]
    mfine = np.logspace(np.log10(m0 * 1e-2), np.log10(m0 * 1e5), 500)

    # mass conservation across the run (cell [1])
    Nm2_all = res["Sigma_cell1"] / B                        # = N * m^2
    dustmass = np.array([np.sum(res["Sigma_cell1"][i]) for i in range(len(res["t"]))])
    massdrift = (dustmass - dustmass[0]) / dustmass[0]

    print(f"\n=== DustPy-alone validation: {kernel} kernel, {res['nb']} bins ===")
    print(f"   dust-mass rel drift over run: max |{np.max(np.abs(massdrift)):.2e}|")
    print(f"   {'dimless_t':>10} {'rel-L1(g) vs exact':>20}")

    fdim = CA.n_constant if kernel == "constant" else CA.n_additive
    for c, dt in zip(colors, res["dimless_times"]):
        # nearest snapshot to this dimensionless time
        i = 1 + int(np.argmin(np.abs(res["t"][1:] * res["rate"] - dt)))
        Nm2 = res["Sigma_cell1"][i] / B                    # DustPy N m^2 at cell 1
        # exact analytic in N m^2: m^2 * (N0/m0) * f_dim(m/m0, dt)
        g_ana_grid = m**2 * (N0 / m0) * (fdim(np.asarray(m / m0), dt)
                                         if kernel == "additive" else fdim(m / m0, dt))
        e = np.sum(np.abs(Nm2 - g_ana_grid)) / np.sum(np.abs(g_ana_grid))
        print(f"   {dt:10.3f} {e:20.4e}")
        g_ana_fine = mfine**2 * (N0 / m0) * (fdim(mfine / m0, dt)
                                             if kernel == "additive" else fdim(mfine / m0, dt))
        ax.loglog(mfine, g_ana_fine, "-", color=c, lw=1.6, label=f"exact t*={dt:g}")
        ax.loglog(m, Nm2, "o", color=c, ms=4, mfc="none")
    ax.set_xlabel("grain mass m [g]"); ax.set_ylabel(r"$N\,m^2$  (= $m^2\,dN/dm$)")
    ax.set_title(f"DustPy (circles) vs exact analytic (lines): {kernel} kernel")
    ax.legend(fontsize=8, frameon=False)
    ymax = (m**2 * (N0 / m0) * (fdim(np.asarray(m / m0), res["dimless_times"][0])
            if kernel == "additive" else fdim(m / m0, res["dimless_times"][0]))).max()
    ax.set_ylim(ymax * 1e-6, ymax * 5)
    out = os.path.join(os.path.dirname(__file__) or ".", "..", "figures",
                       f"dustpy_validate_{kernel}.png")
    out = os.path.abspath(out)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=150); plt.close(fig)
    print(f"   wrote {out}")


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
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    dimless = [float(s) for s in args.times.split(",")]
    outdir = args.outdir or tempfile.mkdtemp(prefix=f"dustpy_{args.kernel}_")
    res = run(args.kernel, args.a_min, args.a_max, args.mass_ratio, args.rho,
              args.m0_radius, args.a_kernel, args.N0, dimless, outdir)
    validate(res)


if __name__ == "__main__":
    main()
