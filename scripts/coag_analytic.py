#!/usr/bin/env python3
"""
scripts/coag_analytic.py -- Phase III Part B (numerical-diffusion validation).

Closed-form Smoluchowski solutions for the CONSTANT and ADDITIVE (Golovin)
kernels with an exponential initial condition n(x,0)=exp(-x), x=m/m0, together
with the exponential-IC generator for NEMO's tabulated grid+IC path. One
definition of (N0, m0) drives BOTH the IC written for NEMO and the analytic
reference, so they agree by construction.

Analytic forms (mathematical facts; the additive form is transcribed from
Lombart & Laibe 2021, arXiv:2011.12298, Eq. 7, which restates Golovin 1963 /
Scott 1968; the constant form is Smoluchowski 1916 / Scott 1968, with its time
normalisation tied to NEMO's Rung-1 decay gate coag_gate.sh). In L&L's Eq. 7 the
symbol "T" denotes 1-exp(-tau); we call that g here to avoid clashing with the
constant-kernel dimensionless time T.

  CONSTANT kernel  K_ij = K0
    dimensionless time  T = K0 * nH * N0 * t
    n(x,T) = (1+T/2)^-2 * exp( -x / (1+T/2) )
    zeroth moment  ∫ n dx = 1/(1+T/2)   (= N(t)/N0; ties to coag_gate.sh)
    first  moment  ∫ x n dx = 1

  ADDITIVE (Golovin) kernel  K_ij = B (m_i+m_j)
    dimensionless time  tau = B * M1_phys * t = B * nH * N0 * m0 * t
    g = 1 - exp(-tau)
    n(x,tau) = (1-g)/(x*sqrt(g)) * I1(2 x sqrt(g)) * exp( -(1+g) x )
    zeroth moment  ∫ n dx = exp(-tau)   (= N(t)/N0)
    first  moment  ∫ x n dx = 1

n is returned in units of N0/m0 (set N0=m0=1); multiply by N0/m0 with x=m/m0 for
physical values. The two self-consistency gates below (first moment constant;
zeroth-moment law) are source-independent: `python3 scripts/coag_analytic.py`
runs them.

Physical <-> dimensionless time (NEMO units: K0/B = constant_kernel_k0,
nH = initial_gas_density [cm^-3], t [s]; N0 [per-H], m0 [g]):
  constant:  T   = K0 * nH * N0 * t
  additive:  tau = B  * nH * N0 * m0 * t     (B has units cm^3 s^-1 g^-1)
The additive map is DERIVED from dN/dt = -B M1 N (standard Golovin number law),
and gate (ii) confirms the analytic zeroth moment is exp(-tau), pinning it.
"""
import argparse
import numpy as np
from scipy.special import ive          # exponentially-scaled I1: ive(1,z)=I1(z)exp(-|z|)
from scipy.integrate import quad

PI = np.pi


# --------------------------------------------------------------------------
# grain grid (replicates dust_grid.f90 exactly)
# --------------------------------------------------------------------------
def mass_of_radius(a, rho):
    return (4.0 / 3.0) * PI * rho * a**3


def radius_of_mass(m, rho):
    return (3.0 * m / (4.0 * PI * rho))**(1.0 / 3.0)


def build_grid(a_min, a_max, mass_ratio, rho):
    """Return (radii, masses, nb) for NEMO's geometric mass grid."""
    m_min = mass_of_radius(a_min, rho)
    m_max = mass_of_radius(a_max, rho)
    nb = int(np.log(m_max / m_min) / np.log(mass_ratio)) + 1
    masses = m_min * mass_ratio**np.arange(nb)
    radii = radius_of_mass(masses, rho)
    return radii, masses, nb


def exp_ic_bins(masses, N0, m0):
    """Exponential IC as number per bin, integrated over the bin CENTERED on each
    grid point m_k (geometric-mean edges [sqrt(m_{k-1}m_k), sqrt(m_k m_{k+1})]).
    This is consistent with representing n_k as grains AT m_k and reconstructing a
    density n_k/Delta m_k on the same edges -- so the IC and the density agree, and
    NEMO matches DustPy / the analytic (both point-evaluate the density at m_k) at
    t0. Returns n_k [number per H per bin]."""
    nb = masses.size
    edges = np.empty(nb + 1)
    edges[1:-1] = np.sqrt(masses[:-1] * masses[1:])
    edges[0] = masses[0] / np.sqrt(masses[1] / masses[0])
    edges[-1] = masses[-1] * np.sqrt(masses[-1] / masses[-2])
    lo = np.exp(-edges / m0)
    return N0 * (lo[:-1] - lo[1:])              # integral over [e_k, e_{k+1}]


# --------------------------------------------------------------------------
# closed-form solutions (dimensionless: x = m/m0, n in units N0/m0)
# --------------------------------------------------------------------------
def n_constant(x, T):
    a = 1.0 + 0.5 * T
    return a**(-2) * np.exp(-x / a)


def n_additive(x, tau):
    x = np.asarray(x, dtype=float)
    g = 1.0 - np.exp(-tau)
    if g <= 0.0:
        return np.exp(-x)                     # tau -> 0 limit: IC
    sg = np.sqrt(g)
    z = 2.0 * x * sg
    # I1(z) = ive(1,z)*exp(z); combine the exponentials for stability:
    #   z - (1+g) x = -x (1 - sqrt(g))^2
    out = np.where(
        x > 0.0,
        (1.0 - g) / (np.where(x > 0.0, x, 1.0) * sg) * ive(1, z)
        * np.exp(-x * (1.0 - sg)**2),
        1.0 - g,                              # x -> 0 limit: n = 1-g = exp(-tau)
    )
    return out


# --------------------------------------------------------------------------
# dimensionless-time maps (physical NEMO time -> T / tau)
# --------------------------------------------------------------------------
def T_of_t(t, K0, nH, N0):
    return K0 * nH * N0 * t


def tau_of_t(t, B, nH, N0, m0):
    return B * nH * N0 * m0 * t


# --------------------------------------------------------------------------
# NEMO tabulated grid+IC table writer (radius, 1/n_k, T_d, T_CR,peak)
# --------------------------------------------------------------------------
def write_grid_table(path, radii, n_k, T_d=10.0, T_CR=15.0):
    """dust_grid_source=tabulated + dust_ic=tabulated read: 4 columns per bin.
    GTODN = 1/n_k is the abundance column NEMO turns into abundances(GRAIN_k).
    T_d/T_CR are irrelevant for constant/additive pure-coag (kernels use gas T /
    are athermal), so a constant placeholder is fine (documented in header)."""
    with open(path, "w") as f:
        f.write("! NEMO tabulated dust grid + IC (generated by coag_analytic.py).\n")
        f.write("! columns: radius[cm]  1/n_k[per-H^-1]  T_d[K]  T_CR_peak[K]\n")
        f.write("! T_d/T_CR are placeholders: constant/additive kernels do not use them.\n")
        for a, nk in zip(radii, n_k):
            gtodn = 1.0e300 if nk <= 1.0e-290 else 1.0 / nk   # empty/underflowed bin (cf. dust_ic_monodisperse)
            f.write(f"{a:.15E}  {gtodn:.15E}  {T_d:.6E}  {T_CR:.6E}\n")


# --------------------------------------------------------------------------
# self-consistency gates (source-independent; the real correctness check)
# --------------------------------------------------------------------------
def _xmax_constant(T):
    return 60.0 * (1.0 + 0.5 * T)


def _xmax_additive(tau):
    g = 1.0 - np.exp(-tau)
    sg = np.sqrt(g)
    L = 1.0 / max((1.0 - sg)**2, 1e-30)       # mass-decay length of the Golovin tail
    return 60.0 * max(L, 5.0)


def _moment(fn, order, arg, xmax):
    if fn is n_additive:
        integrand = lambda x: x**order * n_additive(np.array([x]), arg)[0]
    else:
        integrand = lambda x: x**order * n_constant(x, arg)
    return quad(integrand, 0.0, xmax, limit=800)[0]


def run_gates(verbose=True):
    ok = True
    # test over the regime where the distribution stays on a finite mass grid,
    # which is the only regime the NEMO comparison uses (mean mass < m_max).
    Ts = [0.0, 0.5, 2.0, 10.0, 20.0]
    taus = [0.1, 0.5, 1.0, 2.0]
    rtol_mass, rtol_num = 1e-6, 1e-6

    if verbose:
        print("=== coag_analytic self-consistency gates ===")
        print(" CONSTANT kernel   n(x,T)=(1+T/2)^-2 exp(-x/(1+T/2))")
    for T in Ts:
        xmax = _xmax_constant(T)
        M1 = _moment(n_constant, 1, T, xmax)
        M0 = _moment(n_constant, 0, T, xmax)
        num_expected = 1.0 / (1.0 + 0.5 * T)
        e_mass = abs(M1 - 1.0)
        e_num = abs(M0 - num_expected) / num_expected
        ok &= (e_mass < rtol_mass) and (e_num < rtol_num)
        if verbose:
            print(f"   T={T:6.2f}  mass=∫xn={M1:.10f} (err {e_mass:.1e})   "
                  f"num=∫n={M0:.10f} vs 1/(1+T/2)={num_expected:.10f} (err {e_num:.1e})")

    if verbose:
        print(" ADDITIVE kernel   n(x,tau)=(1-g)/(x√g) I1(2x√g) exp(-(1+g)x),  g=1-e^-tau")
    for tau in taus:
        xmax = _xmax_additive(tau)
        M1 = _moment(n_additive, 1, tau, xmax)
        M0 = _moment(n_additive, 0, tau, xmax)
        num_expected = np.exp(-tau)
        e_mass = abs(M1 - 1.0)
        e_num = abs(M0 - num_expected) / num_expected
        ok &= (e_mass < rtol_mass) and (e_num < rtol_num)
        if verbose:
            print(f"   tau={tau:5.2f}  mass=∫xn={M1:.10f} (err {e_mass:.1e})   "
                  f"num=∫n={M0:.10f} vs e^-tau={num_expected:.10f} (err {e_num:.1e})")

    # IC limit check: n(x,0)=e^-x for both kernels
    xs = np.array([0.1, 0.5, 1.0, 3.0])
    ic_c = np.max(np.abs(n_constant(xs, 0.0) - np.exp(-xs)))
    ic_a = np.max(np.abs(n_additive(xs, 0.0) - np.exp(-xs)))
    ok &= (ic_c < 1e-12) and (ic_a < 1e-12)
    if verbose:
        print(f" IC limit (t=0 -> e^-x): constant err {ic_c:.1e}, additive err {ic_a:.1e}")
        print("RESULT:", "PASS" if ok else "FAIL")
    return ok


# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-ic", metavar="PATH",
                    help="write a NEMO tabulated grid+IC table to PATH")
    ap.add_argument("--a-min", type=float, default=1.0e-6)
    ap.add_argument("--a-max", type=float, default=3.2003e-5)
    ap.add_argument("--mass-ratio", type=float, default=2.0)
    ap.add_argument("--rho", type=float, default=3.0)
    ap.add_argument("--N0", type=float, default=None,
                    help="total grain number per H (default: set by --dtg)")
    ap.add_argument("--dtg", type=float, default=1.0e-2,
                    help="target dust-to-gas: sets N0 so sum(m_k n_k) matches")
    ap.add_argument("--m0-bin", type=int, default=4,
                    help="mean initial mass m0 = mass of this (1-based) bin")
    args = ap.parse_args()

    if not args.write_ic:
        run_gates()
        return

    radii, masses, nb = build_grid(args.a_min, args.a_max, args.mass_ratio, args.rho)
    m0 = masses[args.m0_bin - 1]
    if args.N0 is not None:
        N0 = args.N0
    else:
        # choose N0 so total dust mass per H = dtg * (mean gas mass per H ~ 1.4 amu)
        AMU = 1.66053892e-24
        target_M1 = args.dtg * 1.4 * AMU
        # for the binned exponential, sum(m_k n_k) ~ N0 * m0 * O(1); solve numerically
        n_unit = exp_ic_bins(masses, 1.0, m0)
        N0 = target_M1 / np.sum(masses * n_unit)
    n_k = exp_ic_bins(masses, N0, m0)
    write_grid_table(args.write_ic, radii, n_k)
    M1 = np.sum(masses * n_k)
    print(f"wrote {args.write_ic}: nb={nb} bins, m0={m0:.4e} g (bin {args.m0_bin}), "
          f"N0={N0:.4e} /H, sum(m_k n_k)={M1:.4e} g/H, sum(n_k)={np.sum(n_k):.4e} /H")


if __name__ == "__main__":
    main()
