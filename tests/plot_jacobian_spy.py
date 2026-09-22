#!/usr/bin/env python3
# ===========================================================================
# scripts/plot_jacobian_spy.py -- Phase III (docs/PhaseIII_jacobian_figure_brief.md, sec 3)
#
# Pure post-processing of the provenance dump (tests/fixture_jac_figure/
# jac_pattern.tsv + jac_species.tsv). No code units, no solver. Produces two
# block-structured spy figures of the coupled Jacobian, colored by provenance.
#
#   primary  : full pattern, one marker per nonzero, colored by mechanism with
#              novelty-first PRECEDENCE for multi-bit cells (default).
#   secondary: coupling zoom -- grain-abundance coupling block + coag/ice blocks
#              + e- row, with the gtodn_recip subset drawn as an overlay and the
#              couplings-ride-on-other-physics cells (ice_transport x grain_coupling,
#              surface_net x grain_coupling) marked, charge_electron highlighted.
#
# Reorder (resolved layout): [ gas | grain(bin 1..N, 0 then -) | ice(bin 1..N, per base) ].
# No separate surface band: all surface species in reduced_CHO are ice (J);
# surface_net is a provenance COLOUR, not a species band.
#
# First-look: writes PNGs by default (eyeball legibility before the PDF/styling pass).
# ===========================================================================
import argparse, os, sys
import numpy as np
import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# bit index -> (name); precedence is defined separately (novelty-first)
BIT = dict(diagonal=0, gas_chem=1, surface_net=2, coag_grain=3,
           ice_transport=4, charge_electron=5, grain_coupling=6, gtodn_recip=7)

# novelty-first precedence (highest first) -- must match the caption text below
PRECEDENCE = [
    (r'charge (e$^-$)', 5, "#222222"),  # near-black — rare, must be findable, no clash
    ("reciprocal term",     7, "#882255"),  # wine  — value subset of grain_coupling
    ("grain-abundance coupling",  6, "#CC6677"),  # rose  — same family as gtodn (light vs dark)
    ("ice transport",   4, "#44AA99"),  # teal  — big block, recedes
    ("grain coagulation",      3, "#DDCC77"),  # sand
    ("surface reactions",     2, "#117733"),  # green
    ("gas-phase reactions",   1, "#332288"),  # indigo
    ("diagonal",        0, "#DDDDDD"),  # light grey
]
COLOR = {n: c for (n, b, c) in PRECEDENCE}


def load_species(path):
    rows = []
    with open(path) as f:
        next(f)
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) < 5:
                p += [""] * (5 - len(p))
            rows.append(dict(index=int(p[0]), name=p[1], cls=p[2],
                             bin=int(p[3]) if p[3] else 0, charge=p[4]))
    return rows


def load_pattern(path):
    r, c, m = [], [], []
    with open(path) as f:
        next(f)
        for ln in f:
            a, b, d = ln.split("\t")
            r.append(int(a)); c.append(int(b)); m.append(int(d))
    return np.array(r), np.array(c), np.array(m)


def build_order(species):
    gas = sorted([s for s in species if s["cls"] == "gas"], key=lambda s: s["index"])
    grain = sorted([s for s in species if s["cls"] == "grain"],
                   key=lambda s: (s["bin"], 0 if s["charge"] == "0" else 1))
    ice = sorted([s for s in species if s["cls"] == "ice"],
                 key=lambda s: (s["bin"], s["name"][3:]))
    order = gas + grain + ice
    pos = {s["index"]: i for i, s in enumerate(order)}   # species index -> row/col position
    bounds = dict(n_gas=len(gas), n_grain=len(grain), n_ice=len(ice))
    return order, pos, bounds


def precedence_family(mask):
    for name, bit, _ in PRECEDENCE:
        if mask & (1 << bit):
            return name
    return None


def draw_bands(ax, bounds, N, label_axes=True):
    b1 = bounds["n_gas"] - 0.5
    b2 = bounds["n_gas"] + bounds["n_grain"] - 0.5
    for b in (b1, b2):
        ax.axvline(b, color="0.", lw=1., ls="-", zorder=20)
        ax.axhline(b, color="0.", lw=1., ls="-", zorder=20)
    b_gr = bounds["n_gas"] + 4 - 0.5
    ax.axvline(b_gr, color="0.2", lw=0.8, ls="--", zorder=20)
    ax.axhline(b_gr, color="0.2", lw=0.8, ls="--", zorder=20)
    # b1_ice = bounds["n_gas"] + bounds["n_grain"] + 11 - 0.5
    # b2_ice = bounds["n_gas"] + bounds["n_grain"] + 22 - 0.5
    # b3_ice = bounds["n_gas"] + bounds["n_grain"] + 33 - 0.5
    # for b in (b1_ice, b2_ice, b3_ice):
    #     ax.axvline(b, color="0.4", lw=0.7, ls=":", zorder=5)
    #     ax.axhline(b, color="0.4", lw=0.7, ls=":", zorder=5)
    if not label_axes:
        return
    centers = [bounds["n_gas"] / 2,
               bounds["n_gas"] + bounds["n_grain"] / 2,
               bounds["n_gas"] + bounds["n_grain"] + bounds["n_ice"] / 2]
    names = ["gas", "grain (0,-)", "ice (bin, base)"]          # (1) ice on one line
    for cpos, nm in zip(centers, names):
        ax.text(cpos, -0.02 * N, nm, ha="center", va="bottom",
                fontsize=9, fontweight='bold', color="0.25", clip_on=False)                 # columns: top
        ax.text(N + 1.2, cpos, nm, ha="center", va="center",
                fontsize=9, fontweight='bold', color="0.25", rotation=90, clip_on=False)    # (2) rows: right, vertical


def draw_bands_zoom(ax, bounds, N, label_axes=True):
    b1 = bounds["n_gas"] - 0.5
    b2 = bounds["n_gas"] + bounds["n_grain"] - 0.5
    ax.axvline(b1, color="0.", lw=1., ls="-", zorder=20)
    ax.axhline(b1, color="0.", lw=1., ls="-", zorder=20)
    ax.axvline(b2, color="0.", lw=1., ls="-", zorder=20)
    ax.axhline(b2, color="0.", lw=1., ls="-", zorder=20)
    ax.axhline(b2, color="0.", lw=1., ls="-", zorder=20)
    b_gr = bounds["n_gas"] + 4 - 0.5
    ax.axvline(b_gr, color="0.2", lw=0.8, ls="--", zorder=20)
    ax.axhline(b_gr, color="0.2", lw=0.8, ls="--", zorder=20)
    # b1_ice = bounds["n_gas"] + bounds["n_grain"] + 11 - 0.5
    # b2_ice = bounds["n_gas"] + bounds["n_grain"] + 22 - 0.5
    # b3_ice = bounds["n_gas"] + bounds["n_grain"] + 33 - 0.5
    # for b in (b1_ice, b2_ice, b3_ice):
    #     ax.axvline(b, color="0.4", lw=0.7, ls=":", zorder=5)
    #     ax.axhline(b, color="0.4", lw=0.7, ls=":", zorder=5)
    if not label_axes:
        return
    centers = [bounds["n_gas"] / 2,
               bounds["n_gas"] + bounds["n_grain"] / 2,
               bounds["n_gas"] + bounds["n_grain"] + bounds["n_ice"] / 2]
    centers_col = [bounds["n_gas"] + bounds["n_grain"] / 2,
               bounds["n_gas"] + bounds["n_grain"] + bounds["n_ice"] / 2]
    names = ["gas", "grain (0,-)", "ice (bin, base)"]          # (1) ice on one line
    names_col = ["grain (0,-)", "ice (bin, base)"]          # (1) ice on one line
    for cpos, nm in zip(centers, names):                # columns: top
        ax.text(N + 1.2, cpos, nm, ha="center", va="center",
                fontsize=7, fontweight='bold', color="0.25", rotation=90, clip_on=False)    # (2) rows: right, vertical
    for cpos, nm in zip(centers_col, names_col):  
        ax.text(cpos, -0.02 * N, nm, ha="center", va="bottom",
                fontsize=7, fontweight='bold', color="0.25", clip_on=False) 


def plot_primary(rr, cc, mm, pos, bounds, order, overlap, outpath):
    N = len(order)
    x = np.array([pos[c] for c in cc]); y = np.array([pos[r] for r in rr])
    fig, ax = plt.subplots(figsize=(8.4/1.1, 8.0/1.1))
    ms = 36

    if overlap == "precedence":
        fam = np.array([precedence_family(m) for m in mm])
        for name, bit, col in PRECEDENCE:
            sel = fam == name
            if sel.any():
                ax.scatter(x[sel], y[sel], s=ms, marker="o", c=col,
                           linewidths=0, label=name, zorder=3)
    elif overlap == "mixed":
        pc = np.array([bin(m).count("1") for m in mm])
        single = pc == 1
        for name, bit, col in PRECEDENCE:
            sel = single & ((mm & (1 << bit)) != 0)
            if sel.any():
                ax.scatter(x[sel], y[sel], s=ms, marker="o", c=col, linewidths=0,
                           label=name, zorder=3)
        sel = pc >= 2
        if sel.any():
            ax.scatter(x[sel], y[sel], s=ms, marker="o", c="k", linewidths=0,
                       label="mixed (>=2)", zorder=4)
    else:  # overlay -- translucent overplot, one pass per family
        for name, bit, col in PRECEDENCE:
            sel = (mm & (1 << bit)) != 0
            if sel.any():
                ax.scatter(x[sel], y[sel], s=ms, marker="o", c=col, linewidths=0,
                           alpha=0.45, label=name, zorder=3)
    #ax.set_title("Coupled Jacobian sparsity by provenance (4 bins, reduced network)")
    draw_bands(ax, bounds, N)
    ax.set_xlim(-0.5, N - 0.5); ax.set_ylim(-0.5, N - 0.5)
    ax.invert_yaxis(); ax.set_aspect("equal")
    ax.set_xlabel("Perturbed species (column)", fontsize=15); ax.set_ylabel("Affected species (row)", fontsize=15)
    # ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8,
    #           markerscale=2.0, frameon=False, title="provenance")
    ax.legend(loc="best", fontsize=8, markerscale=1.5, frameon=True,   # (3) inside + framed
              facecolor="#D9D9D9",  framealpha=0.8, edgecolor="0.7")
    # cap = ("Overlap rule: novelty-first precedence  "
    #        "charge_electron > gtodn_recip > grain_coupling > ice_transport > "
    #        "coag_grain > surface_net > gas_chem > diagonal.\n"
    #        "State vector reordered by species class for legibility; production ordering differs. "
    #        "Structure scales to the 20-bin fiducial.")
    #fig.text(0.01, 0.005, cap, fontsize=7, color="0.3", va="bottom")
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(outpath, dpi=15, bbox_inches="tight")
    #plt.close(fig)
    #plt.show()
    plt.close(fig)


def plot_secondary(rr, cc, mm, pos, bounds, order, outpath):
    N = len(order)
    n_gas = bounds["n_gas"]
    # columns of interest: grain + ice (drop gas columns) -> the coupling / coag / ice blocks
    col_lo = n_gas - 0.5
    col_hi = N - 0.5
    xcol = np.array([pos[c] for c in cc]); yrow = np.array([pos[r] for r in rr])
    keep = xcol >= n_gas                       # grain + ice columns
    xk, yk, mk = xcol[keep], yrow[keep], mm[keep]

    fig, ax = plt.subplots(figsize=(4.8/1.1, 8.0/1.1))
    ms = 36

    b6 = (mk & (1 << BIT["grain_coupling"])) != 0
    b7 = (mk & (1 << BIT["gtodn_recip"])) != 0
    b5 = (mk & (1 << BIT["charge_electron"])) != 0
    b4 = (mk & (1 << BIT["ice_transport"])) != 0
    b3 = (mk & (1 << BIT["coag_grain"])) != 0
    b2 = (mk & (1 << BIT["surface_net"])) != 0

    # Marker sizes
    s_bg = 18
    s_main = 18
    s_overlay = 16
    s_recip = 14
    s_charge = 14

    # ------------------------------------------------------------------
    # Background: coagulation / ice transport
    # ------------------------------------------------------------------
    bg = (b3 | b4) & ~b6

    ax.scatter(
        xk[bg], yk[bg],
        s=s_bg,
        marker="o",
        c="0.80",
        linewidths=0,
        zorder=1,
        label="grain coagulation / ice transport",
    )


    # ------------------------------------------------------------------
    # Grain-abundance coupling
    # ------------------------------------------------------------------
    ax.scatter(
        xk[b6], yk[b6],
        s=s_main,
        marker="o",
        c=COLOR["grain-abundance coupling"],
        linewidths=0,
        zorder=2,
        label="grain-abundance coupling",
    )


    # ------------------------------------------------------------------
    # Grain coupling + ice transport
    # ------------------------------------------------------------------
    ov1 = b6 & b4

    ax.scatter(
        xk[ov1], yk[ov1],
        s=s_overlay,
        marker="o",
        facecolors="none",
        edgecolors=COLOR["ice transport"],
        linewidths=1.0,
        zorder=3,
        label="grain coupling × ice transport",
    )


    # ------------------------------------------------------------------
    # Grain coupling + surface reactions
    # ------------------------------------------------------------------
    ov2 = b6 & b2

    ax.scatter(
        xk[ov2], yk[ov2],
        s=s_overlay,
        marker="o",
        facecolors="none",
        edgecolors=COLOR["surface reactions"],
        linewidths=1.0,
        zorder=3,
        label="grain coupling × surface reactions",
    )


    # ------------------------------------------------------------------
    # Reciprocal term
    # ------------------------------------------------------------------
    ax.scatter(
        xk[b7], yk[b7],
        s=s_recip,
        marker="+",
        c=COLOR["reciprocal term"],
        linewidths=0.9,
        zorder=4,
        label="reciprocal term",
    )


    # ------------------------------------------------------------------
    # Charge / electron
    # ------------------------------------------------------------------
    ax.scatter(
        xk[b5], yk[b5],
        s=s_charge,
        marker="D",
        facecolors=COLOR[r'charge (e$^-$)'],
        edgecolors="none",
        zorder=5,
        label=r"charge ($e^-$)",
    )

    draw_bands_zoom(ax, bounds, N)
    # band line between grain and ice columns
    ax.axvline(n_gas + bounds["n_grain"] - 0.5, color="0.0", lw=1, ls="-", zorder=10)
    ax.set_xlim(col_lo, col_hi); ax.set_ylim(-0.5, N - 0.5)
    ax.invert_yaxis(); ax.set_aspect("equal")
    ax.tick_params(labelsize=9)
    ax.set_xlabel("grain (0,-) | ice", fontsize=12)
    ax.set_ylabel("all species (row)", fontsize=12)
    #ax.set_title("Coupling zoom: grain-abundance couplings, gtodn overlay, charge entries")
    #ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8, frameon=False)
    ax.legend(loc="best", fontsize=8, markerscale=1.5, frameon=True,   # (3) inside + framed
              facecolor="#F7F7F7",  framealpha=0.8, edgecolor="0.7")
    # cap = ("gtodn_recip is the value subset of grain_coupling (co-located 160/160), drawn as an "
    #        "overlay. Outlined cells show grain_coupling entries that other physics also populates "
    #        "(ice_transport, surface_net) -- the companion to the analytic-Jacobian verification.")
    # fig.text(0.01, 0.005, cap, fontsize=7, color="0.3", va="bottom")
    #fig.tight_layout(rect=[0, 0.05, 1, 1])
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(outpath, dpi=15, bbox_inches="tight")
    #plt.close(fig)
    #plt.show()
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", default="fixture_jac_figure/jac_pattern.tsv")
    ap.add_argument("--species", default="fixture_jac_figure/jac_species.tsv")
    ap.add_argument("--overlap", choices=["precedence", "mixed", "overlay"],
                    default="precedence")
    ap.add_argument("--outdir", default="figures")
    ap.add_argument("--ext", default="png", choices=["png", "pdf"])
    args = ap.parse_args()

    species = load_species(args.species)
    rr, cc, mm = load_pattern(args.pattern)
    order, pos, bounds = build_order(species)
    os.makedirs(args.outdir, exist_ok=True)

    p1 = os.path.join(args.outdir, f"jacobian_spy_primary.{args.ext}")
    p2 = os.path.join(args.outdir, f"jacobian_coupling_zoom.{args.ext}")
    plot_primary(rr, cc, mm, pos, bounds, order, args.overlap, p1)
    plot_secondary(rr, cc, mm, pos, bounds, order, p2)
    print(f"N={len(order)}  bands: gas={bounds['n_gas']} grain={bounds['n_grain']} "
          f"ice={bounds['n_ice']}")
    print(f"wrote {p1}\nwrote {p2}")


if __name__ == "__main__":
    main()
