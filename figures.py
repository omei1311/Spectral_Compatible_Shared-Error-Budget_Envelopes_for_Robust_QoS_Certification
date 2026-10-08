# -*- coding: utf-8 -*-
"""Figures 1 and 2 (2 panels each), generated only from computed results.

Categorical palette = validated reference instance (dataviz skill), fixed
slot order; every series additionally carries a distinct marker/line style
so identity is never color-alone.  Q/Q_G/L_quad curves are numerical 1-D
approximations and are labelled as such -- certification uses only the
interval lower bounds.
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import utils

PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
GRAY = "#666666"

plt.rcParams.update({
    "font.size": 9, "axes.grid": True, "grid.alpha": 0.25,
    "grid.linewidth": 0.5, "axes.spines.top": False,
    "axes.spines.right": False, "figure.dpi": 150,
    "savefig.bbox": "tight", "legend.frameon": False,
})


def make_fig1(a):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.6))

    # (a) A1: closed-form R vs direct scalar extremum reference
    rand = [r for r in a["a1"] if r["kind"] == "random"]
    ref = np.array([r["ref_max"] for r in rand])
    cl = np.array([r["R_closed"] for r in rand])
    gm = np.array([r["grid_max"] for r in rand])
    lo = 0.9 * min(gm.min(), ref.min())
    hi = 1.02 * max(ref.max(), cl.max())
    ax1.plot([lo, hi], [lo, hi], color=GRAY, lw=0.8, zorder=1)
    ax1.scatter(ref, gm, s=16, color=PAL[1], marker="s", zorder=2,
                label="4097-point grid only")
    ax1.scatter(ref, cl, s=16, color=PAL[0], zorder=3,
                label="closed form $R$")
    ax1.set_xlabel("reference value $\\sqrt{\\max_s f(s)}$")
    ax1.set_ylabel("candidate value")
    ax1.set_title("(a) A1: directional extremum, 30 random sets")
    ax1.legend(loc="upper left")
    txt = (f"max $|R^2-\\mathrm{{ref}}^2|$ = "
           f"{max(r['abs_err_R2'] for r in rand):.1e}\n"
           f"max grid underestimation = "
           f"{max(r['grid_under'] for r in rand):.1e}")
    ax1.text(0.97, 0.05, txt, transform=ax1.transAxes, ha="right",
             va="bottom", fontsize=7.5, color="#444444")

    # (b) A2/A3: construction quality vs eps
    a3 = a["a3"]["construction"]
    eps3 = np.array([r["eps"] for r in a3])
    conv = np.array([abs(r["obj_minus_QG"]) for r in a3])
    a2 = a["a2"]
    eps2 = np.array([r["eps"] for r in a2])
    match = np.array([max(abs(r["match_alpha"]), abs(r["match_beta"]),
                          abs(r["match_beta0"])) for r in a2])
    lam = np.array([max(r["lam_min_M"], 0.0) for r in a2])
    FLOOR = 1e-17  # exact zeros (match errors are 0 to machine precision)
    ax2.loglog(eps3, conv, "o-", color=PAL[0], ms=5, lw=1.2,
               label="A3: $|q(e^*)-Q_G|$ at constructed error")
    ax2.loglog(eps2, np.maximum(match, FLOOR), "s--", color=PAL[2], ms=5,
               lw=1.2, label="A2: max six-scalar match error (=0, floored)")
    ax2.loglog(eps2, np.maximum(lam, FLOOR), "^:", color=PAL[4], ms=5,
               lw=1.2,
               label="A2: $\\lambda_{\\min}(T-\\ell_\\varepsilon\\ell_\\varepsilon^H)$")
    ax2.axhline(1e-16, color=GRAY, lw=0.7, ls=":")
    ax2.text(0.98, 1.3e-16, "machine precision", transform=ax2.get_yaxis_transform(),
             ha="right", va="bottom", fontsize=7, color="#444444")
    ax2.set_xlabel("construction shrinkage $\\varepsilon$")
    ax2.set_ylabel("error / residual")
    ax2.set_title("(b) A2/A3: construction quality vs $\\varepsilon$")
    ax2.legend(loc="lower left")
    fig.tight_layout()
    utils.save_fig(fig, "fig1_direction_construction")
    plt.close(fig)


def make_fig2(rows):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.6))
    x = np.array([r["x_axis"] for r in rows])

    def col(k):
        return np.array([r[k] for r in rows])

    series = [("D", PAL[0], "-", "o", 3.0),
              ("Q_old", PAL[1], "-", "s", 3.0),
              ("Q_G", PAL[2], "-", "D", 3.0),
              ("L_quad", PAL[3], "-", "^", 3.0),
              ("q_star", PAL[4], "-", "v", 3.0)]
    labels = {"D": "$D$", "Q_old": "$Q$ (old, approx.)",
              "Q_G": "$Q_G$ (approx.)", "L_quad": "$L_{\\rm quad}$",
              "q_star": "exact $q^*$"}
    for k, c, ls, mk, ms in series:
        ax1.plot(x, col(k), ls, color=c, marker=mk, ms=ms, lw=1.3,
                 label=labels[k])
    ax1.plot(x, col("Q_G_LB"), "--", color=PAL[2], lw=1.3,
             label="$Q_{G,\\rm LB}$ (certified, $J$=256)")
    ax1.axhline(0.0, color=GRAY, lw=0.9)
    ax1.set_xlabel("$\\beta_0/\\sqrt{\\alpha^2+\\beta^2}$")
    ax1.set_ylabel("worst-case margin")
    ax1.set_title("(a) B: spectral scan margins")
    ax1.legend(loc="lower left", fontsize=7.5, ncol=2)

    ax2.plot(x, col("q_star_minus_Q_G"), "-", color=PAL[4], marker="v",
             ms=3, lw=1.2, label="$q^*-Q_G$")
    ax2.plot(x, col("q_star_minus_QG_LB"), "-", color=PAL[2], marker="D",
             ms=3, lw=1.2, label="$q^*-Q_{G,\\rm LB}$")
    ax2.plot(x, col("Q_G_minus_Q_old"), "-", color=PAL[1], marker="s",
             ms=3, lw=1.2, label="$Q_G-Q$ (old)")
    ax2.plot(x, col("Q_G") - col("L_quad"), "-", color=PAL[3], marker="^",
             ms=3, lw=1.2, label="$Q_G-L_{\\rm quad}$")
    ax2.axhline(0.0, color=GRAY, lw=0.9)
    ax2.set_xlabel("$\\beta_0/\\sqrt{\\alpha^2+\\beta^2}$")
    ax2.set_ylabel("gap")
    ax2.set_title("(b) B: gaps between methods")
    ax2.legend(loc="best", fontsize=7.5)
    fig.tight_layout()
    utils.save_fig(fig, "fig2_spectral_scan")
    plt.close(fig)
