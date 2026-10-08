# -*- coding: utf-8 -*-
"""Experiment B: 21-point spectral parameter scan.

Fixed p=4, a=1.3, I0=0.5, alpha=beta=1, gamma=1, sigma^2=1.5, r=0.5;
beta0 on 21 equispaced points from 1.02 to sqrt(2).  Per point the actual
matrices/channel are constructed (a family of instances, not one physical
system); the channel norm changes along the scan and is recorded.
"""
import os

import numpy as np

import core
import utils

P, A_TGT, I0_TGT = 4.0, 1.3, 0.5
ALPHA = BETA = 1.0
GAMMA, NOISE, R = 1.0, 1.5, 0.5
N_POINTS = 21
B0_GRID = np.linspace(1.02, np.sqrt(2.0), N_POINTS)


def instance_at_beta0(beta0):
    delta = np.sqrt(beta0 ** 2 - 1.0)
    CH = np.array([[1.0, delta, 0, 0, 0],
                   [0.0, 0.0, 1.0, 0, 0],
                   [0.0, 0.0, 0.0, 0, 0]], dtype=complex)
    h = np.array([A_TGT, -1.8 / delta, 0, 0, 0], dtype=complex).reshape(-1, 1)
    b = np.array([2.0, 0, 0, 0, 0], dtype=complex).reshape(-1, 1)
    return CH, h, b, delta


def run_b():
    utils.ensure_dirs()
    rows = []
    for i, beta0 in enumerate(B0_GRID):
        CH, h, b, delta = instance_at_beta0(beta0)
        s = core.spectral_summary(b, h, CH.conj().T, GAMMA, NOISE)
        six_err = max(abs(s["p"] - P), abs(s["a"] - A_TGT),
                      abs(s["I0"] - I0_TGT), abs(s["alpha"] - ALPHA),
                      abs(s["beta"] - BETA), abs(s["beta0"] - beta0))
        D = core.D_bound(s, R, GAMMA, NOISE)
        q_old, _, _ = core.Q_approx(s, R, GAMMA, NOISE, "old")
        qg, _, _ = core.Q_approx(s, R, GAMMA, NOISE, "G")
        lq = core.L_quad(s, R, GAMMA)
        ref = core.q_star_reference(s, R, GAMMA, NOISE)
        q_lb = core.interval_lower_bound(s, R, GAMMA, NOISE, 256, "old")
        qg_lb = core.interval_lower_bound(s, R, GAMMA, NOISE, 256, "G")
        tau = ref["tau"]
        rows.append(dict(
            idx=i, beta0=beta0, beta0_measured=s["beta0"], delta=delta,
            x_axis=beta0 / np.sqrt(ALPHA ** 2 + BETA ** 2),
            h_norm=float(np.linalg.norm(h)), six_scalar_err=six_err,
            D=D, Q_old=q_old, Q_G=qg, L_quad=lq, q_star=ref["q"],
            Q_LB=q_lb, Q_G_LB=qg_lb, tau=tau,
            reference_valid=ref["reference_valid"], ref_mode=ref["mode"],
            Q_G_minus_Q_old=qg - q_old, q_star_minus_Q_G=ref["q"] - qg,
            q_star_minus_QG_LB=ref["q"] - qg_lb))
    fields = ["idx", "beta0", "beta0_measured", "delta", "x_axis", "h_norm",
              "six_scalar_err", "D", "Q_old", "Q_G", "L_quad", "q_star",
              "Q_LB", "Q_G_LB", "tau", "reference_valid", "ref_mode",
              "Q_G_minus_Q_old", "q_star_minus_Q_G", "q_star_minus_QG_LB"]
    utils.write_csv(os.path.join(utils.RESULTS_DIR, "spectral_scan.csv"),
                    fields, rows)

    # summary diagnostics
    strict = [r for r in rows if r["Q_G_minus_Q_old"] > 1e-9]
    print(f"[B] {len(rows)} points; six-scalar max err over scan: "
          f"{max(r['six_scalar_err'] for r in rows):.3e}")
    print(f"[B] strict improvement Q_G > Q_old (tol 1e-9): {len(strict)} "
          f"points, beta0 in [{min(r['beta0'] for r in strict):.4f}, "
          f"{max(r['beta0'] for r in strict):.4f}]"
          if strict else "[B] no strict improvement points")
    neg_gap = [r for r in rows
               if r["q_star_minus_Q_G"] < -1e-9 and r["reference_valid"]]
    print(f"[B] points with q* < Q_G - tol (potential violation): "
          f"{len(neg_gap)}")
    lb_gap = min(r["q_star_minus_QG_LB"] for r in rows if r["reference_valid"])
    print(f"[B] min over points of (q* - Q_G_LB): {lb_gap:.6e} (must be >= ~0)")
    lq_vs = [(r["beta0"], r["Q_G"] - r["L_quad"]) for r in rows]
    n_lq_better = sum(1 for _, d in lq_vs if d < -1e-9)
    print(f"[B] points where L_quad > Q_G: {n_lq_better} "
          f"(no presupposed dominance)")
    flips = [r for r in rows if r["Q_G_LB"] > r["tau"] and r["Q_LB"] <= 0]
    print(f"[B] certification flips (Q_G_LB > tau while Q_LB <= 0): "
          f"{len(flips)} points"
          + (f", beta0 from {min(r['beta0'] for r in flips):.4f} to "
             f"{max(r['beta0'] for r in flips):.4f}" if flips else ""))
    print(f"[B] h norm range along scan: "
          f"{min(r['h_norm'] for r in rows):.2f} .. "
          f"{max(r['h_norm'] for r in rows):.2f}")
    return rows


if __name__ == "__main__":
    run_b()
