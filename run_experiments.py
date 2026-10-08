# -*- coding: utf-8 -*-
"""Run all experiments (A, B, C), produce figures, config and the
post-run verification checklist."""
import os
import platform
import sys

import matplotlib
import numpy as np
import scipy

import experiments_a
import experiments_b
import experiments_c
import figures
import utils


def write_config():
    cfg = dict(
        seed_A1=experiments_a.RNG_SEED_A1,
        seeds_C=experiments_c.SEEDS,
        rho_list=experiments_c.RHO_LIST,
        gamma_C=experiments_c.GAMMA_U,
        noise_power_C=experiments_c.NOISE_U,
        A3_params=experiments_a.A3_PARAMS,
        B=dict(p=experiments_b.P, a=experiments_b.A_TGT,
               I0=experiments_b.I0_TGT, alpha=experiments_b.ALPHA,
               beta=experiments_b.BETA, gamma=experiments_b.GAMMA,
               noise_power=experiments_b.NOISE, r=experiments_b.R,
               beta0_grid="linspace(1.02, sqrt(2), 21)"),
        interval_J_default=256, interval_J_escalation=[256, 1024, 4096],
        tau_rule="1e-8 * max(1, |c|, 2r||g||, r^2||A||_2); accept only if "
                 "LB > tau; near-zero => ambiguous",
        versions=dict(python=sys.version, numpy=np.__version__,
                      scipy=scipy.__version__, matplotlib=matplotlib.__version__,
                      platform=platform.platform()),
    )
    utils.write_json(os.path.join(utils.RESULTS_DIR, "config.json"), cfg)


def post_run_verification(ares, brows, cres):
    print("\n=== post-run verification checklist ===")
    ok = True

    # 1. R upper-bound property + branch coverage
    worst = 0.0
    rng = np.random.default_rng(99)
    for _ in range(2000):
        alpha, beta = rng.uniform(0.2, 1.5, 2)
        beta0 = rng.uniform(max(alpha, beta),
                            np.sqrt(alpha ** 2 + beta ** 2))
        x, y = rng.uniform(0.0, 0.6, 2)
        import core
        worst = max(worst, float(core.R_envelope(x, y, alpha, beta, beta0)
                                 - min(beta0 * np.hypot(x, y),
                                       alpha * x + beta * y)))
    branches = set()
    delta = np.sqrt(1.1 ** 2 - 1)
    nu = np.sqrt(1.1 ** 2 - 1)
    for (x, y) in [(0.2, 0.4), (0.4, 0.3), (0.45, 0.1)]:
        branches.add(1 if delta * x <= nu * y else
                     (3 if delta * x >= 1.0 * y else 2))
    print(f"1a. R <= min(beta0*norm, alpha x + beta y): worst = {worst:.2e}")
    print(f"1b. branch coverage of fixed examples: {sorted(branches)} "
          f"(need {{1,2,3}})")
    ok &= worst <= 1e-12 and branches == {1, 2, 3}

    # 2. six-scalar construction error (A2)
    m2 = max(max(abs(r["match_alpha"]), abs(r["match_beta"]),
                 abs(r["match_beta0"])) for r in ares["a2"])
    print(f"2. A2 six-scalar construction max error: {m2:.2e}")
    ok &= m2 < 1e-10

    # 3. Q_G_LB <= q* (B)
    b_bad = [r for r in brows
             if r["q_star_minus_QG_LB"] < -1e-9 and r["reference_valid"]]
    print(f"3a. B: points with Q_G_LB > q* + tol: {len(b_bad)}")
    ok &= len(b_bad) == 0

    # 4. Q_G_LB <= q* (C) and misaccepts
    c_bad = [r for r in cres["rows"]
             if r["reference_valid"] and r["Q_G"] > r["q_star"] + 1e-9]
    n_mis = sum(v for v in
                [cres["agg"]["stats"][m]["misaccept"] for m in
                 experiments_c.METHODS])
    print(f"3b. C: user problems with Q_G_LB > q* + tol: {len(c_bad)}")
    print(f"4.  C: misaccepts (all methods): {n_mis}")
    ok &= len(c_bad) == 0 and n_mis == 0

    # 5. artifacts exist
    arts = ["direction_checks.csv", "construction_checks.csv",
            "boundary_checks.csv", "spectral_scan.csv",
            "pool_user_details.csv", "pool_summary.csv", "timings.csv",
            "config.json"]
    missing = [f for f in arts
               if not os.path.exists(os.path.join(utils.RESULTS_DIR, f))]
    figs = [f for f in ["fig1_direction_construction.pdf",
                        "fig2_spectral_scan.pdf"]
            if not os.path.exists(os.path.join(utils.FIGURES_DIR, f))]
    print(f"5. artifacts missing: results={missing} figures={figs}")
    ok &= not missing and not figs
    print(f"=== verification {'PASSED' if ok else 'FAILED'} ===")
    return ok


def main():
    utils.ensure_dirs()
    print("--- Experiment A ---")
    ares = experiments_a.run_all()
    print("\n--- Experiment B ---")
    brows = experiments_b.run_b()
    print("\n--- Experiment C ---")
    cres = experiments_c.run_c()
    print("\n--- figures ---")
    figures.make_fig1(ares)
    figures.make_fig2(brows)
    write_config()
    ok = post_run_verification(ares, brows, cres)
    print("\nAll experiments complete. verification_ok =", ok)
    if not ok:
        raise RuntimeError('Post-run verification failed')


if __name__ == "__main__":
    main()
