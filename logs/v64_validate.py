# -*- coding: utf-8 -*-
"""V6.4 validation: strict clients vs core implementations.

1. seven-scalar client == core.L_quad (which is itself scalar-form).
2. compact TR cache == core.q_star_reference across validation seeds and
   the degenerate battery (PSD interior, singular, hard case, r=0,
   repeated eigenvalues).
3. six-scalar client == core.interval_lower_bound.
"""
import os

import numpy as np

import core
import experiments_c as ec
import summary_clients as sc

OUT = "results_v64"
FAILS = []


def check(name, ok, detail=""):
    print(f"[v64:{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        FAILS.append(name)


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(64)
    gamma, noise = 1.0, 0.05

    # --- degenerate battery for compact cache -------------------------
    cases = [
        ("psd_interior", np.eye(3), np.array([[.1], [.2], [0]]), 5.0),
        ("singular_solvable", np.diag([1.0, 0.0]),
         np.array([[1.0], [0.0]]), 5.0),
        ("hard_eq", np.diag([-1.0, 1.0]), np.array([[1.0], [0.0]]), 0.5),
        ("repeated_min", np.diag([-1.0, -1.0, 2.0]),
         np.array([[0.0], [0.0], [1.0]]), 0.7),
        ("negcurv", np.diag([-2.0, 1.0]), np.array([[0.3], [0.4]]), 0.5),
    ]
    worst_cache = 0.0
    for tag, A, g, r in cases:
        lam, U = np.linalg.eigh(np.asarray(A, complex))
        cache = dict(lam=lam,
                     w=np.abs(U.conj().T @ np.asarray(g, complex)) .ravel() ** 2,
                     c=1.0,
                     scale=max(1.0, float(np.abs(A).max()),
                               float(np.linalg.norm(g))))
        val, mode = sc.tr_qstar_from_cache(cache, r)
        ref = core.trust_region_min(A, g, r)
        ref_total = 1.0 + ref["quad_val"]
        worst_cache = max(worst_cache, abs(val - ref_total))
        check(f"cache/{tag}", abs(val - ref_total) < 1e-8,
              f"cache={val:.10f} ref_total={ref_total:.10f} mode={mode}")

    # --- pool instances (seeds 401-402, full grids) --------------------
    worst_lq = worst_qg = 0.0
    n = 0
    for seed in (401, 402):
        cands = ec.all_candidates(seed)
        for gamma_ in (0.5, 1.0, 2.0):
            for cand in cands:
                for u in range(4):
                    s = core.spectral_summary(cand["b"][u], cand["h"][u],
                                              cand["C"][u], gamma_, noise)
                    six = (s["p"], s["a"], s["I0"], s["alpha"], s["beta"],
                           s["beta0"])
                    kappa = float(np.linalg.norm(s["g"]))
                    r = 0.15 * float(np.linalg.norm(cand["h"][u]))
                    # seven-scalar client vs core
                    lq_c = sc.lquad_client(six, kappa, r, gamma_, noise)
                    worst_lq = max(worst_lq, abs(lq_c - core.L_quad(s, r,
                                                                    gamma_)))
                    # six-scalar client vs core
                    qg_c = sc.qg_lb_client(six, r, gamma_, noise, 256)
                    worst_qg = max(worst_qg, abs(
                        qg_c - core.interval_lower_bound(s, r, gamma_, noise,
                                                         256, "G")))
                    # compact cache vs reference
                    cache = sc.build_compact_cache(cand["b"][u],
                                                   cand["h"][u],
                                                   cand["C"][u], gamma_,
                                                   noise)
                    val, _ = sc.tr_qstar_from_cache(cache, r)
                    ref = core.q_star_reference(s, r, gamma_, noise,
                                                mu_rel_tol=1e-12)
                    worst_cache = max(worst_cache, abs(ref["q"] - val))
                    n += 1
    check("client/L_quad_matches_core", worst_lq < 1e-12,
          f"max diff {worst_lq:.2e} over {n}")
    check("client/QG_matches_core", worst_qg < 1e-12,
          f"max diff {worst_qg:.2e}")
    check("cache/matches_reference", worst_cache < 1e-7,
          f"max diff {worst_cache:.2e} over {n + len(cases)} problems")
    print(f"validated {n} pool problems + {len(cases)} degenerate cases")
    return not FAILS


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
