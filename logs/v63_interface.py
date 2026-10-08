# -*- coding: utf-8 -*-
"""V6.3 focused interface-value experiment.

Question: why retain the six-scalar interface instead of adding kappa or
keeping full geometry?  Scenario: cached configuration descriptors, repeated
certification queries at multiple error radii.

Branches (same instances, same radii, paired):
  A six scalars   -> Q_G interval LB (J=256)
  B six + kappa   -> L_quad and max(L_quad, Q_G)
  C full geometry -> tolerance-stopped TR reference with cached eigh

Seeds 401-402 are unused elsewhere.  Config frozen to results_v63/
before any timing.  Summary-only clients never touch A, h_hat, kappa or
eigendecompositions at query time; tau uses a summary-only conservative
scale (see tau_summary_only) so no hidden full-geometry access occurs.
"""
import json
import os
from datetime import datetime
from time import perf_counter

import numpy as np

import core
import experiments_c as ec

OUT = "results_v63"
SEEDS = (401, 402)
GAMMA = 1.0
NOISE = ec.NOISE_U
RHO_GRID = np.round(np.logspace(-2.2, -0.6, 10), 5)
QUERIES = (1, 10, 100)
REF_TOL = 1e-12
J_CELLS = 256


def tau_summary_only(scal, r):
    """Summary-only conservative scale: bound |c_q| <= p a^2 + gamma I0^2
    (triangle inequality on h^H A h = |b^Hh|^2 - gamma||C^Hh||^2),
    2r||g|| <= 2r(pa + gamma I0(beta0 + I0))?  g = A h; ||g|| <= ||bb^Hh||
    + gamma||CC^Hh|| <= p a + gamma I0 beta0.  ||A||_2 <= p + gamma beta0^2.
    All quantities derive from the six summaries only."""
    p, a, I0 = scal["p"], scal["a"], scal["I0"]
    b0 = scal["beta0"]
    c_up = p * a * a + gamma_safe * I0 * I0
    g_up = p * a + gamma_safe * I0 * b0
    A_up = p + gamma_safe * b0 * b0
    return 1e-8 * max(1.0, c_up, 2 * r * g_up, r * r * A_up)


gamma_safe = GAMMA  # module-level for the helper


def main():
    os.makedirs(OUT, exist_ok=True)
    cfg = dict(written_at=datetime.now().isoformat(timespec="seconds"),
               seeds=list(SEEDS), gamma=GAMMA, noise=NOISE,
               rho_grid=[float(x) for x in RHO_GRID],
               queries=list(QUERIES), J=J_CELLS, ref_tol=REF_TOL,
               candidate_rules="identical to main pools (6 phases x "
                               "MRT/RZF, per-BS power 1)",
               tau="summary-only conservative scale (no full-geometry "
                   "access at query time)",
               tasks=["certification-only", "exact-decision-with-fallback"])
    with open(os.path.join(OUT, "interface_config.json"), "w",
              encoding="utf-8") as f:
        json.dump(cfg, f, indent=1)

    per_problem = []
    summary_rows = []
    for seed in SEEDS:
        cands = ec.all_candidates(seed)
        # ---- descriptor build (once per user problem) ----
        t0 = perf_counter()
        descs = []
        for cand in cands:
            for u in range(4):
                descs.append(core.spectral_summary(
                    cand["b"][u], cand["h"][u], cand["C"][u], GAMMA, NOISE))
        t_build6 = (perf_counter() - t0) / len(descs)
        t0 = perf_counter()
        kappas = [float(np.linalg.norm(d["g"])) for d in descs]
        t_build7 = (perf_counter() - t0) / len(descs) + t_build6
        t0 = perf_counter()
        eigs = [np.linalg.eigh(d["A"]) for d in descs]
        hnorms = [float(np.linalg.norm(cands[i // 4]["h"][i % 4]))
                  for i in range(len(descs))]
        t_buildC = (perf_counter() - t0) / len(descs) + t_build6

        bytes6 = 6 * 8
        bytes7 = 7 * 8
        n = cands[0]["b"][0].shape[0]
        mC = cands[0]["C"][0].shape[1]
        bytesC = (n + n + n * mC) * 16 + (n * n * 2) * 8  # data + eig cache

        # ---- per-query timing ----
        tq = {"A": [], "B": [], "C": []}
        for i, d in enumerate(descs):
            for rho in RHO_GRID:
                r = rho * hnorms[i]
                tau_s = tau_summary_only(d, r)
                tau_ref = core.engineering_tau(d, r)[0]
                t0 = perf_counter()
                qa = core.interval_lower_bound(d, r, GAMMA, NOISE, J_CELLS,
                                               "G")
                tq["A"].append(perf_counter() - t0)
                t0 = perf_counter()
                qb = core.L_quad(d, r, GAMMA)
                tq["B"].append(perf_counter() - t0)
                t0 = perf_counter()
                ref = core.q_star_reference(d, r, GAMMA, NOISE,
                                            eig=eigs[i], mu_rel_tol=REF_TOL)
                tq["C"].append(perf_counter() - t0)
                per_problem.append(dict(
                    seed=seed, cand=i // 4, user=i % 4, rho=float(rho),
                    r=r, q_G_LB=qa, L_quad=qb,
                    joint=max(qa, qb), q_star=ref["q"],
                    tr_valid=ref["reference_valid"],
                    status=ec.user_status(ref),
                    tau_summary=tau_s, tau_ref=tau_ref,
                    accept_A=bool(qa > tau_s), accept_B=bool(qb > tau_s),
                    accept_joint=bool(max(qa, qb) > tau_s)))
        med = {k: float(np.median(v)) for k, v in tq.items()}
        # repeated-query totals (certification-only task)
        reps = {q: {k: med[k] * q * len(per_problem[:0] or [1]) * 0
                    + med[k] * q} for q in QUERIES for k in tq}
        # coverage + fallback for the exact-decision task
        n_units = 0
        acc_A = acc_B = acc_J = 0
        fallback_A = fallback_B = fallback_J = 0
        units = {}
        for row in per_problem:
            if row["seed"] != seed:
                continue
            units.setdefault((row["cand"], row["rho"]), []).append(row)
        for key, us in units.items():
            n_units += 1
            aA = all(u["accept_A"] for u in us)
            aB = all(u["accept_B"] for u in us)
            aJ = all(u["accept_joint"] for u in us)
            acc_A += aA
            acc_B += aB
            acc_J += aJ
            fallback_A += not aA
            fallback_B += not aB
            fallback_J += not aJ
        summary_rows.append(dict(
            seed=seed, n_user_problems=len(descs),
            build_ms_A=1e3 * t_build6, build_ms_B=1e3 * t_build7,
            build_ms_C=1e3 * t_buildC, bytes_A=bytes6, bytes_B=bytes7,
            bytes_C=bytesC, query_us_A=1e6 * med["A"],
            query_us_B=1e6 * med["B"], query_us_C=1e6 * med["C"],
            units=n_units, direct_accept_A=acc_A, direct_accept_B=acc_B,
            direct_accept_joint=acc_J,
            fallback_A=fallback_A, fallback_B=fallback_B,
            fallback_J=fallback_J))
        print(f"[v63:info] seed {seed}: query med A {1e6*med['A']:.1f}us "
              f"B {1e6*med['B']:.1f}us C {1e6*med['C']:.1f}us; direct "
              f"accept A {acc_A}/{n_units} B {acc_B} J {acc_J}")
    with open(os.path.join(OUT, "interface_details.csv"), "w",
              newline="", encoding="utf-8") as f:
        import csv
        w = csv.DictWriter(f, fieldnames=list(per_problem[0].keys()))
        w.writeheader()
        w.writerows(per_problem)
    with open(os.path.join(OUT, "interface_summary.csv"), "w",
              newline="", encoding="utf-8") as f:
        import csv
        w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        w.writeheader()
        w.writerows(summary_rows)
    r = summary_rows[0]
    print(f"[v63:info] bytes: six {r['bytes_A']}B, 6+k {r['bytes_B']}B, "
          f"full+cache {r['bytes_C']}B")
    return True


if __name__ == "__main__":
    main()
