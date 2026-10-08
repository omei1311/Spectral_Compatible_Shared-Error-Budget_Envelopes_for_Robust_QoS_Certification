# -*- coding: utf-8 -*-
"""V6.4 repeated-query and fallback experiment (strict clients).

Same instances as V6.3 (seeds 401-402, 10-radius grid, gamma=1).  Batch
lengths 1/10/100 over the SAME deterministic radius sequence (longer
batches cycle the grid); batch totals are MEASURED, not extrapolated.
Warm-up + 3 interleaved repetitions with method rotation; raw reps kept.

Branches (query-side cost only; descriptor build reported separately):
  A six-scalar Q_G_LB client        (48 B payload, float64)
  B seven-scalar L_quad client      (56 B payload)
  C seven-scalar joint rule         (max(L_quad, Q_G_LB); both computed)
  D compact TR cache                (lam, w, c + scale, actual nbytes)
  E full TR reference               (supplementary audit branch)

Tasks: (1) summary-only certification (accept / unresolved);
       (2) summary screening + compact-TR fallback (exact decisions).
Fallback = summary query fails AND compact TR queried for that user.
'candidates_without_direct_cert' counts candidates needing fallback;
'tr_queries' counts actual TR query invocations (early stop on clearly
infeasible user)."""
import csv
import hashlib
import json
import os
from datetime import datetime
from time import perf_counter

import numpy as np

import core
import experiments_c as ec
import summary_clients as sc

OUT = "results_v64"
SEEDS = (401, 402)
GAMMA = 1.0
NOISE = ec.NOISE_U
RHO = np.round(np.logspace(-2.2, -0.6, 10), 5)
BATCH = (1, 10, 100)
REPS = 3
J = 256


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def main():
    os.makedirs(OUT, exist_ok=True)
    cfg = dict(written_at=datetime.now().isoformat(timespec="seconds"),
               seeds=list(SEEDS), gamma=GAMMA, noise=NOISE,
               rho_grid=[float(x) for x in RHO], batch=list(BATCH),
               reps=REPS, J=J, scenario="cached descriptors, repeated "
               "absolute-radius queries; radius sequence = the grid, "
               "cycled for longer batches",
               branches="A six-scalar QG | B L_quad | C joint | "
                        "D compact TR | E full TR (audit)",
               storage="float64 payload bytes; compact cache = actual "
                       "array nbytes (lam, w, c, scale)",
               code_hashes={f: sha(f) for f in
                            ("summary_clients.py", "v64_experiment.py")})
    with open(os.path.join(OUT, "config.json"), "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=1)

    detail, timing, fallback_rows = [], [], []
    for seed in SEEDS:
        cands = ec.all_candidates(seed)
        # ---------- build side ----------
        t0 = perf_counter()
        sixes, kappas, caches = [], [], []
        for cand in cands:
            for u in range(4):
                s = core.spectral_summary(cand["b"][u], cand["h"][u],
                                          cand["C"][u], GAMMA, NOISE)
                sixes.append((s["p"], s["a"], s["I0"], s["alpha"],
                              s["beta"], s["beta0"]))
                kappas.append(float(np.linalg.norm(s["g"])))
                caches.append(sc.build_compact_cache(cand["b"][u],
                                                     cand["h"][u],
                                                     cand["C"][u], GAMMA,
                                                     NOISE))
        build_six = perf_counter() - t0
        t0 = perf_counter()
        full_cache = []
        for cand in cands:
            for u in range(4):
                s = core.spectral_summary(cand["b"][u], cand["h"][u],
                                          cand["C"][u], GAMMA, NOISE)
                full_cache.append((s, np.linalg.eigh(s["A"])))
        build_full = build_six + (perf_counter() - t0)
        per_build6 = 1e3 * build_six / len(sixes)
        per_buildD = 1e3 * (build_six + (perf_counter() - t0) * 0) / 1  # D uses same summary; recompute properly below
        # compact-cache build timing (separate)
        t0 = perf_counter()
        for cand in cands:
            for u in range(4):
                sc.build_compact_cache(cand["b"][u], cand["h"][u],
                                       cand["C"][u], GAMMA, NOISE)
        per_buildD = 1e3 * (perf_counter() - t0) / len(sixes)
        bytes_D = float(np.mean([sc.compact_cache_bytes(c) for c in
                                 caches]))
        n_dim = cands[0]["b"][0].shape[0]
        mC = cands[0]["C"][0].shape[1]
        bytes_E = float((n_dim + n_dim + n_dim * mC) * 16
                        + (n_dim * n_dim * 2) * 8)

        # ---------- per-user radius list (absolute) ----------
        hnorms = [float(np.linalg.norm(cands[i // 4]["h"][i % 4]))
                  for i in range(len(sixes))]
        radii = [[rho * hn for rho in RHO] for hn in hnorms]

        # ---------- batch timing (measured) ----------
        def run_batch(branch, nq):
            t0 = perf_counter()
            for i in range(0, nq, len(RHO)):
                for j, r in enumerate(radii[i % len(sixes)]):
                    if branch == "A":
                        sc.qg_lb_client(sixes[i % len(sixes)], r, GAMMA,
                                        NOISE, J)
                    elif branch == "B":
                        sc.lquad_client(sixes[i % len(sixes)], kappas[
                            i % len(sixes)], r, GAMMA, NOISE)
                    elif branch == "C":
                        sc.joint_client(sixes[i % len(sixes)], kappas[
                            i % len(sixes)], r, GAMMA, NOISE, J)
                    elif branch == "D":
                        sc.tr_qstar_from_cache(caches[i % len(sixes)], r)
            return perf_counter() - t0

        # warm-up
        for br in "ABCD":
            run_batch(br, 10)
        rows_t = []
        for rep in range(REPS):
            order = np.roll(np.array(list("ABCD")), rep)
            for br in order:
                for nq in BATCH:
                    rows_t.append(dict(seed=seed, rep=rep, branch=br,
                                       n_queries=nq,
                                       batch_seconds=run_batch(br, nq)))
        timing.extend(rows_t)

        # ---------- coverage + fallback task ----------
        tau = [sc.tau_summary(sixes[i][0], sixes[i][1], sixes[i][2],
                              sixes[i][5], GAMMA, NOISE, radii[i][0])
               for i in range(len(sixes))]
        units = {}
        for i in range(len(sixes)):
            cand_i, u = i // 4, i % 4
            for j, rho in enumerate(RHO):
                r = radii[i][j]
                t = sc.tau_summary(sixes[i][0], sixes[i][1], sixes[i][2],
                                   sixes[i][5], GAMMA, NOISE, r)
                qa = sc.qg_lb_client(sixes[i], r, GAMMA, NOISE, J)
                qb = sc.lquad_client(sixes[i], kappas[i], r, GAMMA, NOISE)
                qj = max(qa, qb)
                qd, _ = sc.tr_qstar_from_cache(caches[i], r)
                ref = core.q_star_reference(
                    core.spectral_summary(cands[cand_i]["b"][u],
                                          cands[cand_i]["h"][u],
                                          cands[cand_i]["C"][u], GAMMA,
                                          NOISE), r, GAMMA, NOISE,
                    mu_rel_tol=1e-12)
                st = ec.user_status(ref)
                units.setdefault((cand_i, float(rho)), []).append(
                    dict(qa=qa, qb=qb, qj=qj, qd=qd, tau=t, st=st,
                         i=i, r=r))
                detail.append(dict(seed=seed, cand=cand_i, user=u,
                                   rho=float(rho), r=r, tau=t, Q_G=qa,
                                   L_quad=qb, joint=qj, compact=qd,
                                   q_star=ref["q"], status=st))
        for (cand_i, rho), us in units.items():
            accA = all(u["qa"] > u["tau"] for u in us)
            accB = all(u["qb"] > u["tau"] for u in us)
            accJ = all(u["qj"] > u["tau"] for u in us)
            feas = all(u["st"] == "clearly_feasible" for u in us)
            # fallback task for the joint branch
            tr_q = 0
            final = None
            if accJ:
                final = "accepted"
            else:
                for u in us:
                    tr_q += 1
                    if u["qd"] < -u["tau"]:
                        final = "rejected"
                        break
                if final is None:
                    if all(u["qd"] > u["tau"] for u in us):
                        final = "accepted"
                    else:
                        final = "unresolved"
            fallback_rows.append(dict(seed=seed, cand=cand_i, rho=rho,
                                     direct_A=accA, direct_B=accB,
                                     direct_J=accJ, ref_feasible=feas,
                                     tr_queries=tr_q, final=final))
        fa = [f for f in fallback_rows if f["seed"] == seed]
        print(f"[v64:info] seed {seed}: units {len(fa)}; direct A "
              f"{sum(f['direct_A'] for f in fa)} B "
              f"{sum(f['direct_B'] for f in fa)} J "
              f"{sum(f['direct_J'] for f in fa)}; ref-feasible "
              f"{sum(f['ref_feasible'] for f in fa)}; TR queries "
              f"{sum(f['tr_queries'] for f in fa)}")
    with open(os.path.join(OUT, "batch_timings.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(timing[0].keys()))
        w.writeheader()
        w.writerows(timing)
    with open(os.path.join(OUT, "per_problem.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(detail[0].keys()))
        w.writeheader()
        w.writerows(detail)
    with open(os.path.join(OUT, "fallback_units.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(fallback_rows[0].keys()))
        w.writeheader()
        w.writerows(fallback_rows)
    # summary aggregates
    from collections import defaultdict
    agg = defaultdict(list)
    for t in timing:
        agg[(t["branch"], t["n_queries"])].append(t["batch_seconds"])
    print("\nbranch  nq  median_batch_s  per_query_us")
    for (br, nq), v in sorted(agg.items()):
        print(f"  {br}  {nq:3d}  {np.median(v):.6f}  "
              f"{1e6*np.median(v)/(nq*len(sixes)):.2f}")
    fa = fallback_rows
    n = len(fa)
    print(f"\nunits {n}: direct A {sum(f['direct_A'] for f in fa)} "
          f"B {sum(f['direct_B'] for f in fa)} J "
          f"{sum(f['direct_J'] for f in fa)}; ref-feasible "
          f"{sum(f['ref_feasible'] for f in fa)}; TR queries "
          f"{sum(f['tr_queries'] for f in fa)}; unresolved "
          f"{sum(1 for f in fa if f['final']=='unresolved')}")
    print(f"build: six {per_build6:.3f} ms/user, compact "
          f"{per_buildD:.3f} ms/user; bytes: six 48, seven 56, "
          f"compact {bytes_D:.0f}, full+cache {bytes_E:.0f}")
    json.dump(dict(direct_A=sum(f["direct_A"] for f in fa),
                   direct_B=sum(f["direct_B"] for f in fa),
                   direct_J=sum(f["direct_J"] for f in fa),
                   ref_feasible=sum(f["ref_feasible"] for f in fa),
                   tr_queries=sum(f["tr_queries"] for f in fa),
                   unresolved=sum(1 for f in fa
                                  if f["final"] == "unresolved"),
                   units=n,
                   build_ms_six=per_build6, build_ms_compact=per_buildD,
                   bytes_compact=bytes_D, bytes_full=bytes_E),
              open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    return True


if __name__ == "__main__":
    main()
