# -*- coding: utf-8 -*-
"""Optimization round on top of the V4.1 suite.

Joint certification (max of certified lower bounds), staged screening with
configurable order, adaptive interval refinement, fair paired timing with
the validated tolerance-stopped reference, strategy freezing with fresh
confirmation seeds, and small/medium/large scale configurations.

All new results go to results_v41_optimization/; results/, results_review/
and results_v41/ are never modified.  Entry point: run_optimization.py.
"""
import hashlib
import json
import os
from datetime import datetime
from time import perf_counter

import numpy as np

import core
import experiments_c as ec
import utils

OUT_TAG = "results_v41_optimization"
GAMMAS = (0.5, 1.0, 2.0)
RHOS = (0.06, 0.15, 0.30)
EXPLORATION_SEEDS = (11, 22, 33, 44, 55, 66, 77, 88, 99, 110)
VALIDATION_SEEDS = (211, 212, 213, 214, 215)      # diagnostic only now
CONFIRM_SEEDS = (221, 222, 223, 224, 225)         # fresh, frozen before use
SCALE_SEEDS = (231, 232, 233)
ADAPTIVE_J_SEQ = (64, 256, 1024, 4096)
REF_TOL = 1e-12          # validated tolerance-stopped reference everywhere
BOUND_STAGES = ("D", "Q_old", "Q_G", "L_quad")

FAILURES = []


def check(name, ok, detail=""):
    print(f"[opt:{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        FAILURES.append(f"{name}: {detail}")
    return ok


# ==========================================================================
# 1. shared precompute (radius independent), reused by every bound and by
#    the reference across radii
# ==========================================================================
def prepare(cands, gamma, noise=ec.NOISE_U):
    """Per user: spectral summary + one eigendecomposition of A, reused by
    all bound stages, the threshold and every reference call."""
    out = []
    for cand in cands:
        users = []
        for u in range(len(cand["b"])):
            s = core.spectral_summary(cand["b"][u], cand["h"][u],
                                      cand["C"][u], gamma, noise)
            users.append(dict(s=s, eig=np.linalg.eigh(s["A"]),
                              hnorm=float(np.linalg.norm(cand["h"][u]))))
        out.append(users)
    return out


def bound_value(prep_u, r, gamma, kind, noise=ec.NOISE_U, J=256):
    s = prep_u["s"]
    if kind == "D":
        return core.D_bound(s, r, gamma, noise)
    if kind == "L_quad":
        return core.L_quad(s, r, gamma)
    if kind == "Q_old":
        return core.interval_lower_bound(s, r, gamma, noise, J, "old")
    if kind == "Q_G":
        return core.interval_lower_bound(s, r, gamma, noise, J, "G")
    raise ValueError(kind)


# ==========================================================================
# 2. adaptive interval refinement for Q_G
# ==========================================================================
def adaptive_qg(prep_u, r, gamma, tau, noise=ec.NOISE_U, J_seq=ADAPTIVE_J_SEQ,
                margin_rel=1e-9):
    """Certified-LB-driven refinement.

    Returns dict(LB, reason, J_used_list, cells_evaluated) with reason in
    {'accept', 'cannot_cross', 'budget'}.  The stopping upper estimate is
    the value of the one-dimensional objective at feasible cell midpoints
    (an upper bound of the continuous Q_G since every midpoint is feasible
    for the scalar problem), guarded by a rounding margin; stopping never
    declares the true problem infeasible.
    """
    s = prep_u["s"]
    p, a, I0 = s["p"], s["a"], s["I0"]
    alpha, beta, beta0 = s["alpha"], s["beta"], s["beta0"]
    best_lb = -np.inf
    cells = 0
    used = []
    for J in J_seq:
        lb = core.interval_lower_bound(s, r, gamma, noise, J, "G")
        th = np.linspace(0.0, np.pi / 2, J + 1)
        th_m = 0.5 * (th[:-1] + th[1:])
        xm, ym = r * np.sin(th_m), r * np.cos(th_m)
        rho = core.R_envelope(xm, ym, alpha, beta, beta0)
        phi = p * (a - xm) ** 2 - gamma * (I0 + rho) ** 2 - gamma * noise
        ub = float(np.min(phi))
        margin = margin_rel * max(1.0, abs(ub), abs(tau))
        cells += J
        used.append(J)
        best_lb = max(best_lb, lb)
        if best_lb > tau:
            return dict(LB=best_lb, reason="accept", J_used_list=used,
                        cells_evaluated=cells, ub_at_stop=ub)
        if ub < tau - margin:
            return dict(LB=best_lb, reason="cannot_cross", J_used_list=used,
                        cells_evaluated=cells, ub_at_stop=ub)
    return dict(LB=best_lb, reason="budget", J_used_list=used,
                cells_evaluated=cells, ub_at_stop=ub)


# ==========================================================================
# 3. staged screening cascade (the single code path for every pipeline,
#    including the pure-reference baseline: stages = ())
# ==========================================================================
def cascade(cands, prep, gamma, stages, noise=ec.NOISE_U, rhos=RHOS,
            adaptive=False, ref_tol=REF_TOL, collect=False):
    """Same final decision task for all pipelines.

    stages: ordered bound names evaluated on not-yet-certified users; the
    exact reference is always the last stage.  adaptive=True turns the Q_G
    stage into the J-sequence refinement.  Certified users never re-enter a
    stage; a failed bound is never a rejection; a clearly infeasible user
    (valid reference) ends the candidate; invalid references leave the
    candidate unresolved and are reported, never silently decided.
    """
    t0 = perf_counter()
    n_users = len(cands[0]["b"])
    stats = dict(stage_evals={k: 0 for k in stages},
                 stage_pass={k: 0 for k in stages}, ref_calls=0,
                 invalid_refs=0, cells_evaluated=0, adaptive_stops={})
    decisions = []
    detail_rows = []
    for ci, cand in enumerate(cands):
        for rho in rhos:
            certified = [False] * n_users
            passed_by = [None] * n_users
            for st_i, stage in enumerate(stages):
                for u in range(n_users):
                    if certified[u]:
                        continue
                    r = rho * prep[ci][u]["hnorm"]
                    tau, _ = core.engineering_tau(prep[ci][u]["s"], r)
                    if adaptive and stage == "Q_G":
                        ad = adaptive_qg(prep[ci][u], r, gamma, tau,
                                         noise)
                        stats["cells_evaluated"] += ad["cells_evaluated"]
                        stats["adaptive_stops"][ad["reason"]] = \
                            stats["adaptive_stops"].get(ad["reason"], 0) + 1
                        lb = ad["LB"]
                    else:
                        lb = bound_value(prep[ci][u], r, gamma, stage, noise)
                    stats["stage_evals"][stage] += 1
                    if lb > tau:
                        certified[u] = True
                        passed_by[u] = stage if not (adaptive
                                                     and stage == "Q_G") \
                            else "Q_G_adaptive"
                        stats["stage_pass"][stage] += 1
            if all(certified):
                decisions.append("accepted")
                if collect:
                    detail_rows.append(dict(cand=ci, rho=rho,
                                            decision="accepted",
                                            mode="direct",
                                            passed_by=";".join(
                                                str(x) for x in passed_by)))
                continue
            statuses = []
            for u in range(n_users):
                if certified[u]:
                    statuses.append("clearly_feasible")
                    continue
                ref = core.q_star_reference(prep[ci][u]["s"],
                                            rho * prep[ci][u]["hnorm"], gamma,
                                            noise, eig=prep[ci][u]["eig"],
                                            mu_rel_tol=ref_tol)
                stats["ref_calls"] += 1
                if not ref["reference_valid"]:
                    stats["invalid_refs"] += 1
                    statuses.append("invalid")
                    break
                statuses.append(ec.user_status(ref))
                if statuses[-1] == "clearly_infeasible":
                    break
            if "invalid" in statuses:
                decisions.append("unresolved_invalid")
            elif "clearly_infeasible" in statuses:
                decisions.append("rejected")
            elif all(x == "clearly_feasible" for x in statuses):
                decisions.append("accepted")
            else:
                decisions.append("unresolved")
            if collect:
                detail_rows.append(dict(cand=ci, rho=rho,
                                        decision=decisions[-1],
                                        mode="reference",
                                        passed_by=";".join(
                                            str(x) for x in passed_by)))
    return dict(seconds=perf_counter() - t0, decisions=decisions,
                stats=stats, detail_rows=detail_rows)


# ==========================================================================
# 4. complementarity analysis from existing V4.1 details (no recomputation
#    of references; read-only on results_v41/)
# ==========================================================================
def complementarity(cohort, v41_dir="results_v41"):
    import csv
    rows = list(csv.DictReader(
        open(os.path.join(v41_dir, f"{cohort}_details.csv"),
             encoding="utf-8")))
    for r in rows:
        for k in ("gamma", "rho", "q_star", "tau", "D", "L_quad", "Q_old",
                  "Q_G"):
            r[k] = float(r[k])
    units = {}
    for r in rows:
        units.setdefault((r["seed"], r["cand_idx"], r["rho"], r["gamma"]),
                         []).append(r)
    out_rows = []
    cnt = dict(units=len(units), qg=0, lq=0, joint=0, feasible=0,
               only_qg=0, only_lq=0, jadd_neither=0, mis_joint=0,
               viol_joint=0, invalid=0)
    for key, us in sorted(units.items()):
        aq = all(u["Q_G"] > u["tau"] for u in us)
        al = all(u["L_quad"] > u["tau"] for u in us)
        aj = all(max(u["Q_G"], u["L_quad"]) > u["tau"] for u in us)
        feas = all(u["status"] == "clearly_feasible" for u in us)
        infeas = any(u["status"] == "clearly_infeasible" for u in us)
        cnt["qg"] += aq
        cnt["lq"] += al
        cnt["joint"] += aj
        cnt["only_qg"] += aq and not al
        cnt["only_lq"] += al and not aq
        cnt["jadd_neither"] += aj and not aq and not al
        cnt["feasible"] += feas
        cnt["invalid"] += any(u["status"] == "reference_invalid" for u in us)
        cnt["mis_joint"] += aj and infeas
        cnt["viol_joint"] += any(max(u["Q_G"], u["L_quad"]) > u["q_star"] + 1e-8
                                 for u in us)
        # per-user diagnostics for units the joint method newly accepts
        if aj and not al:
            for u in us:
                out_rows.append(dict(
                    seed=key[0], cand_idx=key[1], rho=key[2], gamma=key[3],
                    user=u["user"], note="joint-added over L_quad",
                    Q_G_LB=u["Q_G"], L_quad=u["L_quad"], tau=u["tau"],
                    q_star=u["q_star"], status=u["status"],
                    passes_QG=bool(u["Q_G"] > u["tau"]),
                    passes_Lq=bool(u["L_quad"] > u["tau"]),
                    QG_minus_Lq=u["Q_G"] - u["L_quad"]))
    # user-level partition
    ucnt = dict(n=len(rows), only_qg=0, only_lq=0, both=0, neither=0)
    for r in rows:
        g = r["Q_G"] > r["tau"]
        l = r["L_quad"] > r["tau"]
        ucnt["only_qg"] += g and not l
        ucnt["only_lq"] += l and not g
        ucnt["both"] += g and l
        ucnt["neither"] += not g and not l
    utils.write_csv(os.path.join(OUT_TAG, f"complementarity_{cohort}.csv"),
                    list(out_rows[0].keys()) if out_rows
                    else ["seed"], out_rows)
    res = dict(cohort=cohort, **cnt, user_level=ucnt,
               coverage_qg=cnt["qg"] / cnt["feasible"] if cnt["feasible"]
               else None,
               coverage_lq=cnt["lq"] / cnt["feasible"] if cnt["feasible"]
               else None,
               coverage_joint=cnt["joint"] / cnt["feasible"]
               if cnt["feasible"] else None)
    utils.write_json(os.path.join(OUT_TAG,
                                  f"complementarity_{cohort}.json"), res)
    print(f"[opt:info] complementarity {cohort}: units {cnt['units']} "
          f"feasible {cnt['feasible']} | Q_G {cnt['qg']} L_quad {cnt['lq']} "
          f"joint {cnt['joint']} | onlyQG {cnt['only_qg']} onlyLq "
          f"{cnt['only_lq']} joint-adds-over-neither {cnt['jadd_neither']} "
          f"| violations {cnt['viol_joint']} misaccepts {cnt['mis_joint']}")
    return res


# ==========================================================================
# 5. cascade comparison + fair paired timing on one cohort
# ==========================================================================
PIPELINES = {
    "pure_ref": (),
    "Q_old_ref": ("Q_old",),
    "Q_G_ref": ("Q_G",),
    "L_quad_ref": ("L_quad",),
    "QG_Lq_ref": ("Q_G", "L_quad"),
    "Lq_QG_ref": ("L_quad", "Q_G"),
    "D_QG_Lq_ref": ("D", "Q_G", "L_quad"),
    "D_Lq_QG_ref": ("D", "L_quad", "Q_G"),
    "QG_adaptive_Lq_ref": ("Q_G", "L_quad"),   # adaptive=True variant
}


def run_cascade_study(seeds, cohort, quick=False, pipelines=None,
                      reps=3, label=None):
    """Per (seed, gamma): shared precompute once per pipeline invocation
    (inside the timer), warm-up, interleaved repetitions, decision
    consistency against the pure-reference pipeline."""
    pipelines = pipelines or list(PIPELINES)
    timing_rows = []
    accept_rows = []
    for seed in seeds:
        cands = ec.all_candidates(seed)
        for gamma in GAMMAS:
            runs = {p: [] for p in pipelines}
            for rep in range(reps):
                order = np.roll(np.array(pipelines,
                                         dtype=object), rep)
                for p in order:
                    t0 = perf_counter()
                    prep = prepare(cands, gamma)   # counted work
                    adaptive = p == "QG_adaptive_Lq_ref"
                    out = cascade(cands, prep, gamma,
                                  PIPELINES[p], adaptive=adaptive)
                    runs[p].append((perf_counter() - t0, out))
            base_dec = runs["pure_ref"][0][1]["decisions"]
            for p in pipelines:
                for _, out in runs[p]:
                    if out["decisions"] != base_dec:
                        FAILURES.append(f"{cohort}/{label or ''}{seed}/"
                                        f"{gamma}/{p}: decision mismatch")
                secs = [t for t, _ in runs[p]]
                out0 = runs[p][0][1]
                n = len(base_dec)
                acc = sum(1 for d in out0["decisions"] if d == "accepted")
                timing_rows.append(dict(
                    cohort=cohort, seed=seed, gamma=gamma, pipeline=p,
                    seconds_med=float(np.median(secs)),
                    seconds_min=min(secs), seconds_max=max(secs),
                    ref_calls=out0["stats"]["ref_calls"],
                    invalid_refs=out0["stats"]["invalid_refs"],
                    cells_evaluated=out0["stats"]["cells_evaluated"],
                    units=n, accepted=acc))
                accept_rows.append(dict(
                    cohort=cohort, seed=seed, gamma=gamma, pipeline=p,
                    units=n, accepted=acc,
                    rejected=sum(1 for d in out0["decisions"]
                                 if d == "rejected"),
                    unresolved=sum(1 for d in out0["decisions"]
                                   if d.startswith("unresolved"))))
        print(f"[opt:info] cascade {cohort}: seed {seed} done", flush=True)
    return timing_rows, accept_rows


def summarize_timing(timing_rows, baselines=("pure_ref", "L_quad_ref")):
    """Paired per-(seed,gamma) ratios; explicitly distinguishes per-pair
    speedups from the ratio of medians."""
    by = {}
    for r in timing_rows:
        by.setdefault((r["cohort"], r["seed"], r["gamma"]), {})[r["pipeline"]] = r
    pipelines = sorted({r["pipeline"] for r in timing_rows})
    out = []
    for p in pipelines:
        for base in baselines:
            if p == base:
                continue
            pairs = [v for v in by.values() if p in v and base in v]
            ratios = [v[base]["seconds_med"] / v[p]["seconds_med"]
                      for v in pairs]
            wins = sum(1 for rr in ratios if rr > 1.0)
            med_ratio_of_medians = (np.median([v[base]["seconds_med"]
                                               for v in pairs])
                                    / np.median([v[p]["seconds_med"]
                                                 for v in pairs]))
            out.append(dict(
                pipeline=p, baseline=base, pairs=len(pairs), wins=wins,
                median_pair_ratio=float(np.median(ratios)),
                min_pair_ratio=float(np.min(ratios)),
                max_pair_ratio=float(np.max(ratios)),
                ratio_of_medians=float(med_ratio_of_medians),
                ref_calls_pipeline=sum(v[p]["ref_calls"] for v in pairs),
                ref_calls_baseline=sum(v[base]["ref_calls"] for v in pairs)))
            print(f"[opt:info] {p} vs {base}: wins {wins}/{len(pairs)}, "
                  f"median pair ratio "
                  f"{out[-1]['median_pair_ratio']:.2f}, ratio of medians "
                  f"{out[-1]['ratio_of_medians']:.2f}, refs "
                  f"{out[-1]['ref_calls_pipeline']} vs "
                  f"{out[-1]['ref_calls_baseline']}")
    return out


# ==========================================================================
# 6. adaptive vs fixed resolution (acceptance must be identical to
#    fixed J_max; only cost may differ)
# ==========================================================================
def adaptive_vs_fixed(seeds, cohort, quick=False):
    rows = []
    agree_all = True
    for seed in seeds:
        cands = ec.all_candidates(seed)
        for gamma in GAMMAS:
            prep = prepare(cands, gamma)
            for rho in RHOS:
                # acceptance at fixed resolutions vs adaptive
                acc = {}
                for J in (256, 1024, 4096):
                    ok_units = 0
                    for ci in range(len(cands)):
                        good = True
                        for u in range(len(cands[ci]["b"])):
                            r = rho * prep[ci][u]["hnorm"]
                            tau, _ = core.engineering_tau(prep[ci][u]["s"], r)
                            lb = core.interval_lower_bound(
                                prep[ci][u]["s"], r, gamma, ec.NOISE_U, J,
                                "G")
                            if not lb > tau:
                                good = False
                                break
                        ok_units += good
                    acc[J] = ok_units
                ada_units = 0
                cells_adaptive = 0
                cells_fixed256 = 0
                for ci in range(len(cands)):
                    good = True
                    for u in range(len(cands[ci]["b"])):
                        r = rho * prep[ci][u]["hnorm"]
                        tau, _ = core.engineering_tau(prep[ci][u]["s"], r)
                        ad = adaptive_qg(prep[ci][u], r, gamma, tau)
                        cells_adaptive += ad["cells_evaluated"]
                        cells_fixed256 += 256
                        if ad["reason"] != "accept":
                            good = False
                            break
                    ada_units += good
                rows.append(dict(cohort=cohort, seed=seed, gamma=gamma,
                                 rho=rho,
                                 acc_fixed256=acc[256], acc_fixed1024=acc[1024],
                                 acc_fixed4096=acc[4096], acc_adaptive=ada_units,
                                 cells_adaptive=cells_adaptive,
                                 cells_fixed256=cells_fixed256))
                agree_all &= (ada_units == acc[4096])
    utils.write_csv(os.path.join(OUT_TAG, "adaptive_vs_fixed.csv"),
                    list(rows[0].keys()), rows)
    tot = {k: sum(r[k] for r in rows) for k in
           ("acc_fixed256", "acc_fixed4096", "acc_adaptive", "cells_adaptive",
            "cells_fixed256")}
    print(f"[opt:info] adaptive vs fixed ({cohort}): accepts 256/4096/"
          f"adaptive = {tot['acc_fixed256']}/{tot['acc_fixed4096']}/"
          f"{tot['acc_adaptive']}; cells adaptive {tot['cells_adaptive']} vs "
          f"fixed-256 {tot['cells_fixed256']}")
    check("adaptive acceptance equals fixed J=4096 acceptance", agree_all)
    return rows, agree_all


# ==========================================================================
# 7. strategy freeze + fresh-seed confirmation
# ==========================================================================
def _sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def freeze_strategy(selected, selection_rule, timing_rows, quick):
    used = set(EXPLORATION_SEEDS) | set(VALIDATION_SEEDS) | \
        set(CONFIRM_SEEDS) | set(SCALE_SEEDS) | set(ec.SEEDS)
    assert not (set(CONFIRM_SEEDS) & (used - set(CONFIRM_SEEDS))), \
        "confirmation seeds must be unused"
    here = os.path.dirname(os.path.abspath(__file__))
    cfg = dict(
        written_at=datetime.now().isoformat(timespec="seconds"),
        mode="quick" if quick else "full",
        selected_pipeline=selected,
        selection_rule=selection_rule,
        selection_data="exploration cohort only (validation seeds already "
                       "observed in the V4.1 round; used as diagnostics "
                       "only, not as independent confirmation)",
        confirmation_seeds=list(CONFIRM_SEEDS),
        gammas=list(GAMMAS), rhos=list(RHOS),
        acceptance="certified interval lower bound > tau per user; joint = "
                   "max over staged bounds; nonpositive never a rejection",
        tau_rule="1e-8 * max(1,|c_q|,2r||g||,r^2||A||_2)",
        adaptive_J_seq=list(ADAPTIVE_J_SEQ),
        adaptive_stop_margin_rel=1e-9,
        reference=dict(mu_rel_tol=REF_TOL,
                       note="validated tolerance-stopped multiplier search "
                            "(|dq|<=5.7e-10 vs 200-iteration, validity "
                            "preserved)"),
        timing=dict(repetitions=3, interleaved=True, warmup=True,
                    per_user_eigh_cached=True, skip_certified_users=True,
                    early_stop="first clearly infeasible user",
                    candidate_generation_timed=False),
        exploration_timing_summary=[
            dict(pipeline=r["pipeline"], seed=r["seed"], gamma=r["gamma"],
                 seconds_med=r["seconds_med"], ref_calls=r["ref_calls"])
            for r in timing_rows],
        code_hashes={f: _sha256(os.path.join(here, f)) for f in
                     ("core.py", "experiments_c.py", "optimization_v41.py")},
    )
    utils.write_json(os.path.join(OUT_TAG, "frozen_strategy_config.json"),
                     cfg)
    return cfg


def select_pipeline(timing_rows):
    """Mechanical rule fixed before confirmation: among pipelines whose
    stage set contains both Q_G and L_quad (joint acceptance), pick the
    lowest median seconds over exploration pairs; tie-break on reference
    calls."""
    cand = [p for p in PIPELINES
            if set(PIPELINES[p]) >= {"Q_G", "L_quad"}]
    agg = {}
    for r in timing_rows:
        if r["pipeline"] in cand:
            a = agg.setdefault(r["pipeline"], [0.0, 0, 0])
            a[0] += r["seconds_med"]
            a[1] += r["ref_calls"]
            a[2] += 1
    ranked = sorted(agg.items(), key=lambda kv: (kv[1][0], kv[1][1]))
    selected = ranked[0][0]
    rule = ("argmin over {Q_G, L_quad}-containing pipelines of summed "
            "median seconds on the exploration cohort; tie-break on total "
            "reference calls. Single-bound pipelines are excluded because "
            "their accepted set is strictly smaller than the joint set.")
    print(f"[opt:info] selected pipeline: {selected} "
          f"(totals: {[(k, round(v[0], 3), v[1]) for k, v in ranked]})")
    return selected, rule


# ==========================================================================
# 8. scale experiment (generalized candidate generation)
# ==========================================================================
def _cn(rng, shape, var):
    return (rng.normal(size=shape) + 1j * rng.normal(size=shape)) \
        * np.sqrt(var / 2.0)


def gen_candidates(seed, n_ant, users_per_cell, n_ris=16):
    """Same generation rules as the main pool, generalized: per-BS total
    power 1 split equally over its users (per-stream 1/users_per_cell),
    same/cross cell scales 1/0.5, RZF regularization 0.1, first phase all
    ones + 5 random 2-bit states, MRT and RZF beams, 12 candidates."""
    rng = np.random.default_rng(seed)
    n_bs = 2
    n_users = n_bs * users_per_cell
    G = [_cn(rng, (n_ris, n_ant), 0.25) for _ in range(n_bs)]
    f = [_cn(rng, (n_ris,), 1.0 / 16.0) for _ in range(n_users)]
    d = [[_cn(rng, (n_ant,), 0.25) for _ in range(n_users)]
         for _ in range(n_bs)]
    thetas = [np.ones(n_ris, dtype=complex)]
    for _ in range(5):
        thetas.append(ec.PHASE_VALUES[rng.integers(0, 4, size=n_ris)])
    per_stream = 1.0 / users_per_cell
    own = {j: tuple(range(j * users_per_cell, (j + 1) * users_per_cell))
           for j in range(n_bs)}
    own_bs = {u: j for j in range(n_bs) for u in own[j]}
    cands = []
    for theta in thetas:
        h = [[None] * n_users for _ in range(n_bs)]
        for j in range(n_bs):
            for u in range(n_users):
                base = d[j][u].reshape(-1, 1) \
                    + G[j].conj().T @ (theta * f[u]).reshape(-1, 1)
                h[j][u] = base if own_bs[u] == j else 0.5 * base
        for btype in ("MRT", "RZF"):
            W = {}
            for j in range(n_bs):
                if btype == "MRT":
                    for u in own[j]:
                        w = h[j][u] / np.linalg.norm(h[j][u])
                        W[u] = w * np.sqrt(per_stream)
                else:
                    Hj = np.hstack([h[j][u] for u in own[j]])
                    Wj = Hj @ np.linalg.inv(Hj.conj().T @ Hj
                                            + 0.1 * np.eye(users_per_cell))
                    for i, u in enumerate(own[j]):
                        w = Wj[:, [i]] / np.linalg.norm(Wj[:, [i]])
                        W[u] = w * np.sqrt(per_stream)
            beam_full = {}
            for u in range(n_users):
                col = np.zeros((n_bs * n_ant, 1), dtype=complex)
                col[own_bs[u] * n_ant:(own_bs[u] + 1) * n_ant, :] = W[u]
                beam_full[u] = col
            b_list, C_list, hs = [], [], []
            for u in range(n_users):
                hs.append(np.vstack([h[0][u], h[1][u]]))
                b_list.append(beam_full[u])
                C_list.append(np.hstack([beam_full[v] for v in range(n_users)
                                         if v != u]))
            cands.append(dict(seed=seed, phase_idx=None, btype=btype, h=hs,
                              b=b_list, C=C_list))
    return cands


SCALE_CONFIGS = [
    dict(tag="S", n_ant=2, users_per_cell=2),
    dict(tag="M", n_ant=4, users_per_cell=2),
    dict(tag="L", n_ant=8, users_per_cell=3),
]
SCALE_GAMMA, SCALE_RHO = 1.0, 0.15


def run_scale(quick=False, pipelines=("pure_ref", "Q_G_ref", "L_quad_ref",
                                      "QG_adaptive_Lq_ref")):
    rows = []
    for cfgp in SCALE_CONFIGS:
        n_dim = 2 * cfgp["n_ant"]
        n_users = 2 * cfgp["users_per_cell"]
        per_stream = 1.0 / cfgp["users_per_cell"]
        for seed in SCALE_SEEDS[:2] if quick else SCALE_SEEDS:
            cands = gen_candidates(seed, cfgp["n_ant"],
                                   cfgp["users_per_cell"])
            prep = prepare(cands, SCALE_GAMMA)
            runs = {p: [] for p in pipelines}
            for rep in range(3):
                for p in np.roll(np.array(pipelines), rep):
                    t0 = perf_counter()
                    out = cascade(cands, prep, SCALE_GAMMA, PIPELINES[p],
                                  rhos=(SCALE_RHO,),
                                  adaptive=(p == "QG_adaptive_Lq_ref"))
                    runs[p].append((perf_counter() - t0, out))
            base = runs["pure_ref"][0][1]["decisions"]
            for p in pipelines:
                for _, out in runs[p]:
                    if out["decisions"] != base:
                        FAILURES.append(f"scale {cfgp['tag']}/{seed}/{p}: "
                                        "decision mismatch")
                d = runs[p][0][1]["decisions"]
                rows.append(dict(
                    config=cfgp["tag"], n_ant=cfgp["n_ant"],
                    users_per_cell=cfgp["users_per_cell"], dim=n_dim,
                    n_users=n_users, per_stream_power=per_stream,
                    total_bs_power=1.0, seed=seed, pipeline=p,
                    units=len(d),
                    accepted=sum(1 for x in d if x == "accepted"),
                    rejected=sum(1 for x in d if x == "rejected"),
                    seconds_med=float(np.median([t for t, _ in runs[p]])),
                    ref_calls=runs[p][0][1]["stats"]["ref_calls"]))
        print(f"[opt:info] scale {cfgp['tag']} done", flush=True)
    utils.write_csv(os.path.join(OUT_TAG, "scale_results.csv"),
                    list(rows[0].keys()), rows)
    return rows


# ==========================================================================
# 9. verification suite
# ==========================================================================
def run_verification(confirm_seeds):
    # (a) bound safety on fresh confirmation instances: every implemented
    # lower bound (D, Q_old, Q_G fixed and adaptive) <= fresh exact q*
    worst = 0.0
    n_checked = 0
    for seed in confirm_seeds[:2]:
        cands = ec.all_candidates(seed)
        for gamma in GAMMAS:
            prep = prepare(cands, gamma)
            for rho in RHOS:
                for ci in range(len(cands)):
                    for u in range(len(cands[ci]["b"])):
                        pu = prep[ci][u]
                        r = rho * pu["hnorm"]
                        tau, _ = core.engineering_tau(pu["s"], r)
                        ref = core.q_star_reference(pu["s"], r, gamma,
                                                    ec.NOISE_U,
                                                    eig=pu["eig"],
                                                    mu_rel_tol=REF_TOL)
                        if not ref["reference_valid"]:
                            FAILURES.append(f"verification: invalid ref "
                                            f"{seed}/{gamma}/{rho}/{u}")
                            continue
                        vals = [bound_value(pu, r, gamma, k)
                                for k in BOUND_STAGES]
                        vals.append(adaptive_qg(pu, r, gamma, tau)["LB"])
                        worst = max(worst, max(vals) - ref["q"])
                        n_checked += 1
    check("all implemented lower bounds <= exact reference",
          worst <= 1e-8, f"max excess={worst:.2e} over {n_checked} checks")

    # (b) stage order does not change the accepted set (only cost)
    cands = ec.all_candidates(confirm_seeds[0])
    prep = prepare(cands, 1.0)
    a = cascade(cands, prep, 1.0, ("Q_G", "L_quad"))["decisions"]
    b = cascade(cands, prep, 1.0, ("L_quad", "Q_G"))["decisions"]
    check("stage order invariance of accepted set", a == b)
    d1 = cascade(cands, prep, 1.0, ("Q_G", "L_quad"))["decisions"]
    d2 = cascade(cands, prep, 1.0, ("Q_G", "L_quad"))["decisions"]
    check("determinism of repeated runs", d1 == d2)

    # (c) joint acceptance superset of singles (from complementarity JSONs)
    for cohort in ("exploration", "validation"):
        with open(os.path.join(OUT_TAG, f"complementarity_{cohort}.json"),
                  encoding="utf-8") as fh:
            cj = json.load(fh)
        check(f"joint acceptance contains singles ({cohort})",
              cj["joint"] >= max(cj["qg"], cj["lq"]))
    return worst


# ==========================================================================
# 10. figures from saved data
# ==========================================================================
PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
GRAY = "#666666"


def make_figures():
    import csv

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.grid": True,
                         "grid.alpha": 0.25, "axes.spines.top": False,
                         "axes.spines.right": False, "figure.dpi": 150,
                         "savefig.bbox": "tight", "legend.frameon": False})
    figdir = os.path.join(OUT_TAG, "figures")
    os.makedirs(figdir, exist_ok=True)

    def rd(name):
        with open(os.path.join(OUT_TAG, name), encoding="utf-8") as f:
            return list(csv.DictReader(f))

    # fig 1: complementarity
    comp = {c: json.load(open(os.path.join(OUT_TAG,
                                           f"complementarity_{c}.json"),
                              encoding="utf-8"))
            for c in ("exploration", "validation")}
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.6))
    x = np.arange(2)
    w = 0.26
    for i, (m, lab) in enumerate((("qg", "$Q_G$"), ("lq", "$L_{\\rm quad}$"),
                                  ("joint", "joint (max)"))):
        vals = [comp[c][m] for c in ("exploration", "validation")]
        feas = [comp[c]["feasible"] for c in ("exploration", "validation")]
        ax1.bar(x + (i - 1) * w, vals, w, color=PAL[i], label=lab)
        for xx, v, f_ in zip(x + (i - 1) * w, vals, feas):
            ax1.text(xx, v + 4, f"{v}/{f_}", ha="center", fontsize=7,
                     color="#444444")
    ax1.set_xticks(x)
    ax1.set_xticklabels(["exploration (1080 units)",
                         "validation (540 units)"])
    ax1.set_ylabel("directly accepted units")
    ax1.set_title("(a) unit-level certification by method")
    ax1.legend(fontsize=8)
    ul = [comp[c]["user_level"] for c in ("exploration", "validation")]
    parts = [("both", "certified by both", PAL[2]),
             ("only_qg", "only $Q_G$", PAL[0]),
             ("only_lq", "only $L_{\\rm quad}$", PAL[1]),
             ("neither", "neither", "#bbbbbb")]
    bottom = np.zeros(2)
    for k, lab, col in parts:
        vals = np.array([u[k] for u in ul], dtype=float)
        ax2.bar(x, vals, 0.5, bottom=bottom, color=col, label=lab)
        bottom += vals
    ax2.set_xticks(x)
    ax2.set_xticklabels(["exploration (4320 users)",
                         "validation (2160 users)"])
    ax2.set_ylabel("user-level problems")
    ax2.set_title("(b) per-user complementarity")
    ax2.legend(fontsize=7.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "fig_opt_complementarity.pdf"))
    fig.savefig(os.path.join(figdir, "fig_opt_complementarity.png"), dpi=150)
    plt.close(fig)

    # fig 2: confirmation timing, paired per (seed, gamma)
    tt = rd("confirmation_timings.csv")
    if tt:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.6))
        by = {}
        for r in tt:
            by.setdefault((r["seed"], r["gamma"]), {})[r["pipeline"]] = \
                float(r["seconds_med"])
        pipes = sorted({r["pipeline"] for r in tt} - {"pure_ref"})
        for i, p in enumerate(pipes):
            ratios = [v["pure_ref"] / v[p] for v in by.values() if p in v]
            ax1.scatter(np.full(len(ratios), i) + np.random.default_rng(3)
                        .uniform(-0.12, 0.12, len(ratios)), ratios, s=18,
                        color=PAL[i % len(PAL)], zorder=3)
            ax1.hlines(np.median(ratios), i - 0.25, i + 0.25,
                       color=PAL[i % len(PAL)], lw=2, zorder=4)
        ax1.axhline(1.0, color=GRAY, lw=0.9)
        ax1.set_xticks(range(len(pipes)))
        ax1.set_xticklabels(pipes, rotation=20, ha="right", fontsize=7.5)
        ax1.set_ylabel("paired speedup vs pure reference")
        ax1.set_title("(a) per (seed,$\\gamma$) paired ratios")
        calls = {p: sum(int(r["ref_calls"]) for r in tt if r["pipeline"] == p)
                 for p in pipes + ["pure_ref"]}
        ax2.bar(range(len(calls)), list(calls.values()), 0.6,
                color=PAL[:len(calls)])
        ax2.set_xticks(range(len(calls)))
        ax2.set_xticklabels(list(calls), rotation=20, ha="right", fontsize=7.5)
        ax2.set_ylabel("total reference calls")
        ax2.set_title("(b) reference calls (confirmation batch)")
        fig.tight_layout()
        fig.savefig(os.path.join(figdir, "fig_opt_timing.pdf"))
        fig.savefig(os.path.join(figdir, "fig_opt_timing.png"), dpi=150)
        plt.close(fig)
    print(f"[opt:info] figures written to {figdir}")


# ==========================================================================
# 11. driver
# ==========================================================================
def run_all(quick=False):
    utils.ensure_dirs()
    print("--- optimization: complementarity from existing data ---")
    comp = {c: complementarity(c) for c in ("exploration", "validation")}

    print("--- optimization: cascade comparison + timing (exploration) ---")
    ex_pipelines = list(PIPELINES)
    if quick:
        ex_pipelines = ["pure_ref", "Q_G_ref", "L_quad_ref", "QG_Lq_ref",
                        "QG_adaptive_Lq_ref"]
    ex_tim, ex_acc = run_cascade_study(
        EXPLORATION_SEEDS[:2] if quick else EXPLORATION_SEEDS,
        "exploration", quick, pipelines=ex_pipelines)
    utils.write_csv(os.path.join(OUT_TAG, "exploration_cascade_timings.csv"),
                    list(ex_tim[0].keys()), ex_tim)
    utils.write_csv(os.path.join(OUT_TAG, "exploration_cascade_accepts.csv"),
                    list(ex_acc[0].keys()), ex_acc)
    pair_rows = summarize_timing(ex_tim)
    utils.write_csv(os.path.join(OUT_TAG, "exploration_paired_timing.csv"),
                    list(pair_rows[0].keys()), pair_rows)

    print("--- optimization: adaptive vs fixed resolution ---")
    av_rows, _ = adaptive_vs_fixed(
        EXPLORATION_SEEDS[:2] if quick else EXPLORATION_SEEDS,
        "exploration", quick)

    print("--- optimization: freeze strategy ---")
    selected, rule = select_pipeline(ex_tim)
    frozen = freeze_strategy(selected, rule, ex_tim, quick)

    print("--- optimization: confirmation on fresh seeds ---")
    conf_pipelines = ["pure_ref", "Q_G_ref", "L_quad_ref", selected]
    cf_tim, cf_acc = run_cascade_study(
        CONFIRM_SEEDS[:2] if quick else CONFIRM_SEEDS,
        "confirmation", quick, pipelines=conf_pipelines)
    for name, rr in (("confirmation_timings", cf_tim),
                     ("confirmation_accepts", cf_acc)):
        utils.write_csv(os.path.join(OUT_TAG, name + ".csv"),
                        list(rr[0].keys()), rr)
    cf_pairs = summarize_timing(cf_tim)
    utils.write_csv(os.path.join(OUT_TAG, "confirmation_paired_timing.csv"),
                    list(cf_pairs[0].keys()), cf_pairs)

    print("--- optimization: verification suite ---")
    run_verification(CONFIRM_SEEDS[:2] if quick else CONFIRM_SEEDS)

    print("--- optimization: scale experiment ---")
    sc_rows = run_scale(quick)

    print("--- optimization: figures ---")
    make_figures()

    manifest = dict(
        mode="quick" if quick else "full",
        finished_at=datetime.now().isoformat(timespec="seconds"),
        selected_pipeline=selected,
        failures=list(FAILURES),
        complementarity={k: dict(units=v["units"], feasible=v["feasible"],
                                 qg=v["qg"], lq=v["lq"], joint=v["joint"])
                         for k, v in comp.items()})
    utils.write_json(os.path.join(OUT_TAG, "run_manifest.json"), manifest)
    print(f"\n=== optimization suite {'FAILED' if FAILURES else 'PASSED'} ===")
    for f in FAILURES:
        print("  FAIL:", f)
    return not FAILURES
