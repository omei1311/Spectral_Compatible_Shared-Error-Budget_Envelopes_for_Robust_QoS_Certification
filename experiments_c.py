# -*- coding: utf-8 -*-
"""Experiment C: small fixed RIS candidate pool.

2 BS x 4 antennas, 2 users per cell (4 users), 16-element RIS with 2-bit
phases, 5 seeds {11,22,33,44,55}, 6 phase states x 2 beam types = 12
candidates per seed, relative radii rho in {0.01, 0.03, 0.06}.
Total user-level ball problems: 5*12*3*4 = 720, dimension 8.

Normalized dimensionless links (no geographic deployment, path-loss fitting
or hardware model): G_j ~ CN(0,1/4), f_u ~ CN(0,1/16), d_{j,u} ~ CN(0,1/4);
same-cell scaling s=1, cross-cell s=0.5.  The ball r_u = rho*||h_u|| is a
specified effective-channel ball, not a claim of calibrated statistical
coverage.
"""
import os
from time import perf_counter

import numpy as np

import core
import utils

SEEDS = [11, 22, 33, 44, 55]
RHO_LIST = [0.01, 0.03, 0.06]
GAMMA_U = 0.5
NOISE_U = 0.05
N_BS, N_ANT, N_USER, N_RIS = 2, 4, 4, 16
OWN = {0: (0, 1), 1: (2, 3)}
OWN_BS = {0: 0, 1: 0, 2: 1, 3: 1}
PHASE_VALUES = np.array([1.0, 1.0j, -1.0, -1.0j], dtype=complex)
METHODS = ("D", "Q_old", "Q_G", "L_quad")


# --------------------------------------------------------------------------
# channel / candidate generation
# --------------------------------------------------------------------------
def _cn(rng, shape, var):
    return (rng.normal(size=shape) + 1j * rng.normal(size=shape)) \
        * np.sqrt(var / 2.0)


def generate_seed(seed):
    rng = np.random.default_rng(seed)
    G = [_cn(rng, (N_RIS, N_ANT), 0.25) for _ in range(N_BS)]
    f = [_cn(rng, (N_RIS,), 1.0 / 16.0) for _ in range(N_USER)]
    d = [[_cn(rng, (N_ANT,), 0.25) for _ in range(N_USER)] for _ in range(N_BS)]
    thetas = [np.ones(N_RIS, dtype=complex)]
    for _ in range(5):
        thetas.append(PHASE_VALUES[rng.integers(0, 4, size=N_RIS)])
    return G, f, d, thetas


def effective_channels(G, f, d, theta):
    """h_{j,u} = s_{j,u} [ d_{j,u} + G_j^H diag(theta) f_u ]  (4x1 each)."""
    h = [[None] * N_USER for _ in range(N_BS)]
    for j in range(N_BS):
        for u in range(N_USER):
            base = d[j][u].reshape(-1, 1) \
                + G[j].conj().T @ (theta * f[u]).reshape(-1, 1)
            h[j][u] = base if OWN_BS[u] == j else 0.5 * base
    return h


def beams_for(h, btype):
    """Per-BS beams for its 2 own users; per-stream power 0.5."""
    W = [[None, None] for _ in range(N_BS)]
    for j in range(N_BS):
        users = OWN[j]
        if btype == "MRT":
            for i, u in enumerate(users):
                w = h[j][u] / np.linalg.norm(h[j][u])
                W[j][i] = w * np.sqrt(0.5)
        elif btype == "RZF":
            Hj = np.hstack([h[j][u] for u in users])
            Wj = Hj @ np.linalg.inv(Hj.conj().T @ Hj + 0.1 * np.eye(2))
            for i, u in enumerate(users):
                w = Wj[:, [i]] / np.linalg.norm(Wj[:, [i]])
                W[j][i] = w * np.sqrt(0.5)
        else:
            raise ValueError(btype)
    return W


def build_unit(h, W):
    """Stacked 8x1 channels, per-user beams b_u (8x1) and interferers
    C_u (8x3, stacked beams of the other 3 users)."""
    hs = [np.vstack([h[0][u], h[1][u]]) for u in range(N_USER)]
    beam_full = [None] * N_USER
    for u in range(N_USER):
        j = OWN_BS[u]
        i = OWN[j].index(u)
        col = np.zeros((2 * N_ANT, 1), dtype=complex)
        col[j * N_ANT:(j + 1) * N_ANT, :] = W[j][i]
        beam_full[u] = col
    b_list, C_list = [], []
    for u in range(N_USER):
        b_list.append(beam_full[u])
        C_list.append(np.hstack([beam_full[v] for v in range(N_USER)
                                 if v != u]))
    return hs, b_list, C_list


def all_candidates(seed):
    """12 candidates: 6 phase states x {MRT, RZF}."""
    G, f, d, thetas = generate_seed(seed)
    cands = []
    for p_idx, theta in enumerate(thetas):
        h = effective_channels(G, f, d, theta)
        for btype in ("MRT", "RZF"):
            W = beams_for(h, btype)
            hs, b_list, C_list = build_unit(h, W)
            cands.append(dict(seed=seed, phase_idx=p_idx, btype=btype,
                              h=hs, b=b_list, C=C_list))
    return cands


# --------------------------------------------------------------------------
# per-user evaluation helpers
# --------------------------------------------------------------------------
def user_reference(b, h, C, eig_cache):
    """A, g, c + eigendecomposition once per user; q* for all radii."""
    s = core.spectral_summary(b, h, C, GAMMA_U, NOISE_U)
    eig = np.linalg.eigh(s["A"])
    eig_cache.append(eig)
    return s, eig


def bound_values(s, r):
    return dict(D=core.D_bound(s, r, GAMMA_U, NOISE_U),
                Q_old=core.interval_lower_bound(s, r, GAMMA_U, NOISE_U,
                                                256, "old"),
                Q_G=core.interval_lower_bound(s, r, GAMMA_U, NOISE_U,
                                              256, "G"),
                L_quad=core.L_quad(s, r, GAMMA_U))


def user_status(ref_out):
    if not ref_out["reference_valid"]:
        return "reference_invalid"
    if ref_out["q"] > ref_out["tau"]:
        return "clearly_feasible"
    if ref_out["q"] < -ref_out["tau"]:
        return "clearly_infeasible"
    return "ambiguous"


def unit_status(st_list):
    if any(x == "reference_invalid" for x in st_list):
        return "reference_invalid"
    if any(x == "clearly_infeasible" for x in st_list):
        return "clearly_infeasible"
    if all(x == "clearly_feasible" for x in st_list):
        return "clearly_feasible"
    return "ambiguous"


# --------------------------------------------------------------------------
# truth pass (untimed): everything computed and recorded
# --------------------------------------------------------------------------
def truth_pass(cands):
    rows = []
    for c_idx, cand in enumerate(cands):
        for u in range(N_USER):
            cache = []
            s, eig = user_reference(cand["b"][u], cand["h"][u], cand["C"][u],
                                    cache)
            h_norm = float(np.linalg.norm(cand["h"][u]))
            for rho in RHO_LIST:
                r = rho * h_norm
                bv = bound_values(s, r)
                ref = core.q_star_reference(s, r, GAMMA_U, NOISE_U, eig=eig)
                st = user_status(ref)
                row = dict(seed=cand["seed"], cand_idx=c_idx,
                           phase_idx=cand["phase_idx"], btype=cand["btype"],
                           user=u, rho=rho, r=r, h_norm=h_norm,
                           tau=ref["tau"], status=st,
                           ref_mode=ref["mode"],
                           reference_valid=ref["reference_valid"],
                           q_star=ref["q"])
                row.update(bv)
                rows.append(row)
    return rows


def aggregate(rows):
    """Candidate(=seed,cand,rho)-level metrics per method."""
    units = {}
    for row in rows:
        key = (row["seed"], row["cand_idx"], row["rho"])
        units.setdefault(key, dict(users=[]))
        units[key]["users"].append(row)
    stats = {m: dict(accepted=0, misaccept=0, unknown=0) for m in METHODS}
    n_units = len(units)
    ref_feasible_units = 0
    ref_infeasible_units = 0
    ambiguous_units = 0
    invalid_units = 0
    qg_added = 0
    misaccept_details = []
    for key, unit in units.items():
        sts = [r["status"] for r in unit["users"]]
        us = unit_status(sts)
        if us == "clearly_feasible":
            ref_feasible_units += 1
        elif us == "clearly_infeasible":
            ref_infeasible_units += 1
        elif us == "reference_invalid":
            invalid_units += 1
        else:
            ambiguous_units += 1
        acc = {}
        for m in METHODS:
            ok = all(r[m] > r["tau"] for r in unit["users"])
            acc[m] = ok
            if ok:
                stats[m]["accepted"] += 1
                if us == "clearly_infeasible":
                    stats[m]["misaccept"] += 1
                    misaccept_details.append((key, m))
            else:
                stats[m]["unknown"] += 1
        if acc["Q_G"] and not acc["Q_old"]:
            qg_added += 1
    out = dict(n_units=n_units, stats=stats,
               ref_feasible_units=ref_feasible_units,
               ref_infeasible_units=ref_infeasible_units,
               ambiguous_units=ambiguous_units, invalid_units=invalid_units,
               qg_added=qg_added, misaccept_details=misaccept_details)
    # coverage: accepted AND reference-clearly-feasible
    for m in METHODS:
        cov_count = 0
        for key, unit in units.items():
            if unit_status([r["status"] for r in unit["users"]]) \
                    == "clearly_feasible" \
                    and all(r[m] > r["tau"] for r in unit["users"]):
                cov_count += 1
        out["coverage_" + m] = (cov_count / ref_feasible_units
                                if ref_feasible_units else None)
    return out, units


# --------------------------------------------------------------------------
# timing passes
# --------------------------------------------------------------------------
def time_screening(cands):
    """Summary construction + 4-method bound evaluation, all radii."""
    t0 = perf_counter()
    for cand in cands:
        for u in range(N_USER):
            s = core.spectral_summary(cand["b"][u], cand["h"][u],
                                      cand["C"][u], GAMMA_U, NOISE_U)
            h_norm = float(np.linalg.norm(cand["h"][u]))
            for rho in RHO_LIST:
                r = rho * h_norm
                core.engineering_tau(s, r)
                bound_values(s, r)
    return perf_counter() - t0


def time_full_reference(cands):
    """All 144 user-level reference solves per seed (one eigh reused over
    the 3 radii)."""
    t0 = perf_counter()
    for cand in cands:
        for u in range(N_USER):
            b, h, C = cand["b"][u], cand["h"][u], cand["C"][u]
            A = b @ b.conj().T - GAMMA_U * (C @ C.conj().T)
            g = A @ h
            c = float(np.real((h.conj().T @ A @ h)[0, 0])) - GAMMA_U * NOISE_U
            eig = np.linalg.eigh(A)
            h_norm = float(np.linalg.norm(h))
            for rho in RHO_LIST:
                r = rho * h_norm
                core.trust_region_min(A, g, r, eig=eig)
                # tau for status decisions (scale data)
                s_scale = max(1.0, abs(c),
                              2.0 * r * float(np.linalg.norm(g)),
                              r * r * float(np.abs(eig[0]).max()
                                            if eig[0].size else 1.0))
                _ = 1e-8 * s_scale
    return perf_counter() - t0


def time_end_to_end(cands):
    """Screening + serial per-user reference fallback for unknown
    candidates, with early stop at the first clearly-infeasible user.
    Computes everything fresh (no truth cache)."""
    t0 = perf_counter()
    fallback_calls = 0
    accepted = 0
    for cand in cands:
        summaries = []
        h_norms = []
        eig_cache = {}
        for u in range(N_USER):
            s = core.spectral_summary(cand["b"][u], cand["h"][u],
                                      cand["C"][u], GAMMA_U, NOISE_U)
            summaries.append(s)
            h_norms.append(float(np.linalg.norm(cand["h"][u])))
        for rho in RHO_LIST:
            ok = True
            for u in range(N_USER):
                r = rho * h_norms[u]
                tau, _ = core.engineering_tau(summaries[u], r)
                lb = core.interval_lower_bound(summaries[u], r, GAMMA_U,
                                               NOISE_U, 256, "G")
                if not (lb > tau):
                    ok = False
                    break
            if ok:
                accepted += 1
                continue
            # serial fallback with early stop
            for u in range(N_USER):
                b, h, C = cand["b"][u], cand["h"][u], cand["C"][u]
                A = b @ b.conj().T - GAMMA_U * (C @ C.conj().T)
                g = A @ h
                if u not in eig_cache:
                    eig_cache[u] = np.linalg.eigh(A)
                eig = eig_cache[u]
                r = rho * h_norms[u]
                s = summaries[u]
                tau, _ = core.engineering_tau(s, r)
                ref = core.q_star_reference(s, r, GAMMA_U, NOISE_U, eig=eig)
                fallback_calls += 1
                if user_status(ref) == "clearly_infeasible":
                    break
    return perf_counter() - t0, fallback_calls, accepted


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------
def run_c():
    utils.ensure_dirs()
    all_rows = []
    timing_rows = []
    total_fallback = 0
    for seed in SEEDS:
        cands = all_candidates(seed)
        rows = truth_pass(cands)
        all_rows.extend(rows)
    agg, units = aggregate(all_rows)
    n_user_problems = len(all_rows)

    fields = ["seed", "cand_idx", "phase_idx", "btype", "user", "rho", "r",
              "h_norm", "tau", "D", "Q_old", "Q_G", "L_quad", "q_star",
              "status", "ref_mode", "reference_valid"]
    utils.write_csv(os.path.join(utils.RESULTS_DIR, "pool_user_details.csv"),
                    fields, all_rows)

    # summary CSV (candidate-level over seed x cand x rho units)
    n_u = agg["n_units"]
    sum_rows = []
    for m in METHODS:
        st = agg["stats"][m]
        cov = agg["coverage_" + m]
        sum_rows.append(dict(
            method=m, units=n_u, accepted=st["accepted"],
            accept_rate=st["accepted"] / n_u,
            misaccept=st["misaccept"],
            unknown=st["unknown"], unknown_rate=st["unknown"] / n_u,
            coverage=(cov if cov is not None else "NA"),
            qg_added_vs_Qold=(agg["qg_added"] if m == "Q_G" else "")))
    sum_rows.append(dict(method="_reference", units=n_u, accepted="",
                         accept_rate="", misaccept="",
                         unknown="", unknown_rate="",
                         coverage="",
                         qg_added_vs_Qold=""))
    sum_rows.append(dict(
        method="_counts", units=n_u,
        accepted=f"ref_feasible={agg['ref_feasible_units']}",
        accept_rate=f"ref_infeasible={agg['ref_infeasible_units']}",
        misaccept=f"ambiguous={agg['ambiguous_units']}",
        unknown=f"invalid={agg['invalid_units']}",
        unknown_rate=f"user_problems={n_user_problems}", coverage="",
        qg_added_vs_Qold=""))
    utils.write_csv(os.path.join(utils.RESULTS_DIR, "pool_summary.csv"),
                    ["method", "units", "accepted", "accept_rate", "misaccept",
                     "unknown", "unknown_rate", "coverage",
                     "qg_added_vs_Qold"], sum_rows)

    print(f"[C] user-level problems: {n_user_problems} (cap 720), "
          f"candidate-radius units: {n_u}")
    for m in METHODS:
        st = agg["stats"][m]
        cov = agg["coverage_" + m]
        cov_s = "NA" if cov is None else f"{cov:.4f}"
        print(f"[C] {m:7s} accepted={st['accepted']:3d}/{n_u} "
              f"misaccept={st['misaccept']} unknown={st['unknown']:3d} "
              f"coverage={cov_s}")
    print(f"[C] reference: feasible={agg['ref_feasible_units']} "
          f"infeasible={agg['ref_infeasible_units']} "
          f"ambiguous={agg['ambiguous_units']} "
          f"invalid={agg['invalid_units']}")
    print(f"[C] Q_G-added direct accepts vs Q_old: {agg['qg_added']}")
    if agg["misaccept_details"]:
        print("[C] *** MISACCEPT DETECTED *** -- performance conclusions "
              "suspended; investigate first:")
        for key, m in agg["misaccept_details"]:
            print("      ", key, m)

    # ---- timing (after truth; fresh computation, no cache reads) ----
    warm_cands = all_candidates(SEEDS[0])
    time_screening(warm_cands[:1])
    time_full_reference(warm_cands[:1])
    time_end_to_end(warm_cands[:1])

    for seed in SEEDS:
        cands = all_candidates(seed)  # deterministic regeneration
        t_scr = [time_screening(cands) for _ in range(3)]
        t_ref = [time_full_reference(cands) for _ in range(3)]
        t_e2e, calls, acc = [], [], []
        for _ in range(3):
            dt, nc, na = time_end_to_end(cands)
            t_e2e.append(dt)
            calls.append(nc)
            acc.append(na)
        timing_rows.append(dict(
            seed=seed,
            t_screen_med=float(np.median(t_scr)),
            t_screen_min=min(t_scr), t_screen_max=max(t_scr),
            t_fullref_med=float(np.median(t_ref)),
            t_fullref_min=min(t_ref), t_fullref_max=max(t_ref),
            t_e2e_med=float(np.median(t_e2e)),
            t_e2e_min=min(t_e2e), t_e2e_max=max(t_e2e),
            fallback_calls_med=int(np.median(calls)),
            e2e_accepted_med=int(np.median(acc))))
        total_fallback += int(np.median(calls))
    utils.write_csv(
        os.path.join(utils.RESULTS_DIR, "timings.csv"),
        ["seed", "t_screen_med", "t_screen_min", "t_screen_max",
         "t_fullref_med", "t_fullref_min", "t_fullref_max",
         "t_e2e_med", "t_e2e_min", "t_e2e_max", "fallback_calls_med",
         "e2e_accepted_med"], timing_rows)

    med = lambda k: float(np.median([r[k] for r in timing_rows]))
    rng_ = lambda k: (min(r[k] for r in timing_rows),
                      max(r[k] for r in timing_rows))
    print(f"[C] timing (5 seeds, median [range], seconds): screen "
          f"{med('t_screen_med'):.4f} {rng_('t_screen_med')}, "
          f"full-ref {med('t_fullref_med'):.4f} {rng_('t_fullref_med')}, "
          f"screen+fallback {med('t_e2e_med'):.4f} {rng_('t_e2e_med')}")
    print(f"[C] user-level fallback reference calls (sum of per-seed "
          f"medians): {total_fallback} vs {n_user_problems} full-reference "
          f"user problems")
    return dict(rows=all_rows, agg=agg, timing_rows=timing_rows,
                total_fallback=total_fallback)


if __name__ == "__main__":
    run_c()
