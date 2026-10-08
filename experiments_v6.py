# -*- coding: utf-8 -*-
"""V6 preparation experiments.

F-A: realistic CSI-error interpretation.
     (i)  rho-to-NMSE mapping table for the relative ball r = rho*||h_hat||;
     (ii) a path-loss/Rician mini-pool: 2 BS, 16-element RIS, distance-
          dependent path loss, Rician BS-RIS and RIS-UE links (Rayleigh
          BS-UE), bounded normalized CSI error, same candidate rules,
          gamma in {0.5,1,2}, rho in {0.06,0.15,0.30} = {-24.4,-16.5,
          -10.5} dB NMSE.  Illustrative heterogeneous scaling, NOT a
          calibrated link budget; all cells reported.

F-B: controlled scaling experiment (replaces the V5 scale table).
     Confounders fixed: 4 users (2 per cell), per-BS power 1 (per-stream
     0.5), gamma = 1, rho = 0.15 (relative), noise 0.05, identical
     candidate rules; only the antenna count varies, dim in {4,8,16,32}.
     Pipelines: tolerance-stopped TR pure, Q_G+fallback, L_quad+fallback,
     L_quad->Q_G+fallback.  Measures bound cost vs reference cost vs dim.

Writes results_v6/; earlier result directories untouched.
"""
import json
import os
from datetime import datetime
from time import perf_counter

import numpy as np

import core
import experiments_c as ec
import optimization_v41 as opt
import utils

OUT = "results_v6"
FA_SEEDS = (301, 302, 303, 304, 305)
FB_SEEDS = (311, 312, 313)
FB_DIMS = (2, 4, 8, 16)          # antennas per BS -> stacked dim 4..32
GAMMAS = (0.5, 1.0, 2.0)
RHOS = (0.06, 0.15, 0.30)
REF_TOL = 1e-12

FAILURES = []


def check(name, ok, detail=""):
    print(f"[v6:{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        FAILURES.append(f"{name}: {detail}")
    return ok


# ==========================================================================
# F-A part 1: rho <-> NMSE mapping (relative ball)
# ==========================================================================
def rho_nmse_table():
    rows = []
    for rho in RHOS:
        rows.append(dict(rho=rho, nmse_db=20 * np.log10(rho),
                         interpretation="worst-case relative error ball "
                                        "r = rho*||h_hat||; rho^2 equals "
                                        "the NMSE if the bound is set from "
                                        "E||e||^2/E||h||^2"))
    for db in (-20.0, -15.0, -10.0):
        rows.append(dict(rho=10 ** (db / 20.0), nmse_db=db,
                         interpretation="reverse mapping: dB target -> "
                                        "relative radius"))
    utils.write_csv(os.path.join(OUT, "fa_rho_nmse_map.csv"),
                    ["rho", "nmse_db", "interpretation"], rows)
    for r in rows[:3]:
        print(f"[v6:info] rho={r['rho']:.2f} <-> {r['nmse_db']:.1f} dB "
              f"NMSE")
    return rows


# ==========================================================================
# F-A part 2: path-loss / Rician candidate pool
# ==========================================================================
def _rician(rng, shape, k, scale):
    """sqrt(K/(K+1))*deterministic-ish LoS + sqrt(1/(K+1))*NLoS CN(0,1),
    times amplitude scale.  LoS phases are random per realization."""
    los = np.exp(1j * rng.uniform(0, 2 * np.pi, size=shape))
    nlos = (rng.normal(size=shape) + 1j * rng.normal(size=shape)) \
        * np.sqrt(0.5)
    if k <= 0:
        return scale * nlos
    return scale * (np.sqrt(k / (k + 1.0)) * los
                    + np.sqrt(1.0 / (k + 1.0)) * nlos)


def _db(x):
    return 10.0 ** (x / 20.0)


def fa_candidates(seed, n_ant=4, n_ris=16, users_per_cell=2):
    """Geometry: BS1(0,0), BS2(60,0), RIS(30,20); users uniform in a
    service area, assigned to the nearest BS (2 per cell by order).
    Path loss (illustrative): BS-RIS 12 dB fixed; BS-UE and RIS-UE
    distance-dependent, PL(d) = PL0 + 10*n*log10(d/d0) with PL0 = 10 dB at
    d0 = 10 m and n = 3.  Rician K = 3 on BS-RIS and RIS-UE, Rayleigh
    BS-UE.  Same MRT/RZF candidate rules as the main pool."""
    rng = np.random.default_rng(seed)
    n_users = 2 * users_per_cell
    pos_bs = [np.array([0.0, 0.0]), np.array([60.0, 0.0])]
    pos_ris = np.array([30.0, 20.0])
    users = []
    while len(users) < n_users:
        p = np.array([rng.uniform(5, 55), rng.uniform(-25, 25)])
        cell = 0 if np.linalg.norm(p - pos_bs[0]) \
            <= np.linalg.norm(p - pos_bs[1]) else 1
        users.append((p, cell))
    # balance: exactly users_per_cell per cell (swap farthest if needed)
    while sum(1 for _, c in users if c == 0) != users_per_cell:
        c0 = [i for i, (_, c) in enumerate(users) if c == 0]
        c1 = [i for i, (_, c) in enumerate(users) if c == 1]
        if len(c0) > users_per_cell:
            users[c0[0]] = (users[c0[0]][0], 1)
        else:
            users[c1[0]] = (users[c1[0]][0], 0)

    def pl(dist):
        return _db(10.0 + 30.0 * np.log10(max(dist, 1.0) / 10.0))

    G = [_rician(rng, (n_ris, n_ant), 3.0, _db(12.0))
         for _ in range(2)]
    f, d, dcell = [], [], []
    for u, (p, cell) in enumerate(users):
        f.append(_rician(rng, (n_ris,), 3.0, pl(np.linalg.norm(p - pos_ris))))
        dj = []
        for j in range(2):
            dj.append(_rician(rng, (n_ant,), 0.0,
                              pl(np.linalg.norm(p - pos_bs[j]))))
        d.append(dj)
        dcell.append(cell)
    thetas = [np.ones(n_ris, dtype=complex)]
    for _ in range(5):
        thetas.append(ec.PHASE_VALUES[rng.integers(0, 4, size=n_ris)])
    own = {j: tuple(u for u in range(n_users) if dcell[u] == j)
           for j in range(2)}
    per_stream = 1.0 / users_per_cell
    cands = []
    for theta in thetas:
        h = [[None] * n_users for _ in range(2)]
        for j in range(2):
            for u in range(n_users):
                base = d[u][j].reshape(-1, 1) \
                    + G[j].conj().T @ (theta * f[u]).reshape(-1, 1)
                h[j][u] = base if dcell[u] == j else 0.5 * base
        for btype in ("MRT", "RZF"):
            W = {}
            for j in range(2):
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
            beam = {}
            for u in range(n_users):
                col = np.zeros((2 * n_ant, 1), dtype=complex)
                col[dcell[u] * n_ant:(dcell[u] + 1) * n_ant, :] = W[u]
                beam[u] = col
            b_list, C_list, hs = [], [], []
            for u in range(n_users):
                hs.append(np.vstack([h[0][u], h[1][u]]))
                b_list.append(beam[u])
                C_list.append(np.hstack([beam[v] for v in range(n_users)
                                         if v != u]))
            cands.append(dict(seed=seed, btype=btype, h=hs, b=b_list,
                              C=C_list))
    return cands, [np.linalg.norm(p - pos_bs[c]) for p, c in users]


def run_fa(quick=False):
    utils.ensure_dirs()
    rows = []
    seeds = FA_SEEDS[:2] if quick else FA_SEEDS
    for seed in seeds:
        cands, dists = fa_candidates(seed)
        print(f"[v6:info] F-A seed {seed}: user BS distances "
              f"{[f'{x:.0f}m' for x in dists]}", flush=True)
        for gamma in GAMMAS:
            prep = opt.prepare(cands, gamma)
            for rho in RHOS:
                for ci in range(len(cands)):
                    dec = {}
                    feas = True
                    for u in range(len(cands[ci]["b"])):
                        pu = prep[ci][u]
                        r = rho * pu["hnorm"]
                        tau, _ = core.engineering_tau(pu["s"], r)
                        ref = core.q_star_reference(pu["s"], r, gamma,
                                                    ec.NOISE_U,
                                                    eig=pu["eig"],
                                                    mu_rel_tol=REF_TOL)
                        if not ref["reference_valid"]:
                            FAILURES.append(f"F-A invalid ref {seed}")
                        vals = {m: opt.bound_value(pu, r, gamma, m)
                                for m in opt.BOUND_STAGES}
                        st = ec.user_status(ref)
                        rows.append(dict(
                            seed=seed, gamma=gamma, rho=rho, cand_idx=ci,
                            user=u, r=r, tau=tau, q_star=ref["q"],
                            status=st,
                            h_norm=pu["hnorm"], **vals,
                            joint=max(vals["Q_G"], vals["L_quad"])))
                        if st != "clearly_feasible":
                            feas = False
                        if max(vals.values()) > ref["q"] + 1e-8:
                            FAILURES.append(f"F-A bound above q*: seed "
                                            f"{seed} u{u}")
                    for m in list(opt.BOUND_STAGES) + ["joint"]:
                        dec[m] = True
                    # unit acceptance computed below in aggregation
        print(f"[v6:info] F-A seed {seed} done", flush=True)
    utils.write_csv(os.path.join(OUT, "fa_realistic_details.csv"),
                    list(rows[0].keys()), rows)
    # aggregate
    units = {}
    for r_ in rows:
        units.setdefault((r_["seed"], r_["gamma"], r_["rho"], r_["cand_idx"]),
                         []).append(r_)
    agg = {}
    for key, us in units.items():
        feas = all(x["status"] == "clearly_feasible" for x in us)
        a = {"feasible": feas}
        for m in list(opt.BOUND_STAGES) + ["joint"]:
            a[m] = all(x[m] > x["tau"] for x in us)
        a["flip_QG"] = a["Q_G"] and not a["Q_old"]
        agg[key] = a
    summary = []
    for gamma in GAMMAS:
        for rho in RHOS:
            sel = [a for (s_, g_, r_, _), a in agg.items()
                   if g_ == gamma and r_ == rho]
            n = len(sel)
            if not n:
                continue
            summary.append(dict(
                gamma=gamma, rho=rho,
                nmse_db=round(20 * np.log10(rho), 1), units=n,
                feasible=sum(a["feasible"] for a in sel),
                D=sum(a["D"] for a in sel),
                Q_old=sum(a["Q_old"] for a in sel),
                Q_G=sum(a["Q_G"] for a in sel),
                L_quad=sum(a["L_quad"] for a in sel),
                joint=sum(a["joint"] for a in sel),
                QG_flips=sum(a["flip_QG"] for a in sel)))
    utils.write_csv(os.path.join(OUT, "fa_realistic_summary.csv"),
                    list(summary[0].keys()), summary)
    for s_ in summary:
        print(f"[v6:info] F-A gamma={s_['gamma']} rho={s_['rho']}"
              f"({s_['nmse_db']}dB): feas {s_['feasible']}/{s_['units']} | "
              f"old {s_['Q_old']} Q_G {s_['Q_G']} Lq {s_['L_quad']} "
              f"joint {s_['joint']} flips {s_['QG_flips']}")
    check("F-A: no bound violation", not any("F-A bound" in f
                                             for f in FAILURES))
    check("F-A: all references valid",
          not any("F-A invalid" in f for f in FAILURES))
    return rows, summary


# ==========================================================================
# F-B: controlled scaling
# ==========================================================================
def run_fb(quick=False):
    utils.ensure_dirs()
    rows = []
    seeds = FB_SEEDS[:1] if quick else FB_SEEDS
    for n_ant in FB_DIMS:
        dim = 2 * n_ant
        for seed in seeds:
            cands = opt.gen_candidates(seed, n_ant, users_per_cell=2)
            # per-user-problem cost measurement (fresh, interleaved)
            prep = opt.prepare(cands, 1.0)
            t_lb, t_ref = [], []
            for rep in range(3):
                t0 = perf_counter()
                for ci in range(len(cands)):
                    for u in range(4):
                        r = 0.15 * prep[ci][u]["hnorm"]
                        core.interval_lower_bound(prep[ci][u]["s"], r, 1.0,
                                                  ec.NOISE_U, 256, "G")
                t_lb.append((perf_counter() - t0) / (len(cands) * 4))
                t0 = perf_counter()
                for ci in range(len(cands)):
                    for u in range(4):
                        r = 0.15 * prep[ci][u]["hnorm"]
                        core.q_star_reference(prep[ci][u]["s"], r, 1.0,
                                              ec.NOISE_U,
                                              eig=prep[ci][u]["eig"],
                                              mu_rel_tol=REF_TOL)
                t_ref.append((perf_counter() - t0) / (len(cands) * 4))
            # pipelines on the same decision task
            pipes = ("TR_pure", "Q_G_ref", "L_quad_ref", "Lq_QG_ref")
            for p in pipes:
                opt_pre = opt.prepare(cands[:1], 1.0)
                if p == "TR_pure":
                    opt.cascade(cands[:1], opt_pre, 1.0, ())
                else:
                    opt.cascade(cands[:1], opt_pre, 1.0,
                                opt.PIPELINES[p])          # warm-up
            runs = {p: [] for p in pipes}
            for rep in range(3):
                for p in np.roll(np.array(pipes), rep):
                    t0 = perf_counter()
                    pr = opt.prepare(cands, 1.0)
                    if p == "TR_pure":
                        out = opt.cascade(cands, pr, 1.0, ())
                    else:
                        out = opt.cascade(cands, pr, 1.0, opt.PIPELINES[p])
                    runs[p].append((perf_counter() - t0, out))
            base = runs["TR_pure"][0][1]["decisions"]
            for p in pipes:
                for _, out in runs[p]:
                    if out["decisions"] != base:
                        FAILURES.append(f"F-B {dim}d {p}: decision mismatch")
                d = runs[p][0][1]["decisions"]
                rows.append(dict(
                    dim=dim, n_ant=n_ant, seed=seed, pipeline=p,
                    units=len(d),
                    accepted=sum(1 for x in d if x == "accepted"),
                    seconds_med=float(np.median([t for t, _ in runs[p]])),
                    ref_calls=runs[p][0][1]["stats"].get(
                        "ref_calls",
                        runs[p][0][1]["stats"].get("sdp_calls", 0)),
                    lb_ms_per_user=float(np.median(t_lb) * 1000),
                    ref_ms_per_user=float(np.median(t_ref) * 1000)))
        med_lb = np.median([r_["lb_ms_per_user"] for r_ in rows
                            if r_["dim"] == dim])
        med_ref = np.median([r_["ref_ms_per_user"] for r_ in rows
                             if r_["dim"] == dim])
        print(f"[v6:info] F-B dim={dim}: LB {med_lb:.3f} ms/user-problem, "
              f"TR ref {med_ref:.3f} ms/user-problem", flush=True)
    utils.write_csv(os.path.join(OUT, "fb_scaling.csv"),
                    list(rows[0].keys()), rows)
    check("F-B: pipeline decisions match reference everywhere",
          not any("F-B" in f and "mismatch" in f for f in FAILURES))
    return rows


def run_all(quick=False):
    utils.ensure_dirs()
    print("--- V6 F-A: rho/NMSE mapping ---")
    rho_nmse_table()
    print("--- V6 F-A: realistic pool ---")
    run_fa(quick)
    print("--- V6 F-B: controlled scaling ---")
    run_fb(quick)
    utils.write_json(os.path.join(OUT, "run_manifest.json"),
                     dict(mode="quick" if quick else "full",
                          finished_at=datetime.now().isoformat(
                              timespec="seconds"),
                          failures=list(FAILURES),
                          fa_seeds=list(FA_SEEDS), fb_seeds=list(FB_SEEDS),
                          confounders="F-B: 4 users, per-BS power 1 "
                                      "(per-stream 0.5), gamma=1, rho=0.15, "
                                      "noise 0.05, identical candidate "
                                      "rules; only antennas vary"))
    print(f"\n=== V6 experiments {'FAILED' if FAILURES else 'PASSED'} ===")
    for f in FAILURES:
        print("  FAIL:", f)
    return not FAILURES


if __name__ == "__main__":
    import sys
    sys.exit(0 if run_all("--quick" in sys.argv) else 1)
