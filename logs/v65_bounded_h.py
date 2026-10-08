# -*- coding: utf-8 -*-
"""V6.5 bounded-nominal-channel-norm experiments (numerical, bracketed).

Question: how much does the unrestricted six-scalar sharp envelope lose
when the admissible class is restricted by ||h_hat|| <= H_max?

Bracket:  Q_G(s) = B*(s)  <=  B_H*(s)  <=  best_found_feasible_value.
The upper estimate is a MULTI-START feasible search (random unitary
conjugations fixing v, random orthogonal h-perp directions with exact
I0 matching by solving a scalar quadratic), evaluated by the exact
trust-region solver.  It is NOT a proven restricted-class optimum.

Experiment A: boundary-approach family (Sec. 5) at several beta0 gaps x
H_max grid; watches B_H* descending toward Q_G as H_max grows.
Experiment B: pool instances (seeds 211-215, gamma=1, rho=0.15) with
H_max = c * ||h_actual||, c in {1, 1.25, 1.5, 2, 5}.

Writes results_v65/.  No historical result is modified.
"""
import csv
import json
import os
from datetime import datetime

import numpy as np

import core
import experiments_c as ec

OUT = "results_v65"
GAMMA, NOISE = 1.0, ec.NOISE_U
SEED_A = 65
SEED_B = 66
N_STARTS = 400


def rand_unitary_complement(rng, n):
    """Random unitary Q with Q e1 = e1 (identity on span{e1}): built from
    a random Haar unitary on the last n-1 coordinates."""
    import scipy.stats as st
    m = n - 1
    z = (rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m)))
    q, r = np.linalg.qr(z)
    d = np.diag(r).copy()
    q = q @ np.diag(d / np.abs(d))
    Q = np.eye(n, dtype=complex)
    Q[1:, 1:] = q
    return Q


def rand_output_unitary(rng, m):
    z = (rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m)))
    q, r = np.linalg.qr(z)
    d = np.diag(r).copy()
    return q @ np.diag(d / np.abs(d))


def solve_rho(Ce1, Ct, a, I0):
    """Find rho >= 0 with ||a*Ce1 + rho*Ct|| = I0 for unit direction t.
    Returns rho or None. ||a Ce1 + rho Ct||^2 = rho^2 w + 2 rho Re(g) + v0."""
    g = np.vdot(Ct, Ce1)
    w = float(np.vdot(Ct, Ct).real)
    v0 = float(np.vdot(Ce1, Ce1).real) * a * a
    # v0 + 2 a rho Re(g) + rho^2 w = I0^2
    A = w
    B = 2.0 * a * float(np.real(g))
    Cc = v0 - I0 * I0
    if A <= 1e-15:
        if abs(B) > 1e-15:
            rho = -Cc / B
            return rho if rho >= 0 else None
        return None
    disc = B * B - 4.0 * A * Cc
    if disc < 0:
        return None
    sq = np.sqrt(disc)
    for rho in ((-B + sq) / (2 * A), (-B - sq) / (2 * A)):
        if rho >= 0:
            return float(rho)
    return None


def sample_class_instance(rng, six, r, gamma, noise, H_max, n, m,
                           base_CH=None):
    """Random instance with EXACT six summaries and ||h|| <= H_max.
    Strategy: unitary-conjugate a base C (fixes v=e1 via complement
    rotation; output rotation free), then solve h on a random orthogonal
    direction with exact I0.  Returns (b, h, C, hnorm) or None."""
    p, a, I0, alpha, beta, beta0 = six
    if base_CH is None:
        # Lemma-1 style D-branch base
        delta = np.sqrt(max(beta0 ** 2 - alpha ** 2, 0.0))
        base_CH = np.zeros((m, n), complex)
        if alpha == 0 and beta > 0:
            base_CH[0, 1] = beta
        elif beta0 > alpha > 0 or (alpha > 0 and beta0 > alpha):
            base_CH[0, 0], base_CH[0, 1] = alpha, delta
            base_CH[1, 2] = beta if n > 2 and m > 1 else 0.0
        elif beta0 == alpha:
            base_CH[0, 0], base_CH[1, 1] = alpha, beta
    Qn = rand_unitary_complement(rng, n)
    Qm = rand_output_unitary(rng, m)
    CH = Qm @ base_CH @ Qn
    C = CH.conj().T
    Ce1 = CH[:, 0]
    b = np.zeros((n, 1), complex)
    b[0, 0] = np.sqrt(p)
    # random unit h_perp direction
    t = (rng.normal(size=n) + 1j * rng.normal(size=n))
    t[0] = 0.0
    t /= np.linalg.norm(t)
    Ct = CH @ t.reshape(-1, 1)
    rho = solve_rho(Ce1.reshape(-1, 1), Ct, a, I0)
    if rho is None:
        return None
    hn2 = a * a + rho * rho
    if hn2 > H_max * H_max + 1e-12:
        return None
    ph = rng.uniform(0, 2 * np.pi)
    h = a * np.exp(1j * ph) * np.eye(n, 1, dtype=complex) \
        + rho * np.exp(1j * rng.uniform(0, 2 * np.pi)) \
        * t.reshape(-1, 1)
    return b, h, C, float(np.sqrt(hn2))


def estimate_BH(six, r, gamma, noise, H_max, n, m, rng, n_starts,
                base_CH=None, include=None):
    """Upper estimate of B_H*(s) via multi-start feasible search.
    include: optional extra instance (b,h,C) to seed the search."""
    p, a, I0, alpha, beta, beta0 = six
    best = np.inf
    found = 0
    trials = []
    if include is not None:
        b, h, C = include
        if float(np.linalg.norm(h)) <= H_max + 1e-12:
            s = core.spectral_summary(b, h, C, gamma, noise)
            ok = max(abs(s[k] - v) for k, v in
                     zip(("p", "a", "I0", "alpha", "beta", "beta0"), six))
            if ok < 1e-9:
                ref = core.q_star_reference(s, r, gamma, noise,
                                            mu_rel_tol=1e-12)
                if ref["reference_valid"]:
                    best = min(best, ref["q"])
                    found += 1
    for _ in range(n_starts):
        got = sample_class_instance(rng, six, r, gamma, noise, H_max, n, m,
                                    base_CH=base_CH)
        if got is None:
            continue
        b, h, C, hn = got
        s = core.spectral_summary(b, h, C, gamma, noise)
        ok = max(abs(s[k] - v) for k, v in
                 zip(("p", "a", "I0", "alpha", "beta", "beta0"), six))
        if ok > 1e-9:
            continue
        ref = core.q_star_reference(s, r, gamma, noise, mu_rel_tol=1e-12)
        if not ref["reference_valid"]:
            continue
        best = min(best, ref["q"])
        found += 1
    return best, found


def experiment_a():
    """Boundary-approach family x H_max grid."""
    rows = []
    H_grid = (2.0, 5.0, 10.0, 20.0, 50.0, 100.0)
    for delta in (0.1, 0.05, 0.03, 0.01):
        CH = np.array([[1.0, delta, 0.0], [0.0, 0.0, 0.7]], complex)
        h = np.array([1.0, -2.0 / delta, 0.0]).reshape(-1, 1)
        b = np.array([[np.sqrt(2.0)], [0.0], [0.0]])
        s = core.spectral_summary(b, h, CH.conj().T, GAMMA, NOISE)
        six = (s["p"], s["a"], s["I0"], s["alpha"], s["beta"], s["beta0"])
        qg = core.Q_approx(dict(zip(("p", "a", "I0", "alpha", "beta",
                                     "beta0"), six)), 0.2, GAMMA, NOISE,
                           "G")[0]
        hn_fam = float(np.linalg.norm(h))
        for H in H_grid:
            rng = np.random.default_rng(SEED_A)
            best, found = estimate_BH(six, 0.2, GAMMA, NOISE, H, 3, 2, rng,
                                      N_STARTS, base_CH=CH,
                                      include=(b, h, CH.conj().T))
            gap = (best - qg) if np.isfinite(best) else np.nan
            rows.append(dict(experiment="A_boundary_family", delta=delta,
                             beta0=s["beta0"], h_family=hn_fam, H_max=H,
                             Q_G=qg, BH_upper=best, gap_upper=gap,
                             n_feasible=found,
                             family_member_admissible=bool(hn_fam <= H)))
            print(f"[A] δ={delta:.2f} H={H:5.1f}: Q_G={qg:.6f} "
                  f"BH_up={best if np.isfinite(best) else float('nan'):.6f} "
                  f"gap={gap if gap == gap else float('nan'):.3e} "
                  f"(feas={found})")
    return rows


def experiment_b():
    """Pool instances with H_max = c * ||h_actual||."""
    rows = []
    c_grid = (1.0, 1.25, 1.5, 2.0, 5.0)
    rng_master = np.random.default_rng(SEED_B)
    picked = 0
    for seed in (211, 212, 213):
        cands = ec.all_candidates(seed)
        # fixed rule declared before running: 3 users per seed
        # (users 0..2 of candidate 0), gamma=1, rho=0.15
        for u in range(3):
            cand = cands[0]
            b, h, C = cand["b"][u], cand["h"][u], cand["C"][u]
            n = b.shape[0]
            m = C.shape[1]
            s = core.spectral_summary(b, h, C, GAMMA, NOISE)
            six = (s["p"], s["a"], s["I0"], s["alpha"], s["beta"],
                   s["beta0"])
            r = 0.15 * float(np.linalg.norm(h))
            hn = float(np.linalg.norm(h))
            qg = core.Q_approx(dict(zip(("p", "a", "I0", "alpha", "beta",
                                         "beta0"), six)), r, GAMMA,
                                NOISE, "G")[0]
            ref = core.q_star_reference(s, r, GAMMA, NOISE,
                                        mu_rel_tol=1e-12)
            lq = core.L_quad(s, r, GAMMA)
            for c in c_grid:
                H = c * hn
                rng = np.random.default_rng(rng_master.integers(1 << 30))
                best, found = estimate_BH(six, r, GAMMA, NOISE, H, n, m,
                                          rng, 120, include=(b, h, C))
                gap = (best - qg) if np.isfinite(best) else np.nan
                rows.append(dict(experiment="B_pool", seed=seed, user=u,
                                 c=c, H_max=H, h_actual=hn, Q_G=qg,
                                 q_star_actual=ref["q"], L_quad=lq,
                                 BH_upper=best, gap_upper=gap,
                                 n_feasible=found))
            picked += 1
            print(f"[B] seed {seed} u{u}: h={hn:.2f} Q_G={qg:.4f} "
                  f"q*={ref['q']:.4f} Lq={lq:.4f}")
    return rows


def main():
    os.makedirs(OUT, exist_ok=True)
    cfg = dict(written_at=datetime.now().isoformat(timespec="seconds"),
               gamma=GAMMA, noise=NOISE,
               A=dict(family="Sec.5 boundary family", deltas=[0.1, 0.05,
                     0.03, 0.01], r=0.2, H_grid=[2, 5, 10, 20, 50, 100],
                     n_starts=N_STARTS, seed=SEED_A),
               B=dict(seeds=[211, 212, 213], rule="users 0-2 of candidate "
                     "0, rho=0.15", c_grid=[1, 1.25, 1.5, 2, 5],
                     n_starts=120, seed=SEED_B),
               semantics="BH_upper = best found feasible instance value "
                         "= UPPER estimate of B_H*(s); lower bound is "
                         "Q_G(s)=B*(s); no global optimality claimed")
    with open(os.path.join(OUT, "config.json"), "w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=1)
    rows_a = experiment_a()
    rows_b = experiment_b()
    with open(os.path.join(OUT, "bounded_h_details.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_a[0].keys())
                          if rows_a else ["experiment"])
        w.writeheader()
        w.writerows(rows_a)
    with open(os.path.join(OUT, "bounded_h_pool.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_b[0].keys()))
        w.writeheader()
        w.writerows(rows_b)
    # summary statistics for experiment B
    import numpy as np
    gaps = [r_["gap_upper"] for r_ in rows_b if r_["gap_upper"] == r_["gap_upper"]]
    near0 = sum(1 for g in gaps if g is not None and g < 1e-3)
    print(f"\n[B summary] restricted-gap upper estimates: n={len(gaps)}, "
          f"median={np.median(gaps):.4e}, p90={np.percentile(gaps,90):.4e}, "
          f"max={np.max(gaps):.4e}, near-zero(<1e-3)={near0}/{len(gaps)}")
    return True


if __name__ == "__main__":
    main()
