# -*- coding: utf-8 -*-
"""Experiment A: closed-form / construction / boundary checks.

A1  directional extremum: closed-form R vs direct scalar maximization
A2  three-branch direction construction (input dim 5, output dim 3)
A3  deterministic certification flip + Q_G construction instance
A4  boundary-mechanism short check (beta0 -> alpha)
"""
import os

import numpy as np
from scipy.optimize import minimize_scalar

import core
import utils

RNG_SEED_A1 = 20261006


# --------------------------------------------------------------------------
# A1 helpers
# --------------------------------------------------------------------------
def scalar_extremum_ref(alpha, beta, beta0, x, y, n_grid=4097):
    """max_{0<=s<=beta} {alpha^2 x^2 + s^2 y^2 + 2xy min(alpha s, delta
    sqrt(beta0^2-s^2))}: 4097-point grid, bounded scalar refinement around
    every grid local peak, plus the special points s=0, delta, beta."""
    delta = np.sqrt(max(beta0 ** 2 - alpha ** 2, 0.0))

    def f(s):
        s = np.asarray(s, dtype=float)
        return (alpha ** 2 * x ** 2 + s ** 2 * y ** 2
                + 2.0 * x * y * np.minimum(alpha * s,
                                           delta * np.sqrt(np.maximum(beta0 ** 2 - s ** 2, 0.0))))

    s_grid = np.linspace(0.0, beta, n_grid)
    v = f(s_grid)
    grid_max = float(v.max())
    h = beta / (n_grid - 1)
    cand = {0.0, beta, min(delta, beta), float(s_grid[int(np.argmax(v))])}
    interior = np.where((v[1:-1] >= v[:-2]) & (v[1:-1] >= v[2:]))[0] + 1
    for i in interior:
        cand.add(float(s_grid[i]))
    best = max(grid_max, *(float(f(c)) for c in cand))
    for c in cand:
        lo, hi = max(0.0, c - 2 * h), min(beta, c + 2 * h)
        if hi <= lo:
            continue
        res = minimize_scalar(lambda s: -float(f(s)), bounds=(lo, hi),
                              method="bounded",
                              options=dict(xatol=1e-13 * max(1.0, beta)))
        best = max(best, -float(res.fun))
    return np.sqrt(best), np.sqrt(grid_max)


def run_a1():
    rng = np.random.default_rng(RNG_SEED_A1)
    rows = []

    # --- 30 random compatible parameter sets
    for i in range(30):
        alpha = float(rng.uniform(0.2, 1.5))
        beta = float(rng.uniform(0.2, 1.5))
        lo = max(alpha, beta)
        hi = np.sqrt(alpha ** 2 + beta ** 2)
        beta0 = float(rng.uniform(lo, hi))
        x = float(rng.uniform(0.05, 0.6))
        y = float(rng.uniform(0.05, 0.6))
        R = core.R_envelope(x, y, alpha, beta, beta0)
        ref_max, grid_max = scalar_extremum_ref(alpha, beta, beta0, x, y)
        rows.append(dict(kind="random", idx=i, alpha=alpha, beta=beta,
                         beta0=beta0, x=x, y=y, R_closed=R, ref_max=ref_max,
                         grid_max=grid_max,
                         abs_err_R2=abs(R ** 2 - ref_max ** 2),
                         grid_under=ref_max - grid_max, note=""))

    # --- fixed three-branch examples
    for j, (x, y) in enumerate([(0.2, 0.4), (0.4, 0.3), (0.45, 0.1)]):
        alpha = beta = 1.0
        beta0 = 1.1
        delta = np.sqrt(beta0 ** 2 - alpha ** 2)
        nu = np.sqrt(beta0 ** 2 - beta ** 2)
        if delta * x <= nu * y:
            br = 1
        elif delta * x >= alpha * y:
            br = 3
        else:
            br = 2
        R = core.R_envelope(x, y, alpha, beta, beta0)
        ref_max, grid_max = scalar_extremum_ref(alpha, beta, beta0, x, y)
        rows.append(dict(kind="fixed_branch", idx=j, alpha=alpha, beta=beta,
                         beta0=beta0, x=x, y=y, R_closed=R, ref_max=ref_max,
                         grid_max=grid_max, abs_err_R2=abs(R ** 2 - ref_max ** 2),
                         grid_under=ref_max - grid_max, note=f"branch={br}"))

    # --- branch continuity at the two internal boundaries
    alpha = beta = 1.0
    beta0 = 1.1
    delta = np.sqrt(beta0 ** 2 - alpha ** 2)
    nu = np.sqrt(beta0 ** 2 - beta ** 2)
    for t in [0.15, 0.3, 0.5]:
        # branch1/2 boundary: delta x = nu y  ->  (x, y) = (nu t, delta t)
        xb, yb = nu * t, delta * t
        vals = [core.R_envelope(xb, yb, alpha, beta, beta0),
                core.R_envelope(xb * 1.000001, yb, alpha, beta, beta0),
                core.R_envelope(xb * 0.999999, yb, alpha, beta, beta0),
                core.R_envelope(xb, yb * 1.000001, alpha, beta, beta0),
                core.R_envelope(xb, yb * 0.999999, alpha, beta, beta0)]
        jump = float(np.max(vals) - np.min(vals))
        rows.append(dict(kind="continuity_12", idx=int(t * 100), alpha=alpha,
                         beta=beta, beta0=beta0, x=xb, y=yb,
                         R_closed=vals[0], ref_max=np.nan, grid_max=np.nan,
                         abs_err_R2=np.nan, grid_under=np.nan,
                         note=f"max_jump={jump:.3e}"))
        # branch2/3 boundary: delta x = alpha y  ->  (x, y) = (alpha t, delta t)
        xb, yb = alpha * t, delta * t
        vals = [core.R_envelope(xb, yb, alpha, beta, beta0),
                core.R_envelope(xb * 1.000001, yb, alpha, beta, beta0),
                core.R_envelope(xb * 0.999999, yb, alpha, beta, beta0),
                core.R_envelope(xb, yb * 1.000001, alpha, beta, beta0),
                core.R_envelope(xb, yb * 0.999999, alpha, beta, beta0)]
        jump = float(np.max(vals) - np.min(vals))
        rows.append(dict(kind="continuity_23", idx=int(t * 100), alpha=alpha,
                         beta=beta, beta0=beta0, x=xb, y=yb,
                         R_closed=vals[0], ref_max=np.nan, grid_max=np.nan,
                         abs_err_R2=np.nan, grid_under=np.nan,
                         note=f"max_jump={jump:.3e}"))

    # --- global inequality R <= min(beta0*norm, alpha x + beta y)
    worst = 0.0
    for _ in range(1000):
        x = float(rng.uniform(0.001, 0.6))
        y = float(rng.uniform(0.001, 0.6))
        alpha = float(rng.uniform(0.2, 1.5))
        beta = float(rng.uniform(0.2, 1.5))
        beta0 = float(rng.uniform(max(alpha, beta),
                                  np.sqrt(alpha ** 2 + beta ** 2)))
        R = core.R_envelope(x, y, alpha, beta, beta0)
        cap = min(beta0 * np.hypot(x, y), alpha * x + beta * y)
        worst = max(worst, float(R - cap))
    rows.append(dict(kind="inequality", idx=0, alpha=np.nan, beta=np.nan,
                     beta0=np.nan, x=np.nan, y=np.nan, R_closed=np.nan,
                     ref_max=np.nan, grid_max=np.nan, abs_err_R2=np.nan,
                     grid_under=np.nan,
                     note=f"max(R-min(cap,linear)) over 1000 draws = {worst:.3e}"))

    # --- hand endpoints
    endpoints = []
    cases = [
        ("alpha=0", dict(alpha=0.0, beta=0.9, beta0=0.9, x=0.3, y=0.5), 0.9 * 0.5),
        ("beta=0", dict(alpha=1.2, beta=0.0, beta0=1.2, x=0.4, y=0.2), 1.2 * 0.4),
        ("x=0", dict(alpha=1.0, beta=0.8, beta0=1.1, x=0.0, y=0.5), 0.8 * 0.5),
        ("y=0", dict(alpha=1.0, beta=0.8, beta0=1.1, x=0.4, y=0.0), 1.0 * 0.4),
    ]
    for j, (name, kw, expected) in enumerate(cases):
        R = core.R_envelope(kw["x"], kw["y"], kw["alpha"], kw["beta"], kw["beta0"])
        ref_max, grid_max = scalar_extremum_ref(kw["alpha"], kw["beta"],
                                                kw["beta0"], kw["x"], kw["y"])
        endpoints.append(dict(kind="endpoint", idx=j, **kw, R_closed=R,
                              ref_max=ref_max, grid_max=grid_max,
                              abs_err_R2=abs(R ** 2 - expected ** 2),
                              grid_under=ref_max - grid_max,
                              note=f"{name}: expected R={expected:.6f}"))

    rows.extend(endpoints)
    fields = ["kind", "idx", "alpha", "beta", "beta0", "x", "y", "R_closed",
              "ref_max", "grid_max", "abs_err_R2", "grid_under", "note"]
    utils.write_csv(os.path.join(utils.RESULTS_DIR, "direction_checks.csv"),
                    fields, rows)

    rand_rows = [r for r in rows if r["kind"] == "random"]
    max_err = max(r["abs_err_R2"] for r in rand_rows)
    max_grid = max(r["grid_under"] for r in rand_rows)
    print(f"[A1] 30 random sets: max |R^2 - ref^2| = {max_err:.3e}, "
          f"max pure-grid underestimation = {max_grid:.3e}")
    print(f"[A1] inequality worst violation: {worst:.3e}")
    return rows


# --------------------------------------------------------------------------
# A2
# --------------------------------------------------------------------------
def run_a2():
    alpha = beta = 1.0
    beta0 = 1.1
    rows = []
    for (x, y) in [(0.2, 0.4), (0.4, 0.3), (0.45, 0.1)]:
        for eps in [1e-2, 1e-3, 1e-4]:
            con = core.build_direction_construction(alpha, beta, beta0, x, y, eps)
            C_H = con["C_H"]
            n = C_H.shape[1]
            v = np.zeros((n, 1), dtype=complex)
            v[0, 0] = 1.0
            P = np.eye(n, dtype=complex) - v @ v.conj().T
            a_chk = float(np.linalg.norm(C_H @ v))
            b0_chk = float(np.sqrt(np.linalg.eigvalsh(C_H @ C_H.conj().T)[-1]))
            b_chk = float(np.sqrt(np.linalg.eigvalsh(C_H @ P @ C_H.conj().T)[-1]))
            e_test = np.zeros((n, 1), dtype=complex)
            e_test[0, 0], e_test[1, 0] = x, y
            dir_eps = float(np.linalg.norm(C_H @ e_test))
            rows.append(dict(
                x=x, y=y, eps=eps, t_star=con["t"],
                bracket_ok=con["bracket_ok"],
                match_alpha=a_chk - alpha, match_beta=b_chk - beta,
                match_beta0=b0_chk - beta0,
                exact_dir_err=con["R_val"] - core.R_envelope(x, y, alpha, beta, beta0),
                lam_min_M=con["lam_min_M"],
                eps_dir_deficit=core.R_envelope(x, y, alpha, beta, beta0) - dir_eps,
                f0=con["f0"], f1=con["f1"]))
    fields = ["x", "y", "eps", "t_star", "bracket_ok", "match_alpha",
              "match_beta", "match_beta0", "exact_dir_err", "lam_min_M",
              "eps_dir_deficit", "f0", "f1"]
    utils.write_csv(os.path.join(utils.RESULTS_DIR, "construction_checks.csv"),
                    fields, rows)
    print(f"[A2] {len(rows)} constructions; worst |alpha match| = "
          f"{max(abs(r['match_alpha']) for r in rows):.3e}, "
          f"|beta match| = {max(abs(r['match_beta']) for r in rows):.3e}, "
          f"|beta0 match| = {max(abs(r['match_beta0']) for r in rows):.3e}, "
          f"min lam_min(T-ll^H) = {min(r['lam_min_M'] for r in rows):.3e}")
    return rows


# --------------------------------------------------------------------------
# A3
# --------------------------------------------------------------------------
A3_PARAMS = dict(p=4.0, a=1.3, I0=0.5, alpha=1.0, beta=1.0, beta0=1.1,
                 gamma=1.0, noise_power=1.5, r=0.5)


def run_a3():
    P = A3_PARAMS
    scal = dict(p=P["p"], a=P["a"], I0=P["I0"], alpha=P["alpha"],
                beta=P["beta"], beta0=P["beta0"])
    qg_val, x_star, y_star = core.Q_approx(scal, P["r"], P["gamma"], P["noise_power"], "G")
    print(f"[A3] Q_G 1-D optimum: x* = {x_star:.8f}, y* = {y_star:.8f}, "
          f"Q_G = {qg_val:.8f}")
    out = {"params": P, "x_star": x_star, "y_star": y_star,
           "Q_G_approx": qg_val, "construction": []}
    for eps in [1e-2, 1e-3, 1e-4]:
        inst = core.build_A3_instance(P["p"], P["a"], P["I0"], P["alpha"],
                                      P["beta"], P["beta0"], x_star, y_star,
                                      eps, P["gamma"], P["noise_power"])
        s = core.spectral_summary(inst["b"], inst["h"], inst["C"],
                                  P["gamma"], P["noise_power"])
        ref = core.q_star_reference(s, P["r"], P["gamma"], P["noise_power"])
        n = 5
        # worst-case direction: C^H e = -(x* d + y* ell_e) aligns the
        # interference with C^H h (norm I0 + R_eps) while the signal drops
        # to p (a - x*)^2
        e_test = np.zeros((n, 1), dtype=complex)
        e_test[0, 0], e_test[1, 0] = -x_star, -y_star
        obj_e = core.objective_at(inst["b"], inst["h"], inst["C"], P["gamma"],
                                  P["noise_power"], e_test)
        rec = dict(
            eps=eps,
            p_err=s["p"] - P["p"], a_err=s["a"] - P["a"],
            I0_err=s["I0"] - P["I0"], alpha_err=s["alpha"] - P["alpha"],
            beta_err=s["beta"] - P["beta"], beta0_err=s["beta0"] - P["beta0"],
            q_star=ref["q"], Q_G_approx=qg_val, obj_at_e=obj_e,
            q_star_minus_QG=ref["q"] - qg_val,
            obj_minus_QG=obj_e - qg_val,
            reference_valid=ref["reference_valid"], mode=ref["mode"])
        out["construction"].append(rec)
        print(f"[A3] eps={eps:.0e}: six-scalar max err="
              f"{max(abs(v) for k, v in rec.items() if k.endswith('_err')):.3e}, "
              f"q*={ref['q']:.8f}, obj(e*)={obj_e:.8f}, Q_G={qg_val:.8f}")

    # deterministic flip example
    CH = np.array([[1, np.sqrt(0.21), 0, 0, 0],
                   [0, 0, 1, 0, 0],
                   [0, 0, 0, 0, 0]], dtype=complex)
    h = np.array([1.3, -1.8 / np.sqrt(0.21), 0, 0, 0],
                 dtype=complex).reshape(-1, 1)
    b = np.array([2, 0, 0, 0, 0], dtype=complex).reshape(-1, 1)
    gam, npw, r = 1.0, 1.5, 0.5
    s = core.spectral_summary(b, h, CH.conj().T, gam, npw)
    D = core.D_bound(s, r, gam, npw)
    q_old, _, _ = core.Q_approx(s, r, gam, npw, "old")
    qg, _, _ = core.Q_approx(s, r, gam, npw, "G")
    q_lb = core.interval_lower_bound(s, r, gam, npw, 256, "old")
    qg_lb = core.interval_lower_bound(s, r, gam, npw, 256, "G")
    qg_lb_fine = core.interval_lower_bound(s, r, gam, npw, 65536, "G")
    ref = core.q_star_reference(s, r, gam, npw)
    flip = dict(D=D, Q_old=q_old, Q_G=qg, Q_LB_256=q_lb, Q_G_LB_256=qg_lb,
                Q_G_LB_65536=qg_lb_fine, q_star=ref["q"],
                reference_valid=ref["reference_valid"], mode=ref["mode"],
                tau=ref["tau"])
    out["flip"] = flip
    print("[A3] flip: D={D:.6f} Q_old={Q_old:.6f} Q_G={Q_G:.6f} "
          "Q_LB={Q_LB_256:.6f} Q_G_LB={Q_G_LB_256:.6f} q*={q_star:.6f} "
          "valid={reference_valid}".format(**flip))
    utils.write_json(os.path.join(utils.RESULTS_DIR, "a3_flip.json"), out)
    return out


# --------------------------------------------------------------------------
# A4
# --------------------------------------------------------------------------
def run_a4():
    gam, npw, r = 1.0, 0.05, 0.2
    rows = []
    # beta0 = alpha boundary case (six scalars only)
    scal = dict(p=2.0, a=1.0, I0=1.0, alpha=1.0, beta=0.7, beta0=1.0)
    qg, _, _ = core.Q_approx(scal, r, gam, npw, "G")
    qg_lb = core.interval_lower_bound(scal, r, gam, npw, 256, "G")
    qpartial, _ = core.Q_boundary_approx(scal, r, gam, npw)
    CH_boundary = np.array([[1., 0., 0.], [0., 0., 0.7]], dtype=complex)
    h_boundary = np.array([[1.], [0.], [0.]], dtype=complex)
    b_boundary = np.array([[np.sqrt(2.)], [0.], [0.]], dtype=complex)
    sb = core.spectral_summary(b_boundary, h_boundary, CH_boundary.conj().T, gam, npw)
    rb = core.q_star_reference(sb, r, gam, npw)
    rows.append(dict(delta=0.0, beta0_target=1.0, beta0_measured=1.0,
                     h_norm=1.0, Q_G=qg, Q_G_LB=qg_lb, Q_partial=qpartial,
                     obj_at_e=core.objective_at(b_boundary, h_boundary, CH_boundary.conj().T,
                                              gam, npw, np.array([[-r], [0.], [0.]])),
                     q_star=rb['q'], note="separate boundary Q_partial; general Q_G retained"))
    for delta in [0.1, 0.03, 0.01]:
        CH = np.array([[1.0, delta, 0.0], [0.0, 0.0, 0.7]], dtype=complex)
        h = np.array([1.0, -2.0 / delta, 0.0], dtype=complex).reshape(-1, 1)
        b = np.array([np.sqrt(2.0), 0.0, 0.0], dtype=complex).reshape(-1, 1)
        s = core.spectral_summary(b, h, CH.conj().T, gam, npw)
        qg, _, _ = core.Q_approx(s, r, gam, npw, "G")
        qg_lb = core.interval_lower_bound(s, r, gam, npw, 256, "G")
        e_eval = np.array([-0.2, 0.0, 0.0], dtype=complex).reshape(-1, 1)
        obj = core.objective_at(b, h, CH.conj().T, gam, npw, e_eval)
        ref = core.q_star_reference(s, r, gam, npw)
        rows.append(dict(delta=delta, beta0_target=np.sqrt(1 + delta ** 2),
                         beta0_measured=s["beta0"],
                         h_norm=float(np.linalg.norm(h)), Q_G=qg,
                         Q_G_LB=qg_lb, obj_at_e=obj, q_star=ref["q"],
                         note="general branch"))
    fields = ["delta", "beta0_target", "beta0_measured", "h_norm", "Q_G",
              "Q_G_LB", "Q_partial", "obj_at_e", "q_star", "note"]
    utils.write_csv(os.path.join(utils.RESULTS_DIR, "boundary_checks.csv"),
                    fields, rows)
    for row in rows:
        hnorm = row["h_norm"]
        hnorm_s = "n/a" if hnorm != hnorm else f"{hnorm:.2f}"
        print(f"[A4] delta={row['delta']:.2f} beta0={row['beta0_measured']:.6f} "
              f"Q_G={row['Q_G']:.6f} Q_partial={row.get('Q_partial', 'n/a')} "
              f"|h|={hnorm_s} obj(e)={row['obj_at_e']}")
    return rows


def run_all():
    utils.ensure_dirs()
    return dict(a1=run_a1(), a2=run_a2(), a3=run_a3(), a4=run_a4())


if __name__ == "__main__":
    run_all()
