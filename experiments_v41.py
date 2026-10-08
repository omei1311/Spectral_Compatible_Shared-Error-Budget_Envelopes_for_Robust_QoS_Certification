# -*- coding: utf-8 -*-
"""V4.1 formal verification suite.

Integrates the V4.1 analytical constructions into the experiment flow and
adds the pre-registered validation pool, resolution ablation and paired
timing.  All outputs go to results_v41/; results/ and results_review/ are
never touched.

Sections:
  theory   : Theorem 2 construction battery (output-space PSD comparison +
             planar 2x2 spectral completion), Theorem 3 boundary battery
             (four corrected instances + random boundary + nondegenerate
             family), degenerate/endpoints battery
  scan     : experiment B under the V4.1 diagnostics (violation counts,
             max deviations)
  pools    : baseline easy pool + exploration grid + pre-registered
             validation seeds, overlap/dedup accounting
  resolution : J = 256/1024/4096 ablation + theoretical-vs-numerical loss
  timing   : paired interleaved same-decision pipelines; tolerance-stopped
             trust-region variant validation and comparison
"""
import hashlib
import json
import os
from datetime import datetime

import numpy as np

import audit_and_stress as az
import core
import experiments_b
import experiments_c as ec
import utils

GAMMAS = (0.5, 1.0, 2.0)
RHOS = (0.06, 0.15, 0.30)
EXPLORATION_SEEDS = (11, 22, 33, 44, 55, 66, 77, 88, 99, 110)
VALIDATION_SEEDS = (211, 212, 213, 214, 215)
THEORY_SEED = 202610061

FAILURES = []  # collected across sections; non-empty => nonzero exit


def check(name, ok, detail=""):
    status = "PASS" if ok else "FAIL"
    print(f"[v41:{status}] {name} {detail}")
    if not ok:
        FAILURES.append(f"{name}: {detail}")
    return ok


# ==========================================================================
# 1. theory: Theorem 2 construction battery
# ==========================================================================
def construction_row(tag, p, a, I0, alpha, beta, beta0, r, gamma,
                     noise_power, x, y, eps):
    scal = dict(p=p, a=a, I0=I0, alpha=alpha, beta=beta, beta0=beta0)
    qg, _, _ = core.Q_approx(scal, r, gamma, noise_power, "G")
    try:
        con = core.construct_v41(alpha, beta, beta0, x, y, eps, p, a, I0)
    except Exception as exc:  # bracket/PD failures must be reported
        return dict(tag=tag, p=p, a=a, I0=I0, alpha=alpha, beta=beta,
                    beta0=beta0, r=r, x=x, y=y, eps=eps,
                    construct_ok=False, error=str(exc), Q_G=qg)
    summ = core.spectral_summary(con["b"], con["h"], con["C"], gamma,
                                 noise_power)
    six = max(abs(summ[k] - scal[k]) for k in scal)
    obj = core.objective_at(con["b"], con["h"], con["C"], gamma, noise_power,
                            con["e"])
    R = core.R_envelope(x, y, alpha, beta, beta0)
    bound22 = gamma * (2 * I0 + 2 * R) * eps * y * beta
    # envelope objective at the same (x, y): the pure epsilon-loss identity
    # q(e_eps) - envelope_obj(x, y) = gamma (R - R_eps)(2 I0 + R + R_eps)
    # holds for ANY (x, y); the paper's bound (22) against Q_G additionally
    # requires (x, y) to minimize (9).
    obj_at_xy = (p * (a - x) ** 2 - gamma * (I0 + R) ** 2
                 - gamma * noise_power)
    at_min = bool(abs(obj_at_xy - qg) <= 1e-6 * max(1.0, abs(qg)))
    ref = core.q_star_reference(summ, r, gamma, noise_power)
    return dict(
        tag=tag, p=p, a=a, I0=I0, alpha=alpha, beta=beta, beta0=beta0, r=r,
        x=x, y=y, eps=eps, construct_ok=True,
        t_eps=con["t_eps"], k=con["k"],
        S0_minus_b0sq=con["S0_minus_b0sq"], S1_minus_b0sq=con["S1_minus_b0sq"],
        six_scalar_err=six,
        spec_err_alpha=summ["alpha"] - alpha,
        spec_err_beta=summ["beta"] - beta,
        spec_err_beta0=summ["beta0"] - beta0,
        lam_min_M=con["lam_min_M"],
        e_x=-x, e_y=-y, e_norm=float(np.linalg.norm(con["e"])),
        R=R, R_eps=con["R_eps"], R_minus_R_eps=R - con["R_eps"],
        eps_y_ell=eps * y * float(np.linalg.norm(con["ell"])),
        obj_at_e=obj, Q_G=qg, obj_minus_QG=obj - qg,
        obj_at_xy_envelope=obj_at_xy, at_minimizer=at_min,
        eps_loss=obj - obj_at_xy,
        bound22=bound22,
        eps_loss_bound_holds=bool(obj - obj_at_xy <= bound22 + 1e-9),
        bound22_holds=bool(not at_min
                           or obj - qg <= bound22 + 1e-9),
        q_star=ref["q"], reference_valid=ref["reference_valid"],
        ref_mode=ref["mode"],
        q_star_minus_obj=ref["q"] - obj,
        q_star_ge_QG=bool(ref["q"] >= qg - 1e-9))


def run_theory(quick=False):
    utils.ensure_dirs()
    rng = np.random.default_rng(THEORY_SEED)
    rows = []

    # (a) flip-family branch points at the Q_G minimizer
    scal0 = dict(p=4.0, a=1.3, I0=0.5, alpha=1.0, beta=1.0, beta0=1.1)
    qg0, x0, y0 = core.Q_approx(scal0, 0.5, 1.0, 1.5, "G")
    for eps in (1e-2, 1e-3, 1e-4):
        rows.append(construction_row("flip_family", 4.0, 1.3, 0.5, 1.0, 1.0,
                                     1.1, 0.5, 1.0, 1.5, x0, y0, eps))

    # (b) random tuples incl. asymmetric spectra, beta=beta0, upper endpoint
    n_rand = 6 if quick else 40
    tuples = []
    for _ in range(n_rand):
        al, be = rng.uniform(0.2, 1.8, 2)
        b0 = rng.uniform(max(al, be) + 1e-4, np.hypot(al, be))
        tuples.append(("random", al, be, b0))
    tuples.append(("beta_eq_beta0", 1.0, 1.2, 1.2))
    tuples.append(("beta_eq_beta0", 0.7, 1.5, 1.5))
    tuples.append(("upper_endpoint", 1.0, 0.7, np.hypot(1.0, 0.7)))
    tuples.append(("upper_endpoint", 1.3, 0.6, np.hypot(1.3, 0.6)))
    for tag, al, be, b0 in tuples:
        p = float(rng.uniform(0.5, 3.0))
        a = float(rng.uniform(0.0, 2.0))
        I0 = float(rng.uniform(0.0, 2.0))
        r, gam, noise = 0.4, 1.0, 0.05
        sc = dict(p=p, a=a, I0=I0, alpha=al, beta=be, beta0=b0)
        _, x, y = core.Q_approx(sc, r, gam, noise, "G")
        for eps in (1e-2, 1e-3):
            rows.append(construction_row(tag, p, a, I0, al, be, b0, r, gam,
                                         noise, x, y, eps))

    # (c) direction endpoints and I0 = 0 (compatible by construction)
    for (x, y) in [(0.4, 0.0), (0.0, 0.4)]:
        for eps in (1e-2, 1e-3):
            rows.append(construction_row("direction_endpoint", 2.0, 1.0, 0.0,
                                         1.0, 1.0, 1.1, 0.4, 1.0, 0.05, x, y,
                                         eps))
    for eps in (1e-2, 1e-3):
        rows.append(construction_row("I0_zero", 2.0, 1.0, 0.0, 1.0, 1.0, 1.1,
                                     0.4, 1.0, 0.05, 0.25,
                                     np.sqrt(0.4 ** 2 - 0.25 ** 2), eps))

    fields = list(rows[0].keys())
    utils.write_csv(os.path.join(utils.RESULTS_DIR, "theory_constructions.csv"),
                    fields, rows)
    ok_rows = [r for r in rows if r.get("construct_ok")]
    max_six = max(r["six_scalar_err"] for r in ok_rows)
    max_spec = max(max(abs(r["spec_err_alpha"]), abs(r["spec_err_beta"]),
                       abs(r["spec_err_beta0"])) for r in ok_rows)
    min_pd = min(r["lam_min_M"] for r in ok_rows)
    b22 = all(r["bound22_holds"] for r in ok_rows)
    eps_loss_ok = all(r["eps_loss_bound_holds"] for r in ok_rows)
    obj_ge = all(r["obj_minus_QG"] >= -1e-9 for r in ok_rows)
    safe = all(r["q_star_ge_QG"] for r in ok_rows)
    valid = all(r["reference_valid"] for r in ok_rows)
    conv = True
    for tag in ("flip_family",):
        tag_rows = [r for r in ok_rows if r["tag"] == tag]
        conv = conv and all(tag_rows[i]["obj_minus_QG"]
                            > tag_rows[i + 1]["obj_minus_QG"]
                            for i in range(len(tag_rows) - 1))
    check("theory constructions built", len(rows) == len(ok_rows),
          f"({len(ok_rows)}/{len(rows)})")
    check("six-scalar recovery", max_six < 1e-10, f"max={max_six:.2e}")
    check("spectral summary errors", max_spec < 1e-10, f"max={max_spec:.2e}")
    check("M_epsilon positive definite", min_pd > 0, f"min={min_pd:.2e}")
    check("epsilon-loss identity within bound", eps_loss_ok)
    check("paper bound (22) holds at minimizers", b22)
    check("constructed objective >= Q_G", obj_ge)
    check("Q_G <= q* on constructions", safe)
    check("references valid", valid)
    check("objective converges in eps", conv)
    summary = dict(n_rows=len(rows), max_six_scalar_err=max_six,
                   max_spectral_err=max_spec, min_lambda_min_M=min_pd,
                   bound22_all_hold=bool(b22), q_star_ge_QG_all=bool(safe),
                   references_valid_all=bool(valid))
    utils.write_json(os.path.join(utils.RESULTS_DIR, "theory_summary.json"),
                     summary)
    return summary


# ==========================================================================
# 2. theory: Theorem 3 boundary battery + nondegenerate family
# ==========================================================================
BOUNDARY_FOUR = [  # (p, I0); a = alpha = gamma = 1, beta = 0.7 fixed
    (2.0, 1.0), (0.2, 1.0), (2.0, 1.3), (0.2, 1.3)]


def run_boundary(quick=False):
    utils.ensure_dirs()
    rng = np.random.default_rng(THEORY_SEED + 1)
    rows = []
    for p, I0 in BOUNDARY_FOUR:
        bc = core.construct_boundary_v41(p, 1.0, I0, 1.0, 0.7, 1.0, 0.05, 0.2)
        summ = core.spectral_summary(bc["b"], bc["h"], bc["C"], 1.0, 0.05)
        obj = core.objective_at(bc["b"], bc["h"], bc["C"], 1.0, 0.05, bc["e"])
        ref = core.q_star_reference(summ, 0.2, 1.0, 0.05)
        qg, _, _ = core.Q_approx(bc["scal"], 0.2, 1.0, 0.05, "G")
        rows.append(dict(tag="fixed_four", p=p, a=1.0, I0=I0, alpha=1.0,
                         beta=0.7, r=0.2, Jperp=bc["Jperp"],
                         c_parallel=bc["c_parallel"],
                         Q_partial=bc["Q_partial"], obj_at_e=obj,
                         q_star=ref["q"], Q_G_general=qg,
                         reference_valid=ref["reference_valid"],
                         six_scalar_err=max(abs(summ[k] - bc["scal"][k])
                                            for k in bc["scal"])))
    n_rand = 6 if quick else 30
    for _ in range(n_rand):
        al = float(rng.uniform(0.2, 2.0))
        be = float(rng.uniform(0.1, al))
        a = float(rng.uniform(0.0, 2.0))
        J = float(rng.uniform(0.0, 2.0))
        p = float(rng.uniform(0.1, 4.0))
        gam = float(rng.uniform(0.2, 2.0))
        r = float(rng.uniform(0.1, 1.0))
        I0 = float(np.hypot(al * a, J))
        bc = core.construct_boundary_v41(p, a, I0, al, be, gam, 0.05, r)
        summ = core.spectral_summary(bc["b"], bc["h"], bc["C"], gam, 0.05)
        obj = core.objective_at(bc["b"], bc["h"], bc["C"], gam, 0.05, bc["e"])
        ref = core.q_star_reference(summ, r, gam, 0.05)
        rows.append(dict(tag="random_boundary", p=p, a=a, I0=I0, alpha=al,
                         beta=be, r=r, Jperp=bc["Jperp"],
                         c_parallel=bc["c_parallel"],
                         Q_partial=bc["Q_partial"], obj_at_e=obj,
                         q_star=ref["q"], Q_G_general=np.nan,
                         reference_valid=ref["reference_valid"],
                         six_scalar_err=max(abs(summ[k] - bc["scal"][k])
                                            for k in bc["scal"])))
    fields = list(rows[0].keys())
    utils.write_csv(os.path.join(utils.RESULTS_DIR, "boundary_instances.csv"),
                    fields, rows)
    worst_att = max(abs(r["Q_partial"] - r["q_star"]) for r in rows)
    worst_obj = max(abs(r["obj_at_e"] - r["Q_partial"]) for r in rows)
    both_signs = {r["c_parallel"] >= 0 for r in rows} == {True, False}
    check("boundary attainment Q_partial == q*", worst_att < 1e-8,
          f"max dev={worst_att:.2e}")
    check("boundary objective == Q_partial", worst_obj < 1e-8,
          f"max dev={worst_obj:.2e}")
    check("both c_parallel signs covered", both_signs)
    fixed = [r for r in rows if r["tag"] == "fixed_four"]
    for r in fixed:
        print(f"[v41:info] boundary four: p={r['p']}, I0={r['I0']} -> "
              f"Q_partial={r['Q_partial']:.10f} (general Q_G="
              f"{r['Q_G_general']:.6f})")

    # nondegenerate family approaching the boundary + boundary row
    fam = []
    scal_b = dict(p=2.0, a=1.0, I0=1.0, alpha=1.0, beta=0.7, beta0=1.0)
    qg_b, _, _ = core.Q_approx(scal_b, 0.2, 1.0, 0.05, "G")
    fam.append(dict(delta=0.0, beta0=1.0, Q_G=qg_b, Q_partial=0.59,
                    h_norm=1.0, obj_at_e=np.nan, note="boundary instance"))
    for delta in (0.1, 0.03, 0.01):
        CH = np.array([[1.0, delta, 0.0], [0.0, 0.0, 0.7]], dtype=complex)
        h = np.array([1.0, -2.0 / delta, 0.0], dtype=complex).reshape(-1, 1)
        b = np.array([np.sqrt(2.0), 0.0, 0.0], dtype=complex).reshape(-1, 1)
        s = core.spectral_summary(b, h, CH.conj().T, 1.0, 0.05)
        qg, _, _ = core.Q_approx(s, 0.2, 1.0, 0.05, "G")
        ref = core.q_star_reference(s, 0.2, 1.0, 0.05)
        e_eval = np.array([-0.2, 0.0, 0.0], dtype=complex).reshape(-1, 1)
        fam.append(dict(delta=delta, beta0=s["beta0"], Q_G=qg,
                        Q_partial=np.nan,
                        h_norm=float(np.linalg.norm(h)),
                        obj_at_e=core.objective_at(b, h, CH.conj().T, 1.0,
                                                   0.05, e_eval),
                        q_star=ref["q"], note="nondegenerate family"))
    utils.write_csv(os.path.join(utils.RESULTS_DIR, "boundary_family.csv"),
                    list(fam[0].keys()), fam)
    last = fam[-1]
    check("nondegenerate family -> Q_G -> -0.21",
          abs(last["Q_G"] + 0.21) < 5e-3 and last["h_norm"] > 100,
          f"Q_G={last['Q_G']:.6f}, |h|={last['h_norm']:.1f}")
    return rows, fam


# ==========================================================================
# 3. theory: degenerate and endpoint battery
# ==========================================================================
def run_degenerate(quick=False):
    utils.ensure_dirs()
    rng = np.random.default_rng(THEORY_SEED + 2)
    rows = []

    def rec(tag, assumption, detail, value, expected, tol=1e-9):
        dev = abs(value - expected) if value == value else np.inf
        rows.append(dict(tag=tag, assumption=assumption, detail=detail,
                         value=value, expected=expected, deviation=dev,
                         ok=bool(dev <= tol)))

    # alpha = 0 (forces beta0 = beta): explicit V4.1 attainment instance
    for _ in range(3):
        be = float(rng.uniform(0.3, 1.5))
        a = float(rng.uniform(0.1, 2.0))
        I0 = float(rng.uniform(0.0, 1.5))
        p = float(rng.uniform(0.5, 3.0))
        gam, noise, r = 1.0, 0.05, 0.4
        CH = np.array([[0.0, be]], dtype=complex)
        h = np.array([a, I0 / be], dtype=complex).reshape(-1, 1)
        b = np.array([np.sqrt(p), 0.0], dtype=complex).reshape(-1, 1)
        s = core.spectral_summary(b, h, CH.conj().T, gam, noise)
        sc = dict(p=p, a=a, I0=I0, alpha=0.0, beta=be, beta0=be)
        qg, _, _ = core.Q_approx(sc, r, gam, noise, "G")
        ref = core.q_star_reference(s, r, gam, noise)
        rec("alpha_zero", "alpha=0 => beta0=beta", f"beta={be:.3f}",
            ref["q"], qg, 1e-8)
        six = max(abs(s[k] - sc[k]) for k in sc)
        rows[-1].update(six_scalar_err=six)

    # r = 0: every compatible instance has q* = p a^2 - gamma I0^2 - gamma s^2
    for _ in range(3):
        n = 5
        b = np.zeros((n, 1), complex)
        b[0, 0] = np.sqrt(2.0)
        C = (rng.normal(size=(n, 3)) + 1j * rng.normal(size=(n, 3))) * 0.3
        h = (rng.normal(size=(n, 1)) + 1j * rng.normal(size=(n, 1)))
        gam, noise = 0.8, 0.05
        s = core.spectral_summary(b, h, C, gam, noise)
        ident = s["p"] * s["a"] ** 2 - gam * s["I0"] ** 2 - gam * noise
        ref = core.q_star_reference(s, 0.0, gam, noise)
        sc = dict(p=s["p"], a=s["a"], I0=s["I0"], alpha=s["alpha"],
                  beta=s["beta"], beta0=s["beta0"])
        qg, _, _ = core.Q_approx(sc, 0.0, gam, noise, "G")
        rec("r_zero", "r=0: singleton ball", "random instance", ref["q"],
            ident, 1e-10)
        rows[-1].update(Q_G_zero=qg)

    # beta = 0 (forces beta0 = alpha and I0 = alpha a): parallel attainment
    for _ in range(3):
        al = float(rng.uniform(0.4, 1.5))
        a = float(rng.uniform(0.1, 2.0))
        p = float(rng.uniform(0.5, 3.0))
        gam, noise, r = 1.0, 0.05, 0.4
        CH = np.array([[al, 0.0]], dtype=complex)
        h = np.array([a, 0.0], dtype=complex).reshape(-1, 1)
        b = np.array([np.sqrt(p), 0.0], dtype=complex).reshape(-1, 1)
        s = core.spectral_summary(b, h, CH.conj().T, gam, noise)
        sc = dict(p=p, a=a, I0=al * a, alpha=al, beta=0.0, beta0=al)
        qp, _ = core.Q_boundary_approx(sc, r, gam, noise)
        ref = core.q_star_reference(s, r, gam, noise)
        rec("beta_zero", "beta=0 => beta0=alpha, I0=alpha a",
            f"alpha={al:.3f}", ref["q"], qp, 1e-8)
        qg, _, _ = core.Q_approx(sc, r, gam, noise, "G")
        rows[-1].update(Q_G=qg, Q_G_le_qstar=bool(qg <= ref["q"] + 1e-9))

    # branch connection points of the scalar extremum
    from experiments_a import scalar_extremum_ref
    for _ in range(6):
        al, be = rng.uniform(0.2, 1.5, 2)
        b0 = rng.uniform(max(al, be), np.hypot(al, be))
        delta = np.sqrt(b0 ** 2 - al ** 2)
        nu = np.sqrt(b0 ** 2 - be ** 2)
        for t in (0.15, 0.4):
            for name, (xb, yb) in (("conn_12", (nu * t, delta * t)),
                                   ("conn_23", (al * t, delta * t))):
                R = core.R_envelope(xb, yb, al, be, b0)
                refm, _ = scalar_extremum_ref(al, be, b0, xb, yb)
                rec(name, "delta x = nu y / delta x = alpha y",
                    f"t={t}", R ** 2, refm ** 2, 1e-10)

    # trust region: hard case equality, repeated min eigenvalue, r equality
    A = np.diag([-1.0, 1.0]).astype(complex)
    g = np.array([[1.0], [0.0]])
    out = core.trust_region_min(A, g, 0.5)
    rec("tr_hard_equality", "phi_pi = r^2 exactly (pinv on ball)",
        f"mode={out['mode']}", out["quad_val"], -1.25, 1e-10)
    rows[-1].update(ball_viol=out["checks"]["ball_viol"],
                    stat_res=out["checks"]["stat_res"])
    A = np.diag([-1.0, -1.0, 2.0]).astype(complex)
    g = np.array([[0.0], [0.0], [1.0]])
    out = core.trust_region_min(A, g, 0.7)
    # minimizer: e3 component t = -1/3 (interior), remaining norm in the
    # doubled minimum eigenspace: value = -0.49 - 1/3
    rec("tr_repeated_min", "lambda_min multiplicity 2, g in eigenspace",
        f"mode={out['mode']}", out["quad_val"], -0.49 - 1.0 / 3.0, 1e-10)
    A = np.diag([1.0, 0.0]).astype(complex)
    g = np.array([[1.0], [0.0]])
    out = core.trust_region_min(A, g, 1.0)
    rec("tr_r_equality", "PSD pinv norm equals r (boundary interior)",
        f"mode={out['mode']}", out["quad_val"], -1.0, 1e-10)
    rows[-1].update(norm_e=float(np.linalg.norm(out["e"])))

    fields = ["tag", "assumption", "detail", "value", "expected", "deviation",
              "ok"]
    extra = sorted({k for r in rows for k in r} - set(fields))
    utils.write_csv(os.path.join(utils.RESULTS_DIR, "degenerate_endpoints.csv"),
                    fields + extra, rows)
    bad = [r for r in rows if not r["ok"]]
    check("degenerate/endpoints battery", not bad,
          f"{len(rows) - len(bad)}/{len(rows)} ok"
          + (f"; failing: {[r['tag'] for r in bad]}" if bad else ""))
    return rows


# ==========================================================================
# 4. scan (experiment B) under V4.1 diagnostics
# ==========================================================================
def run_scan():
    brows = experiments_b.run_b()
    tol = 1e-9
    viol_lb = [r for r in brows
               if r["reference_valid"] and r["q_star_minus_QG_LB"] < -tol]
    viol_ub = [r for r in brows
               if r["reference_valid"] and r["q_star_minus_Q_G"] < -tol]
    lq_above = sum(1 for r in brows if r["Q_G"] - r["L_quad"] < -tol)
    summary = dict(
        n_points=len(brows),
        reference_valid_all=bool(all(r["reference_valid"] for r in brows)),
        tol=tol,
        lb_violations=len(viol_lb),
        lb_violation_max_dev=(max(-r["q_star_minus_QG_LB"] for r in viol_lb)
                              if viol_lb else 0.0),
        approx_violations=len(viol_ub),
        approx_violation_max_dev=(max(-r["q_star_minus_Q_G"] for r in viol_ub)
                                  if viol_ub else 0.0),
        min_gap_qstar_to_lb=min(r["q_star_minus_QG_LB"] for r in brows),
        l_quad_above_QG=lq_above,
        note="approximate-minimum inversions below tol are roundoff, "
             "never an acceptance criterion")
    utils.write_json(os.path.join(utils.RESULTS_DIR, "b_scan_summary.json"),
                     summary)
    check("B scan: certified LB never above reference", len(viol_lb) == 0)
    check("B scan: approx minima above reference only within roundoff",
          len(viol_ub) == 0
          or summary["approx_violation_max_dev"] < 1e-12,
          f"max dev={summary['approx_violation_max_dev']:.2e}")
    return brows, summary


# ==========================================================================
# 5. pools
# ==========================================================================
def run_baseline_pool():
    """Original easy pool (gamma=0.5, rho in {0.01,0.03,0.06}); re-run into
    results_v41 for a self-contained record.  This pool saturates and cannot
    demonstrate fallback economics."""
    cres = ec.run_c()
    n_fb = cres["total_fallback"]
    check("baseline pool: no misaccept",
          all(cres["agg"]["stats"][m]["misaccept"] == 0
              for m in ec.METHODS))
    check("baseline pool easy-regime caveat recorded", n_fb == 0,
          f"fallback calls={n_fb} (saturated; no fallback economics)")
    return cres


def run_grid(seeds, cohort, prefix, quick=False):
    """Fixed-grid study for one cohort: details + per-cell aggregates +
    paired interleaved timing of the same final-decision pipelines."""
    details, summaries, timings = [], [], []
    for seed in seeds:
        cands = ec.all_candidates(seed)
        for gamma in GAMMAS:
            ss, norms = az.summarize(cands, gamma)
            rows = []
            for i, cand in enumerate(cands):
                for u, s in enumerate(ss[i]):
                    eig = np.linalg.eigh(s["A"])
                    for rho in RHOS:
                        r = rho * norms[i][u]
                        ref = core.q_star_reference(s, r, gamma, ec.NOISE_U,
                                                    eig=eig)
                        row = dict(
                            seed=seed, gamma=gamma, cand_idx=i, user=u,
                            rho=rho, btype=cand["btype"], q_star=ref["q"],
                            tau=ref["tau"],
                            reference_valid=ref["reference_valid"],
                            status=ec.user_status(ref),
                            D=core.D_bound(s, r, gamma, ec.NOISE_U),
                            L_quad=core.L_quad(s, r, gamma),
                            Q_old=core.interval_lower_bound(
                                s, r, gamma, ec.NOISE_U, 256, "old"),
                            Q_G=core.interval_lower_bound(
                                s, r, gamma, ec.NOISE_U, 256, "G"))
                        if not ref["reference_valid"]:
                            FAILURES.append(f"{prefix}: invalid reference at "
                                            f"seed={seed} u={u} rho={rho}")
                        if max(row[m] for m in ec.METHODS) > ref["q"] + 1e-8:
                            FAILURES.append(f"{prefix}: bound above q* at "
                                            f"seed={seed} u={u} rho={rho}")
                        rows.append(row)
            details.extend(rows)
            for rho in RHOS:
                agg, _ = ec.aggregate([row for row in rows
                                       if row["rho"] == rho])
                n_units = agg["n_units"]
                for method in ec.METHODS:
                    st = agg["stats"][method]
                    cov = agg["coverage_" + method]
                    summaries.append(dict(
                        seed=seed, cohort=cohort, gamma=gamma, rho=rho,
                        method=method, units=n_units,
                        feasible=agg["ref_feasible_units"],
                        infeasible=agg["ref_infeasible_units"],
                        ambiguous=agg["ambiguous_units"],
                        invalid=agg["invalid_units"],
                        accepted=st["accepted"],
                        accept_rate_all=st["accepted"] / n_units,
                        coverage_feasible=("NA" if cov is None else cov),
                        misaccept=st["misaccept"], unknown=st["unknown"],
                        qg_added=agg["qg_added"] if method == "Q_G" else 0))
            methods = ("reference",) + ec.METHODS
            for method in methods:
                az.pipeline(cands[:1], gamma, method)      # warm-up
            runs = {m: [] for m in methods}
            for rep in range(3):
                for method in np.roll(np.array(methods), rep):
                    runs[method].append(az.pipeline(cands, gamma, method))
            base = runs["reference"][0]
            for method in methods:
                for rr in runs[method]:
                    if rr["decisions"] != base["decisions"]:
                        FAILURES.append(f"{prefix}: decision mismatch "
                                        f"seed={seed} gamma={gamma} "
                                        f"method={method}")
                timings.append(dict(
                    seed=seed, gamma=gamma, method=method,
                    time_median=float(np.median([x["seconds"]
                                                 for x in runs[method]])),
                    time_min=min(x["seconds"] for x in runs[method]),
                    time_max=max(x["seconds"] for x in runs[method]),
                    calls=runs[method][0]["calls"],
                    direct=runs[method][0]["direct"],
                    accepted=runs[method][0]["accepted"],
                    rejected=runs[method][0]["rejected"],
                    unresolved=runs[method][0]["unresolved"]))
        print(f"[v41:info] grid {prefix}: seed {seed} done", flush=True)
    for name, rr in ((prefix + "_details", details),
                     (prefix + "_summary", summaries),
                     (prefix + "_timings", timings)):
        utils.write_csv(os.path.join(utils.RESULTS_DIR, name + ".csv"),
                        list(rr[0].keys()), rr)
    mis = sum(r["misaccept"] for r in summaries)
    check(f"{prefix}: zero misacceptance (unit level)", mis == 0)
    check(f"{prefix}: references valid",
          all(r["reference_valid"] for r in details))
    return details, summaries, timings


def _sha256(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def write_validation_config(mode):
    """Pre-registration: written BEFORE the validation seeds are run."""
    used = set(ec.SEEDS) | set(EXPLORATION_SEEDS)
    assert not (set(VALIDATION_SEEDS) & used), \
        "validation seeds must be unused"
    here = os.path.dirname(os.path.abspath(__file__))
    cfg = dict(
        mode=mode,
        written_at=datetime.now().isoformat(timespec="seconds"),
        purpose="Pre-registered validation batch: fixed protocol declared "
                "before any validation-seed result was computed",
        validation_seeds=list(VALIDATION_SEEDS),
        previously_used_seeds=sorted(used),
        gammas=list(GAMMAS), rhos=list(RHOS),
        candidate_generation=dict(
            n_bs=2, n_antennas=4, n_users=4, n_ris=16, phases_2bit=True,
            phase_states_per_seed=6, first_phase_all_ones=True,
            beam_types=["MRT", "RZF"], rzf_regularization=0.1,
            per_stream_power=0.5, same_cell_scale=1.0, cross_cell_scale=0.5,
            var_G="CN(0,1/4)", var_f="CN(0,1/16)", var_d="CN(0,1/4)",
            effective_link="h_{j,u}=s[d + G_j^H diag(theta) f_u]"),
        acceptance=dict(
            rule="interval lower bound > tau for all four users",
            J_cells=256, kind="Q_G",
            tau_rule="1e-8 * max(1,|c_q|,2r||g||,r^2||A||_2)",
            nonpositive="unknown only; never a rejection"),
        reference_validity=("KKT checks in q_star_reference: ball "
                            "violation, full stationarity residual, "
                            "lambda_min(A+mu I), complementarity"),
        timing=dict(repetitions=3, interleaved=True, warmup=1,
                    per_user_eigh_cached=True, skip_certified_users=True,
                    early_stop="first clearly infeasible user",
                    candidate_generation_timed=False,
                    pipelines=["reference", "D", "Q_old", "Q_G", "L_quad"]),
        resolution_ablation_J=[256, 1024, 4096],
        code_hashes={f: _sha256(os.path.join(here, f)) for f in
                     ("core.py", "experiments_c.py", "experiments_v41.py",
                      "audit_and_stress.py")},
    )
    utils.write_json(os.path.join(utils.RESULTS_DIR, "validation_config.json"),
                     cfg)
    return cfg


def run_resolution(validation_details, quick=False):
    """J = 256/1024/4096 ablation (exploration + validation cohorts) and the
    theoretical-vs-numerical loss decomposition on validation rows."""
    rows = []
    cohorts = [("exploration", EXPLORATION_SEEDS),
               ("validation", VALIDATION_SEEDS)]
    if quick:
        cohorts = [(name, seeds[:2]) for name, seeds in cohorts]
    for cohort, seeds in cohorts:
        counts = {J: dict(old=0, G=0, added=0, feasible=0, units=0, mis=0)
                  for J in (256, 1024, 4096)}
        for seed in seeds:
            cands = ec.all_candidates(seed)
            for gamma in GAMMAS:
                ss, norms = az.summarize(cands, gamma)
                for rho in RHOS:
                    refs = []
                    for i in range(len(cands)):
                        for u, s in enumerate(ss[i]):
                            eig = np.linalg.eigh(s["A"])
                            ref = core.q_star_reference(
                                s, rho * norms[i][u], gamma, ec.NOISE_U,
                                eig=eig)
                            refs.append((s, rho * norms[i][u],
                                         ref["q"], ref["tau"]))
                    for cand_i in range(len(cands)):
                        users = refs[cand_i * 4:(cand_i + 1) * 4]
                        feas = all(q > tau for _, _, q, tau in users)
                        for J in (256, 1024, 4096):
                            acc = {}
                            for kind in ("old", "G"):
                                lbs = [core.interval_lower_bound(
                                    s, r, gamma, ec.NOISE_U, J, kind)
                                    for s, r, _, _ in users]
                                acc[kind] = all(lb > tau
                                                for lb, (_, _, _, tau)
                                                in zip(lbs, users))
                                if acc[kind] and not feas:
                                    counts[J]["mis"] += 1
                            c = counts[J]
                            c["units"] += 1
                            c["old"] += acc["old"]
                            c["G"] += acc["G"]
                            c["added"] += acc["G"] and not acc["old"]
                            c["feasible"] += feas
        for J in (256, 1024, 4096):
            c = counts[J]
            rows.append(dict(cohort=cohort, J=J, units=c["units"],
                             feasible=c["feasible"], old_accepts=c["old"],
                             G_accepts=c["G"], added=c["added"],
                             misaccept=c["mis"],
                             added_pct_of_feasible=(c["added"] / c["feasible"]
                                                    if c["feasible"] else
                                                    "NA")))
    utils.write_csv(os.path.join(utils.RESULTS_DIR, "resolution.csv"),
                    list(rows[0].keys()), rows)
    for r in rows:
        print(f"[v41:info] resolution {r['cohort']} J={r['J']}: "
              f"old={r['old_accepts']} G={r['G_accepts']} "
              f"added={r['added']} feasible={r['feasible']}")
    check("resolution: no misacceptance at any J",
          all(r["misaccept"] == 0 for r in rows))

    # loss decomposition on validation rows: upper estimate vs certified LB
    # vs exact reference (theoretical envelope loss vs interval loss)
    loss = []
    for seed in VALIDATION_SEEDS[:2] if quick else VALIDATION_SEEDS:
        cands = ec.all_candidates(seed)
        for gamma in GAMMAS:
            ss, norms = az.summarize(cands, gamma)
            for i, cand in enumerate(cands):
                for u, s in enumerate(ss[i]):
                    for rho in RHOS:
                        r = rho * norms[i][u]
                        qg_up, _, _ = core.Q_approx(s, r, gamma, ec.NOISE_U,
                                                    "G")
                        lb = core.interval_lower_bound(s, r, gamma,
                                                       ec.NOISE_U, 256, "G")
                        eig = np.linalg.eigh(s["A"])
                        ref = core.q_star_reference(s, r, gamma, ec.NOISE_U,
                                                    eig=eig)
                        loss.append(dict(seed=seed, gamma=gamma, cand_idx=i,
                                         user=u, rho=rho,
                                         Q_G_upper_estimate=qg_up,
                                         Q_G_LB_256=lb, q_star=ref["q"],
                                         envelope_gap=ref["q"] - qg_up,
                                         interval_gap=qg_up - lb))
    utils.write_csv(os.path.join(utils.RESULTS_DIR,
                                 "resolution_loss_decomposition.csv"),
                    list(loss[0].keys()), loss)
    return rows, loss


def run_overlap(baseline_rows, expl_details, val_details):
    def keys(rows, gamma=None):
        out = set()
        for r in rows:
            g = r.get("gamma", gamma if gamma else 0.5)
            out.add((r["seed"], r["cand_idx"], r["user"], r["rho"], g))
        return out

    kb = keys(baseline_rows, 0.5)
    ke = keys(expl_details)
    kv = keys(val_details)
    out = dict(
        baseline_user_problems=len(kb),
        exploration_user_problems=len(ke),
        validation_user_problems=len(kv),
        baseline_in_exploration=len(kb & ke),
        baseline_or_exploration_union=len(kb | ke),
        validation_overlap_with_previous=len(kv & (kb | ke)),
        all_three_union=len(kb | ke | kv),
        note="user-level problem = (seed, candidate, user, rho, gamma); "
             "multiple radii/thresholds on one candidate are correlated")
    check("validation seeds disjoint from previous pools",
          out["validation_overlap_with_previous"] == 0)
    utils.write_json(os.path.join(utils.RESULTS_DIR, "pools_overlap.json"),
                     out)
    return out


# ==========================================================================
# 6. trust-region multiplier search: fixed 200 iterations vs tolerance stop
# ==========================================================================
def run_tr_tolerance(quick=False):
    cases = [([-1., 1.], [0., 0.], 1.), ([-1., 1.], [0., 2.], 1.),
             ([-1., -1., 2.], [0., 0., 1.], 1.), ([0., 1.], [0., 0.5], 1.),
             ([0., 1.], [1., 0.5], 1.), ([1., 2.], [0.1, 0.2], 1.)]
    rng = np.random.default_rng(20261007)
    instances = [(np.diag(l).astype(complex), np.array(g, complex
                                                       ).reshape(-1, 1), r)
                 for l, g, r in cases]
    for _ in range(100 if not quick else 20):
        Z = rng.normal(size=(5, 5)) + 1j * rng.normal(size=(5, 5))
        A = (Z + Z.conj().T) / 2
        g = (rng.normal(size=5) + 1j * rng.normal(size=5)).reshape(-1, 1)
        instances.append((A, g, rng.uniform(0.1, 2.0)))
    rows = []
    for idx, (A, g, r) in enumerate(instances):
        base = core.trust_region_min(A, g, r)
        row = dict(idx=idx, mode=base["mode"], iters_200=base["bisect_used"])
        for tol in (1e-10, 1e-12):
            out = core.trust_region_min(A, g, r, mu_rel_tol=tol)
            sc = max(1.0, abs(base["quad_val"]))
            row[f"iters_tol{tol:.0e}"] = out["bisect_used"]
            row[f"dq_tol{tol:.0e}"] = abs(out["quad_val"] - base["quad_val"])
            row[f"valid_tol{tol:.0e}"] = bool(
                out["checks"]["ball_viol"] <= 1e-8 * max(1.0, r)
                and out["checks"]["stat_res"] <= 1e-6 * sc
                and out["checks"]["lam_min_Amu"] >= -1e-8 * sc)
        rows.append(row)
    utils.write_csv(os.path.join(utils.RESULTS_DIR, "tr_tolerance_audit.csv"),
                    list(rows[0].keys()), rows)
    usable = []
    for tol in (1e-10, 1e-12):
        max_dq = max(r[f"dq_tol{tol:.0e}"] for r in rows)
        all_valid = all(r[f"valid_tol{tol:.0e}"] for r in rows)
        med_it = float(np.median([r[f"iters_tol{tol:.0e}"] for r in rows]))
        meets = bool(max_dq < 1e-8 and all_valid)
        if meets:
            usable.append(tol)
        print(f"[v41:{'PASS' if meets else 'OBSERVATION'}] TR tol-stop "
              f"{tol:.0e}: max|dq|={max_dq:.2e} (bar 1e-8), valid="
              f"{all_valid}, median iters {med_it:.0f} (vs fixed 200)")
    check("TR: a tolerance variant matches the 200-iteration values "
          "within the 1e-8 bar with preserved validity", bool(usable),
          f"usable tolerances={usable}")

    # same final-decision task: reference-200 vs reference-tol vs Q_G+fb
    # (200-iter fallback) vs Q_G+fb with tol-stopped fallback
    def _run_pipe(cands, gamma, label):
        if label == "reference_tol":
            return az.pipeline(cands, gamma, "reference", mu_rel_tol=1e-12)
        if label == "Q_G_tolfb":
            return az.pipeline(cands, gamma, "Q_G", mu_rel_tol=1e-12)
        return az.pipeline(cands, gamma, label)

    methods4 = ("reference", "reference_tol", "Q_G", "Q_G_tolfb")
    trows = []
    seeds = VALIDATION_SEEDS[:2] if quick else VALIDATION_SEEDS
    for seed in seeds:
        cands = ec.all_candidates(seed)
        for gamma in GAMMAS:
            for m in methods4:
                _run_pipe(cands[:1], gamma, m)          # warm-up
            runs = {m: [] for m in methods4}
            for rep in range(3):
                order = np.roll(np.array(methods4), rep)
                for m in order:
                    runs[m].append(_run_pipe(cands, gamma, m))
            for m in runs:
                base = runs["reference"][0]
                for rr in runs[m]:
                    if rr["decisions"] != base["decisions"]:
                        FAILURES.append(f"tr-tol timing decision mismatch "
                                        f"seed={seed} gamma={gamma} m={m}")
                trows.append(dict(
                    seed=seed, gamma=gamma, method=m,
                    time_median=float(np.median([x["seconds"]
                                                 for x in runs[m]])),
                    time_min=min(x["seconds"] for x in runs[m]),
                    time_max=max(x["seconds"] for x in runs[m]),
                    calls=runs[m][0]["calls"]))
    utils.write_csv(os.path.join(utils.RESULTS_DIR, "tr_tolerance_timing.csv"),
                    list(trows[0].keys()), trows)
    for m in methods4:
        med = float(np.median([r["time_median"] for r in trows
                               if r["method"] == m]))
        print(f"[v41:info] same-task timing {m}: median {med:.4f} s")
    return rows, trows


# ==========================================================================
# 7. figures from saved data
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

    def rd(name):
        with open(os.path.join(utils.RESULTS_DIR, name), encoding="utf-8") as f:
            return list(csv.DictReader(f))

    figdir = os.path.join(utils.RESULTS_DIR, "figures")
    os.makedirs(figdir, exist_ok=True)

    # fig 1: theory constructions
    tc = rd("theory_constructions.csv")
    fam = [r for r in tc if r["tag"] == "flip_family"]
    if fam:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.6))
        eps = np.array([float(r["eps"]) for r in fam])
        loss = np.array([float(r["obj_minus_QG"]) for r in fam])
        bnd = np.array([float(r["bound22"]) for r in fam])
        rr = np.array([float(r["R_minus_R_eps"]) for r in fam])
        ax1.loglog(eps, loss, "o-", color=PAL[0], ms=5, lw=1.2,
                   label="$q(e_\\varepsilon)-Q_G$")
        ax1.loglog(eps, rr, "s--", color=PAL[2], ms=5, lw=1.2,
                   label="$R-R_\\varepsilon$")
        ax1.loglog(eps, bnd, "^:", color=PAL[1], ms=5, lw=1.2,
                   label="paper bound (22)")
        ax1.set_xlabel("contraction $\\varepsilon$")
        ax1.set_ylabel("deficit")
        ax1.set_title("(a) V4.1 construction: deficit vs $\\varepsilon$")
        ax1.legend(loc="upper left", fontsize=7.5)
        six = np.maximum([float(r["six_scalar_err"]) for r in fam], 1e-17)
        pd_ = [float(r["lam_min_M"]) for r in fam]
        ax2.loglog(eps, six, "s--", color=PAL[2], ms=5, lw=1.2,
                   label="six-scalar error (=0, floored)")
        ax2.loglog(eps, pd_, "o-", color=PAL[0], ms=5, lw=1.2,
                   label="$\\lambda_{\\min}(M_\\varepsilon)$")
        ax2.axhline(1e-16, color=GRAY, lw=0.7, ls=":")
        ax2.set_xlabel("contraction $\\varepsilon$")
        ax2.set_ylabel("error / residual")
        ax2.set_title("(b) exactness and positive definiteness")
        ax2.legend(loc="center left", fontsize=7.5)
        fig.tight_layout()
        fig.savefig(os.path.join(figdir, "fig_theory_v41.pdf"))
        fig.savefig(os.path.join(figdir, "fig_theory_v41.png"), dpi=150)
        plt.close(fig)

    # fig 2: spectral scan from saved CSV
    sc = rd("spectral_scan.csv")
    if sc:
        x = np.array([float(r["x_axis"]) for r in sc])

        def col(k):
            return np.array([float(r[k]) for r in sc])

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.6))
        series = [("D", PAL[0], "-", "o", "$D$"),
                  ("Q_old", PAL[1], "-", "s", "$Q$ (old, approx.)"),
                  ("Q_G", PAL[2], "-", "D", "$Q_G$ (approx.)"),
                  ("L_quad", PAL[3], "-", "^", "$L_{\\rm quad}$"),
                  ("q_star", PAL[4], "-", "v", "exact $q^*$")]
        for k, c, ls, mk, lab in series:
            ax1.plot(x, col(k), ls, color=c, marker=mk, ms=3, lw=1.2,
                     label=lab)
        ax1.plot(x, col("Q_G_LB"), "--", color=PAL[2], lw=1.2,
                 label="$Q_{G,\\rm LB}$ (certified, $N_\\theta$=256)")
        ax1.axhline(0.0, color=GRAY, lw=0.9)
        ax1.set_xlabel("$\\beta_0/\\sqrt{\\alpha^2+\\beta^2}$")
        ax1.set_ylabel("worst-case margin")
        ax1.set_title("(a) constructed scan (nondegenerate family)")
        ax1.legend(loc="lower left", fontsize=7.5, ncol=2)
        ax2.plot(x, col("q_star_minus_Q_G"), "-", color=PAL[4], marker="v",
                 ms=3, lw=1.2, label="$q^*-Q_G$")
        ax2.plot(x, col("q_star_minus_QG_LB"), "-", color=PAL[2], marker="D",
                 ms=3, lw=1.2, label="$q^*-Q_{G,\\rm LB}$")
        ax2.plot(x, col("Q_G_minus_Q_old"), "-", color=PAL[1], marker="s",
                 ms=3, lw=1.2, label="$Q_G-Q$ (old)")
        ax2.plot(x, col("Q_G") - col("L_quad"), "-", color=PAL[3],
                 marker="^", ms=3, lw=1.2, label="$Q_G-L_{\\rm quad}$")
        ax2.axhline(0.0, color=GRAY, lw=0.9)
        ax2.set_xlabel("$\\beta_0/\\sqrt{\\alpha^2+\\beta^2}$")
        ax2.set_ylabel("gap")
        ax2.set_title("(b) gaps between methods")
        ax2.legend(loc="best", fontsize=7.5)
        fig.tight_layout()
        fig.savefig(os.path.join(figdir, "fig_scan_v41.pdf"))
        fig.savefig(os.path.join(figdir, "fig_scan_v41.png"), dpi=150)
        plt.close(fig)
    print(f"[v41:info] figures written to {figdir}")


# ==========================================================================
# 8. driver
# ==========================================================================
def run_all(quick=False):
    utils.ensure_dirs()
    sections = {}
    print("--- V4.1 theory ---")
    sections["theory"] = run_theory(quick)
    sections["boundary"], fam = run_boundary(quick)
    sections["degenerate"] = run_degenerate(quick)
    print("--- V4.1 scan ---")
    _, scan_summary = run_scan()
    sections["scan"] = scan_summary
    print("--- V4.1 baseline pool ---")
    sections["baseline"] = run_baseline_pool()
    print("--- V4.1 exploration grid ---")
    expl_details, expl_sum, expl_tim = run_grid(EXPLORATION_SEEDS,
                                                "exploration", "exploration",
                                                quick)
    print("--- V4.1 validation (pre-registered) ---")
    write_validation_config("quick" if quick else "full")
    val_details, val_sum, val_tim = run_grid(VALIDATION_SEEDS, "validation",
                                             "validation", quick)
    print("--- V4.1 resolution + loss ---")
    sections["resolution"], loss = run_resolution(val_details, quick)
    print("--- V4.1 overlap ---")
    sections["overlap"] = run_overlap(sections["baseline"]["rows"],
                                      expl_details, val_details)
    print("--- V4.1 TR tolerance ---")
    sections["tr_tol_audit"], sections["tr_tol_timing"] = run_tr_tolerance(
        quick)
    print("--- V4.1 figures ---")
    make_figures()
    manifest = dict(
        mode="quick" if quick else "full",
        finished_at=datetime.now().isoformat(timespec="seconds"),
        failures=list(FAILURES),
        sections={k: (v if isinstance(v, (list, str)) else "ok")
                  for k, v in sections.items()},
        gamma_grid=list(GAMMAS), rho_grid=list(RHOS),
        exploration_seeds=list(EXPLORATION_SEEDS),
        validation_seeds=list(VALIDATION_SEEDS))
    utils.write_json(os.path.join(utils.RESULTS_DIR, "run_manifest.json"),
                     manifest)
    print(f"\n=== V4.1 suite {'FAILED' if FAILURES else 'PASSED'} ===")
    if FAILURES:
        for f in FAILURES:
            print("  FAIL:", f)
    return not FAILURES

