# -*- coding: utf-8 -*-
"""Corrected S-procedure-SDP baseline (v2).

Fixes the v1 decision-logic error: a certified LOWER bound t_cert < -tau
was wrongly classified as infeasible.  A negative lower bound proves
nothing about q* < 0.  The corrected per-user decision uses explicit
bounds:

    L = t_cert (residual-corrected dual-feasible lower bound):
        L > tau            -> certified_feasible
    U = q(e_witness) for a feasible error vector recovered from the
        (corrected) dual multiplier:
        U < -tau           -> certified_infeasible
    otherwise              -> unresolved   (never a rejection)

Witness recovery: from the corrected multiplier lam' the Lagrangian
minimizer solves (A + lam' I) e = -g  (least squares if singular); if
||e|| > r the vector is projected onto the ball (e_w = r e/||e||) and the
TRUE objective U = e_w^H A e_w + 2 Re(e_w^H g) + c is re-evaluated.  U is
a valid upper bound of q* because e_w is ball-feasible by construction;
no optimality of e_w is claimed.  The raw SDP value t is never used as an
upper bound.

Residual correction (retained from v1, restructured for consistency):
the shift (t,lam) -> (t-eta(1+r^2), lam+eta) adds blockdiag(eta*I, eta)
to the LMI matrix, so Weyl's inequality restores positive
semidefiniteness up to eigensolver rounding.  The retry loop is
restructured so the SAVED eta, t_cert, lam' and lam_min_corrected always
correspond to the SAME final re-checked matrix; 'psd_margin_achieved'
distinguishes a strictly positive margin from approximate PSD within
tolerance.  This is an empirical rounding allowance (ordinary eigvalsh),
NOT a formal error bound; no interval-arithmetic certificate is claimed.

r = 0 returns the nominal objective c directly; r < 0 raises.
All results go to results_v41_sdp_corrected/ (v1 results retained as
history in results_v41_sdp/).
"""
import hashlib
import json
import os
from datetime import datetime
from time import perf_counter

import numpy as np

import core
import experiments_c as ec
import optimization_v41 as opt
import utils

OUT_TAG = "results_v41_sdp_corrected"
GAMMAS = (0.5, 1.0, 2.0)
RHOS = (0.06, 0.15, 0.30)
VALIDATION_SEEDS = (211, 212, 213, 214, 215)
SOLVER = "CLARABEL"
ROUNDING_ALLOWANCE_K = 10.0

FAILURES = []


def check(name, ok, detail=""):
    print(f"[sdp2:{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        FAILURES.append(f"{name}: {detail}")
    return ok


# ---------------------------------------------------------------------------
# solver with witness recovery
# ---------------------------------------------------------------------------
def real_embed(A, g):
    n = A.shape[0]
    AR = np.block([[A.real, -A.imag], [A.imag, A.real]])
    gR = np.concatenate([g.real.ravel(), g.imag.ravel()]).reshape(-1, 1)
    return AR, gR


def recover_witness(A, g, c, r, lam_corr):
    """Feasible error vector from the corrected multiplier and its TRUE
    objective value (an upper bound of q*).  Returns dict or None."""
    n = A.shape[0]
    M = np.asarray(A, dtype=complex) + lam_corr * np.eye(n, dtype=complex)
    if not np.all(np.isfinite(M)) or not np.all(np.isfinite(g)):
        return None
    e, *_ = np.linalg.lstsq(M, -np.asarray(g, dtype=complex).reshape(-1, 1),
                            rcond=None)
    if not np.all(np.isfinite(e)):
        return None
    norm_raw = float(np.linalg.norm(e))
    if norm_raw == 0.0:
        e_w = e
        projected = False
    elif norm_raw > r:
        e_w = e * (r / norm_raw)
        projected = True
    else:
        e_w = e
        projected = False
    norm_w = float(np.linalg.norm(e_w))
    ball_resid = max(0.0, norm_w - r)
    U = (float(np.real((e_w.conj().T @ np.asarray(A, dtype=complex) @ e_w)
                       [0, 0]))
         + 2.0 * float(np.real((e_w.conj().T @ np.asarray(g, dtype=complex)
                                .reshape(-1, 1))[0, 0])) + float(c))
    return dict(e_norm_raw=norm_raw, e_norm=norm_w, projected=projected,
                ball_residual=ball_resid, U=U)


def sdp_qstar(A, g, c, r, solver=SOLVER, want_witness=True):
    """Solve the S-procedure SDP, apply the residual correction and recover
    a feasible witness.  Returns dict with L (=t_cert), U, statuses and
    diagnostics; never uses the raw t as a decision quantity."""
    import cvxpy as cp

    A = np.asarray(A, dtype=complex)
    g = np.asarray(g, dtype=complex).reshape(-1, 1)
    c = float(c)
    if r < 0:
        raise ValueError("negative radius is invalid")
    n = A.shape[0]
    out = dict(status=None, t=None, lam=None, eta=np.nan, L=None, U=None,
               lam_min_lmi=None, lam_min_corrected=None,
               psd_margin_achieved=None, witness=None, build_s=np.nan,
               solve_s=np.nan, wall_s=np.nan, solver=solver, valid=False,
               note="")
    herm = float(np.max(np.abs(A - A.conj().T)))
    out["herm_dev"] = herm
    if not (np.isfinite(A).all() and np.isfinite(g).all()
            and np.isfinite(c)):
        out.update(status="invalid_input", note="non-finite data")
        return out
    if herm > 1e-10 * max(1.0, float(np.linalg.norm(A, 2))):
        out.update(status="invalid_input", note="A not Hermitian")
        return out
    if r == 0.0:
        out.update(status="r_zero", t=c, lam=0.0, eta=0.0, L=c,
                   psd_margin_achieved=True, valid=True,
                   note="singleton ball: nominal objective; L=U=c")
        if want_witness:
            out["U"] = c
        return out

    t0 = perf_counter()
    AR, gR = real_embed(A, g)
    m = 2 * n
    lam = cp.Variable(nonneg=True)
    tvar = cp.Variable()
    br = cp.reshape(c - tvar - lam * r ** 2, (1, 1), order="F")
    M = cp.bmat([[AR + lam * np.eye(m), gR], [gR.T, br]])
    prob = cp.Problem(cp.Maximize(tvar), [M >> 0])
    out["build_s"] = perf_counter() - t0
    t1 = perf_counter()
    try:
        prob.solve(solver=solver)
    except Exception as exc:
        out.update(status=f"error:{exc}", wall_s=perf_counter() - t1,
                   note="solver exception")
        return out
    out["wall_s"] = perf_counter() - t1
    ss = prob.solver_stats
    out["status"] = prob.status
    out["solve_s"] = (float(ss.solve_time)
                      if getattr(ss, "solve_time", None) is not None
                      else out["wall_s"])
    if prob.status not in ("optimal", "optimal_inaccurate") \
            or tvar.value is None or lam.value is None \
            or not (np.isfinite(tvar.value) and np.isfinite(lam.value)):
        out["note"] = f"unusable solver output ({prob.status})"
        return out
    tv, lv = float(tvar.value), float(lam.value)
    out["t"], out["lam"] = tv, lv

    def lmi(tt, ll):
        Mv = np.zeros((m + 1, m + 1))
        Mv[:m, :m] = AR + ll * np.eye(m)
        Mv[:m, m:] = gR
        Mv[m:, :m] = gR.T
        Mv[m, m] = c - tt - ll * r ** 2
        return float(np.linalg.eigvalsh(Mv)[0]), float(np.linalg.norm(Mv, 2))

    lmin, mnorm = lmi(tv, lv)
    out["lam_min_lmi"] = lmin
    out["lmi_norm"] = mnorm
    tol_psd = ROUNDING_ALLOWANCE_K * (m + 1) * max(mnorm, 1.0) * 2.0 ** -52
    eta = max(0.0, -lmin) + tol_psd
    lam_corr = lv + eta
    t_cert = tv - eta * (1.0 + r * r)
    for _ in range(4):
        lmin_c, _ = lmi(t_cert, lam_corr)
        if lmin_c > 0.0:
            break
        eta *= 2.0
        lam_corr = lv + eta
        t_cert = tv - eta * (1.0 + r * r)
    else:
        lmin_c, _ = lmi(t_cert, lam_corr)
    # saved values now all correspond to the final re-checked matrix
    out.update(eta=eta, L=t_cert, lam_min_corrected=lmin_c,
               psd_margin_achieved=bool(lmin_c > 0.0),
               correction=eta * (1.0 + r * r))
    sc = max(1.0, abs(c), float(np.linalg.norm(g)), mnorm)
    lam_ok = lam_corr >= -1e-12 * sc        # corrected multiplier
    finite_ok = bool(np.isfinite(t_cert) and np.isfinite(lam_corr)
                     and np.isfinite(lmin_c))
    psd_ok = bool(lmin_c >= -tol_psd * 10)
    out["valid"] = bool(finite_ok and lam_ok and psd_ok
                        and prob.status == "optimal")
    if not out["valid"]:
        out["note"] = (f"validity failed: finite={finite_ok}, "
                       f"lam_corr>=0={lam_ok}, psd={psd_ok}, "
                       f"status={prob.status}")
    if want_witness and out["valid"]:
        wit = recover_witness(A, g, c, r, lam_corr)
        out["witness"] = wit
        out["U"] = wit["U"] if wit else None
    return out


def decide(L, U, tau, valid):
    """Corrected decision rule.  A negative LOWER bound never rejects."""
    if not valid or L is None or L != L:
        return "unresolved_invalid"
    if L > tau:
        return "certified_feasible"
    if U is not None and U == U and U < -tau:
        return "certified_infeasible"
    return "unresolved"


def sdp_user_status(res, tau):
    return decide(res.get("L"), res.get("U"), tau, res["valid"])


# ---------------------------------------------------------------------------
# test battery exposing the v1 error modes
# ---------------------------------------------------------------------------
def flip_instance():
    scal0 = dict(p=4.0, a=1.3, I0=0.5, alpha=1.0, beta=1.0, beta0=1.1)
    qg, x0, y0 = core.Q_approx(scal0, 0.5, 1.0, 1.5, "G")
    con = core.construct_v41(1.0, 1.0, 1.1, x0, y0, 1e-3, 4.0, 1.3, 0.5)
    s = core.spectral_summary(con["b"], con["h"], con["C"], 1.0, 1.5)
    return s, 0.5, 1.0, 1.5


def infeasible_instance():
    cands = ec.all_candidates(211)
    cand = cands[0]
    s = core.spectral_summary(cand["b"][0], cand["h"][0], cand["C"][0],
                              2.0, ec.NOISE_U)
    r = 0.30 * float(np.linalg.norm(cand["h"][0]))
    return s, r, 2.0, ec.NOISE_U


def run_tests():
    """Tests designed to catch the v1 'negative lower bound = infeasible'
    error and related failure modes.  Results to tests_sdp_v2.csv."""
    utils.ensure_dirs()
    rows = []

    def rec(name, expected, got, ok, note=""):
        rows.append(dict(test=name, expected=expected, got=got, ok=bool(ok),
                         note=note))
        print(f"[sdp2:{'PASS' if ok else 'FAIL'}] test {name}: {got} "
              f"(expected {expected}) {note}")
        if not ok:
            FAILURES.append(f"test {name}: got {got}, expected {expected}")

    # 1. valid negative lower bound, truly feasible problem: must NOT reject
    s, r, gamma, noise = flip_instance()
    res = sdp_qstar(s["A"], s["g"], s["c"], r)
    tau, _ = core.engineering_tau(s, r)
    qref = core.q_star_reference(s, r, gamma, noise, mu_rel_tol=opt.REF_TOL)
    loose_L = min(res.get("L") if res.get("L") is not None else 0.0, -0.5)
    st = decide(loose_L, None, tau, True)
    rec("neg_lb_feasible_not_rejected", "unresolved", st,
        st == "unresolved",
        note=f"q*={qref['q']:.4f}; v1 logic would have rejected on "
             f"L={loose_L}")
    st_u = decide(loose_L, res.get("U"), tau, True)
    rec("neg_lb_with_witness_no_false_reject", "not infeasible", st_u,
        st_u != "certified_infeasible")

    # 2. feasible witness with negative objective confirms infeasible
    s2, r2, g2, n2 = infeasible_instance()
    res2 = sdp_qstar(s2["A"], s2["g"], s2["c"], r2)
    tau2, _ = core.engineering_tau(s2, r2)
    q2 = core.q_star_reference(s2, r2, g2, n2, mu_rel_tol=opt.REF_TOL)
    st2 = sdp_user_status(res2, tau2)
    rec("witness_confirms_infeasible", "certified_infeasible", st2,
        st2 == "certified_infeasible" and q2["q"] < -tau2,
        note=f"U={res2.get('U')}, q*={q2['q']}")

    # 3. straddling bounds -> unresolved
    rec("straddling_unresolved", "unresolved", decide(-0.1, 0.1, 1e-8, True),
        decide(-0.1, 0.1, 1e-8, True) == "unresolved")

    # 4. near-threshold: search for the closest real instance, then verify
    # logical consistency; synthetic threshold checks cover the boundary
    best, bestd = None, np.inf
    cands = ec.all_candidates(212)
    for ci in range(3):
        for u in range(4):
            s4 = core.spectral_summary(cands[ci]["b"][u], cands[ci]["h"][u],
                                       cands[ci]["C"][u], 2.0, ec.NOISE_U)
            r4 = 0.30 * float(np.linalg.norm(cands[ci]["h"][u]))
            q4 = core.q_star_reference(s4, r4, 2.0, ec.NOISE_U,
                                       mu_rel_tol=opt.REF_TOL)
            d = abs(q4["q"]) / max(q4["tau"], 1e-30)
            if d < bestd:
                best, bestd = (s4, r4, q4), d
    s4, r4, q4 = best
    res4 = sdp_qstar(s4["A"], s4["g"], s4["c"], r4)
    tau4 = q4["tau"]
    st4 = sdp_user_status(res4, tau4)
    ok4 = ((st4 == "certified_feasible" and q4["q"] > tau4)
           or (st4 == "certified_infeasible" and q4["q"] < -tau4)
           or st4 == "unresolved")
    thr_ok = (decide(2 * tau4, None, tau4, True) == "certified_feasible"
              and decide(0.5 * tau4, -0.5 * tau4, tau4, True) == "unresolved"
              and decide(0.5 * tau4, -2 * tau4, tau4, True)
              == "certified_infeasible")
    rec("near_threshold_consistent", "consistent + threshold logic", st4,
        ok4 and thr_ok, note=f"|q*|/tau={bestd:.1f}; synthetic threshold "
                             f"checks {'ok' if thr_ok else 'FAILED'}")

    # 5. witness outside the ball must be projected and re-evaluated
    A5 = np.eye(3, dtype=complex)
    g5 = np.array([[3.0], [0.0], [0.0]], dtype=complex)
    c5 = -1.0
    r5 = 0.5
    wit = recover_witness(A5, g5, c5, r5, 0.0)
    e = np.array([[-0.5], [0.0], [0.0]], dtype=complex)  # projected point
    U_check = (float(np.real(e.conj().T @ A5 @ e)[0, 0])
               + 2 * float(np.real(e.conj().T @ g5)[0, 0]) + c5)
    rec("witness_projected", "projected, ball residual 0",
        f"projected={wit['projected']}, resid={wit['ball_residual']:.1e}",
        wit["projected"] and wit["ball_residual"] <= 1e-12
        and abs(wit["U"] - U_check) < 1e-12,
        note=f"raw norm {wit['e_norm_raw']:.2f} > r={r5}")

    # 6. invalid inputs
    Anan = np.eye(2, dtype=complex)
    Anan[0, 0] = np.nan
    r6 = sdp_qstar(Anan, np.zeros((2, 1), complex), 1.0, 0.5)
    rec("nonfinite_input_invalid", "invalid_input", r6["status"],
        r6["status"] == "invalid_input")
    Anon = np.array([[1.0, 2.0], [0.0, 1.0]], dtype=complex)
    r6b = sdp_qstar(Anon, np.zeros((2, 1), complex), 1.0, 0.5)
    rec("non_hermitian_invalid", "invalid_input", r6b["status"],
        r6b["status"] == "invalid_input")
    try:
        sdp_qstar(np.eye(2), np.zeros((2, 1)), 1.0, -0.1)
        rec("negative_radius_raises", "ValueError", "no error", False)
    except ValueError:
        rec("negative_radius_raises", "ValueError", "ValueError", True)
    rec("invalid_flag_unresolved", "unresolved_invalid",
        decide(None, None, 1e-8, False),
        decide(None, None, 1e-8, False) == "unresolved_invalid")

    # 7. scaling invariance of decisions (tau is scale-relative)
    s7, r7, g7, n7 = infeasible_instance()
    decisions = []
    for scale in (1e-3, 1.0, 1e3):
        tauS = 1e-8 * max(1.0, abs(scale * s7["c"]),
                           2 * r7 * float(np.linalg.norm(scale * s7["g"])),
                           r7 * r7 * float(np.linalg.norm(scale * s7["A"],
                                                          2)))
        resS = sdp_qstar(scale * s7["A"], scale * s7["g"], scale * s7["c"],
                         r7)
        decisions.append(sdp_user_status(resS, tauS))
    rec("scaling_invariant_decisions", "identical across 1e-3..1e3",
        "/".join(decisions), len(set(decisions)) == 1,
        note=f"scales tested 1e-3, 1, 1e3 -> {decisions}")

    utils.write_csv(os.path.join(OUT_TAG, "tests_sdp_v2.csv"),
                    list(rows[0].keys()), rows)
    return rows


def cand_h(cands, u):
    return cands[0]["h"][u]


def cand_C(cands, u):
    return cands[0]["C"][u]


# ---------------------------------------------------------------------------
# smallcheck (same pre-fixed selection as v1) and full batch
# ---------------------------------------------------------------------------
def smallcheck_instances():
    import experiments_v41 as ev41
    import sdp_baseline as v1

    inst = []
    scal0 = dict(p=4.0, a=1.3, I0=0.5, alpha=1.0, beta=1.0, beta0=1.1)
    qg, x0, y0 = core.Q_approx(scal0, 0.5, 1.0, 1.5, "G")
    for eps in (1e-2, 1e-3, 1e-4):
        con = core.construct_v41(1.0, 1.0, 1.1, x0, y0, eps, 4.0, 1.3, 0.5)
        s = core.spectral_summary(con["b"], con["h"], con["C"], 1.0, 1.5)
        inst.append(dict(tag=f"thm2_eps{eps:.0e}", s=s, r=0.5, gamma=1.0,
                         noise=1.5))
    for p, I0 in ev41.BOUNDARY_FOUR:
        bc = core.construct_boundary_v41(p, 1.0, I0, 1.0, 0.7, 1.0, 0.05, 0.2)
        s = core.spectral_summary(bc["b"], bc["h"], bc["C"], 1.0, 0.05)
        inst.append(dict(tag=f"thm3_p{p}_I0{I0}", s=s, r=0.2, gamma=1.0,
                         noise=0.05))
    inst += [it for it in v1.smallcheck_instances()
             if it["tag"].startswith(("alpha_zero", "beta_zero", "r_zero",
                                      "tr_"))]
    inst += [it for it in v1.smallcheck_instances()
             if it["tag"].startswith("pool_")]
    return inst


def run_smallcheck():
    utils.ensure_dirs()
    inst = smallcheck_instances()
    rows = []
    for it in inst:
        s, r, gamma, noise = it["s"], it["r"], it["gamma"], it["noise"]
        tau, _ = core.engineering_tau(s, r)
        ref = core.q_star_reference(s, r, gamma, noise, mu_rel_tol=opt.REF_TOL)
        res = sdp_qstar(s["A"], s["g"], s["c"], r)
        sc = max(1.0, abs(ref["q"]))
        rows.append(dict(
            tag=it["tag"], r=r, tau=tau, q_star_tr=ref["q"],
            tr_valid=ref["reference_valid"],
            sdp_status=res["status"], t=res["t"], L=res["L"], U=res["U"],
            eta=res["eta"], lam_min_lmi=res["lam_min_lmi"],
            lam_min_corrected=res["lam_min_corrected"],
            psd_margin=res["psd_margin_achieved"], valid=res["valid"],
            decision=sdp_user_status(res, tau),
            witness_projected=(res["witness"] or {}).get("projected"),
            witness_norm=(res["witness"] or {}).get("e_norm"),
            abs_diff=(res["t"] - ref["q"]) if res["t"] is not None
            else np.nan,
            rel_diff=((res["t"] - ref["q"]) / sc
                      if res["t"] is not None else np.nan),
            L_minus_qstar=(res["L"] - ref["q"]) if res["L"] is not None
            else np.nan,
            U_minus_qstar=(res["U"] - ref["q"]) if res["U"] is not None
            else np.nan))
    utils.write_csv(os.path.join(OUT_TAG, "smallcheck_sdp_v2.csv"),
                    list(rows[0].keys()), rows)
    solved = [r for r in rows if r["t"] is not None and r["tr_valid"]]
    max_abs = max(abs(r["abs_diff"]) for r in solved)
    max_L_ex = max(r["L_minus_qstar"] for r in solved
                   if r["L_minus_qstar"] == r["L_minus_qstar"])
    max_U_below = max(-r["U_minus_qstar"] for r in solved
                      if r["U_minus_qstar"] == r["U_minus_qstar"])
    check("smallcheck: SDP t matches reference", max_abs < 1e-4,
          f"max abs={max_abs:.3e}")
    check("smallcheck: L never above reference", max_L_ex <= 1e-8,
          f"max excess={max_L_ex:.3e}")
    check("smallcheck: U never below reference", max_U_below <= 1e-7,
          f"max deficit={max_U_below:.3e}")
    return rows


def run_batch():
    utils.ensure_dirs()
    rows = []
    for seed in VALIDATION_SEEDS:
        cands = ec.all_candidates(seed)
        for gamma in GAMMAS:
            prep = opt.prepare(cands, gamma)
            for rho in RHOS:
                for ci in range(len(cands)):
                    for u in range(len(cands[ci]["b"])):
                        pu = prep[ci][u]
                        r = rho * pu["hnorm"]
                        tau, _ = core.engineering_tau(pu["s"], r)
                        ref = core.q_star_reference(
                            pu["s"], r, gamma, ec.NOISE_U, eig=pu["eig"],
                            mu_rel_tol=opt.REF_TOL)
                        res = sdp_qstar(pu["s"]["A"], pu["s"]["g"],
                                        pu["s"]["c"], r)
                        sc = max(1.0, abs(ref["q"]))
                        rows.append(dict(
                            seed=seed, gamma=gamma, rho=rho, cand_idx=ci,
                            user=u, r=r, tau=tau, q_star_tr=ref["q"],
                            tr_valid=ref["reference_valid"],
                            tr_status=ec.user_status(ref),
                            sdp_status=res["status"], t=res["t"],
                            L=res["L"], U=res["U"], eta=res["eta"],
                            lam=res["lam"],
                            lam_min_corrected=res["lam_min_corrected"],
                            psd_margin=res["psd_margin_achieved"],
                            sdp_valid=res["valid"],
                            decision=sdp_user_status(res, tau),
                            witness_norm=(res["witness"] or {}).get("e_norm"),
                            witness_projected=(res["witness"] or {})
                            .get("projected"),
                            ball_residual=(res["witness"] or {})
                            .get("ball_residual"),
                            build_s=res["build_s"], solve_s=res["solve_s"],
                            wall_s=res["wall_s"],
                            abs_diff=(res["t"] - ref["q"]
                                      if res["t"] is not None else np.nan),
                            rel_diff=(res["t"] - ref["q"]) / sc
                            if res["t"] is not None else np.nan,
                            L_minus_qstar=(res["L"] - ref["q"]
                                           if res["L"] is not None
                                           else np.nan),
                            U_minus_qstar=(res["U"] - ref["q"]
                                           if res["U"] is not None
                                           else np.nan)))
        print(f"[sdp2:info] batch seed {seed} done", flush=True)
    utils.write_csv(os.path.join(OUT_TAG, "batch_sdp_v2_details.csv"),
                    list(rows[0].keys()), rows)
    okrows = [r for r in rows if r["tr_valid"]]
    max_abs = max(abs(r["abs_diff"]) for r in okrows
                  if r["abs_diff"] == r["abs_diff"])
    max_L = max(r["L_minus_qstar"] for r in okrows
                if r["L_minus_qstar"] == r["L_minus_qstar"])
    max_U_below = max(-r["U_minus_qstar"] for r in okrows
                      if r["U_minus_qstar"] == r["U_minus_qstar"])
    misaccept = sum(1 for r in okrows
                    if r["decision"] == "certified_feasible"
                    and r["tr_status"] != "clearly_feasible")
    misreject = sum(1 for r in okrows
                    if r["decision"] == "certified_infeasible"
                    and r["tr_status"] != "clearly_infeasible")
    bound_order_bad = sum(1 for r in okrows
                          if r["L"] is not None and r["U"] is not None
                          and r["L"] > r["U"] + 1e-7 * max(1.0, abs(r["L"])))
    ball_bad = sum(1 for r in rows if (r["ball_residual"] or 0) > 1e-9)
    dec = {}
    for r in rows:
        dec[r["decision"]] = dec.get(r["decision"], 0) + 1
    check("batch: L (certified lower bound) never above reference",
          max_L <= 1e-8, f"max excess={max_L:.3e}")
    check("batch: U (witness upper bound) never below reference",
          max_U_below <= 1e-6, f"max deficit={max_U_below:.3e}")
    check("batch: bound ordering L <= U", bound_order_bad == 0,
          f"violations={bound_order_bad}")
    check("batch: zero false acceptance", misaccept == 0)
    check("batch: zero false rejection", misreject == 0)
    check("batch: witnesses inside the ball", ball_bad == 0)
    print(f"[sdp2:info] batch n={len(rows)} decisions={dec}; max|t-q*|="
          f"{max_abs:.3e}; max L-q*={max_L:.3e}; max q*-U={max_U_below:.3e}")
    summary = dict(n=len(rows), decisions=dec, max_abs_diff=max_abs,
                   max_L_excess=max_L, max_U_deficit=max_U_below,
                   misaccept=misaccept, misreject=misreject,
                   bound_order_violations=bound_order_bad,
                   med_wall=float(np.nanmedian([r["wall_s"] for r in rows])))
    utils.write_json(os.path.join(OUT_TAG, "batch_sdp_v2_summary.json"),
                     summary)
    return rows, summary


def compare_with_v1(rows):
    """Instances whose v1 decision basis changed: v1 classified negative
    lower bounds as infeasible; v2 must re-decide them from the witness."""
    import csv

    path = os.path.join("results_v41_sdp", "batch_sdp_details.csv")
    if not os.path.exists(path):
        return []
    v1 = {(r["seed"], r["gamma"], r["rho"], r["cand_idx"], r["user"]): r
          for r in csv.DictReader(open(path, encoding="utf-8"))}

    def k5(r):
        # normalize to the string forms used in the v1 CSV
        return (f"{int(float(r['seed']))}", f"{float(r['gamma']):g}",
                f"{float(r['rho']):g}", f"{int(float(r['cand_idx']))}",
                f"{int(float(r['user']))}")

    out = []
    for r in rows:
        key = k5(r)
        old = v1.get(key)
        if old is None or old["sdp_status_user"] != "clearly_infeasible":
            continue
        out.append(dict(**{k: r[k] for k in ("seed", "gamma", "rho",
                                             "cand_idx", "user")},
                        v1_basis="negative lower bound (invalid rule)",
                        v1_status="clearly_infeasible",
                        v2_status=r["decision"],
                        v2_basis=("witness U=%.6g" % r["U"])
                        if r["U"] is not None else "no witness",
                        q_star_tr=r["q_star_tr"],
                        changed=r["decision"] != "certified_infeasible"))
    utils.write_csv(os.path.join(OUT_TAG, "decision_changes_v1_to_v2.csv"),
                    list(out[0].keys()) if out else ["seed"], out)
    n_conf = sum(1 for o in out if not o["changed"])
    print(f"[sdp2:info] v1 negative-LB 'infeasible' instances: {len(out)}; "
          f"v2 re-confirmed via witness: {n_conf}; changed to unresolved: "
          f"{len(out) - n_conf}")
    return out


# ---------------------------------------------------------------------------
# pipelines: SDP-only (unresolved allowed) and explicitly named
# SDP + TR-completion variants for the full-decision comparison
# ---------------------------------------------------------------------------
def cascade_sdp_v2(cands, prep, gamma, stages, rhos=RHOS,
                   tr_completion=False):
    """Screening as in optimization_v41.cascade; fallback = corrected SDP
    (L/U decision).  With tr_completion=True, users left unresolved by the
    SDP are finished by the trust-region reference; those TR calls are
    counted separately in stats['tr_completion_calls'] and the pipeline is
    explicitly named '..._TR' — never presented as pure SDP."""
    t0 = perf_counter()
    n_users = len(cands[0]["b"])
    stats = dict(stage_evals={k: 0 for k in stages if k != "__none__"},
                 stage_pass={k: 0 for k in stages if k != "__none__"},
                 sdp_calls=0, tr_completion_calls=0, invalid=0,
                 build_s=0.0, solve_s=0.0, wall_s=0.0)
    decisions = []
    for ci, cand in enumerate(cands):
        for rho in rhos:
            certified = [False] * n_users
            for stage in stages:
                if stage == "__none__":
                    continue
                for u in range(n_users):
                    if certified[u]:
                        continue
                    r = rho * prep[ci][u]["hnorm"]
                    tau, _ = core.engineering_tau(prep[ci][u]["s"], r)
                    lb = opt.bound_value(prep[ci][u], r, gamma, stage)
                    stats["stage_evals"][stage] += 1
                    if lb > tau:
                        certified[u] = True
                        stats["stage_pass"][stage] += 1
            if all(certified):
                decisions.append("accepted")
                continue
            statuses = []
            for u in range(n_users):
                if certified[u]:
                    statuses.append("certified_feasible")
                    continue
                pu = prep[ci][u]
                r = rho * pu["hnorm"]
                res = sdp_qstar(pu["s"]["A"], pu["s"]["g"], pu["s"]["c"], r)
                stats["sdp_calls"] += 1
                for key in ("build_s", "solve_s", "wall_s"):
                    v = res[key]
                    stats[key] += v if v == v else 0.0
                tau, _ = core.engineering_tau(pu["s"], r)
                st = sdp_user_status(res, tau)
                if st == "unresolved_invalid":
                    stats["invalid"] += 1
                if st == "unresolved" and tr_completion:
                    ref = core.q_star_reference(pu["s"], r, gamma,
                                                ec.NOISE_U,
                                                eig=pu["eig"],
                                                mu_rel_tol=opt.REF_TOL)
                    stats["tr_completion_calls"] += 1
                    st = ec.user_status(ref)
                statuses.append(st)
                if st == "certified_infeasible" \
                        or st == "clearly_infeasible":
                    break
            if any(x in ("unresolved_invalid", "invalid") for x in statuses):
                decisions.append("unresolved_invalid")
            elif "certified_infeasible" in statuses \
                    or "clearly_infeasible" in statuses:
                decisions.append("rejected")
            elif all(x in ("certified_feasible", "clearly_feasible")
                     for x in statuses):
                decisions.append("accepted")
            else:
                decisions.append("unresolved")
    return dict(seconds=perf_counter() - t0, decisions=decisions,
                stats=stats)


V2_PIPELINES = {
    "TR_pure": dict(stages=(), sdp=False, trcomp=False),
    "SDP_only": dict(stages=("__none__",), sdp=True, trcomp=False),
    "SDP_plus_TR": dict(stages=("__none__",), sdp=True, trcomp=True),
    "Q_old_SDP": dict(stages=("Q_old",), sdp=True, trcomp=False),
    "Q_G_SDP": dict(stages=("Q_G",), sdp=True, trcomp=False),
    "Lq_QG_SDP": dict(stages=("L_quad", "Q_G"), sdp=True, trcomp=False),
    "Lq_QG_SDP_TR": dict(stages=("L_quad", "Q_G"), sdp=True, trcomp=True),
}


def run_pipelines(reps=3):
    timing_rows = []
    for seed in VALIDATION_SEEDS:
        cands = ec.all_candidates(seed)
        for gamma in GAMMAS:
            for p, cfgp in V2_PIPELINES.items():
                pw = opt.prepare(cands[:1], gamma)
                if cfgp["sdp"]:
                    cascade_sdp_v2(cands[:1], pw, gamma, cfgp["stages"],
                                   tr_completion=cfgp["trcomp"])
                else:
                    opt.cascade(cands[:1], pw, gamma, ())     # warm-up TR
            runs = {p: [] for p in V2_PIPELINES}
            for rep in range(reps):
                order = np.roll(np.array(list(V2_PIPELINES)), rep)
                for p in order:
                    cfgp = V2_PIPELINES[p]
                    t0 = perf_counter()
                    prep = opt.prepare(cands, gamma)
                    if cfgp["sdp"]:
                        out = cascade_sdp_v2(cands, prep, gamma,
                                             cfgp["stages"],
                                             tr_completion=cfgp["trcomp"])
                    else:
                        out = opt.cascade(cands, prep, gamma, ())
                    runs[p].append((perf_counter() - t0, out))
            base = runs["TR_pure"][0][1]["decisions"]
            for p in V2_PIPELINES:
                secs = [t for t, _ in runs[p]]
                o0 = runs[p][0][1]
                st = o0["stats"]
                dec = o0["decisions"]
                agree = sum(1 for a, b in zip(dec, base) if a == b)
                complete = all(d in ("accepted", "rejected")
                               for d in dec)
                if complete and p != "TR_pure":
                    if dec != base:
                        FAILURES.append(f"{p} seed {seed} gamma {gamma}: "
                                        "complete-decision mismatch vs TR")
                timing_rows.append(dict(
                    seed=seed, gamma=gamma, pipeline=p,
                    seconds_med=float(np.median(secs)),
                    seconds_min=min(secs), seconds_max=max(secs),
                    sdp_calls=st.get("sdp_calls",
                                     st.get("ref_calls", 0)),
                    tr_completion_calls=st.get("tr_completion_calls", 0),
                    invalid=st.get("invalid", st.get("invalid_refs", 0)),
                    sdp_build_s=st.get("build_s", np.nan),
                    sdp_solve_s=st.get("solve_s", np.nan),
                    unresolved=sum(1 for d in dec
                                   if d.startswith("unresolved")),
                    complete_decision=bool(complete),
                    agreement_vs_TR=f"{agree}/{len(base)}",
                    accepted=sum(1 for d in dec if d == "accepted")))
        print(f"[sdp2:info] pipelines seed {seed} done", flush=True)
    utils.write_csv(os.path.join(OUT_TAG, "pipeline_timings_v2.csv"),
                    list(timing_rows[0].keys()), timing_rows)
    by = {}
    for r in timing_rows:
        by.setdefault(r["pipeline"], []).append(r)
    for p, rs in sorted(by.items()):
        med = float(np.median([r["seconds_med"] for r in rs]))
        unres = sum(r["unresolved"] for r in rs)
        trc = sum(r["tr_completion_calls"] for r in rs)
        print(f"[sdp2:info] {p}: median {med:.4f} s; unresolved units "
              f"{unres}; TR-completion calls {trc}; SDP calls "
              f"{sum(r['sdp_calls'] for r in rs)}")
    return timing_rows


# ---------------------------------------------------------------------------
# figures + driver
# ---------------------------------------------------------------------------
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

    rows = rd("batch_sdp_v2_details.csv")
    qs = np.array([float(r["q_star_tr"]) for r in rows])
    Ls = np.array([float(r["L"]) for r in rows])
    Us = np.array([float(r["U"]) for r in rows])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.6))
    lo, hi = min(qs.min(), Ls.min()), max(qs.max(), Us.max())
    m = 0.05 * (hi - lo + 1e-3)
    ax1.plot([lo - m, hi + m], [lo - m, hi + m], color=GRAY, lw=0.8)
    ax1.scatter(qs, Us, s=10, color=PAL[1], marker="s", zorder=2,
                label="witness upper bound $U$")
    ax1.scatter(qs, Ls, s=10, color=PAL[0], zorder=3,
                label="certified lower bound $L$")
    ax1.set_xlabel("trust-region reference $q^*$")
    ax1.set_ylabel("SDP bounds")
    ax1.set_title("(a) two-sided bounds vs reference")
    ax1.legend(loc="upper left", fontsize=8)
    gap = Us - Ls
    ax2.hist(np.log10(np.maximum(gap, 1e-18)), bins=40, color=PAL[2])
    ax2.set_xlabel(r"$\log_{10}(U-L)$")
    ax2.set_ylabel("user problems")
    ax2.set_title(f"(b) bound gap (median {np.median(gap):.1e})")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "fig_sdp_v2_bounds.pdf"))
    fig.savefig(os.path.join(figdir, "fig_sdp_v2_bounds.png"), dpi=150)
    plt.close(fig)
    print(f"[sdp2:info] figures written to {figdir}")


def run_all(quick=False):
    utils.ensure_dirs()
    here = os.path.dirname(os.path.abspath(__file__))

    def sha(f):
        with open(os.path.join(here, f), "rb") as fh:
            return hashlib.sha256(fh.read()).hexdigest()

    import cvxpy as cp
    cfg = dict(
        written_at=datetime.now().isoformat(timespec="seconds"),
        correction_note="v2: negative lower bounds never reject; "
                        "infeasibility requires witness upper bound; "
                        "v1 results retained in results_v41_sdp/ as history",
        solver=SOLVER, rounding_allowance_k=ROUNDING_ALLOWANCE_K,
        cvxpy_version=cp.__version__,
        selection_rules="identical to v1 sdp_config.json (pre-fixed)",
        code_hashes={f: sha(f) for f in ("core.py", "sdp_baseline.py",
                                         "sdp_baseline_v2.py",
                                         "optimization_v41.py")})
    utils.write_json(os.path.join(OUT_TAG, "sdp_config_v2.json"), cfg)

    print("--- v2 tests exposing the old error modes ---")
    run_tests()
    print("--- v2 smallcheck ---")
    sc = run_smallcheck()
    if quick:
        print("[sdp2:info] quick mode stops after tests + smallcheck")
        return not FAILURES
    print("--- v2 full batch (2160 user problems) ---")
    rows, summary = run_batch()
    print("--- v1 -> v2 decision-basis comparison ---")
    changes = compare_with_v1(rows)
    print("--- v2 pipelines (SDP-only and SDP+TR completion) ---")
    tim = run_pipelines()
    print("--- v2 figures ---")
    make_figures()
    manifest = dict(finished_at=datetime.now().isoformat(timespec="seconds"),
                    solver=SOLVER, failures=list(FAILURES),
                    batch_summary=summary,
                    v1_negative_lb_instances=len(changes))
    utils.write_json(os.path.join(OUT_TAG, "run_manifest.json"), manifest)
    print(f"\n=== SDP v2 round {'FAILED' if FAILURES else 'PASSED'} ===")
    for f in FAILURES:
        print("  FAIL:", f)
    return not FAILURES
