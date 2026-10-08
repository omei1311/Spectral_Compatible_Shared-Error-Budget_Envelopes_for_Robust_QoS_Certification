# -*- coding: utf-8 -*-
"""S-procedure-SDP baseline for the V4.1 robust QoS certification round.

Baseline provenance: the standard single-sphere S-procedure LMI adapted to
the fixed-candidate certification problem of this project (shared error
ball, one QCQP per user).  This is a classical technique, NOT a
reproduction of the full energy-efficient robust beamforming algorithm of
Zheng et al., IEEE Trans. Commun. 73(8):6306-6321, 2025
(DOI 10.1109/TCOMM.2024.3516491, arXiv:2412.08032), and it is not a "new
2025 S-procedure".  That paper is cited as representative context for
S-procedure usage in RIS robust design.

For one user problem (identical data to the existing pipelines):
    q(e) = e^H A e + 2 Re(e^H g) + c,   ||e||_2 <= r,
    A = b b^H - gamma C C^H,  g = A h_hat,  c = h^H A h - gamma sigma^2,
the SDP is
    maximize t  s.t.  lambda >= 0,
    [[A + lambda I, g], [g^H, c - t - lambda r^2]] >= 0,
whose optimum equals q*(r) by the S-lemma (single sphere).  Solved on the
real embedding [[Re A, -Im A],[Im A, Re A]] with CLARABEL (fallback SCS),
no warm start.

Residual-corrected certified lower bound (derivation): let M(t,lam)
denote the LMI matrix and
    eta = max(0, -lambda_min(M)) + K*(2n+1)*||M||_2*2^-52
(K = 10 rounding allowance).  The shifted point
    (t', lam') = (t - eta*(1+r^2), lam + eta)
satisfies exactly
    M(t', lam') = M(t, lam) + blockdiag(eta*I, eta)
(top-left gains eta*I from lam+eta; the (2,2) entry changes by
+eta(1+r^2) - eta*r^2 = +eta), so by Weyl's inequality
    lambda_min(M') >= lambda_min(M) + eta >= K*(2n+1)*||M||*2^-52 > 0
up to the eigensolver's own rounding, which is re-checked numerically
after the correction (one retry with doubled eta if needed).  The value
    t_cert = t - eta*(1+r^2)
is a dual-feasible certified lower bound in floating point (residual
re-checked), not an interval-arithmetic formal certificate.  Uncorrected
t is only a numerical estimate and is never used for acceptance.

All results go to results_v41_sdp/; earlier result directories untouched.
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

OUT_TAG = "results_v41_sdp"
GAMMAS = (0.5, 1.0, 2.0)
RHOS = (0.06, 0.15, 0.30)
VALIDATION_SEEDS = (211, 212, 213, 214, 215)
SOLVER = "CLARABEL"
ROUNDING_ALLOWANCE_K = 10.0

FAILURES = []


def check(name, ok, detail=""):
    print(f"[sdp:{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        FAILURES.append(f"{name}: {detail}")
    return ok


# ---------------------------------------------------------------------------
# solver wrapper
# ---------------------------------------------------------------------------
def real_embed(A, g):
    """Complex Hermitian A, complex vector g -> real symmetric embedding."""
    n = A.shape[0]
    AR = np.block([[A.real, -A.imag], [A.imag, A.real]])
    gR = np.concatenate([g.real.ravel(), g.imag.ravel()]).reshape(-1, 1)
    return AR, gR


def sdp_qstar(A, g, c, r, solver=SOLVER):
    """Solve the S-procedure SDP for one user problem.

    Returns dict(status, t, lam, t_cert, lam_min_lmi, lam_min_corrected,
    build_s, solve_s, wall_s, solver, dual_gap_or_nan, valid, note).
    r == 0 short-circuits to the nominal objective value c.
    """
    import cvxpy as cp

    n = A.shape[0]
    out = dict(status=None, t=None, lam=None, t_cert=None,
               lam_min_lmi=None, lam_min_corrected=None, build_s=np.nan,
               solve_s=np.nan, wall_s=np.nan, solver=solver,
               dual_gap=np.nan, valid=False, note="")
    if r <= 0.0:
        out.update(status="r_zero", t=c, lam=0.0, t_cert=c,
                   lam_min_lmi=np.nan, lam_min_corrected=np.nan, valid=True,
                   note="singleton ball: nominal objective")
        return out

    t0 = perf_counter()
    AR, gR = real_embed(np.asarray(A), np.asarray(g))
    m = 2 * n
    lam = cp.Variable(nonneg=True)
    tvar = cp.Variable()
    br = cp.reshape(c - tvar - lam * r ** 2, (1, 1), order="F")
    M = cp.bmat([[AR + lam * np.eye(m), gR],
                 [gR.T, br]])
    prob = cp.Problem(cp.Maximize(tvar), [M >> 0])
    out["build_s"] = perf_counter() - t0
    t1 = perf_counter()
    try:
        prob.solve(solver=solver)
    except Exception as exc:
        out.update(status=f"error:{exc}", wall_s=perf_counter() - t1)
        out["note"] = "solver exception"
        return out
    out["wall_s"] = perf_counter() - t1
    ss = prob.solver_stats
    out["status"] = prob.status
    if hasattr(ss, "solve_time") and ss.solve_time is not None:
        out["solve_s"] = float(ss.solve_time)
    else:
        out["solve_s"] = out["wall_s"]
    if prob.status not in ("optimal", "optimal_inaccurate") \
            or tvar.value is None or lam.value is None:
        out["note"] = f"solver status {prob.status}"
        return out
    tv, lv = float(tvar.value), float(lam.value)

    def lmi_min(tt, ll):
        Mv = np.zeros((m + 1, m + 1))
        Mv[:m, :m] = AR + ll * np.eye(m)
        Mv[:m, m:] = gR
        Mv[m:, :m] = gR.T
        Mv[m, m] = c - tt - ll * r ** 2
        return float(np.linalg.eigvalsh(Mv)[0]), float(np.linalg.norm(Mv, 2))

    lmin, mnorm = lmi_min(tv, lv)
    out.update(t=tv, lam=lv, lam_min_lmi=lmin)
    # residual-corrected certified bound (derivation in module docstring)
    eta = max(0.0, -lmin) + ROUNDING_ALLOWANCE_K * (m + 1) * mnorm * 2.0 ** -52
    for _ in range(2):
        lmin_c, _ = lmi_min(tv - eta * (1 + r * r), lv + eta)
        if lmin_c >= -ROUNDING_ALLOWANCE_K * (m + 1) * max(mnorm, 1.0) \
                * 2.0 ** -52 * 10:
            break
        eta *= 2.0
    out["t_cert"] = tv - eta * (1.0 + r * r)
    out["lam_min_corrected"] = lmin_c
    out["correction"] = eta * (1.0 + r * r)
    sc = max(1.0, abs(c), float(np.linalg.norm(g)), mnorm)
    out["valid"] = bool(lv >= -1e-8 * sc
                        and lmin_c >= -1e-9 * max(1.0, mnorm)
                        and prob.status == "optimal")
    if prob.status == "optimal_inaccurate":
        out["note"] = "optimal_inaccurate"
    return out


def sdp_user_status(res, tau):
    if not res["valid"]:
        return "sdp_invalid"
    if res["t_cert"] > tau:
        return "clearly_feasible"
    if res["t_cert"] < -tau:
        return "clearly_infeasible"
    return "ambiguous"


# ---------------------------------------------------------------------------
# pre-fixed small-scale verification set
# ---------------------------------------------------------------------------
def smallcheck_instances():
    """Selection rules fixed independently of SDP behaviour."""
    import experiments_v41 as ev41

    inst = []
    # Theorem-2 construction family (flip parameters), three epsilon levels
    scal0 = dict(p=4.0, a=1.3, I0=0.5, alpha=1.0, beta=1.0, beta0=1.1)
    qg, x0, y0 = core.Q_approx(scal0, 0.5, 1.0, 1.5, "G")
    for eps in (1e-2, 1e-3, 1e-4):
        con = core.construct_v41(1.0, 1.0, 1.1, x0, y0, eps, 4.0, 1.3, 0.5)
        s = core.spectral_summary(con["b"], con["h"], con["C"], 1.0, 1.5)
        inst.append(dict(tag=f"thm2_eps{eps:.0e}", s=s, r=0.5, gamma=1.0,
                         noise=1.5))
    # Theorem-3 boundary four instances (both c_parallel signs)
    for p, I0 in ev41.BOUNDARY_FOUR:
        bc = core.construct_boundary_v41(p, 1.0, I0, 1.0, 0.7, 1.0, 0.05, 0.2)
        s = core.spectral_summary(bc["b"], bc["h"], bc["C"], 1.0, 0.05)
        inst.append(dict(tag=f"thm3_p{p}_I0{I0}", s=s, r=0.2, gamma=1.0,
                         noise=0.05))
    # degenerate: alpha=0, beta=0, r=0
    CH = np.array([[0.0, 0.9]], dtype=complex)
    h = np.array([1.2, 0.5], dtype=complex).reshape(-1, 1)
    b = np.array([np.sqrt(1.7), 0.0], dtype=complex).reshape(-1, 1)
    inst.append(dict(tag="alpha_zero", s=core.spectral_summary(
        b, h, CH.conj().T, 1.0, 0.05), r=0.4, gamma=1.0, noise=0.05))
    CH = np.array([[1.1, 0.0]], dtype=complex)
    h = np.array([0.9, 0.0], dtype=complex).reshape(-1, 1)
    b = np.array([np.sqrt(2.0), 0.0], dtype=complex).reshape(-1, 1)
    inst.append(dict(tag="beta_zero", s=core.spectral_summary(
        b, h, CH.conj().T, 1.0, 0.05), r=0.4, gamma=1.0, noise=0.05))
    inst.append(dict(tag="r_zero", s=inst[-1]["s"], r=0.0, gamma=1.0,
                     noise=0.05))
    # trust-region special structures
    special = [
        ("tr_psd_interior", np.eye(3), np.array([[0.1], [0.2], [0.0]]), 5.0),
        ("tr_singular_solvable", np.diag([1.0, 0.0]),
         np.array([[1.0], [0.0]]), 5.0),
        ("tr_hard_equality", np.diag([-1.0, 1.0]),
         np.array([[1.0], [0.0]]), 0.5),
        ("tr_repeated_min", np.diag([-1.0, -1.0, 2.0]),
         np.array([[0.0], [0.0], [1.0]]), 0.7),
    ]
    for tag, A, g, r in special:
        s = dict(A=np.asarray(A, dtype=complex), g=np.asarray(g, complex),
                 c=0.0)
        inst.append(dict(tag=tag, s=s, r=r, gamma=1.0, noise=0.0))
    # pool instances: first candidate of each validation seed, rho=0.15,
    # all gamma (fixed rule, spans feasible/critical/infeasible)
    for seed in VALIDATION_SEEDS:
        cands = ec.all_candidates(seed)
        cand = cands[0]
        for gamma in GAMMAS:
            for u in range(4):
                s = core.spectral_summary(cand["b"][u], cand["h"][u],
                                          cand["C"][u], gamma, ec.NOISE_U)
                r = 0.15 * float(np.linalg.norm(cand["h"][u]))
                inst.append(dict(tag=f"pool_s{seed}_g{gamma}_u{u}", s=s, r=r,
                                 gamma=gamma, noise=ec.NOISE_U))
    return inst


def run_smallcheck(write_config_first=True):
    utils.ensure_dirs()
    inst = smallcheck_instances()
    if write_config_first:
        here = os.path.dirname(os.path.abspath(__file__))

        def sha(f):
            with open(os.path.join(here, f), "rb") as fh:
                return hashlib.sha256(fh.read()).hexdigest()
        import cvxpy as cp
        cfg = dict(
            written_at=datetime.now().isoformat(timespec="seconds"),
            purpose="small-scale SDP verification set; selection rules fixed "
                    "independently of SDP behaviour",
            selection_rules=dict(
                theorem2="flip family at the Q_G minimizer, eps in "
                         "{1e-2,1e-3,1e-4}",
                theorem3="four fixed (p,I0) boundary instances",
                degenerate="alpha=0, beta=0, r=0 singletons",
                trust_region="PSD interior, singular, hard-case equality, "
                             "repeated minimum eigenvalue",
                pool="first candidate of each validation seed 211-215, "
                     "rho=0.15, all gamma, all users"),
            solver=SOLVER, rounding_allowance_k=ROUNDING_ALLOWANCE_K,
            cvxpy_version=cp.__version__,
            code_hashes={f: sha(f) for f in
                         ("core.py", "sdp_baseline.py",
                          "optimization_v41.py")},
        )
        utils.write_json(os.path.join(OUT_TAG, "sdp_config.json"), cfg)
    rows = []
    for it in inst:
        s, r, gamma, noise = it["s"], it["r"], it["gamma"], it["noise"]
        ref = core.q_star_reference(s, r, gamma, noise,
                                    mu_rel_tol=opt.REF_TOL)
        res = sdp_qstar(s["A"], s["g"], s["c"], r)
        sc = max(1.0, abs(ref["q"]))
        rows.append(dict(
            tag=it["tag"], r=r, q_star_tr=ref["q"], tr_valid=ref[
                "reference_valid"],
            sdp_status=res["status"], t=res["t"], lam=res["lam"],
            t_cert=res["t_cert"], correction=res.get("correction", np.nan),
            lam_min_lmi=res["lam_min_lmi"],
            lam_min_corrected=res["lam_min_corrected"],
            valid=res["valid"], build_s=res["build_s"],
            solve_s=res["solve_s"], wall_s=res["wall_s"],
            abs_diff=(res["t"] - ref["q"]) if res["t"] is not None
            else np.nan,
            rel_diff=(res["t"] - ref["q"]) / sc if res["t"] is not None
            else np.nan,
            cert_minus_qstar=(res["t_cert"] - ref["q"])
            if res["t_cert"] is not None else np.nan))
    utils.write_csv(os.path.join(OUT_TAG, "smallcheck_sdp.csv"),
                    list(rows[0].keys()), rows)
    solved = [r for r in rows if r["t"] is not None and r["tr_valid"]]
    max_abs = max(abs(r["abs_diff"]) for r in solved)
    max_rel = max(abs(r["rel_diff"]) for r in solved)
    worst_cert = max(r["cert_minus_qstar"] for r in solved)
    med_wall = float(np.nanmedian([r["wall_s"] for r in rows]))
    check("smallcheck: all instances solved or r_zero",
          all(r["t"] is not None for r in rows))
    check("smallcheck: SDP t matches trust-region reference",
          max_abs < 1e-4 * max(1.0, max(abs(r["q_star_tr"])
                                        for r in solved)),
          f"max abs diff={max_abs:.3e}, max rel={max_rel:.3e}")
    check("smallcheck: corrected bound never exceeds reference",
          worst_cert <= 1e-8 * max(1.0, max(abs(r["q_star_tr"])
                                            for r in solved)),
          f"max excess={worst_cert:.3e}")
    print(f"[sdp:info] smallcheck n={len(rows)}, median wall {med_wall:.4f} s"
          f"/solve; projected 2160-user batch ~{med_wall * 2160:.0f} s "
          f"per pass")
    return rows, med_wall


# ---------------------------------------------------------------------------
# full validation batch: user-level SDP vs trust-region reference
# ---------------------------------------------------------------------------
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
                            user=u, r=r, tau=tau,
                            q_star_tr=ref["q"], tr_valid=ref["reference_valid"],
                            tr_status=ec.user_status(ref),
                            sdp_status=res["status"], t=res["t"],
                            lam=res["lam"], t_cert=res["t_cert"],
                            correction=res.get("correction", np.nan),
                            lam_min_lmi=res["lam_min_lmi"],
                            lam_min_corrected=res["lam_min_corrected"],
                            sdp_valid=res["valid"],
                            sdp_status_user=sdp_user_status(res, tau),
                            build_s=res["build_s"], solve_s=res["solve_s"],
                            wall_s=res["wall_s"],
                            abs_diff=(res["t"] - ref["q"]
                                      if res["t"] is not None else np.nan),
                            rel_diff=(res["t"] - ref["q"]) / sc
                            if res["t"] is not None else np.nan))
        print(f"[sdp:info] batch seed {seed} done", flush=True)
    utils.write_csv(os.path.join(OUT_TAG, "batch_sdp_details.csv"),
                    list(rows[0].keys()), rows)
    solved = [r for r in rows if r["t"] is not None and r["tr_valid"]]
    max_abs = max(abs(r["abs_diff"]) for r in solved)
    max_rel = max(abs(r["rel_diff"]) for r in solved)
    worst_cert = max(r["t_cert"] - r["q_star_tr"] for r in solved
                     if r["t_cert"] is not None)
    agree = sum(1 for r in rows
                if r["sdp_status_user"] == r["tr_status"])
    check("batch: all SDP solves valid",
          all(r["sdp_valid"] for r in rows))
    check("batch: SDP t matches trust-region reference",
          max_abs < 1e-5, f"max abs={max_abs:.3e}, max rel={max_rel:.3e}")
    check("batch: corrected bound never exceeds reference",
          worst_cert <= 1e-7, f"max excess={worst_cert:.3e}")
    print(f"[sdp:info] batch n={len(rows)}; status agreement with TR: "
          f"{agree}/{len(rows)} ({100 * agree / len(rows):.2f}%); median "
          f"wall {np.nanmedian([r['wall_s'] for r in rows]):.4f} s")
    summary = dict(n=len(rows), max_abs_diff=max_abs, max_rel_diff=max_rel,
                   worst_cert_excess=worst_cert,
                   status_agreement=agree,
                   med_wall=float(np.nanmedian([r['wall_s'] for r in rows])),
                   med_build=float(np.nanmedian([r['build_s']
                                                 for r in rows])),
                   med_solve_reported=float(np.nanmedian(
                       [r['solve_s'] for r in rows])))
    utils.write_json(os.path.join(OUT_TAG, "batch_sdp_summary.json"),
                     summary)
    return rows, summary


# ---------------------------------------------------------------------------
# pipelines with SDP fallback (same final decision task; caching identical)
# ---------------------------------------------------------------------------
def cascade_sdp(cands, prep, gamma, stages, rhos=RHOS):
    """Screening stages as in optimization_v41.cascade, with the
    S-procedure-SDP (corrected certified bound) as the fallback stage."""
    t0 = perf_counter()
    n_users = len(cands[0]["b"])
    stats = dict(stage_evals={k: 0 for k in stages},
                 stage_pass={k: 0 for k in stages}, sdp_calls=0,
                 invalid=0, build_s=0.0, solve_s=0.0, wall_s=0.0)
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
                    statuses.append("clearly_feasible")
                    continue
                pu = prep[ci][u]
                res = sdp_qstar(pu["s"]["A"], pu["s"]["g"], pu["s"]["c"],
                                rho * pu["hnorm"])
                stats["sdp_calls"] += 1
                for key in ("build_s", "wall_s"):
                    v = res[key]
                    stats[key] += v if v == v else 0.0
                v = res["solve_s"] if res["solve_s"] == res["solve_s"] \
                    else res["wall_s"]
                stats["solve_s"] += v if v == v else 0.0
                tau, _ = core.engineering_tau(pu["s"], rho * pu["hnorm"])
                st = sdp_user_status(res, tau)
                if st == "sdp_invalid":
                    stats["invalid"] += 1
                statuses.append(st)
                if st == "clearly_infeasible":
                    break
            if "sdp_invalid" in statuses:
                decisions.append("unresolved_invalid")
            elif "clearly_infeasible" in statuses:
                decisions.append("rejected")
            elif all(x == "clearly_feasible" for x in statuses):
                decisions.append("accepted")
            else:
                decisions.append("unresolved")
    return dict(seconds=perf_counter() - t0, decisions=decisions,
                stats=stats)


SDP_PIPELINES = {
    "TR_pure": (),
    "SDP_pure": ("__none__",),          # every user goes to the SDP stage
    "Q_old_SDP": ("Q_old",),
    "Q_G_SDP": ("Q_G",),
    "Lq_QG_SDP": ("L_quad", "Q_G"),
}


def run_pipelines(reps=3):
    timing_rows = []
    for seed in VALIDATION_SEEDS:
        cands = ec.all_candidates(seed)
        for gamma in GAMMAS:
            for p in SDP_PIPELINES:
                prep_w = opt.prepare(cands[:1], gamma)
                if p == "TR_pure":
                    opt.cascade(cands[:1], prep_w, gamma, ())      # warm-up
                else:
                    cascade_sdp(cands[:1], prep_w, gamma,
                                SDP_PIPELINES[p])                  # warm-up
            runs = {p: [] for p in SDP_PIPELINES}
            for rep in range(reps):
                order = np.roll(np.array(list(SDP_PIPELINES)), rep)
                for p in order:
                    t0 = perf_counter()
                    prep = opt.prepare(cands, gamma)
                    if p == "TR_pure":
                        out = opt.cascade(cands, prep, gamma, ())
                    else:
                        out = cascade_sdp(cands, prep, gamma,
                                          SDP_PIPELINES[p])
                    runs[p].append((perf_counter() - t0, out))
            base = runs["TR_pure"][0][1]["decisions"]
            for p in SDP_PIPELINES:
                for _, out in runs[p]:
                    pass  # agreement reported, near-threshold divergence
                    # between SDP and TR numerics is expected and counted
                secs = [t for t, _ in runs[p]]
                o0 = runs[p][0][1]
                agree = sum(1 for a, b in zip(o0["decisions"], base)
                            if a == b)
                st = o0["stats"]
                timing_rows.append(dict(
                    seed=seed, gamma=gamma, pipeline=p,
                    seconds_med=float(np.median(secs)),
                    seconds_min=min(secs), seconds_max=max(secs),
                    sdp_calls=st.get("sdp_calls",
                                     st.get("ref_calls", 0)),
                    invalid=st.get("invalid", st.get("invalid_refs", 0)),
                    sdp_build_s=st.get("build_s", np.nan),
                    sdp_solve_s=st.get("solve_s", np.nan),
                    agreement_vs_TR=f"{agree}/{len(base)}",
                    accepted=sum(1 for d in o0["decisions"]
                                 if d == "accepted")))
        print(f"[sdp:info] pipelines seed {seed} done", flush=True)
    utils.write_csv(os.path.join(OUT_TAG, "pipeline_timings.csv"),
                    list(timing_rows[0].keys()), timing_rows)
    by = {}
    for r in timing_rows:
        by.setdefault(r["pipeline"], []).append(r)
    for p, rs in sorted(by.items()):
        med = float(np.median([r["seconds_med"] for r in rs]))
        print(f"[sdp:info] {p}: median {med:.4f} s per (seed,gamma); "
              f"SDP calls {sum(r['sdp_calls'] for r in rs)}; agreement "
              f"min {min(int(r['agreement_vs_TR'].split('/')[0]) for r in rs)}"
              f"/{int(rs[0]['agreement_vs_TR'].split('/')[1])}")
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

    rows = rd("batch_sdp_details.csv")
    qs = np.array([float(r["q_star_tr"]) for r in rows])
    tt = np.array([float(r["t"]) for r in rows])
    tc = np.array([float(r["t_cert"]) for r in rows])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.6))
    lo, hi = min(qs.min(), tt.min()), max(qs.max(), tt.max())
    m = 0.05 * (hi - lo + 1e-3)
    ax1.plot([lo - m, hi + m], [lo - m, hi + m], color=GRAY, lw=0.8)
    ax1.scatter(qs, tt, s=10, color=PAL[0], zorder=3, label="SDP $t$")
    ax1.scatter(qs, tc, s=10, color=PAL[1], marker="s", zorder=2,
                label="corrected $t_{cert}$")
    ax1.set_xlabel("trust-region reference $q^*$")
    ax1.set_ylabel("SDP value")
    ax1.set_title(f"(a) consistency (max |t-q*| = "
                  f"{np.max(np.abs(tt - qs)):.1e})")
    ax1.legend(loc="upper left", fontsize=8)
    ax2.hist(np.log10(np.maximum(np.abs(tt - qs), 1e-18)), bins=40,
             color=PAL[2])
    ax2.set_xlabel(r"$\log_{10}|t-q^*|$")
    ax2.set_ylabel("user problems")
    ax2.set_title("(b) error distribution (2160 user problems)")
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "fig_sdp_consistency.pdf"))
    fig.savefig(os.path.join(figdir, "fig_sdp_consistency.png"), dpi=150)
    plt.close(fig)

    tim = rd("pipeline_timings.csv")
    by = {}
    for r in tim:
        by.setdefault(r["pipeline"], []).append(float(r["seconds_med"]))
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    names = sorted(by)
    meds = [float(np.median(by[p])) for p in names]
    order = np.argsort(meds)
    names = [names[i] for i in order]
    meds = [meds[i] for i in order]
    ax.bar(range(len(names)), meds, 0.6, color=PAL[:len(names)])
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("median seconds per (seed,$\\gamma$)")
    ax.set_title("same final decision task (validation batch)")
    for i, v in enumerate(meds):
        ax.text(i, v, f"{v:.3f}", ha="center", va="bottom", fontsize=7.5)
    fig.tight_layout()
    fig.savefig(os.path.join(figdir, "fig_sdp_timing.pdf"))
    fig.savefig(os.path.join(figdir, "fig_sdp_timing.png"), dpi=150)
    plt.close(fig)
    print(f"[sdp:info] figures written to {figdir}")


def run_all(quick=False):
    utils.ensure_dirs()
    print("--- SDP smallcheck (pre-fixed selection rules) ---")
    sc_rows, med_wall = run_smallcheck()
    if quick:
        print("[sdp:info] quick mode stops after smallcheck")
        return not FAILURES
    if med_wall * 2160 > 1200:
        print(f"[sdp:info] projected batch cost {med_wall * 2160:.0f} s "
              "exceeds budget; run full batch anyway is NOT started -- "
              "adjust and rerun (no silent subsetting)")
        return False
    print("--- SDP full validation batch ---")
    rows, summary = run_batch()
    print("--- SDP pipeline comparison ---")
    tim = run_pipelines()
    print("--- SDP figures ---")
    make_figures()
    manifest = dict(
        finished_at=datetime.now().isoformat(timespec="seconds"),
        solver=SOLVER, failures=list(FAILURES),
        batch_summary=summary)
    utils.write_json(os.path.join(OUT_TAG, "run_manifest.json"), manifest)
    print(f"\n=== SDP round {'FAILED' if FAILURES else 'PASSED'} ===")
    for f in FAILURES:
        print("  FAIL:", f)
    return not FAILURES
