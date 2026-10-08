# -*- coding: utf-8 -*-
"""V6.3 targeted checks.

1. Norm-bounded construction at the upper spectral endpoint (counter-
   example to the blanket divergence claim): p=a=I0=alpha=beta=gamma=1,
   beta0=sqrt(2), r=0.2, sigma^2=0.05.
2. Candidate finite-attainment instance at that endpoint.
3. Rectangular dimensions n=5/m=3, n=7/m=4, m>n zero-padding.
4. Degenerate spectra, branch switches, Lemma 1/2 recovery.
Writes results_v63/.  Checks verify implementations, not general proofs.
"""
import csv
import json
import os

import numpy as np

import core
import experiments_c as ec

OUT = "results_v63"
FAILS = []


def rec(rows, name, ok, detail):
    rows.append(dict(case=name, ok=bool(ok), detail=detail))
    print(f"[v63:{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        FAILS.append(name)


def eps_at_endpoint(rows):
    """Counter-example: upper spectral endpoint keeps ||h|| bounded."""
    alpha = beta = gamma = 1.0
    beta0 = np.sqrt(2.0)
    a = I0 = 1.0
    p = 2.0
    r, noise = 0.2, 0.05
    scal = dict(p=p, a=a, I0=I0, alpha=alpha, beta=beta, beta0=beta0)
    qg, x, y = core.Q_approx(scal, r, gamma, noise, "G")
    rec(rows, "endpoint/QG", abs(qg - (-0.3119466960)) < 1e-6,
        f"Q_G={qg:.10f} (expected ~ -0.3119466960)")
    rec(rows, "endpoint/argmin", abs(x - 0.1833242752) < 1e-5
        and abs(y - 0.0799512984) < 1e-5, f"x={x:.10f} y={y:.10f}")
    norms = []
    gaps = []
    for eps in (1e-1, 1e-2, 1e-3, 1e-5):
        con = core.construct_v41(alpha, beta, beta0, x, y, eps, p, a, I0)
        s = core.spectral_summary(con["b"], con["h"], con["C"], gamma, noise)
        sixerr = max(abs(s[k] - scal[k]) for k in scal)
        obj = core.objective_at(con["b"], con["h"], con["C"], gamma, noise,
                                con["e"])
        hn = float(np.linalg.norm(con["h"]))
        norms.append(hn)
        gaps.append(obj - qg)
        rec(rows, f"endpoint/eps{eps:.0e}", sixerr < 1e-12,
            f"six-err {sixerr:.1e}, ||h||={hn:.10f}, obj-Q_G={obj-qg:.2e}, "
            f"t_eps={con['t_eps']:.6f}")
    rec(rows, "endpoint/hnorm_bounded",
        all(abs(hn - np.sqrt(5.0)) < 1e-9 for hn in norms),
        f"all ||h|| ~ sqrt(5)={np.sqrt(5):.10f}")
    rec(rows, "endpoint/obj_converges",
        all(g >= -1e-12 for g in gaps) and gaps[-1] < 1e-5,
        f"final gap {gaps[-1]:.2e} = O(eps) as proven (bound 22)")


def finite_attainment(rows):
    """Candidate finite-norm attainment at the upper endpoint (n=5, m=3)."""
    alpha = beta = gamma = 1.0
    beta0 = np.sqrt(2.0)
    a = I0 = 1.0
    p = 2.0
    r, noise = 0.2, 0.05
    CH = np.array([[1, 1, 0, 0, 0],
                   [0, 0, 0, 1, 0],
                   [0, 0, 0, 0, 1]], dtype=float)
    h = np.array([1.0, -2.0, 0, 0, 0]).reshape(-1, 1)
    b = np.array([np.sqrt(2.0), 0, 0, 0, 0]).reshape(-1, 1)
    s = core.spectral_summary(b, h, CH.conj().T, gamma, noise)
    want = dict(p=p, a=a, I0=I0, alpha=alpha, beta=beta, beta0=beta0)
    sixerr = max(abs(s[k] - want[k]) for k in want)
    rec(rows, "attain/six", sixerr < 1e-12, f"six-err {sixerr:.2e}")
    qg, x, y = core.Q_approx(want, r, gamma, noise, "G")
    e = np.array([-x, -y, 0, 0, 0]).reshape(-1, 1)
    rec(rows, "attain/ball", np.linalg.norm(e) <= r + 1e-12,
        f"||e||={np.linalg.norm(e):.10f} vs r={r}")
    obj = core.objective_at(b, h, CH.conj().T, gamma, noise, e)
    ref = core.q_star_reference(s, r, gamma, noise)
    rec(rows, "attain/equation", abs(obj - qg) < 1e-10,
        f"obj={obj:.10f} vs Q_G={qg:.10f} (diff {obj-qg:.2e})")
    rec(rows, "attain/reference", abs(ref["q"] - qg) < 1e-8,
        f"q*={ref['q']:.10f}")


def rectangular(rows):
    """Lemma 2 / Lemma 1 constructions on rectangular dims."""
    rng = np.random.default_rng(63)

    def lemma2_C(alpha, beta, beta0, x, y, n, m):
        # outputs live in C^m; choose d, l, l2 in C^m
        if alpha == 0:
            d = np.zeros(m, complex)
            ell = np.zeros(m, complex)
            ell[0] = beta
            ell2 = np.zeros(m, complex)
            ell2[1] = beta
        else:
            delta = np.sqrt(beta0 ** 2 - alpha ** 2)
            nu = np.sqrt(beta0 ** 2 - beta ** 2)
            k = alpha if y == 0 else float(np.clip(delta * x / y, nu,
                                                   alpha))
            d = np.zeros(m, complex)
            d[0] = alpha
            ell = np.zeros(m, complex)
            ell[0] = delta * k / alpha
            ell[1] = beta0 * np.sqrt(max(1 - (k / alpha) ** 2, 0.0))
            ell2 = np.zeros(m, complex)
            ell2[2] = beta
        VU = np.zeros((n, 3), complex)
        VU[0, 0], VU[1, 1], VU[2, 2] = 1, 1, 1     # [v, u, u']
        return np.column_stack((d, ell, ell2)) @ VU.conj().T  # m x n

    for (n, m, tag) in [(5, 3, "5x3"), (7, 4, "7x4"), (3, 5, "m>n")]:
        for trial, (al, be) in enumerate([(1.0, 0.7), (0.8, 1.2)]):
            b0 = np.hypot(al, be)          # upper endpoint exercises k=alpha
            x, y = 0.3, 0.2
            CH = lemma2_C(al, be, b0, x, y, n, m)
            C = CH.conj().T
            bvec = np.zeros((n, 1), complex)
            bvec[0, 0] = 1.0
            hv = np.zeros((n, 1), complex)
            s = core.spectral_summary(bvec, hv, C, 1.0, 0.05)
            R = core.R_envelope(x, y, al, be, b0)
            e = np.zeros((n, 1), complex)
            e[0, 0], e[1, 0] = x, y
            out = float(np.linalg.norm(C.conj().T @ e))
            errs = max(abs(s["alpha"] - al), abs(s["beta"] - be),
                       abs(s["beta0"] - b0))
            rec(rows, f"rect/{tag}/t{trial}",
                errs < 1e-12 and abs(out - R) < 1e-12,
                f"six {errs:.1e}, |out-R| {abs(out-R):.1e}")

    # zero-padding instance with m>n handled above; also beta0=alpha boundary
    for (n, m) in [(5, 3), (7, 4)]:
        al, be = 1.0, 0.6
        CH = lemma2_C(al, be, al, 0.25, 0.15, n, m)   # boundary branch: k=alpha
        C = CH.conj().T
        bvec = np.zeros((n, 1), complex)
        bvec[0, 0] = 1.0
        hv = np.zeros((n, 1), complex)
        s = core.spectral_summary(bvec, hv, C, 1.0, 0.05)
        errs = max(abs(s["alpha"] - al), abs(s["beta"] - be),
                   abs(s["beta0"] - al))
        rec(rows, f"rect/boundary_{n}x{m}", errs < 1e-11,
            f"six {errs:.1e}")


def degenerate_and_switches(rows):
    """Zero spectra, beta=0, alpha=0, x/y=0, branch switch points."""
    # Lemma 2 direction attainment at switch points of Theorem 1
    for (x, y, tag) in [(0.0, 0.4, "x=0"), (0.4, 0.0, "y=0")]:
        al, be, b0 = 1.0, 0.8, 1.1
        # reuse 3x3 construction from v62 checker logic
        delta = np.sqrt(b0 ** 2 - al ** 2)
        nu = np.sqrt(b0 ** 2 - be ** 2)
        k = al if y == 0 else float(np.clip(delta * x / y, nu, al))
        d = np.array([al, 0, 0], complex)
        ell = np.array([delta * k / al,
                        b0 * np.sqrt(max(1 - (k / al) ** 2, 0)), 0],
                       complex)
        ell2 = np.array([0, 0, be], complex)
        VU = np.eye(3, dtype=complex)
        CH = np.column_stack((d, ell, ell2)) @ VU.conj().T
        C = CH.conj().T
        bvec = np.zeros((3, 1), complex)
        bvec[0, 0] = 1.0
        hv = np.zeros((3, 1), complex)
        s = core.spectral_summary(bvec, hv, C, 1.0, 0.05)
        R = core.R_envelope(x, y, al, be, b0)
        e = np.array([x, y, 0]).reshape(-1, 1)
        out = float(np.linalg.norm(C.conj().T @ e))
        errs = max(abs(s["alpha"] - al), abs(s["beta"] - be),
                   abs(s["beta0"] - b0))
        rec(rows, f"switch/{tag}", errs < 1e-12 and abs(out - R) < 1e-12,
            f"six {errs:.1e}, |out-R| {abs(out-R):.1e}")
    # exact branch-switch curve: delta*x = nu*y and delta*x = alpha*y
    al, be, b0 = 1.0, 0.8, 1.1
    delta = np.sqrt(b0 ** 2 - al ** 2)
    nu = np.sqrt(b0 ** 2 - be ** 2)
    for tag, (x, y) in [("b12", (nu * 0.3, delta * 0.3)),
                        ("b23", (al * 0.3, delta * 0.3))]:
        k = float(np.clip(delta * x / y, nu, al))
        d = np.array([al, 0, 0], complex)
        ell = np.array([delta * k / al,
                        b0 * np.sqrt(max(1 - (k / al) ** 2, 0)), 0],
                       complex)
        ell2 = np.array([0, 0, be], complex)
        CH = np.column_stack((d, ell, ell2))
        C = CH.conj().T
        bvec = np.zeros((3, 1), complex)
        bvec[0, 0] = 1.0
        hv = np.zeros((3, 1), complex)
        s = core.spectral_summary(bvec, hv, C, 1.0, 0.05)
        R = core.R_envelope(x, y, al, be, b0)
        e = np.array([x, y, 0]).reshape(-1, 1)
        out = float(np.linalg.norm(C.conj().T @ e))
        rec(rows, f"switch/{tag}", abs(out - R) < 1e-12,
            f"|out-R| {abs(out-R):.1e}")
    # Lemma 1 recovery on degenerate families
    def l1(p, a, I0, al, be, b0, tag):
        # build minimal instance per branch, n=5, m=3 via explicit rows
        CH = np.zeros((3, 5), complex)
        hv = np.zeros((5, 1), complex)
        hv[0, 0] = a
        if al == 0 and be == 0:
            pass                                     # C = 0
        elif be == 0:
            CH[0, 0] = al
            hv[0, 0] = a                              # I0 = alpha a
        elif al == 0:
            CH[0, 0], CH[0, 1] = 0.0, be
            hv[1, 0] = I0 / be
        elif b0 > al:
            dl = np.sqrt(b0 ** 2 - al ** 2)
            CH[0, 0], CH[0, 1] = al, dl
            CH[1, 2] = be
            hv[1, 0] = (I0 - al * a) / dl
        else:
            Jp = np.sqrt(max(I0 ** 2 - al ** 2 * a ** 2, 0.0))
            CH[0, 0], CH[1, 1] = al, be
            hv[1, 0] = Jp / be
        bv = np.zeros((5, 1), complex)
        bv[0, 0] = np.sqrt(p)
        s = core.spectral_summary(bv, hv, CH.conj().T, 1.0, 0.05)
        want = dict(p=p, a=a, I0=I0, alpha=al, beta=be, beta0=b0)
        err = max(abs(s[k] - want[k]) for k in want)
        rec(rows, f"lemma1/{tag}", err < 1e-12, f"six {err:.1e}")

    l1(2.0, 1.3, 0.0, 0, 0, 0, "zero_op")
    l1(2.0, 1.3, 1.95, 1.5, 0, 1.5, "beta0_eq_alpha_beta0")
    l1(2.0, 0.5, 0.9, 0.0, 0.9, 0.9, "alpha0")
    l1(2.0, 1.0, 0.4, 1.0, 0.7, np.hypot(1, 0.7), "upper_endpt")
    l1(2.0, 1.0, 1.3, 1.0, 0.7, 1.0, "boundary_Jpos")
    l1(2.0, 1.0, 1.0, 1.0, 0.7, 1.0, "boundary_J0")


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = []
    eps_at_endpoint(rows)
    finite_attainment(rows)
    rectangular(rows)
    degenerate_and_switches(rows)
    with open(os.path.join(OUT, "construction_checks_v63.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["case", "ok", "detail"])
        w.writeheader()
        w.writerows(rows)
    print(f"\n{len(rows)} checks; failures: {FAILS if FAILS else 'none'}")
    return not FAILS


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
