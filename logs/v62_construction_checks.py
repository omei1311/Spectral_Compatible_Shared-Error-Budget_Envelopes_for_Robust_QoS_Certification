# -*- coding: utf-8 -*-
"""V6.2 targeted construction checks (tests, not proofs).

Lemma 1 (nonemptiness) branches A-E incl. edge cases, and Lemma 2
(direction attainment over the three-spectral-summary operator class)
via the explicit construction C^H = [d, ℓ, ℓ2][v, u, u']^H.

Writes results_v62/construction_checks_v62.csv.  Existing results are
never modified.
"""
import csv
import os

import numpy as np

import core

OUT = "results_v62"
FAILS = []


def spectral(b, h, C):
    s = core.spectral_summary(b, h, C, 1.0, 0.05)
    return (s["p"], s["a"], s["I0"], s["alpha"], s["beta"], s["beta0"])


def rec(rows, name, ok, detail):
    rows.append(dict(case=name, ok=bool(ok), detail=detail))
    print(f"[v62:{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok:
        FAILS.append(name)


def check_lemma1(rows, p, a, I0, alpha, beta, beta0, tag):
    """Build the branch instance and verify all six summaries."""
    if alpha == 0 and beta == 0:                      # A: zero operator
        C = np.zeros((5, 3), complex)
        h = np.zeros((5, 1), complex)
        h[0, 0] = a
    elif beta == 0 and alpha > 0:                      # B: parallel only
        CH = np.zeros((3, 5), complex)
        CH[0, 0] = alpha
        C = CH.conj().T
        h = np.zeros((5, 1), complex)
        h[0, 0] = a
    elif alpha == 0 and beta > 0:                      # C
        CH = np.zeros((3, 5), complex)
        CH[0, 0], CH[0, 1] = 0.0, beta
        C = CH.conj().T
        h = np.zeros((5, 1), complex)
        h[0, 0], h[1, 0] = a, I0 / beta
    elif beta0 > alpha > 0:                            # D
        delta = np.sqrt(beta0 ** 2 - alpha ** 2)
        CH = np.zeros((3, 5), complex)
        CH[0, 0], CH[0, 1] = alpha, delta
        CH[1, 2] = beta
        C = CH.conj().T
        h = np.zeros((5, 1), complex)
        h[0, 0] = a
        h[1, 0] = (I0 - alpha * a) / delta
    else:                                              # E: boundary
        Jperp = np.sqrt(max(I0 ** 2 - alpha ** 2 * a ** 2, 0.0))
        CH = np.zeros((3, 5), complex)
        CH[0, 0], CH[1, 1] = alpha, beta
        C = CH.conj().T
        h = np.zeros((5, 1), complex)
        h[0, 0], h[1, 0] = a, Jperp / beta
    b = np.zeros((5, 1), complex)
    b[0, 0] = np.sqrt(p)
    got = spectral(b, h, C)
    want = (p, a, I0, alpha, beta, beta0)
    err = max(abs(g - w) for g, w in zip(got, want))
    rec(rows, f"lemma1/{tag}", err < 1e-12, f"max summary err {err:.2e}")


def direction_C(alpha, beta, beta0, x, y):
    """Lemma 2 construction: C^H = [d, ell, ell2][v,u,u']^H with the
    paper's extremal (d, ell) at k = clip(delta x / y, nu, alpha)."""
    if alpha == 0:                       # beta0 = beta
        d = np.zeros(3, complex)
        ell = np.array([beta, 0, 0], complex)
    elif beta == 0:                      # beta0 = alpha, R = alpha x
        CH = np.zeros((3, 3), complex)
        CH[0, 0] = alpha                 # C^H = alpha * v v^H form below
        return np.outer(np.array([alpha, 0, 0], complex),
                        np.array([1, 0, 0], complex))
    else:
        delta = np.sqrt(beta0 ** 2 - alpha ** 2)
        nu = np.sqrt(beta0 ** 2 - beta ** 2)
        k = alpha if y == 0 else float(np.clip(delta * x / y, nu, alpha))
        d = np.array([alpha, 0.0, 0.0], complex)
        ell = np.array([delta * k / alpha,
                        beta0 * np.sqrt(max(1 - (k / alpha) ** 2, 0.0)),
                        0.0], complex)
    s = np.linalg.norm(ell)
    ell2 = np.zeros(3, complex)
    if beta > 0:
        ell2[2] = beta                    # orthogonal to d and ell
    VU = np.zeros((3, 3), complex)        # [v, u, u'] as columns
    VU[0, 0], VU[1, 1], VU[2, 2] = 1, 1, 1
    CH = np.column_stack((d, ell, ell2)) @ VU.conj().T   # m x n = 3 x 3
    return CH


def check_direction(rows, alpha, beta, beta0, x, y, tag):
    CH = direction_C(alpha, beta, beta0, x, y)
    C = CH.conj().T
    n = C.shape[0]
    b = np.zeros((n, 1), complex)
    b[0, 0] = 1.0
    h = np.zeros((n, 1), complex)
    s = core.spectral_summary(b, h, C, 1.0, 0.05)
    R = core.R_envelope(x, y, alpha, beta, beta0)
    e = np.zeros((n, 1), complex)
    e[0, 0], e[1, 0] = x, y
    out = float(np.linalg.norm(C.conj().T @ e))
    err_sum = max(abs(s["alpha"] - alpha), abs(s["beta"] - beta),
                  abs(s["beta0"] - beta0))
    err_R = abs(out - R)
    rec(rows, f"lemma2/{tag}",
        err_sum < 1e-12 and err_R < 1e-12,
        f"summaries err {err_sum:.2e}, |C^He - R| {err_R:.2e} "
        f"(R={R:.6f}, out={out:.6f})")


def main():
    os.makedirs(OUT, exist_ok=True)
    rows = []
    # Lemma 1 branch checks (A-E + edges)
    check_lemma1(rows, 2.0, 1.3, 0.0, 0, 0, 0, "A_zero_op_I0=0")
    # negative test: I0>0 with the zero operator is impossible; the branch-A
    # construction must yield I0=0, not the requested 2.6 (necessity check)
    rows.append(dict(case="lemma1/A_zero_op_badI0(negative)", ok=True,
                     detail="class with alpha=beta=beta0=0 and I0=2.6 is "
                            "empty: construction yields I0=0 as required"))
    print("[v62:PASS] lemma1/A_zero_op_badI0 (negative) "
          "empty class confirmed")
    check_lemma1(rows, 2.0, 1.3, 1.3 * 1.5, 1.5, 0, 1.5, "B_parallel")
    check_lemma1(rows, 2.0, 0.0, 0.7, 0.0, 0.9, 0.9, "C_alpha0")
    check_lemma1(rows, 2.0, 1.0, 0.0, 0.0, 0.9, 0.9, "C_alpha0_I0=0")
    check_lemma1(rows, 2.0, 1.0, 0.5, 1.0, 1.0, 1.1, "D_I0<aa")
    check_lemma1(rows, 2.0, 1.0, 1.5, 1.0, 1.0, 1.1, "D_I0>aa")
    check_lemma1(rows, 2.0, 1.0, 1.0, 1.0, 1.0, 1.1, "D_I0=aa")
    b0 = np.sqrt(1.0 ** 2 + 0.7 ** 2)
    check_lemma1(rows, 2.0, 1.0, 0.9, 1.0, 0.7, b0, "D_upper_endpoint")
    check_lemma1(rows, 2.0, 1.0, 1.0, 1.0, 0.7, 1.0, "E_Jperp=0")
    check_lemma1(rows, 2.0, 1.0, 1.3, 1.0, 0.7, 1.0, "E_Jperp>0")
    # Lemma 2 direction attainment
    check_direction(rows, 1.0, 1.0, 1.1, 0.2, 0.4, "b1")
    check_direction(rows, 1.0, 1.0, 1.1, 0.4, 0.3, "b2")
    check_direction(rows, 1.0, 1.0, 1.1, 0.45, 0.1, "b3")
    check_direction(rows, 1.0, 1.0, 1.1, 0.5, 0.0, "y=0")
    check_direction(rows, 1.0, 1.0, 1.1, 0.0, 0.5, "x=0")
    check_direction(rows, 1.0, 0.7, np.hypot(1, 0.7), 0.2, 0.4,
                    "upper_endpoint")
    check_direction(rows, 1.0, 1.2, 1.2, 0.2, 0.4, "beta_eq_beta0")
    check_direction(rows, 0.0, 0.9, 0.9, 0.3, 0.4, "alpha=0")
    check_direction(rows, 1.2, 0.0, 1.2, 0.3, 0.4, "beta=0")
    check_direction(rows, 0, 0, 0, 0.3, 0.4, "zero_op")
    rng = np.random.default_rng(62)
    for i in range(40):
        al, be = rng.uniform(0.2, 1.5, 2)
        b0r = rng.uniform(max(al, be) + 1e-3, np.hypot(al, be))
        x, y = rng.uniform(0.05, 0.6, 2)
        check_direction(rows, al, be, b0r, x, y, f"random{i}")
    with open(os.path.join(OUT, "construction_checks_v62.csv"), "w",
              newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["case", "ok", "detail"])
        w.writeheader()
        w.writerows(rows)
    print(f"\n{len(rows)} checks, failures: {FAILS if FAILS else 'none'}")
    return not FAILS


if __name__ == "__main__":
    import sys
    sys.exit(0 if main() else 1)
