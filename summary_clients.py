# -*- coding: utf-8 -*-
"""Strict summary clients + compact TR cache for the V6.4 interface study.

Clients receive ONLY the permitted scalars (plus external r, gamma, sigma2,
J).  They never touch b, C, h_hat, A, g, eigenvectors or descriptor dicts.
The compact TR client receives only (lam, w, c) plus scale metadata.

tau_summary: conservative threshold computable from permitted scalars.
Derivation (six-scalar quantities only, gamma fixed):
  c_q = h^H A h - gamma sigma^2 = |b^H h|^2 - gamma ||C^H h||^2 - gamma sigma^2
      = p a^2 - gamma I0^2 - gamma sigma^2
  |c_q| <= p a^2 + gamma I0^2 + gamma sigma^2                       (T1)
  ||g|| = ||A h|| <= ||b b^H h|| + gamma ||C C^H h||
        <= p a + gamma beta0 I0                                       (T2)
  ||A||_2 <= ||b b^H|| + gamma ||C||^2 = p + gamma beta0^2             (T3)
  tau = 1e-8 * max(1, |c_q|_up, 2 r ||g||_up, r^2 ||A||_up)
Noise enters via the sigma^2 term in (T1); all scales derive from the
permitted summaries.  The seven-scalar client additionally knows kappa,
which only tightens (T2) to kappa itself; we keep ONE common tau from the
six-scalar bound so all branches share an identical decision rule.
"""
from __future__ import annotations

import numpy as np


# --------------------------------------------------------------------------
# summary-only threshold
# --------------------------------------------------------------------------
def tau_summary(p, a, I0, beta0, gamma, sigma2, r):
    c_up = p * a * a + gamma * I0 * I0 + gamma * sigma2
    g_up = p * a + gamma * beta0 * I0
    A_up = p + gamma * beta0 * beta0
    return 1e-8 * max(1.0, c_up, 2.0 * r * g_up, r * r * A_up)


# --------------------------------------------------------------------------
# A. six-scalar client
# --------------------------------------------------------------------------
def qg_lb_client(six, r, gamma, sigma2, J=256):
    """six: tuple (p, a, I0, alpha, beta, beta0) of floats ONLY."""
    p, a, I0, alpha, beta, beta0 = six
    th = np.linspace(0.0, np.pi / 2, J + 1)
    th_lo, th_hi = th[:-1], th[1:]
    th_m = 0.5 * (th_lo + th_hi)
    dth = th_hi - th_lo
    x_m, y_m = r * np.sin(th_m), r * np.cos(th_m)
    # R envelope (three-branch closed form, scalar inputs only)
    delta = np.sqrt(max(beta0 ** 2 - alpha ** 2, 0.0))
    nu = np.sqrt(max(beta0 ** 2 - beta ** 2, 0.0))
    s1 = np.sqrt(np.maximum(alpha ** 2 * x_m ** 2 + beta ** 2 * y_m ** 2
                            + 2.0 * delta * nu * x_m * y_m, 0.0))
    s2 = beta0 * np.sqrt(x_m ** 2 + y_m ** 2)
    s3 = alpha * x_m + delta * y_m
    b1 = delta * x_m <= nu * y_m
    b3 = delta * x_m >= alpha * y_m
    Rm = np.where(b1, s1, np.where(b3, s3, s2))
    R_U = np.minimum(beta0 * r, Rm + r * (alpha + beta) * dth / 2.0)
    x_lo, x_hi = r * np.sin(th_lo), r * np.sin(th_hi)
    dist = np.maximum(np.maximum(x_lo - a, a - x_hi), 0.0)
    cell = p * dist * dist - gamma * (I0 + R_U) ** 2 - gamma * sigma2
    return float(np.min(cell))


# --------------------------------------------------------------------------
# B. seven-scalar client (L_quad from scalars ONLY)
# --------------------------------------------------------------------------
def lquad_client(six, kappa, r, gamma, sigma2):
    p, a, I0, alpha, beta, beta0 = six
    c = p * a * a - gamma * I0 * I0 - gamma * sigma2
    return c - 2.0 * r * kappa - gamma * beta0 * beta0 * r * r


def joint_client(six, kappa, r, gamma, sigma2, J=256):
    return max(lquad_client(six, kappa, r, gamma, sigma2),
               qg_lb_client(six, r, gamma, sigma2, J))


# --------------------------------------------------------------------------
# D. compact TR cache client
# --------------------------------------------------------------------------
def build_compact_cache(b, h, C, gamma, sigma2):
    """Builder side: A = bb^H - gamma CC^H; keep ONLY lam, w=|U^H g|^2, c."""
    A = b @ b.conj().T - gamma * (C @ C.conj().T)
    g = A @ h
    c = float(np.real((h.conj().T @ A @ h)[0, 0])) - gamma * sigma2
    lam, U = np.linalg.eigh(A)
    w = np.abs(U.conj().T @ g).ravel() ** 2
    scale = max(1.0, float(np.abs(lam).max()), float(np.linalg.norm(g)))
    return dict(lam=lam, w=w, c=c, scale=scale)


def compact_cache_bytes(cache):
    return (int(cache["lam"].nbytes) + int(cache["w"].nbytes)
            + 8 + 8)  # c (float64) + scale metadata


def tr_qstar_from_cache(cache, r, bisect_iters=200):
    """Exact optimal value query from (lam, w, c) only.

    Handles: r=0; PSD interior (mu=0); boundary via bracketed bisection on
    phi(mu)=sum w_i/(lam_i+mu)^2 = r^2 (pseudo-inverse convention at
    singular denominators); hard case via norm completion in the min-eig
    subspace (value is direction-independent).  No eigengap perturbation.
    """
    lam = cache["lam"]
    w = cache["w"]
    c = cache["c"]
    n = lam.size
    if r <= 0.0:
        return c, "r_zero"
    tol = 1e-12 * cache["scale"]
    lam_min = float(lam[0])
    psd = lam_min >= -tol

    def phi_safe(mu):
        den = lam + mu
        good = den > tol
        bad = ~good
        if bad.any() and (w[bad] > 1e-10 * cache["scale"]).any():
            return np.inf
        return float(np.sum(w[good] / den[good] ** 2))

    mu_min = max(0.0, -lam_min)
    if psd:
        keep = lam > tol
        e2 = float(np.sum(w[keep] / lam[keep] ** 2)) if keep.any() else 0.0
        gnull = w[~keep].sum() if (~keep).any() else 0.0
        if gnull <= 1e-10 * cache["scale"] and e2 <= r * r:
            return c - phi_sum(w, lam, 0.0), "interior_psd"   # mu = 0
        mu = 0.0
    else:
        phi_pi = phi_safe(mu_min)
        gmin = w[lam <= lam_min + tol].sum()
        if gmin <= 1e-10 * cache["scale"] and phi_pi < r * r:
            # hard case: the dual value at mu = mu_min already equals the
            # optimum; norm completion is only needed to recover e itself
            return (c - mu_min * r ** 2 - phi_sum(w, lam, mu_min),
                    "hard_case")
        mu = mu_min
    mu_hi = max(1.0, cache["scale"]) + mu_min
    for _ in range(1000):
        if phi_safe(mu_hi) <= r * r:
            break
        mu_hi *= 2.0
    for _ in range(bisect_iters):
        mu_mid = 0.5 * (mu + mu_hi)
        if phi_safe(mu_mid) > r * r:
            mu = mu_mid
        else:
            mu_hi = mu_mid
    return c - mu_hi * r * r - phi_sum(w, lam, mu_hi), "boundary"


def phi_sum(w, lam, mu):
    """sum w_i/(lam_i+mu) with pseudo-inverse convention (g_i=0 entries
    dropped); called only when the convention is well-defined."""
    den = lam + mu
    good = den > 1e-12 * max(1.0, float(np.abs(lam).max()))
    return float(np.sum(w[good] / den[good]))
