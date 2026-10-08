# -*- coding: utf-8 -*-
"""Core numerics for: Spectral-Compatible Shared-Error-Budget Envelopes for
Robust QoS Certification.

Conventions (matching the experiment specification):
  * columns vectors, H = conjugate transpose, CN(0,s^2) has independent
    real/imag parts each of variance s^2/2
  * noise parameter is a variance; variable name ``noise_power`` (never
    squared again)
  * signal term of Q and Q_G is the UNtruncated p*(a-x)^2; only D uses
    max(a-r,0)
"""
from __future__ import annotations

import numpy as np

EIG_TOL_REL = 1e-10  # relative eigenvalue tolerance for compatibility checks


# ----------------------------------------------------------------------------
# 1. spectral summary
# ----------------------------------------------------------------------------
def spectral_summary(b, h, C, gamma, noise_power):
    """Six scalars plus quadratic-form data for one (b, h, C) instance.

    Returns dict with keys
      p, a, I0, alpha, beta, beta0  (six scalars, beta via Hermitian eig)
      v, P, A, g, c                 (A = b b^H - gamma C C^H, g = A h,
                                     c = h^H A h - gamma*noise_power)
      compat_max                    (max(alpha,beta) - beta0, must be <= ~0)
    """
    b = np.asarray(b, dtype=complex).reshape(-1, 1)
    h = np.asarray(h, dtype=complex).reshape(-1, 1)
    C = np.asarray(C, dtype=complex)
    n = b.shape[0]

    nb = float(np.linalg.norm(b))
    if nb == 0:
        raise ValueError("b must be nonzero")
    v = b / nb
    P = np.eye(n, dtype=complex) - v @ v.conj().T

    p = nb * nb
    a = float(np.abs(v.conj().T @ h)[0, 0])
    I0 = float(np.linalg.norm(C.conj().T @ h))
    alpha = float(np.linalg.norm(C.conj().T @ v))

    CHC = C.conj().T @ C                      # m x m Hermitian
    CHPC = C.conj().T @ P @ C                 # m x m Hermitian
    w0 = np.linalg.eigvalsh(CHC)
    wP = np.linalg.eigvalsh(CHPC)
    lam0 = float(w0[-1])
    lamP = float(wP[-1])
    scale0 = max(1.0, float(np.abs(w0).max()))
    if lam0 < -EIG_TOL_REL * scale0:
        raise ValueError(f"C^H C not PSD: lambda_max={lam0:.3e}")
    if lamP < -EIG_TOL_REL * scale0:
        raise ValueError(f"C^H P C not PSD: lambda_max={lamP:.3e}")
    beta0 = float(np.sqrt(max(lam0, 0.0)))
    beta = float(np.sqrt(max(lamP, 0.0)))

    viol_up = max(alpha, beta) - beta0
    viol_dn = beta0 - np.sqrt(alpha * alpha + beta * beta)
    tol_c = EIG_TOL_REL * max(1.0, beta0, alpha, beta)
    if viol_up > 1e-8 * max(1.0, beta0):
        raise ValueError(
            f"incompatible scalars: max(alpha,beta)={max(alpha,beta):.6f} "
            f"> beta0={beta0:.6f}")
    if viol_dn > 1e-8 * max(1.0, beta0):
        raise ValueError(
            f"incompatible scalars: beta0={beta0:.6f} > "
            f"sqrt(alpha^2+beta^2)={np.sqrt(alpha**2+beta**2):.6f}")

    A = b @ b.conj().T - gamma * (C @ C.conj().T)
    g = A @ h
    c = float(np.real((h.conj().T @ A @ h)[0, 0])) - gamma * noise_power

    return dict(p=p, a=a, I0=I0, alpha=alpha, beta=beta, beta0=beta0,
                v=v, P=P, A=A, g=g, c=c,
                compat_max=float(viol_up), compat_min=float(viol_dn),
                tol_clip=float(tol_c))


# ----------------------------------------------------------------------------
# 2. directional envelope R and the old R_old
# ----------------------------------------------------------------------------
def R_envelope(x, y, alpha, beta, beta0):
    """Closed-form R(x,y) (vectorized).  Branch order: 1 first, then 3."""
    alpha, beta, beta0 = float(alpha), float(beta), float(beta0)
    delta = np.sqrt(max(beta0 * beta0 - alpha * alpha, 0.0))
    nu = np.sqrt(max(beta0 * beta0 - beta * beta, 0.0))
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    s1 = np.sqrt(np.maximum(alpha**2 * x**2 + beta**2 * y**2
                            + 2.0 * delta * nu * x * y, 0.0))
    s2 = beta0 * np.sqrt(x**2 + y**2)
    s3 = alpha * x + delta * y
    b1 = delta * x <= nu * y
    b3 = delta * x >= alpha * y
    out = np.where(b1, s1, np.where(b3, s3, s2))
    return out if out.ndim else float(out)


def R_old(x, y, alpha, beta, beta0, r=None):
    """Old directional envelope min{beta0*r, alpha*x + beta*y}."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    t = alpha * x + beta * y
    cap = beta0 * (np.sqrt(x**2 + y**2) if r is None else r)
    out = np.minimum(cap, t)
    return out if out.ndim else float(out)


# ----------------------------------------------------------------------------
# 3. comparison bounds D, L_quad
# ----------------------------------------------------------------------------
def D_bound(scal, r, gamma, noise_power):
    a_r = max(scal['a'] - r, 0.0)
    return (scal['p'] * a_r * a_r
            - gamma * (scal['I0'] + scal['beta0'] * r) ** 2
            - gamma * noise_power)


def L_quad(scal, r, gamma):
    """c - 2 r ||g|| - gamma beta0^2 r^2  (uses data beyond the six scalars)."""
    gnorm = float(np.linalg.norm(scal['g']))
    return scal['c'] - 2.0 * r * gnorm - gamma * scal['beta0']**2 * r * r


# ----------------------------------------------------------------------------
# 4. certified interval lower bounds for Q and Q_G
# ----------------------------------------------------------------------------
def interval_lower_bound(scal, r, gamma, noise_power, J=256, kind='G'):
    """Certified lower bound on min_{||e||<=r} of the QoS margin.

    kind='G' : new envelope R
    kind='old' : R_old = min{beta0 r, alpha x + beta y}
    Partition theta in [0, pi/2] into J cells; per cell use
      R_U = min{beta0 r, f(theta_m) + r (alpha+beta) dtheta/2}
      signal LB = p * dist(a, [r sin th_lo, r sin th_hi])^2   (untruncated)
    """
    p, a, I0 = scal['p'], scal['a'], scal['I0']
    alpha, beta, beta0 = scal['alpha'], scal['beta'], scal['beta0']
    th = np.linspace(0.0, np.pi / 2, J + 1)
    th_lo, th_hi = th[:-1], th[1:]
    th_m = 0.5 * (th_lo + th_hi)
    dth = th_hi - th_lo
    x_m, y_m = r * np.sin(th_m), r * np.cos(th_m)
    if kind == 'G':
        f_m = R_envelope(x_m, y_m, alpha, beta, beta0)
    elif kind == 'old':
        f_m = R_old(x_m, y_m, alpha, beta, beta0)
    else:
        raise ValueError(kind)
    R_U = np.minimum(beta0 * r, f_m + r * (alpha + beta) * dth / 2.0)
    x_lo, x_hi = r * np.sin(th_lo), r * np.sin(th_hi)
    dist = np.maximum(np.maximum(x_lo - a, a - x_hi), 0.0)
    sig_lb = p * dist * dist
    cell = sig_lb - gamma * (I0 + R_U) ** 2 - gamma * noise_power
    return float(np.min(cell))


def lb_with_escalation(scal, r, gamma, noise_power, kind='G',
                       J_list=(256, 1024, 4096), ub_gap_tol=None):
    """Refine until acceptance or requested numerical upper/lower gap.

    An approximate upper value is used only as a stopping diagnostic.
    Acceptance always uses the interval lower bound.
    """
    if not J_list or any(J <= 0 for J in J_list):
        raise ValueError('J_list must contain positive cell counts')
    upper = Q_approx(scal, r, gamma, noise_power, kind)[0] if ub_gap_tol is not None else None
    tau = engineering_tau(scal, r)[0]
    LB = -np.inf
    for J in J_list:
        LB = max(LB, interval_lower_bound(scal, r, gamma, noise_power, J, kind))
        if (upper is None and LB > tau) or (upper is not None and upper - LB <= ub_gap_tol):
            break
    return float(LB), J


# ----------------------------------------------------------------------------
# 5. approximate 1-D minimizations (upper estimates -- NOT certificates)
# ----------------------------------------------------------------------------
def _min_1d(fn, r, n_grid=8193, n_refine=60):
    """Dense grid + golden-section refinement of min_{0<=x<=r} fn(x)."""
    xs = np.linspace(0.0, r, n_grid)
    vals = fn(xs)
    i = int(np.argmin(vals))
    lo = xs[max(i - 1, 0)]
    hi = xs[min(i + 1, n_grid - 1)]
    gr = (np.sqrt(5.0) - 1.0) / 2.0
    x1, x2 = hi - gr * (hi - lo), lo + gr * (hi - lo)
    f1, f2 = fn(x1), fn(x2)
    for _ in range(n_refine):
        if f1 <= f2:
            hi, x2, f2 = x2, x1, f1
            x1 = hi - gr * (hi - lo)
            f1 = fn(x1)
        else:
            lo, x1, f1 = x1, x2, f2
            x2 = lo + gr * (hi - lo)
            f2 = fn(x2)
    xm = 0.5 * (lo + hi)
    candidates = [(float(fn(xm)), float(xm)),
                  (float(vals[i]), float(xs[i])),
                  (float(vals[0]), 0.0), (float(vals[-1]), float(r))]
    return min(candidates)


def Q_boundary_approx(scal, r, gamma, noise_power):
    """Paper Theorem 3: separate phase-compatible boundary Q_partial.

    This does not replace the general safe Q_G; it is applicable only
    to beta0=alpha>0 and boundary-compatible nominal summaries.
    """
    alpha, beta = scal['alpha'], scal['beta']
    if alpha <= 0 or not np.isclose(scal['beta0'], alpha, rtol=1e-12, atol=0):
        raise ValueError('Q_partial requires beta0=alpha>0')
    J2 = scal['I0']**2 - alpha**2 * scal['a']**2
    if J2 < -1e-12 * max(1.0, scal['I0']**2) or (beta == 0 and abs(J2) > 1e-12):
        raise ValueError('incompatible boundary nominal summary')
    J = np.sqrt(max(J2, 0.0))
    c = scal['p'] - gamma * alpha**2
    def obj(x):
        y = np.sqrt(np.maximum(r*r - x*x, 0.0))
        signal = (scal['a'] - x)**2 if c >= 0 else (scal['a'] + x)**2
        return c * signal - gamma * (J + beta*y)**2 - gamma * noise_power
    return _min_1d(obj, r)


def Q_approx(scal, r, gamma, noise_power, kind='G'):
    """Numerical min of the 1-D objective (upper estimate of the true min).
    Signal term UNtruncated: p (a-x)^2."""
    p, a, I0 = scal['p'], scal['a'], scal['I0']
    alpha, beta, beta0 = scal['alpha'], scal['beta'], scal['beta0']
    delta = np.sqrt(max(beta0**2 - alpha**2, 0.0))
    nu = np.sqrt(max(beta0**2 - beta**2, 0.0))

    def obj(x):
        x = np.asarray(x, dtype=float)
        y = np.sqrt(np.maximum(r * r - x * x, 0.0))
        if kind == 'G':
            rho = R_envelope(x, y, alpha, beta, beta0)
        else:
            rho = np.minimum(beta0 * r, alpha * x + beta * y)
        return p * (a - x) ** 2 - gamma * (I0 + rho) ** 2 - gamma * noise_power

    val, xm = _min_1d(obj, r)
    ym = float(np.sqrt(max(r * r - xm * xm, 0.0)))
    return val, xm, ym


def Q_grid_value(scal, r, gamma, noise_power, kind='G', n=4097):
    """Plain grid min + endpoint/midpoint evaluation (upper estimate)."""
    p, a, I0 = scal['p'], scal['a'], scal['I0']
    alpha, beta, beta0 = scal['alpha'], scal['beta'], scal['beta0']
    xs = np.linspace(0.0, r, n)
    y = np.sqrt(np.maximum(r * r - xs * xs, 0.0))
    if kind == 'G':
        rho = R_envelope(xs, y, alpha, beta, beta0)
    else:
        rho = np.minimum(beta0 * r, alpha * xs + beta * y)
    vals = p * (a - xs) ** 2 - gamma * (I0 + rho) ** 2 - gamma * noise_power
    return float(np.min(vals))


# ----------------------------------------------------------------------------
# 6. engineering tolerance
# ----------------------------------------------------------------------------
def engineering_tau(scal, r):
    s = max(1.0, abs(scal['c']), 2.0 * r * float(np.linalg.norm(scal['g'])),
            r * r * float(np.linalg.norm(scal['A'], 2)))
    return 1e-8 * s, s


# ----------------------------------------------------------------------------
# 7. exact reference: sphere-constrained trust-region quadratic
# ----------------------------------------------------------------------------
def trust_region_min(A, g, r, bisect_iters=200, eig=None, mu_rel_tol=None):
    """Global min of  c0 + e^H A e + 2 Re(e^H g)  over ||e|| <= r
    (c0 is added by the caller).

    ``eig`` = optional precomputed (lam, U) from np.linalg.eigh(A), so that
    repeated solves on the same instance with different radii reuse one
    factorization (fair caching for the reference method too).

    Returns dict(e, mu, mode, quad_val, checks{...}).
    Hermitian eigendecomposition + 1-D multiplier search; all solves are
    performed in the eigenbasis (gU = U^H g).  Handles PSD interior,
    singular PSD (solvable/unsolvable), negative-curvature boundary, and the
    trust-region hard case -- without adding epsilon to lambda_min
    denominators.
    """
    A = np.asarray(A, dtype=complex)
    g = np.asarray(g, dtype=complex).reshape(-1, 1)
    n = A.shape[0]
    if eig is None:
        lam, U = np.linalg.eigh(A)          # ascending
    else:
        lam, U = eig
    gU = (U.conj().T @ g).ravel()           # g in the eigenbasis
    gl = np.abs(gU)
    scale = max(1.0, float(np.abs(lam).max()) if lam.size else 1.0,
                float(np.linalg.norm(g)))
    tol = 1e-12 * scale
    lam_min = float(lam[0])
    checks = dict(ball_viol=np.nan, stat_res=np.nan, lam_min_Amu=np.nan,
                  comp=np.nan)

    if r <= 0.0:
        # the ball is the single point e = 0
        checks.update(ball_viol=0.0, stat_res=0.0, lam_min_Amu=lam_min,
                      comp=0.0)
        return dict(e=np.zeros((n, 1), dtype=complex), mu=0.0, mode='r_zero',
                    quad_val=0.0, checks=checks)

    mu_min = max(0.0, -lam_min)
    psd = lam_min >= -tol

    def phi_safe(mu):
        """phi(mu) = sum |gU_i|^2 / (lam_i+mu)^2 with pseudo-inverse
        convention at vanishing denominators (+inf if gU_i is nonzero
        there)."""
        den = lam + mu
        good = den > tol
        bad = ~good
        if bad.any() and (gl[bad] > 1e-10 * scale).any():
            return np.inf
        return float(np.sum(gl[good] ** 2 / den[good] ** 2))

    def e_of_mu(mu):
        return U @ (-(gU / (lam + mu)).reshape(-1, 1))

    if psd:
        mask_pi = lam > tol                    # range(A) coordinates
        eU = np.zeros(n, dtype=complex)
        eU[mask_pi] = -gU[mask_pi] / lam[mask_pi]
        e_pi = U @ eU.reshape(-1, 1)
        null_mask = ~mask_pi
        g_in_range = (not null_mask.any()
                      or (gl[null_mask] <= 1e-10 * scale).all())
        if g_in_range and float(np.linalg.norm(e_pi)) <= r * (1.0 + 1e-12):
            e, mu, mode = e_pi, 0.0, 'interior_psd'
        else:
            e, mu, mode = None, 0.0, 'boundary'
    else:
        mu_min = -lam_min
        mask_pi = lam > lam_min + tol          # non-minimum coordinates
        eU = np.zeros(n, dtype=complex)
        eU[mask_pi] = -gU[mask_pi] / (lam[mask_pi] - lam_min)
        e_reg = U @ eU.reshape(-1, 1)
        min_mask = ~mask_pi
        g_orth = (not min_mask.any()
                  or (gl[min_mask] <= 1e-10 * scale).all())
        phi_pi = float(np.sum(gl[mask_pi] ** 2 / (lam[mask_pi] - lam_min) ** 2))
        if g_orth and phi_pi <= r * r:
            # hard case (incl. the equality phi_pi = r^2): g ⟂
            # min-eigenspace, pseudo-inverse solution inside or on the
            # ball; complete any remaining norm inside the min-eigenspace
            fill = np.sqrt(max(r * r - float(np.linalg.norm(e_reg)) ** 2, 0.0))
            i0 = int(np.argmax(min_mask))
            e = e_reg + fill * U[:, [i0]]
            mu, mode = mu_min, 'hard_case_indefinite'
        else:
            e, mu, mode = None, mu_min, 'boundary'

    if mode == 'boundary':
        mu_lo, mu_hi = mu_min, max(1.0, scale) + mu_min
        for _ in range(1000):
            if phi_safe(mu_hi) <= r * r:
                break
            mu_hi *= 2.0
        bisect_used = 0
        for _ in range(bisect_iters):
            if (mu_rel_tol is not None
                    and (mu_hi - mu_lo) <= mu_rel_tol * max(1.0, mu_hi)):
                break
            mu_mid = 0.5 * (mu_lo + mu_hi)
            if phi_safe(mu_mid) > r * r:
                mu_lo = mu_mid
            else:
                mu_hi = mu_mid
            bisect_used += 1
        mu = mu_hi
        e = e_of_mu(mu)
        mode = 'boundary'

    quad = float(np.real(e.conj().T @ A @ e)[0, 0] + 2.0 * np.real(e.conj().T @ g)[0, 0])
    en = float(np.linalg.norm(e))
    Amu_e_g = (A + mu * np.eye(n)) @ e + g
    # Check the full residual, including the minimum eigenspace.
    stat = float(np.linalg.norm(Amu_e_g))
    checks.update(ball_viol=max(0.0, en - r),
                  stat_res=stat,
                  lam_min_Amu=lam_min + mu,
                  comp=abs(mu * (en - r)))
    return dict(e=e, mu=mu, mode=mode, quad_val=quad, checks=checks,
                bisect_used=bisect_used if mode == 'boundary' else 0)


def q_star_reference(scal, r, gamma, noise_power, eig=None, mu_rel_tol=None):
    """Exact q*(r) for one instance; also validity checks and tau."""
    A, g, c = scal['A'], scal['g'], scal['c']
    out = trust_region_min(A, g, r, eig=eig, mu_rel_tol=mu_rel_tol)
    q = c + out['quad_val']
    tau, scale_s = engineering_tau(scal, r)
    sc = max(1.0, abs(c), float(np.linalg.norm(g)), float(np.linalg.norm(A, 2)))
    ok = (out['mode'] == 'r_zero'
          or (out['checks']['ball_viol'] <= 1e-8 * max(1.0, r)
              and out['checks']['stat_res'] <= 1e-6 * sc
              and out['checks']['lam_min_Amu'] >= -1e-8 * sc
              and out['checks']['comp'] <= 1e-6 * sc))
    out.update(q=q, tau=tau, scale=scale_s, reference_valid=bool(ok))
    return out


def objective_at(b, h, C, gamma, noise_power, e):
    """|b^H(h+e)|^2 - gamma ||C^H(h+e)||^2 - gamma*noise_power at a given e."""
    b = np.asarray(b, complex).reshape(-1, 1)
    h = np.asarray(h, complex).reshape(-1, 1)
    C = np.asarray(C, complex)
    e = np.asarray(e, complex).reshape(-1, 1)
    sig = abs(b.conj().T @ (h + e))[0, 0] ** 2
    itf = float(np.linalg.norm(C.conj().T @ (h + e)) ** 2)
    return float(sig) - gamma * itf - gamma * noise_power


# ----------------------------------------------------------------------------
# 8. constructions (experiment A2 / A3)
# ----------------------------------------------------------------------------
def _psd_sqrt(M):
    w, V = np.linalg.eigh(M)
    w = np.maximum(w, 0.0)
    return (V * np.sqrt(w)) @ V.conj().T


def build_direction_construction(alpha, beta, beta0, x, y, eps,
                                 t_iters=200):
    """A2 completion: returns dict with d, ell, ell_e, T, K, B, C_H, t, R_val."""
    alpha, beta, beta0 = float(alpha), float(beta), float(beta0)
    delta = np.sqrt(max(beta0**2 - alpha**2, 0.0))
    nu = np.sqrt(max(beta0**2 - beta**2, 0.0))
    d = np.array([alpha, 0.0, 0.0], dtype=complex)
    k = float(np.clip(delta * x / y, nu, alpha))
    rem = beta0**2 - k**2 - delta**2 * k**2 / alpha**2
    ell = np.array([delta * k / alpha, np.sqrt(max(rem, 0.0)), 0.0],
                   dtype=complex)
    R_val = float(np.linalg.norm(x * d + y * ell))
    ell_e = (1.0 - eps) * ell
    L2 = np.outer(ell_e[:2], ell_e[:2].conj())      # plane restriction

    dd = np.outer(d, d.conj())

    def T_of_t(t):
        T = np.zeros((3, 3), dtype=complex)
        T[:2, :2] = (1.0 - t) * L2 + t * beta**2 * np.eye(2, dtype=complex)
        T[2, 2] = beta**2
        return T

    def f_t(t):
        return float(np.linalg.eigvalsh(dd + T_of_t(t))[-1]) - beta0**2

    f0, f1 = f_t(0.0), f_t(1.0)
    if f0 > 1e-12 * max(1.0, beta0**2) or f1 < -1e-12 * max(1.0, beta0**2):
        # bisection bracket fails; clamp to the end that is feasible
        t_star = 0.0 if f0 <= 0 else 1.0
        bracket_ok = False
    else:
        lo, hi = 0.0, 1.0
        for _ in range(t_iters):
            mid = 0.5 * (lo + hi)
            if f_t(mid) > 0.0:
                hi = mid
            else:
                lo = mid
        t_star = 0.5 * (lo + hi)
        bracket_ok = True
    T = T_of_t(t_star)
    M = T - np.outer(ell_e, ell_e.conj())
    lam_min_M = float(np.linalg.eigvalsh(M)[0])
    K = _psd_sqrt(M)
    B = np.hstack([ell_e.reshape(3, 1), K])
    C_H = np.hstack([d.reshape(3, 1), B])
    return dict(d=d, ell=ell, ell_e=ell_e, T=T, K=K, B=B, C_H=C_H,
                t=t_star, R_val=R_val, delta=delta, nu=nu, k=k,
                lam_min_M=lam_min_M, bracket_ok=bracket_ok,
                f0=f0, f1=f1)


def build_A3_instance(p, a, I0, alpha, beta, beta0, x_star, y_star, eps,
                      gamma, noise_power):
    """A3: full (b, h, C) from the A2 construction at (x*, y*)."""
    con = build_direction_construction(alpha, beta, beta0, x_star, y_star, eps)
    d, ell_e, B, C_H = con['d'], con['ell_e'], con['B'], con['C_H']
    n = C_H.shape[1]
    b = np.zeros((n, 1), dtype=complex)
    b[0, 0] = np.sqrt(p)
    u = -x_star * d - y_star * ell_e
    z = I0 * u / float(np.linalg.norm(u))
    rhs = z - a * d
    h_perp, *_ = np.linalg.lstsq(B, rhs, rcond=None)
    h = np.vstack([[a], h_perp.reshape(-1, 1)])
    C = C_H.conj().T
    return dict(b=b, h=h, C=C, con=con, u=u)


# ----------------------------------------------------------------------------
# 9. V4.1 constructions: output-space PSD comparison + planar spectral
#    completion (paper V4.1, Theorems 2 and 3)
# ----------------------------------------------------------------------------
def construct_v41(alpha, beta, beta0, x, y, eps, p, a, I0,
                  bisect_iters=100):
    """Theorem 2 construction exactly as proven in V4.1.

    Steps (paper EQ 13-20):
      d = (alpha,0,0), ell = (delta k/alpha, beta0 sqrt(1-k^2/alpha^2), 0)
      with k = alpha if y = 0 else clip(delta x / y, nu, alpha);
      output-space PSD ordering  dd^H + l_e l_e^H <= dd^H + l l^H <= b0^2 I
      (NOT a Gram-matrix monotonicity, which is false in general);
      planar completion S(t) = d2 d2^H + (1-t) l_e2 l_e2^H + t beta^2 I2,
      bisected for lambda_max(S(t*)) = beta0^2 -- the PLANAR eigenvalue is
      essential when beta = beta0, where the full 3x3 maximum stays constant
      over an initial t interval;
      T_e = diag((1-t*) l_e2 l_e2^H + t* beta^2 I2, beta^2);
      M_e = T_e - l_e l_e^H = diag(t*(beta^2 I2 - l_e2 l_e2^H), beta^2) > 0;
      K_e = M_e^{1/2} (invertible), B_e = [l_e, K_e], C^H = [d, B_e];
      b = sqrt(p) e1, z = I0 u/||u|| with u = -x d - y l_e, B_e h_perp =
      z - a d, h = (a, h_perp), e = (-x, -y, 0, 0, 0).

    Requires alpha > 0, beta0 > alpha, beta > 0, 0 < eps < 1, (x, y) != 0.
    Returns a dict of parts and diagnostics; caller evaluates the instance.
    """
    if alpha <= 0:
        raise ValueError("construct_v41 requires alpha > 0 (see V4.1 "
                         "zero-alpha construction)")
    if not beta0 > alpha:
        raise ValueError("construct_v41 requires beta0 > alpha "
                         "(boundary instances use construct_boundary_v41)")
    if beta <= 0:
        raise ValueError("construct_v41 requires beta > 0")
    if not (0.0 < eps < 1.0):
        raise ValueError("construct_v41 requires 0 < eps < 1")
    if x < 0 or y < 0 or (x == 0 and y == 0):
        raise ValueError("construct_v41 requires x, y >= 0 not both zero")
    delta = np.sqrt(beta0 ** 2 - alpha ** 2)
    nu = np.sqrt(max(beta0 ** 2 - beta ** 2, 0.0))
    k = alpha if y == 0 else float(np.clip(delta * x / y, nu, alpha))
    d = np.array([alpha, 0.0, 0.0], dtype=complex)
    ell = np.array([delta * k / alpha,
                    beta0 * np.sqrt(max(1.0 - (k / alpha) ** 2, 0.0)),
                    0.0], dtype=complex)
    le = (1.0 - eps) * ell
    d2, le2 = d[:2], le[:2]

    def lam_max_S(t):
        S = (np.outer(d2, d2.conj()) + (1.0 - t) * np.outer(le2, le2.conj())
             + t * beta ** 2 * np.eye(2, dtype=complex))
        return float(np.linalg.eigvalsh(S)[-1])

    S0, S1 = lam_max_S(0.0), lam_max_S(1.0)
    if not S0 < beta0 ** 2:
        raise RuntimeError(f"planar completion bracket failed: S(0)={S0!r} "
                           f">= beta0^2={beta0**2!r}")
    if not S1 >= beta0 ** 2 - 1e-12 * max(1.0, beta0 ** 2):
        raise RuntimeError(f"planar completion upper endpoint failed: "
                           f"S(1)={S1!r}")
    lo, hi = 0.0, 1.0
    for _ in range(bisect_iters):
        mid = 0.5 * (lo + hi)
        if lam_max_S(mid) > beta0 ** 2:
            hi = mid
        else:
            lo = mid
    t_eps = 0.5 * (lo + hi)
    T = np.zeros((3, 3), dtype=complex)
    T[:2, :2] = ((1.0 - t_eps) * np.outer(le2, le2.conj())
                 + t_eps * beta ** 2 * np.eye(2, dtype=complex))
    T[2, 2] = beta ** 2
    M = T - np.outer(le, le.conj())
    w_M, V_M = np.linalg.eigh(M)
    lam_min_M = float(w_M[0])
    if lam_min_M <= 0:
        raise RuntimeError(f"M_epsilon not positive definite: "
                           f"lambda_min={lam_min_M!r}")
    K = (V_M * np.sqrt(w_M)) @ V_M.conj().T
    B = np.column_stack((le, K))
    C_H = np.column_stack((d, B))
    n = C_H.shape[1]
    b = np.zeros((n, 1), dtype=complex)
    b[0, 0] = np.sqrt(p)
    u = -x * d - y * le
    u_norm = float(np.linalg.norm(u))
    if u_norm == 0:
        raise RuntimeError("constructed u_epsilon has zero norm")
    z = I0 * u / u_norm
    h_perp, *_ = np.linalg.lstsq(B, z - a * d, rcond=None)
    h = np.vstack([[a], h_perp.reshape(-1, 1)])
    e_vec = np.array([-x, -y, 0.0, 0.0, 0.0], dtype=complex).reshape(-1, 1)
    return dict(b=b, h=h, C=C_H.conj().T, e=e_vec, d=d, ell=ell, ell_e=le,
                T=T, K=K, B=B, C_H=C_H, t_eps=t_eps, k=k, delta=delta,
                nu=nu, lam_min_M=lam_min_M, S0_minus_b0sq=S0 - beta0 ** 2,
                S1_minus_b0sq=S1 - beta0 ** 2, u=u, R_eps=u_norm)


def construct_boundary_v41(p, a, I0, alpha, beta, gamma, noise_power, r):
    """Theorem 3 boundary construction (paper V4.1 EQ 27).

    Requires beta0 = alpha > 0, 0 < beta <= alpha, I0 >= alpha*a.
    Returns the instance, the attaining error (phase chosen by the sign of
    c_parallel = p - gamma alpha^2) and the boundary bound data.
    """
    if alpha <= 0 or beta <= 0 or beta > alpha:
        raise ValueError("boundary construction requires 0 < beta <= alpha")
    if I0 < alpha * a - 1e-12 * max(1.0, I0):
        raise ValueError("boundary compatibility requires I0 >= alpha*a")
    Jperp = float(np.sqrt(max(I0 ** 2 - alpha ** 2 * a ** 2, 0.0)))
    c_par = p - gamma * alpha ** 2
    C = np.diag([alpha, beta]).astype(complex)
    h = np.array([[a], [Jperp / beta]], dtype=complex)
    b = np.array([[np.sqrt(p)], [0.0]], dtype=complex)
    scal = dict(p=p, a=a, I0=I0, alpha=alpha, beta=beta, beta0=alpha)
    q_partial, x_star = Q_boundary_approx(scal, r, gamma, noise_power)
    y_star = float(np.sqrt(max(r * r - x_star * x_star, 0.0)))
    x_e = -x_star if c_par >= 0 else x_star
    e = np.array([[x_e], [y_star]], dtype=complex)
    return dict(b=b, h=h, C=C, e=e, scal=scal, Jperp=Jperp, c_parallel=c_par,
                Q_partial=q_partial, x_star=x_star, y_star=y_star)

