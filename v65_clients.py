"""V6.5 compact TR value bounds. V6.4 is retained unchanged.

No nonzero spectral weight is rounded to zero. A dual lower value and a
feasible spectral-coordinate upper value guard decisions near hard cases.
Bounds are floating-point diagnostics, not interval-arithmetic proofs.
"""
import numpy as np
from summary_clients import (qg_lb_client, lquad_client, joint_client,
                             tau_summary, build_compact_cache, compact_cache_bytes)


def tr_bounds(cache, r, rel_tol=1e-12, max_iters=200):
    lam = np.asarray(cache['lam'], float)
    w = np.asarray(cache['w'], float)
    c = float(cache['c'])
    if r < 0 or np.any(w < 0) or not np.all(np.isfinite(lam)) or not np.all(np.isfinite(w)):
        raise ValueError('Invalid radius or spectral cache')
    scale = max(1., abs(c), r*r*float(np.max(np.abs(lam))),
                2*r*float(np.sqrt(np.sum(w))))
    guard = 64*np.finfo(float).eps*scale
    def finish(lower, upper, mode, iterations):
        lower, upper = float(lower)-guard, float(upper)+guard
        valid = np.isfinite(lower) and np.isfinite(upper) and lower <= upper
        return dict(lower=lower, upper=upper, q=(lower+upper)/2,
                    valid=bool(valid), mode=mode, iterations=iterations,
                    gap=upper-lower)
    if r == 0:
        return finish(c, c, 'r_zero', 0)
    active = w > 0.  # Exact zero only; squared weights are not compared to amplitude tolerances.
    amin, aw = lam[active], w[active]
    mu0 = max(0., -float(lam[0]))
    def norm2(mu):
        den = amin+mu
        if np.any(den <= 0):
            return np.inf
        with np.errstate(over='ignore', divide='ignore'):
            return float(np.sum(aw/den**2))
    def values(mu, complete=False):
        den = amin+mu
        z2 = aw/den**2
        dual = c-mu*r*r-float(np.sum(aw/den))
        upper = c+float(np.sum(amin*z2-2*aw/den))
        if complete:
            upper += float(lam[0])*max(0., r*r-float(np.sum(z2)))
        return dual, upper
    if norm2(mu0) <= r*r:
        lower, upper = values(mu0, complete=mu0 > 0)
        return finish(lower, upper, 'hard_case' if mu0 > 0 else 'interior_psd', 0)
    lo = mu0
    hi = mu0+max(1., float(cache['scale']), float(np.sqrt(np.sum(w)))/r)
    for _ in range(1000):
        if norm2(hi) <= r*r:
            break
        hi = mu0+2*(hi-mu0)
    else:
        raise ArithmeticError('Unable to bracket TR multiplier')
    for iterations in range(1, max_iters+1):
        lower, upper = values(hi)
        if upper-lower <= rel_tol*scale:
            break
        mid = lo+(hi-lo)/2
        if mid == lo or mid == hi:
            break
        if norm2(mid) > r*r:
            lo = mid
        else:
            hi = mid
    lower, upper = values(hi)
    return finish(lower, upper, 'boundary', iterations)


def decision(bounds, tau):
    if not bounds['valid']:
        return 'unresolved'
    if bounds['lower'] > tau:
        return 'accepted'
    if bounds['upper'] < -tau:
        return 'rejected'
    return 'unresolved'
