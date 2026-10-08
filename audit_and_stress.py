"""Independent reference checks and a disclosed fixed-grid sensitivity study.

Original results/ files are preserved. No candidate is selected by Q_G gain.
Run: python audit_and_stress.py
"""
import os
from time import perf_counter
import numpy as np
import core
import experiments_c as ec
import utils

OUT = os.path.join(os.path.dirname(__file__), 'results_review')
GAMMAS = (0.5, 1.0, 2.0)
RHOS = (0.06, 0.15, 0.30)
SEEDS = (11, 22, 33, 44, 55, 66, 77, 88, 99, 110)


def audit_reference():
    # Includes singular PSD, repeated minimum eigenvalues, and the exact
    # hard-case equality that otherwise evaluates a 0/0 denominator.
    cases = [([-1., 1.], [0., 0.], 1.),
             ([-1., 1.], [0., 2.], 1.),
             ([-1., -1., 2.], [0., 0., 1.], 1.),
             ([0., 1.], [0., 0.5], 1.),
             ([0., 1.], [1., 0.5], 1.),
             ([1., 2.], [0.1, 0.2], 1.)]
    rng = np.random.default_rng(20261007)
    instances = [(np.diag(l).astype(complex), np.array(g, complex).reshape(-1, 1), r)
                 for l, g, r in cases]
    for _ in range(100):
        Z = rng.normal(size=(5, 5)) + 1j*rng.normal(size=(5, 5))
        A = (Z + Z.conj().T)/2
        g = (rng.normal(size=5) + 1j*rng.normal(size=5)).reshape(-1, 1)
        instances.append((A, g, rng.uniform(0.1, 2.)))
    rows = []
    for idx, (A, g, r) in enumerate(instances):
        out = core.trust_region_min(A, g, r)
        lam, U = np.linalg.eigh(A)
        z = (U.conj().T @ g).ravel()
        den = lam + out['mu']
        keep = den > 1e-11
        dual = -out['mu']*r*r - np.sum(np.abs(z[keep])**2/den[keep])
        gap = out['quad_val'] - dual
        assert np.isfinite(out['quad_val'])
        assert abs(gap) < 1e-8
        assert out['checks']['stat_res'] < 1e-8
        assert np.linalg.norm(out['e']) <= r + 1e-9
        rows.append(dict(idx=idx, mode=out['mode'], primal=out['quad_val'],
                         dual=dual, gap=gap, **out['checks']))
    utils.write_csv(os.path.join(OUT, 'reference_audit.csv'), list(rows[0]), rows)
    print('Reference audits:', len(rows), 'max dual gap', max(abs(x['gap']) for x in rows))


def summarize(cands, gamma):
    summaries = [[core.spectral_summary(c['b'][u], c['h'][u], c['C'][u], gamma, ec.NOISE_U)
                  for u in range(ec.N_USER)] for c in cands]
    norms = [[float(np.linalg.norm(h)) for h in c['h']] for c in cands]
    return summaries, norms


def pipeline(cands, gamma, method, mu_rel_tol=None):
    """Same final decision task for every method; reuse eigh across radii.

    Full-reference baseline also uses user early stopping. Filtered runs
    skip certified users during fallback and validate each reference solve.
    Summary construction is inside the timer; candidate generation is not.
    ``mu_rel_tol`` (optional) forwards a tolerance stop for the trust-region
    multiplier bisection (V4.1 efficiency review); default keeps the fixed
    200-iteration behavior.
    """
    start = perf_counter()
    calls = accepted = rejected = unresolved = direct = 0
    bisect_total = 0
    decisions = []
    for cand in cands:
        cache = {}
        summaries = []
        norms = []
        for u in range(ec.N_USER):
            b, h, C = cand['b'][u], cand['h'][u], cand['C'][u]
            if method == 'reference':
                A = b @ b.conj().T - gamma*(C @ C.conj().T)
                summaries.append(dict(A=A, g=A@h,
                                      c=float(np.real((h.conj().T@A@h)[0, 0])) - gamma*ec.NOISE_U))
            else:
                summaries.append(core.spectral_summary(b, h, C, gamma, ec.NOISE_U))
            norms.append(float(np.linalg.norm(h)))
        for rho in RHOS:
            certified = [False]*ec.N_USER
            if method != 'reference':
                for u, s in enumerate(summaries):
                    r = rho*norms[u]
                    tau = core.engineering_tau(s, r)[0]
                    if method == 'D':
                        lb = core.D_bound(s, r, gamma, ec.NOISE_U)
                    elif method == 'L_quad':
                        lb = core.L_quad(s, r, gamma)
                    else:
                        lb = core.interval_lower_bound(s, r, gamma, ec.NOISE_U,
                                                       256, 'G' if method == 'Q_G' else 'old')
                    certified[u] = lb > tau
            if all(certified):
                direct += 1
                accepted += 1
                decisions.append('accepted')
                continue
            statuses = []
            for u, s in enumerate(summaries):
                if certified[u]:
                    statuses.append('clearly_feasible')
                    continue
                if u not in cache:
                    cache[u] = np.linalg.eigh(s['A'])
                out = core.q_star_reference(s, rho*norms[u], gamma, ec.NOISE_U,
                                            eig=cache[u], mu_rel_tol=mu_rel_tol)
                calls += 1
                bisect_total += out.get('bisect_used', 0)
                statuses.append(ec.user_status(out))
                if statuses[-1] == 'clearly_infeasible':
                    break
            if 'clearly_infeasible' in statuses:
                rejected += 1
                decisions.append('rejected')
            elif all(x == 'clearly_feasible' for x in statuses):
                accepted += 1
                decisions.append('accepted')
            else:
                unresolved += 1
                decisions.append('unresolved')
    return dict(seconds=perf_counter()-start, calls=calls, accepted=accepted,
                rejected=rejected, unresolved=unresolved, direct=direct,
                bisect_total=bisect_total, decisions=decisions)


def run_stress():
    details, summaries, timings = [], [], []
    for seed in SEEDS:
        cands = ec.all_candidates(seed)
        for gamma in GAMMAS:
            ss, norms = summarize(cands, gamma)
            rows = []
            for i, cand in enumerate(cands):
                for u, s in enumerate(ss[i]):
                    eig = np.linalg.eigh(s['A'])
                    for rho in RHOS:
                        r = rho*norms[i][u]
                        ref = core.q_star_reference(s, r, gamma, ec.NOISE_U, eig=eig)
                        row = dict(seed=seed, gamma=gamma, cand_idx=i, user=u, rho=rho,
                                   btype=cand['btype'], q_star=ref['q'], tau=ref['tau'],
                                   reference_valid=ref['reference_valid'], status=ec.user_status(ref),
                                   D=core.D_bound(s, r, gamma, ec.NOISE_U),
                                   L_quad=core.L_quad(s, r, gamma),
                                   Q_old=core.interval_lower_bound(s, r, gamma, ec.NOISE_U, 256, 'old'),
                                   Q_G=core.interval_lower_bound(s, r, gamma, ec.NOISE_U, 256, 'G'))
                        assert ref['reference_valid'], row
                        assert max(row[m] for m in ec.METHODS) <= ref['q'] + 1e-8, row
                        rows.append(row)
            details.extend(rows)
            for rho in RHOS:
                agg, _ = ec.aggregate([row for row in rows if row['rho'] == rho])
                for method in ec.METHODS:
                    summaries.append(dict(seed=seed, cohort='original' if seed <= 55 else 'holdout',
                                          gamma=gamma, rho=rho, method=method,
                                          feasible=agg['ref_feasible_units'], infeasible=agg['ref_infeasible_units'],
                                          ambiguous=agg['ambiguous_units'], invalid=agg['invalid_units'],
                                          qg_added=agg['qg_added'], **agg['stats'][method]))
            # Warm up and interleave methods; time the actual final decisions.
            methods = ('reference',) + ec.METHODS
            for method in methods:
                pipeline(cands[:1], gamma, method)
            runs = {m: [] for m in methods}
            for rep in range(3):
                for method in np.roll(np.array(methods), rep):
                    runs[method].append(pipeline(cands, gamma, method))
            baseline = runs['reference'][0]
            for method in methods:
                rs = runs[method]
                assert all(x['decisions'] == baseline['decisions'] for x in rs)
                timings.append(dict(seed=seed, gamma=gamma, method=method,
                                    time_median=np.median([x['seconds'] for x in rs]),
                                    time_min=min(x['seconds'] for x in rs),
                                    time_max=max(x['seconds'] for x in rs),
                                    calls=rs[0]['calls'], direct=rs[0]['direct'],
                                    accepted=rs[0]['accepted'], rejected=rs[0]['rejected'],
                                    unresolved=rs[0]['unresolved']))
        print('Stress completed seed', seed, flush=True)
    for name, rows in [('stress_details', details), ('stress_summary', summaries), ('stress_timings', timings)]:
        utils.write_csv(os.path.join(OUT, name+'.csv'), list(rows[0]), rows)
    utils.write_json(os.path.join(OUT, 'stress_config.json'),
                     dict(gammas=GAMMAS, rhos=RHOS, seeds=SEEDS, J=256, repetitions=3,
                          purpose='Fixed exploratory sensitivity grid; report all cells; holdout seeds 66..110',
                          timing='End-to-end final decision, cached eigh, certified-user skipping, infeasible early stop'))
    print('Stress problems:', len(details))


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    audit_reference()
    run_stress()
