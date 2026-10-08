"""Regression checks for boundary branches, refinement and final decisions."""
import numpy as np
import core
import experiments_c as ec
import audit_and_stress as study


def verify():
    for p, I0 in [(2., 1.), (.2, 1.), (2., 1.3), (.2, 1.3)]:
        C = np.diag([1., .7, 0.]).astype(complex)
        h = np.array([[1.], [np.sqrt(I0**2-1.)/.7], [0.]], complex)
        b = np.array([[np.sqrt(p)], [0.], [0.]], complex)
        s = core.spectral_summary(b, h, C, 1., .05)
        qb = core.Q_boundary_approx(s, .2, 1., .05)[0]
        ref = core.q_star_reference(s, .2, 1., .05)
        assert ref['reference_valid'] and abs(qb-ref['q']) < 1e-9
        lb, J = core.lb_with_escalation(s, .2, 1., .05,
                                        J_list=(16, 64, 256), ub_gap_tol=1e-5)
        assert lb <= core.Q_approx(s, .2, 1., .05)[0] + 1e-10
        assert J in (16, 64, 256)
        general = dict(s, beta0=1.1)
        _, refined_J = core.lb_with_escalation(general, .2, 1., .05,
                                               J_list=(16, 64, 256), ub_gap_tol=1e-12)
        assert refined_J == 256
        print('Boundary branch verified:', p, I0, qb, flush=True)
    for seed in study.SEEDS:
        cands = ec.all_candidates(seed)
        for gamma in study.GAMMAS:
            truth = study.pipeline(cands, gamma, 'reference')['decisions']
            for method in ec.METHODS:
                assert study.pipeline(cands, gamma, method)['decisions'] == truth
        print('Final decisions verified for seed', seed, flush=True)
    print('PASS: four boundary cases; refinement; 1080 per-unit decisions for each method')


if __name__ == '__main__':
    verify()
