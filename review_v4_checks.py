"""Read-only numerical checks for the v4 review. Original experiments untouched."""
import csv
from collections import Counter, defaultdict
import numpy as np
import core
import experiments_c as ec


def main():
    con = core.build_direction_construction(1., 1., 1.1, .45, .1, .01)
    d, ell, ell_e = con['d'], con['ell'], con['ell_e']
    Z = np.column_stack((d, ell))
    Ze = np.column_stack((d, ell_e))
    print('eig(G-G_epsilon):', np.linalg.eigvalsh(Z.conj().T@Z-Ze.conj().T@Ze))
    print('eig(output rank-one difference):',
          np.linalg.eigvalsh(np.outer(ell, ell.conj())-np.outer(ell_e, ell_e.conj())))
    with open('results_review/reproduction/pool_user_details.csv', encoding='utf-8') as f:
        original = list(csv.DictReader(f))
    with open('results_review/stress_details.csv', encoding='utf-8') as f:
        stress = list(csv.DictReader(f))
    lookup = {(int(r['seed']), int(r['cand_idx']), int(r['user']),
               float(r['rho']), float(r['gamma'])): r for r in stress}
    old_keys = {(int(r['seed']), int(r['cand_idx']), int(r['user']), float(r['rho']), .5)
                for r in original}
    print('C overlap:', len(old_keys & lookup.keys()),
          'unique combined C problems:', len(old_keys | lookup.keys()))
    counts = defaultdict(Counter)
    for seed in (11,22,33,44,55,66,77,88,99,110):
        for idx, cand in enumerate(ec.all_candidates(seed)):
            for gamma in (.5,1.,2.):
                ss = [core.spectral_summary(cand['b'][u],cand['h'][u],cand['C'][u],gamma,.05)
                      for u in range(4)]
                norms = [np.linalg.norm(h) for h in cand['h']]
                for rho in (.06,.15,.30):
                    refs = [lookup[(seed,idx,u,rho,gamma)] for u in range(4)]
                    for J in (256,1024,4096):
                        acc = {}
                        for kind in ('old','G'):
                            lbs = [core.interval_lower_bound(s,rho*norms[u],gamma,.05,J,kind)
                                   for u,s in enumerate(ss)]
                            assert all(lb<=float(refs[u]['q_star'])+1e-8 for u,lb in enumerate(lbs))
                            acc[kind] = all(lb>float(refs[u]['tau']) for u,lb in enumerate(lbs))
                            counts[J][kind] += acc[kind]
                        counts[J]['added'] += acc['G'] and not acc['old']
                        counts[J]['feasible'] += all(r['status']=='clearly_feasible' for r in refs)
    print('Resolution ablation:', dict(counts))
    n = 0
    example = None
    for seed in ec.SEEDS:
        for idx,cand in enumerate(ec.all_candidates(seed)):
            for u in range(4):
                s = core.spectral_summary(cand['b'][u],cand['h'][u],cand['C'][u],.5,.05)
                for rho in ec.RHO_LIST:
                    r = rho*np.linalg.norm(cand['h'][u])
                    lq = core.L_quad(s,r,.5)
                    q,x,y = core.Q_approx(s,r,.5,.05)
                    if lq > q+1e-9:
                        n += 1
                        if example is None:
                            example = dict(seed=seed,candidate=idx,user=u,rho=rho,
                                           L_quad=float(lq),evaluated_QG_objective=q,x=x,y=y)
    print('L_quad above Q_G evaluated upper estimate:',n,'/720; witness:',example)


if __name__ == '__main__':
    main()
