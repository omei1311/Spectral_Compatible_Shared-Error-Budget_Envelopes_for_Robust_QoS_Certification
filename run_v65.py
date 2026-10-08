"""Reproduce V6.5 validation, measured query batches and actual fallback.
Run from the project root. Historical result directories are not modified.
"""
import csv, json, hashlib, time, contextlib, io
from pathlib import Path
from collections import defaultdict
import numpy as np
import core
import experiments_c as ec
import v65_clients as sc

OUT = Path('results_v65')
GAMMA, NOISE, J = 1., .05, 256
RHO = np.round(np.logspace(-2.2, -.6, 10), 5)
REPS = 5

def write_csv(name, rows):
    with (OUT/name).open('w', newline='', encoding='utf-8') as f:
        writer=csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

def six_of(s):
    return tuple(s[k] for k in ('p','a','I0','alpha','beta','beta0'))

def validate():
    rows=[]
    def record(name, ok, **details):
        rows.append(dict(case=name, ok=bool(ok), details=json.dumps(details)))
        if not ok: raise AssertionError((name,details))
    # Analytic family includes small nonzero linear terms around the hard case.
    for eps in (0.,1e-12,1e-9,1e-6,1e-5):
        n=8; b=np.zeros((n,1),complex); b[1]=1
        C=np.zeros((n,3),complex); C[0,0]=1
        h=np.zeros((n,1),complex); h[0]=-eps; h[1]=1
        cache=sc.build_compact_cache(b,h,C,1.,.1399999)
        bounds=sc.tr_bounds(cache,.6)
        e=np.zeros((n,1),complex); e[0]=-np.sqrt(.11); e[1]=-.5
        witness=core.objective_at(b,h,C,1.,.1399999,e)
        record(f'near_hard_{eps}',bounds['valid'] and bounds['lower'] <= witness+1e-13,
               lower=bounds['lower'],upper=bounds['upper'],witness=witness,mode=bounds['mode'])
        if eps==1e-6:
            record('no_false_accept',sc.decision(bounds,1.2000012e-8)=='rejected',**bounds)
    # Scalar diagonal cases: values independently known from elementary minimization.
    for tag,lam,w,c,r,exact in [
        ('psd_interior',[1,1,1],[.01,.04,0],1,5,.95),
        ('singular',[0,1],[0,1],1,5,0),
        ('hard_zero',[-1,1],[0,0],1,2,-3),
        ('hard_equality',[-1,1],[0,1],1,.5,.25),
        ('repeated_min',[-1,-1,2],[0,0,1],1,.7,1-.49-1/3),
        ('r_zero',[-2,1],[1,1],3,0,3),
        ('tiny_nonzero_psd',[0,1],[1e-12,0],1,.6,1-1.2e-6)]:
        ca=dict(lam=np.array(lam,float),w=np.array(w,float),c=c,scale=2.)
        z=sc.tr_bounds(ca,r)
        record(tag,z['valid'] and z['lower']-1e-12<=exact<=z['upper']+1e-12,exact=exact,**z)
    worst=0.; count=0
    for seed in (401,402):
        for gamma in (.5,1.,2.):
            for cand in ec.all_candidates(seed):
                for u in range(4):
                    b,h,C=cand['b'][u],cand['h'][u],cand['C'][u]
                    s=core.spectral_summary(b,h,C,gamma,NOISE); six=six_of(s)
                    cache=sc.build_compact_cache(b,h,C,gamma,NOISE)
                    for rho in (.00631,.15,.25119):
                        r=rho*np.linalg.norm(h); z=sc.tr_bounds(cache,r)
                        ref=core.q_star_reference(s,r,gamma,NOISE,mu_rel_tol=1e-12)
                        diff=abs(z['q']-ref['q']);worst=max(worst,diff);count+=1
                        record(f'pool_{seed}_{gamma}_{count}',z['valid'] and ref['reference_valid'] and diff<1e-7,
                               diff=diff,gap=z['gap'])
                        assert abs(sc.lquad_client(six,np.linalg.norm(s['g']),r,gamma,NOISE)-core.L_quad(s,r,gamma))<1e-11
                        assert abs(sc.qg_lb_client(six,r,gamma,NOISE,J)-core.interval_lower_bound(s,r,gamma,NOISE,J,'G'))<1e-11
    write_csv('validation.csv',rows)
    return dict(checks=len(rows),pool_queries=count,max_reference_difference=worst)

def main():
    OUT.mkdir(exist_ok=True)
    config=dict(seeds=[401,402],gamma=GAMMA,noise=NOISE,J=J,rho=list(RHO),
                batch_lengths=[1,10,100],repetitions=REPS,
                timestamp=time.strftime('%Y-%m-%dT%H:%M:%S'),
                hashes={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in ('run_v65.py','v65_clients.py')},
                task='cached descriptor queries; actual user-level fallback; no network latency',
                timing='Each batch visits all 48 users; actual queries = users * batch_length. Shared tau excluded from query timings; included in pipelines.')
    (OUT/'config.json').write_text(json.dumps(config,indent=2),encoding='utf-8')
    validation=validate()
    builds=[]; timings=[]; details=[]; pipe_times=[]; units_rows=[]
    for seed in (401,402):
        cands=ec.all_candidates(seed)
        inputs=[(c['b'][u],c['h'][u],c['C'][u]) for c in cands for u in range(4)]
        # Matched preprocessing stages, independently timed from identical raw data.
        bt=defaultdict(list)
        def build(branch):
            if branch=='D': return [sc.build_compact_cache(b,h,C,GAMMA,NOISE) for b,h,C in inputs]
            out=[]
            for b,h,C in inputs:
                s=core.spectral_summary(b,h,C,GAMMA,NOISE)
                # Existing common builder also forms A and g; explicitly disclosed, not a minimal six-scalar builder.
                out.append((six_of(s),float(np.linalg.norm(s['g'])) if branch=='B' else None))
            return out
        for branch in ('A','B','D'):build(branch)
        for rep in range(REPS):
            for branch in np.roll(np.array(['A','B','D']),rep):
                start=time.perf_counter();built=build(branch);bt[branch].append(time.perf_counter()-start)
        a=build('A');seven=build('B');caches=build('D');sixes=[x[0] for x in a]; kappas=[x[1] for x in seven]
        for branch in bt:
            for rep,secs in enumerate(bt[branch]):builds.append(dict(seed=seed,branch=branch,rep=rep,seconds=secs,users=len(inputs)))
        radii=[[float(rho*np.linalg.norm(h)) for rho in RHO] for b,h,C in inputs]
        def query(branch,i,r):
            if branch=='A':return sc.qg_lb_client(sixes[i],r,GAMMA,NOISE,J)
            if branch=='B':return sc.lquad_client(sixes[i],kappas[i],r,GAMMA,NOISE)
            if branch=='C':return sc.joint_client(sixes[i],kappas[i],r,GAMMA,NOISE,J)
            return sc.tr_bounds(caches[i],r)
        def batch(branch,nq):
            calls=0;start=time.perf_counter()
            for i in range(len(inputs)):
                for j in range(nq): query(branch,i,radii[i][j%len(RHO)]);calls+=1
            secs=time.perf_counter()-start
            assert calls==len(inputs)*nq
            return secs,calls
        for branch in 'ABCD':batch(branch,10)
        for rep in range(REPS):
            for branch in np.roll(np.array(list('ABCD')),rep):
                for nq in (1,10,100):
                    secs,calls=batch(branch,nq)
                    timings.append(dict(seed=seed,rep=rep,branch=branch,batch_length=nq,actual_queries=calls,seconds=secs))
        expected={}
        for ci in range(len(cands)):
            for j,rho in enumerate(RHO):
                statuses=[]
                for u in range(4):
                    i=ci*4+u;r=radii[i][j]
                    s=core.spectral_summary(*inputs[i],GAMMA,NOISE)
                    ref=core.q_star_reference(s,r,GAMMA,NOISE,mu_rel_tol=1e-12)
                    assert ref['reference_valid']
                    tau=sc.tau_summary(*[sixes[i][k] for k in (0,1,2,5)],GAMMA,NOISE,r)
                    z=query('D',i,r);st=sc.decision(z,tau)
                    rst='accepted' if ref['q']>tau else 'rejected' if ref['q'] < -tau else 'unresolved'
                    assert st==rst
                    statuses.append(st)
                    details.append(dict(seed=seed,cand=ci,user=u,rho=rho,r=r,tau=tau,
                                        Q_G=query('A',i,r),L_quad=query('B',i,r),joint=query('C',i,r),
                                        compact_lower=z['lower'],compact_upper=z['upper'],q_star=ref['q'],status=st))
                expected[ci,j]='rejected' if 'rejected' in statuses else 'accepted' if all(s=='accepted' for s in statuses) else 'unresolved'
        def pipeline(branch):
            rows=[];start=time.perf_counter()
            for ci in range(len(cands)):
                for j,rho in enumerate(RHO):
                    ntr=ns=0;unresolved=False;direct=True;final='accepted'
                    for u in range(4):
                        i=ci*4+u;r=radii[i][j]
                        tau=sc.tau_summary(*[sixes[i][k] for k in (0,1,2,5)],GAMMA,NOISE,r)
                        if branch!='D':
                            bound=query(branch,i,r);ns+=1
                            if bound>tau:continue
                            direct=False
                        # Actual fallback invocation, never read from the coverage table.
                        z=query('D',i,r);ntr+=1;st=sc.decision(z,tau)
                        if st=='rejected':final='rejected';break
                        if st=='unresolved':unresolved=True
                    if final=='accepted' and unresolved:final='unresolved'
                    rows.append(dict(seed=seed,branch=branch,cand=ci,rho=rho,tr_calls=ntr,
                                     summary_calls=ns,direct=direct if branch!='D' else False,final=final))
            seconds=time.perf_counter()-start
            assert all(r['final']==expected[r['cand'],list(RHO).index(r['rho'])] for r in rows)
            return seconds,rows
        for branch in 'ABCD':pipeline(branch)
        for rep in range(REPS):
            for branch in np.roll(np.array(list('ABCD')),rep):
                secs,rows=pipeline(branch)
                pipe_times.append(dict(seed=seed,rep=rep,branch=branch,seconds=secs,
                                       actual_tr_calls=sum(r['tr_calls'] for r in rows),actual_summary_calls=sum(r['summary_calls'] for r in rows)))
                if rep==0:units_rows.extend(rows)
        print('completed seed',seed,flush=True)
    for name,rows in [('build_timings.csv',builds),('batch_timings.csv',timings),('per_problem.csv',details),('pipeline_timings.csv',pipe_times),('fallback_units.csv',units_rows)]:write_csv(name,rows)
    summary=dict(validation=validation,units=240,reference_feasible=sum(r['final']=='accepted' for r in units_rows if r['branch']=='D'),bytes=dict(A=48,B=56,C=56,D=sc.compact_cache_bytes(caches[0])),coverage={},query_us={},pipeline={},build_ms={})
    for branch in 'ABCD':
        urs=[r for r in units_rows if r['branch']==branch]
        summary['coverage'][branch]=sum(r['direct'] for r in urs)
        summary['query_us'][branch]={str(n):float(np.median([1e6*t['seconds']/t['actual_queries'] for t in timings if t['branch']==branch and t['batch_length']==n])) for n in (1,10,100)}
        summary['pipeline'][branch]=dict(median_seconds_per_seed=float(np.median([r['seconds'] for r in pipe_times if r['branch']==branch])),tr_calls_total=sum(r['tr_calls'] for r in urs),unresolved=sum(r['final']=='unresolved' for r in urs))
    for branch in ('A','B','D'):summary['build_ms'][branch]=float(np.median([1e3*r['seconds']/r['users'] for r in builds if r['branch']==branch]))
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
