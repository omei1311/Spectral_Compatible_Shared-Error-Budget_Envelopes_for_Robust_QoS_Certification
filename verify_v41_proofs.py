"""Numerical checks of the analytical constructions used in V4.1."""
import json
from pathlib import Path
import numpy as np
import core


def construct(p,a,I0,alpha,beta,beta0,x,y,eps):
    delta=np.sqrt(max(beta0**2-alpha**2,0))
    nu=np.sqrt(max(beta0**2-beta**2,0))
    k=alpha if y==0 else float(np.clip(delta*x/y,nu,alpha))
    d=np.array([alpha,0.,0.],complex)
    ell=np.array([delta*k/alpha,beta0*np.sqrt(max(1-(k/alpha)**2,0)),0.],complex)
    le=(1-eps)*ell
    def F(t):
        S=np.outer(d[:2],d[:2].conj())+(1-t)*np.outer(le[:2],le[:2].conj())+t*beta**2*np.eye(2)
        return np.linalg.eigvalsh(S)[-1]
    assert F(0)<beta0**2+1e-12
    assert F(1)>=beta0**2-1e-12
    lo,hi=0.,1.
    for _ in range(100):
        mid=(lo+hi)/2
        if F(mid)>beta0**2: hi=mid
        else: lo=mid
    t=(lo+hi)/2
    T=np.zeros((3,3),complex)
    T[:2,:2]=(1-t)*np.outer(le[:2],le[:2].conj())+t*beta**2*np.eye(2)
    T[2,2]=beta**2
    M=T-np.outer(le,le.conj())
    w,V=np.linalg.eigh(M)
    assert w[0]>0
    K=(V*np.sqrt(w))@V.conj().T
    B=np.column_stack((le,K))
    CH=np.column_stack((d,B))
    b=np.zeros((5,1),complex); b[0]=np.sqrt(p)
    u=-x*d-y*le
    z=I0*u/np.linalg.norm(u)
    hp=np.linalg.lstsq(B,z-a*d,rcond=None)[0]
    h=np.r_[a,hp].reshape(-1,1)
    e=np.array([-x,-y,0,0,0],complex).reshape(-1,1)
    return b,h,CH.conj().T,e,t,float(w[0])


def main():
    rng=np.random.default_rng(202610061)
    params=[]
    for _ in range(80):
        al,be=rng.uniform(.2,1.8,2)
        b0=rng.uniform(max(al,be)+1e-4,np.hypot(al,be))
        params.append((al,be,b0))
    params += [(.5,1.,1.),(1.,.5,np.hypot(1.,.5)),(1.,1.,np.sqrt(2.))]
    max_match=max_obj_err=max_kkt=0.
    min_t=1.; count=0
    for al,be,b0 in params:
        p,a,I0=rng.uniform(.5,3),rng.uniform(0,2),rng.uniform(0,2)
        gam,noise,r=1.,.05,.4
        s=dict(p=p,a=a,I0=I0,alpha=al,beta=be,beta0=b0)
        q,x,y=core.Q_approx(s,r,gam,noise)
        for eps in (.01,.001,.0001):
            b,h,C,e,t,mineig=construct(p,a,I0,al,be,b0,x,y,eps)
            got=core.spectral_summary(b,h,C,gam,noise)
            match=max(abs(got[k]-s[k]) for k in s)
            max_match=max(max_match,match)
            val=core.objective_at(b,h,C,gam,noise,e)
            expected=p*(a-x)**2-gam*(I0+np.linalg.norm(C.conj().T@e))**2-gam*noise
            max_obj_err=max(max_obj_err,abs(val-expected))
            ref=core.q_star_reference(got,r,gam,noise)
            assert ref['reference_valid']
            assert ref['q']>=q-1e-7
            assert val>=ref['q']-1e-7
            assert val-q<=gam*(2*I0+2*core.R_envelope(x,y,al,be,b0))*eps*y*be+1e-7
            max_kkt=max(max_kkt,ref['checks']['stat_res'])
            min_t=min(min_t,t);count+=1
    # Direction endpoints independent of the minimizer.
    for x,y in [(0.,.4),(.4,0.)]:
        b,h,C,e,t,_=construct(2.,1.,0.,1.,1.,1.1,x,y,.01)
        assert np.all(np.isfinite(h))
        count+=1
    boundary=0
    for _ in range(60):
        al=rng.uniform(.2,2);be=rng.uniform(.1,al);a=rng.uniform(0,2);J=rng.uniform(0,2)
        p=rng.uniform(.1,4);gam=rng.uniform(.2,2);r=rng.uniform(.1,1)
        C=np.diag([al,be,0.]).astype(complex)
        h=np.array([[a],[J/be],[0.]],complex);b=np.array([[np.sqrt(p)],[0.],[0.]],complex)
        s=core.spectral_summary(b,h,C,gam,.05)
        qb,x=core.Q_boundary_approx(s,r,gam,.05)
        y=np.sqrt(max(r*r-x*x,0))
        e=np.array([[-x if p-gam*al*al>=0 else x],[y],[0.]],complex)
        ref=core.q_star_reference(s,r,gam,.05)
        assert abs(core.objective_at(b,h,C,gam,.05,e)-qb)<1e-8
        assert abs(ref['q']-qb)<1e-8
        boundary+=1
    out=dict(non_degenerate_constructions=count,boundary_instances=boundary,
             max_six_scalar_error=max_match,max_objective_identity_error=max_obj_err,
             max_stationarity_residual=max_kkt,min_interpolation_t=min_t,
             seed=202610061,note='Numerical consistency checks; proofs are analytical in V4.1.')
    dest=Path('results_review/v41_proof_checks.json')
    dest.write_text(json.dumps(out,indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))


if __name__=='__main__': main()
