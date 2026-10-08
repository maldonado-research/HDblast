"""Fresh second-order Einstein-scalar shooting benchmark.

This is a floating-point diagnostic of a conditional local expansion. It is
not an interval existence proof, a stability calculation, or a novelty claim.
It does not import or execute the historical HDBLAST numerical producer.
"""
from pathlib import Path
import sys, json, math, hashlib, platform
sys.path.insert(0, '/private/tmp/aps-research-deps')
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import root

OUT = Path(__file__).resolve().parent
C = 2/1.0357712571566784 - 4/3
CASES = [
    dict(name='registered_HDBLAST', k=1/9, w=2., v=2., f0=1+C, f1=C, f2=0.),
    dict(name='changed_cubic_and_tension', k=1/9, w=2., v=-3., f0=1.3, f1=.4, f2=1.2),
    dict(name='negative_scalar_coupling', k=.2, w=1.3, v=.5, f0=.8, f1=-.4, f2=.6),
    dict(name='small_bulk_curvature', k=.05, w=.7, v=-1., f0=1.1, f1=.8, f2=-.5),
    dict(name='zero_scalar_coupling_control', k=1/9, w=2., v=2., f0=1.4, f1=0., f2=.3),
    dict(name='stronger_tension_curvature', k=.15, w=1., v=4., f0=1., f1=.5, f2=2.),
]
DELTAS = [.002, .001, .0005]

def functions(e, p):
    k,w,v = p['k'],p['w'],p['v']
    W=3*k+w*e*e/2+v*e**3/6
    Wp=w*e+v*e*e/2
    Wpp=w+v*e
    U=Wp*Wp/2-2*W*W/3
    Up=Wp*(Wpp-4*W/3)
    Upp=Wpp*(Wpp-4*W/3)+Wp*(v-4*Wp/3)
    f=p['f0']+p['f1']*e+p['f2']*e*e/2
    fp=p['f1']+p['f2']*e
    return W,Wp,U,Up,Upp,f,fp

def solve_case(p,delta,refined):
    k,w=p['k'],p['w']
    tol=2e-12 if refined else 3e-11
    t0=5e-6 if refined else 1e-5
    hmetric=k*p['f0']*delta/3+p['f0']**2*delta**2/36
    tguess=math.asinh(k/math.sqrt(hmetric))
    m2=w*(w-4*k)/k**2
    def linear_rhs(t,z):
        L,logf=z
        return [m2-4*L/math.tanh(t)-L*L,L]
    lin=solve_ivp(linear_rhs,(t0,tguess),[m2*t0/5,m2*t0*t0/10],
                  rtol=tol,atol=1e-14,max_step=.01,method='DOP853')
    L,logf=lin.y[:,-1]
    eb=-delta*p['f1']/(2*(k*L+w))
    sign=-1 if eb<0 else 1
    def integrate(pars, dense=False):
        logabs,tb=pars
        if tb<=t0 or tb>2*tguess+2 or logabs>0:
            raise ValueError('shooting parameters left local branch domain')
        eh=0. if p['f1']==0 else sign*math.exp(logabs)
        W,Wp,U,Up,Upp,f,fp=functions(eh,p)
        a=-U/(36*k*k)
        b=(U*U/4320-Up*Up/750)/k**4
        q2=Up/(10*k*k)
        q4=Up*(Upp/280+U/630)/k**4
        initial=[t0+a*t0**3+b*t0**5,1+3*a*t0*t0+5*b*t0**4,
                 eh+q2*t0*t0+q4*t0**4,2*q2*t0+4*q4*t0**3]
        def rhs(t,z):
            r,rp,e,ep=z
            W,Wp,U,Up,Upp,f,fp=functions(e,p)
            return [rp,-r*(ep*ep/4+U/(6*k*k)),ep,Up/(k*k)-4*rp*ep/r]
        sol=solve_ivp(rhs,(t0,tb),initial,method='DOP853',rtol=tol,
            atol=[1e-13,1e-13,1e-50,1e-50],max_step=.005 if refined else .01,
            dense_output=dense)
        if not sol.success: raise RuntimeError(sol.message)
        return sol
    def residual(pars):
        sol=integrate(pars)
        r,rp,e,ep=sol.y[:,-1]
        W,Wp,U,Up,Upp,f,fp=functions(e,p)
        return np.array([k*rp/r-(2*W+delta*f)/6,k*ep+(2*Wp+delta*fp)/2])/delta
    initial=[math.log(abs(eb))-logf if eb else -100.,tguess]
    if p['f1']==0:
        params=initial
        root_success=True
    else:
        fit=root(residual,initial,tol=1e-9,options={'eps':1e-8})
        params=fit.x
        root_success=bool(fit.success)
    sol=integrate(params,True)
    res=residual(params)*delta
    r,rp,e,ep=sol.y[:,-1]
    W,Wp,U,Up,Upp,f,fp=functions(e,p)
    h2=k*k/(r*r)
    h2_reduced=delta*(W*f/9-Wp*fp/12)+delta*delta*(f*f/36-fp*fp/48)
    prediction=hmetric-k*p['f1']**2*delta**2/(24*(w-2*k))
    sampled=sol.sol(np.linspace(t0,params[1],257))
    constraints=[]
    for rr,rrp,ee,eep in sampled.T:
        UU=functions(ee,p)[2]
        constraints.append(abs(rrp*rrp-1-rr*rr*(eep*eep/12-UU/(6*k*k)))/max(1,rrp*rrp))
    ok=bool(np.max(abs(res))<2e-11 and max(constraints)<2e-10 and
            abs(h2-h2_reduced)<2e-11 and np.isfinite(h2) and h2>0)
    return dict(delta=delta,refined=refined,solver_reported_success=root_success,
      numerical_checks_pass=ok,rtol=tol,t0=t0,t_b=float(params[1]),
      eta_b=float(e),H2=float(h2),H2_boundary_identity=float(h2_reduced),
      second_order_prediction=float(prediction),metric_only_H2=float(hmetric),
      H2_minus_metric_over_delta2=float((h2-hmetric)/delta**2),
      predicted_scalar_coefficient=float(-k*p['f1']**2/(24*(w-2*k))),
      remainder_over_delta3=float((h2-prediction)/delta**3),
      scalar_shift_over_delta=float(e/delta),
      predicted_scalar_shift=float(-p['f1']/(4*(w-2*k))),
      junction_residual=res.tolist(),max_sampled_constraint_residual=max(constraints),
      H2_identity_difference=float(h2-h2_reduced))

def main():
    rows=[]
    for p in CASES:
        case=dict(parameters=p,rows=[])
        for delta in DELTAS:
            pair=[solve_case(p,delta,False),solve_case(p,delta,True)]
            case['rows'].append(dict(delta=delta,standard=pair[0],refined=pair[1],
              H2_refinement_difference=abs(pair[0]['H2']-pair[1]['H2'])))
            print(p['name'],delta,pair[1]['numerical_checks_pass'],flush=True)
        rows.append(case)
        (OUT/'family-benchmark-results.json').write_text(json.dumps(dict(complete=False,cases=rows),indent=2)+'\n')
    result=dict(complete=True,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        cases=rows,scope='Fresh floating-point second-order Einstein-scalar shooting; no interval bounds, existence proof, stability or external novelty claim.',
        selection='Six analyst-selected model families and three positive detunings; exploratory, not publicly preregistered.',
        pass_all=all(r[t]['numerical_checks_pass'] for c in rows for r in c['rows'] for t in ['standard','refined']))
    (OUT/'family-benchmark-results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(pass_all=result['pass_all'],solves=36)))

if __name__=='__main__': main()
