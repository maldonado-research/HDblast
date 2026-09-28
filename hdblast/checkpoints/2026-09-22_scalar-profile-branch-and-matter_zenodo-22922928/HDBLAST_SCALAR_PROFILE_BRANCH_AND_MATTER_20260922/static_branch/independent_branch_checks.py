"""Independent review: four-variable Einstein-scalar IVP, not the producer RHS.

Uses frozen shooting parameters, does not refit them, and checks both junctions.
NumPy + SciPy required. No nonlinear time evolution is performed.
"""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import json
import numpy as np
from scipy.integrate import solve_ivp

HERE=Path(__file__).resolve().parent
C=2/1.0357712571566784-4/3


def potentials(e):
    # Expanded polynomial independently derived from W=1/3+e²+e³/3.
    u=-2/27+e*e*(14/9+e*(50/27+e*(-1/6+e*(-4/9-2*e/27))))
    up=e*(28/9+e*(50/9+e*(-2/3+e*(-20/9-4*e/9))))
    upp=28/9+e*(100/9+e*(-2+e*(-80/9-20*e/9)))
    return u,up,upp


def independent_integrate(row,y_start):
    e0=row['eta_h'];end=row['y_b'];u,up,upp=potentials(e0)
    a=-u/36;b=u*u/4320-up*up/750;c=up/10;d=up*(upp/280+u/630)
    initial=[y_start+a*y_start**3+b*y_start**5,
             1+3*a*y_start**2+5*b*y_start**4,
             e0+c*y_start**2+d*y_start**4,2*c*y_start+4*d*y_start**3]
    def rhs(y,v):
        rho,ry,e,p=v;u,up,_=potentials(e)
        # Independent second-order Einstein equation. No sqrt first integral.
        return [ry,-rho*(p*p/4+u/6),p,up-4*ry*p/rho]
    sol=solve_ivp(rhs,(y_start,end),initial,method='DOP853',rtol=3e-13,
                  atol=[1e-14,1e-14,1e-43,1e-43],max_step=.06,dense_output=True)
    assert sol.success
    r,ry,e,p=sol.y[:,-1];u,_,_=potentials(e);delta=row['delta']
    sig=2/3+2*e*e+2*e**3/3+delta*(1+C+C*e)
    sig1=4*e+2*e*e+delta*C
    yy=np.linspace(y_start,end,8001);rr,ryy,ee,pp=sol.sol(yy)
    uu=potentials(ee)[0]
    constraint=ryy*ryy-1-rr*rr*(pp*pp/12-uu/6)
    denom=1+ryy*ryy+rr*rr*(pp*pp/12+np.abs(uu)/6)
    return sol,{'y_start':y_start,'metric_junction':float(ry/r-sig/6),
        'scalar_junction':float(p+sig1/2),
        'rho_b_relative_difference':float((r-row['rho_b'])/row['rho_b']),
        'eta_b_relative_difference':float((e-row['eta_b'])/row['eta_b']),
        'H2_relative_difference':float(((1/r**2)-row['H2'])/row['H2']),
        'H2_identity_residual':float(1/r**2-(sig*sig/36-sig1*sig1/48+u/6)),
        'constraint_max_absolute':float(np.max(np.abs(constraint))),
        'constraint_max_normalized':float(np.max(np.abs(constraint)/denom)),
        'min_rho':float(rr.min()),'min_rho_y':float(ryy.min())}


def main():
    raw=json.loads((HERE/'PLUS_BRANCH_RESULTS.json').read_text())
    profiles=np.load(HERE/'PLUS_BRANCH_PROFILES.npz')
    checks={};rows=[]
    checks['producer_hash_matches']=hashlib.sha256((HERE/'solve_plus_branch.py').read_bytes()).hexdigest()==raw['script_sha256']
    # Exact rational checks: the coefficient of c² in the second-order term.
    a=F(-9,64)
    second_coefficient=(2*a*a+a)/27+F(7,27)*a*a-(4*a+1)**2/48
    checks['exact_H2_second_order_c_squared_coefficient']=second_coefficient==F(-1,384)
    checks['exact_regular_AdS_exponent']=F(14,9)**2+F(4,9)*F(14,9)-F(28,9)==0
    checks['exact_leading_scalar_junction']=F(14,9)*a+2*a+F(1,2)==0
    for item in raw['rows']:
        row=item['refined'];d=row['delta'];results=[]
        ys,rhos,etas,ps=profiles[str(d)]
        for start in [2e-4,5e-5,1e-5]:
            sol,r=independent_integrate(row,start)
            mask=ys>=start
            iv=sol.sol(ys[mask])
            r['sampled_profile_rho_max_relative_difference']=float(np.max(np.abs(iv[0]-rhos[mask])/rhos[mask]))
            # Relative to brane amplitude so tiny cone amplitudes don't dominate.
            r['sampled_profile_eta_max_brane_normalized_difference']=float(np.max(np.abs(iv[2]-etas[mask]))/abs(row['eta_b']))
            r['sampled_profile_py_max_brane_normalized_difference']=float(np.max(np.abs(iv[3]-ps[mask]))/abs(row['phi_y_b']))
            results.append(r)
        checks[f'delta_{d}_both_junctions']=max(abs(r[k]) for r in results for k in ['metric_junction','scalar_junction'])<1e-10
        checks[f'delta_{d}_independent_constraint']=max(r['constraint_max_normalized'] for r in results)<1e-11
        checks[f'delta_{d}_profile_agreement']=max(r[k] for r in results for k in ['sampled_profile_rho_max_relative_difference','sampled_profile_eta_max_brane_normalized_difference','sampled_profile_py_max_brane_normalized_difference'])<1e-9
        checks[f'delta_{d}_cone_start_H2_stability']=max(abs(r['H2_relative_difference']) for r in results)<1e-10
        rows.append({'delta':d,'runs':results,
                     'eta_O_delta2_remainder_coefficient':(row['eta_b']-row['eta_first_order'])/d**2,
                     'H2_O_delta3_remainder_coefficient':(row['H2']-row['H2_second_order_expansion'])/d**3})
    out={'passed':all(checks.values()),'checks':checks,'n_checks':len(checks),
         'independent_integrations':sum(len(r['runs']) for r in rows),
         'method':'Four-variable second-order Einstein-scalar IVP using archived shooting parameters; no import of producer or use of its first-integral RHS.',
         'scope':'Floating consistency and regular-cone-start sensitivity, not interval existence or dynamic stability.',
         'results_sha256':hashlib.sha256((HERE/'PLUS_BRANCH_RESULTS.json').read_bytes()).hexdigest(),
         'profiles_sha256':hashlib.sha256((HERE/'PLUS_BRANCH_PROFILES.npz').read_bytes()).hexdigest(),
         'rows':rows}
    (HERE/'INDEPENDENT_BRANCH_CHECKS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'passed':out['passed'],'n_checks':out['n_checks'], 'independent_integrations':out['independent_integrations'],
                      'failures':[k for k,v in checks.items() if not v]},indent=2))
    assert out['passed']


if __name__=='__main__':
    main()
