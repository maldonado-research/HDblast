#!/usr/bin/env python3
"""Portable exact formula audit; no actual-field enclosure or sign claim."""
import argparse
from fractions import Fraction as Q
import hashlib
import importlib.util
from itertools import product
import json
from pathlib import Path
import subprocess
import sys
import types
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ANGULAR = 'research/HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026/angular/inputs/angular_exact.py'
TIME = 'research/HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026/theory/screen_coefficient.py'
REVIEW = 'research/HDBLAST_POSTCONTACT_SCREEN_CURVATURE_20261031/independent/screen-review/check_screen.py'
REDUCTION = 'research/HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026/theory/check_reduction.py'


def need(ok, label):
    if not ok:
        raise ValueError(label)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def module_from_bytes(name, raw, path):
    m = types.ModuleType(name); m.__file__ = str(path)
    exec(compile(raw, str(path), 'exec'), m.__dict__)
    return m


def authenticate(root):
    pins = json.loads((HERE/'SOURCE_PINS.json').read_text())
    content = {}
    for rel, pin in pins.items():
        p = root/rel
        need(not p.is_symlink() and p.is_file(), 'SOURCE_TYPE:'+rel)
        need(all(not q.is_symlink() for q in p.absolute().parents), 'SOURCE_ANCESTOR_SYMLINK')
        raw = p.read_bytes()
        need(len(raw) == pin['bytes'] and digest(raw) == pin['sha256'], 'SOURCE_HASH:'+rel)
        content[rel] = raw
    return pins, content


def symbolic(api, old, review):
    import sympy as s
    keys = 'A B C D E L M At Bt Ct Dt Et Lt Mt Ay By Cy Dy Ey Ly My Mtt Myt Myy'
    m = dict(zip(keys.split(), s.symbols(keys)))
    rho, cc, theta, a, ell = s.symbols('rho cc theta a ell', nonzero=True)
    b = dict(rho=rho,c=cc,theta=theta,a=a,ell=ell)
    R=rho*cc; H=theta/rho
    checks=[]
    def eq(left,right,label):
        need(s.cancel(left-right)==0,label); checks.append(label)
    eq(api.time_quadratic_F(m,b),old.quadratic_F(m,b),'retained_time_implementation')
    eq(api.time_quadratic_F(m,b),review.expanded_F(m,b),'independent_14_term_time_expression')
    eq(sum(review.four_contributions(m,b).values()),review.expanded_F(m,b),'retained_four_time_contributions')
    mr=m['My']/2-m['C']/R; mt=m['Mt']/(2*rho)-m['E']/R
    lr=m['Ly']/2-a*m['A']/2+H*m['B']
    eta=R*(m['Ct']-theta*m['C'])/(2*rho)-R*m['Ey']/2-m['B']/2
    eta_r=-eta-m['B']
    radial_terms=[-m['A']*(mr+a*m['M'])/2+m['B']*(mt+H*m['M'])-m['C']*(m['L']-2*m['M'])/R,
                  -(m['L']/2-m['M'])*mr,m['M']*(-lr+mr),mt*eta_r]
    eq(sum(radial_terms),api.radial_quadratic_F(m,b),'four_radial_contributions')
    mt_t=m['Mtt']/(2*rho)-(m['Et']-theta*m['E'])/R
    mt_y=(m['Myt']-a*m['Mt'])/(2*rho)-(m['Ey']-a*m['E'])/R
    mr_y=m['Myy']/2-(m['Cy']-a*m['C'])/R
    g=-m['L']/(2*R**2)+H*mt-a*mr
    et=-mt_t/rho-2*H*mt-m['D']/(2*R**2)+a*mr
    ey=-mr_y-2*a*mr-m['A']/(2*R**2)
    k=eta/R**2-mt_y-a*mt-H*mr
    inv=api.first_weyl_invariants(m,b)
    eq(g+et,inv['At'],'first_Weyl_At_from_Gauss_Ricci')
    eq(ey-g,inv['Ay'],'first_Weyl_Ay_from_Gauss_Ricci')
    eq(k,inv['K'],'first_Weyl_K_from_Codazzi')
    # Rational parametrization proves the boost identity for every rapidity.
    z,xt,xy=s.symbols('z xi_t xi_y',nonzero=True)
    c=(1+z*z)/(1-z*z); sh=2*z/(1-z*z)
    p=(rho*c*xt-sh*xy)/R; q=(-rho*sh*xt+c*xy)/R
    G=c*c*inv['At']+sh*sh*inv['Ay']+2*c*sh*inv['K']
    K=c*sh*(inv['At']+inv['Ay'])+(c*c+sh*sh)*inv['K']
    eq(p*G+q*K,(rho*xt*(c*inv['At']+sh*inv['K'])+xy*(sh*inv['Ay']+c*inv['K']))/R,'position_only_tilt_reduction')
    ux,ur=s.symbols('ux ur')
    result=api.coefficients(m,b,ux,ur,c,sh,xt,xy)
    eq(result['linear_TT'],3*((c/rho)*(ux/ell)+sh*(-ur/(ell*rho)))/(2*R),'TT_physical_and_scaled_clock_factors')
    eq(result['geodesic_tilt'],-4*(p*G+q*K),'Lorentz_second_derivative_and_TT_factors')
    # Derive the first two geodesic forcing components directly from delta
    # Gamma of the 2D base metric, independently of the displayed ODE.
    g0=s.diag(-rho**2,1); gi=g0.inv()
    perturb=s.Matrix([[rho**2*m['D'],rho*m['B']],[rho*m['B'],m['A']]])
    dgi=-gi*perturb*gi
    dg0=[s.zeros(2),s.diag(-2*a*rho**2,0)]
    dq=[s.Matrix([[rho**2*m['Dt'],rho*m['Bt']],[rho*m['Bt'],m['At']]]),
        s.Matrix([[rho**2*(m['Dy']+2*a*m['D']),rho*(m['By']+a*m['B'])],
                  [rho*(m['By']+a*m['B']),m['Ay']]])]
    dGamma={(u,i,j):sum((gi[u,v]*(dq[i][j,v]+dq[j][i,v]-dq[v][i,j])
                         +dgi[u,v]*(dg0[i][j,v]+dg0[j][i,v]-dg0[v][i,j]))/2
                        for v in range(2)) for u,i,j in product(range(2),repeat=3)}
    V=[c/rho,sh]
    force=[sum(dGamma[u,i,j]*V[i]*V[j] for i,j in product(range(2),repeat=2)) for u in range(2)]
    SP=a*m['B']-(c*c*m['Dt']-sh*sh*m['At'])/(2*rho)-c*sh*m['Dy']-sh*sh*m['By']
    SQ=c*c*(m['Bt']/rho-a*(m['A']+m['D'])-m['Dy']/2)+c*sh*(m['At']/rho-2*a*m['B'])+sh*sh*m['Ay']/2
    eq(rho*force[0],SP,'first_position_time_forcing_from_delta_Christoffel')
    eq(force[1],SQ,'first_position_radial_forcing_from_delta_Christoffel')
    return checks, radial_terms, m,b


def angular_checks(r):
    I=range(4); P=r.P; v=r.grad(r.Y); h=r.hess(r.Y)
    eps={(a,c,d):r.add(*(r.scale(r.EPS.get((j,a,c,d),0),r.X[j]) for j in I)) for a,c,d in product(I,repeat=3)}
    def sym(t):
        return {(a,b):r.scale(Q(1,2),r.add(t[a,b],t[b,a])) for a,b in product(I,repeat=2)}
    axial=sym({(a,b):r.add(*(r.mul(eps[a,c,d],r.mul(v[c,],h[b,d])) for c,d in product(I,repeat=2))) for a,b in product(I,repeat=2)})
    checks=[]
    def eq(left,right,label):
        need(r.equal(left,right),label);checks.append(label)
    angular={(b,e,c,d):r.add(r.mul(P[b,c],h[e,d]),r.mul(P[e,d],h[b,c]),r.scale(-1,r.mul(P[b,d],h[e,c])),r.scale(-1,r.mul(P[e,c],h[b,d]))) for b,e,c,d in product(I,repeat=4)}
    observer=sym({(a,b):r.scale(Q(1,2),r.add(*(r.mul(eps[a,c,d],r.mul(v[e,],angular[b,e,c,d])) for c,d,e in product(I,repeat=3)))) for a,b in product(I,repeat=2)})
    eq(observer,r.tscale(-1,axial),'all_point_observer_angular_Weyl_generator')
    screen=sym({(a,b):r.scale(Q(1,2),r.add(*(r.mul(eps[a,c,d],r.add(r.scale(-1,r.mul(v[c,],h[b,d])),r.mul(v[d,],h[b,c]))) for c,d in product(I,repeat=2)))) for a,b in product(I,repeat=2)})
    eq(screen,r.tscale(-1,axial),'all_point_two_electric_or_mixed_screen_slots')
    isotropic=sym({(a,b):r.add(*(r.mul(eps[a,c,d],r.mul(v[c,],P[b,d])) for c,d in product(I,repeat=2))) for a,b in product(I,repeat=2)})
    eq(isotropic,r.zero_tensor(2),'all_point_isotropic_tilt_cancellation')
    eq(axial,r.tadd(r.tscale(2,r.TT[0]),r.tscale(-2,r.TT[1])),'all_point_fixed_TT_pair_normalization')
    need(r.inner(r.TT[0],r.TT[0])==Q(1,6) and r.inner(r.TT[1],r.TT[1])==Q(1,6) and r.inner(r.TT[0],r.TT[1])==0,'probability_sphere_norms')
    checks.append('probability_sphere_norms')
    return checks


def controls(api):
    """Execute wrong formulas on an explicit nonzero rational fixture."""
    names='A B C D E L M At Bt Ct Dt Et Lt Mt Ay By Cy Dy Ey Ly My Mtt Myt Myy'.split()
    m={k:Q((7*i+3)%19-9,i+2) for i,k in enumerate(names)}
    b=dict(rho=Q(7,5),c=Q(5,4),theta=Q(3,5),a=Q(2,7),ell=Q(1,10**6))
    c=Q(5,4); sh=Q(3,4); xt=Q(2,7); xy=Q(-3,11); ux=Q(7,17); ur=Q(11,19)
    out=api.coefficients(m,b,ux,ur,c,sh,xt,xy); R=b['rho']*b['c']
    p=(b['rho']*c*xt-sh*xy)/R; q=(-b['rho']*sh*xt+c*xy)/R
    G=c*c*out['At']+sh*sh*out['Ay']+2*c*sh*out['K']
    K=c*sh*(out['At']+out['Ay'])+(c*c+sh*sh)*out['K']
    mutants={
      'omit_radial_quadratic':out['linear_TT']+4*c*out['Ft']/R+out['geodesic_tilt'],
      'omit_geodesic_tilt':out['linear_TT']+out['boosted_quadratic'],
      'wrong_radial_TT_sign':3*(c*ux+sh*ur)/(2*b['ell']*b['rho']**2*b['c'])+out['boosted_quadratic']+out['geodesic_tilt'],
      'half_second_TT_amplitude':out['linear_TT']/2+out['boosted_quadratic']+out['geodesic_tilt'],
      'half_frame_cross_term':out['linear_TT']+out['boosted_quadratic']+out['geodesic_tilt']/2,
      'wrong_radial_screen_generator':out['linear_TT']+out['boosted_quadratic']-4*(p*G-q*K),
      'omit_mixed_Weyl_tilt':out['linear_TT']+out['boosted_quadratic']-4*(b['rho']*xt*c*out['At']+xy*sh*out['Ay'])/R,
      'reverse_screen_orientation':-out['plus'],
    }
    records=[]
    for name,got in mutants.items():
        need(got != out['plus'],'WRONG_FORMULA_NOT_REJECTED:'+name)
        records.append(dict(name=name,rejected=True,correct=str(out['plus']),wrong=str(got)))
    return records


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--repo-root',required=True,type=Path)
    ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--run-controls',action='store_true')
    ap.add_argument('--direct',action='store_true',help='also execute exact coordinate full-curvature fixture')
    args=ap.parse_args()
    need(not args.output.exists() and not args.output.is_symlink(),'FRESH_OUTPUT')
    need(args.output.resolve()!=Path(__file__).resolve() and args.output.parent.is_dir(),'OUTPUT_PATH')
    pins,raw=authenticate(args.repo_root)
    own_raw=(HERE/'proper_clock_coefficient.py').read_bytes()
    api=module_from_bytes('formula_api',own_raw,HERE/'proper_clock_coefficient.py')
    old=module_from_bytes('retained_time',raw[TIME],args.repo_root/TIME)
    review=module_from_bytes('retained_time_review',raw[REVIEW],args.repo_root/REVIEW)
    r=module_from_bytes('frozen_angular',raw[ANGULAR],args.repo_root/ANGULAR)
    checks,_,_,_=symbolic(api,old,review)
    allpoint=angular_checks(r)
    reduction=module_from_bytes('source_reduction',raw[REDUCTION],args.repo_root/REDUCTION)
    old_angular=reduction.run(args.repo_root/ANGULAR)
    report=dict(status='PASS_PROPER_CLOCK_FORMULA_AND_ANGULAR_REDUCTION',
                source_pins=pins,coefficient_api_sha256=digest(own_raw),
                symbolic_checks=checks,all_point_new_angular_checks=allpoint,
                retained_all_point_angular_checks=old_angular,
                wrong_formula_controls=controls(api) if args.run_controls else [],
                scope={'actual_proper_time_numeric_coefficient_evaluated':False,
                       'old_sign_transferred':False,'conditional_exact_family_premises_retained':True,
                       'Gaussian_congruence_margins_supplied_elsewhere':True})
    if args.direct:
        cmd=[sys.executable,'-I','-B']+(['-O'] if sys.flags.optimize else [])+[str(HERE/'radial-review/check_direct_radial.py'),'--combined']
        proc=subprocess.run(cmd,text=True,capture_output=True)
        need(proc.returncode==0,'DIRECT_FIXTURE_FAILURE:'+proc.stderr)
        report['direct_curvature_fixture']=json.loads(proc.stdout)
    authenticate(args.repo_root)
    need((HERE/'proper_clock_coefficient.py').read_bytes()==own_raw,'API_CHANGED_DURING_RUN')
    with args.output.open('x') as f:
        f.write(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(report['status'])


if __name__=='__main__':
    main()
