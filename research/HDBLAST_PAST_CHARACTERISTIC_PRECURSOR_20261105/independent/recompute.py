"""Independent rational past-master margins; no producer imports or PDE queries."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

PINS={
 'research/HDBLAST_INTERIOR_MASTER_JETS_20261015/inputs/reference/PHYSICAL_SHELL_RESPONSE.md':'a78f2239893fdd0e5ff653d4d5ce5cd9d3244c57c794ed2d3bb43bd4f5e1484b',
 'research/HDBLAST_INTERIOR_MASTER_JETS_20261015/potential/ACTUAL_POTENTIAL_WZZ_DERIVATION.md':'208dfc2107499b3c0cda6c23396173408ba27653e49e584b23e102b823d0fe0a',
 'research/HDBLAST_INTERIOR_MASTER_JETS_20261015/potential/inputs/BACKGROUND_DERIVATIVE_BOUNDS.json':'fd6fd277176033703a4e75a2f66ab94894318aac49c83264ef979e9a27f0a02c',
 'research/HDBLAST_INTERIOR_MASTER_JETS_20261015/potential/inputs/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json':'8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8',
 'research/HDBLAST_INCOMING_PREPARATION_LEDGER_20261104/sources/REGISTERED_PACKET_AUDIT.md':'ec3ad7ef010a381e898cf9dae7a94dee6d8bd35880c02b0a8e03387f3442f37f',
}

def need(ok,message):
    if not ok:raise ValueError(message)

def plus(a,b):return a[0]+b[0],a[1]+b[1]
def scale(c,a):return (c*a[0],c*a[1]) if c>=0 else (c*a[1],c*a[0])
def times(a,b):
    x=[u*v for u in a for v in b];return min(x),max(x)

def polynomial(coefficients,x):
    return sum((c*x**k for k,c in coefficients.items()),F(0))

def encode(value):
    if isinstance(value,F):return str(value)
    if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [encode(v) for v in value]
    return value

def calculate(repo):
    inputs={}
    for name,pin in PINS.items():
        path=repo/name
        need(path.is_file() and all(not q.is_symlink() for q in [path,*path.parents]),'REGULAR_PINNED_SOURCE')
        data=path.read_bytes();need(hashlib.sha256(data).hexdigest()==pin,'SOURCE_PIN:'+name);inputs[name]=data
    suffix='research/HDBLAST_INTERIOR_MASTER_JETS_20261015/potential/inputs/'
    shell=json.loads(inputs[suffix+'PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json'])
    wider=json.loads(inputs[suffix+'BACKGROUND_DERIVATIVE_BOUNDS.json'])
    ell=F(1,10**6);Z=8*ell
    need(F(wider['candidate_state_exact']['Z'])==Z and F(wider['candidate_state_exact']['ell'])==ell,'WIDE_COLLAR_LABEL')
    outer={'rho':(F(129),F(131)),'eta':(-F(9,100000),-F(8,100000)),
           'p':(-F(7,50000),-F(3,25000)),'a':(F(111,1000),F(14,125))}
    rho,eta,p,a=(outer[k] for k in ['rho','eta','p','a'])
    u1={1:F(28,9),2:F(50,9),3:-F(2,3),4:-F(20,9),5:-F(4,9)}
    u2={0:F(28,9),1:F(100,9),2:-F(2),3:-F(80,9),4:-F(20,9)}
    # U3>0 makes U2 increasing; then U2>0 makes U1 increasing.
    u3_lower=F(100,9)-F(80,3)*eta[0]**2
    need(u3_lower>0,'U2_MONOTONICITY')
    u2range=tuple(polynomial(u2,x) for x in eta)
    need(F(31,10)<u2range[0]<u2range[1]<F(25,8),'U2_ENDPOINT_RANGE')
    u1range=tuple(polynomial(u1,x) for x in eta)
    need(u1range[0]<u1range[1]<0,'U1_ENDPOINT_RANGE')
    rhs={
      'rho':scale(-1,times((rho[0]**2,rho[1]**2),a)),
      'eta':scale(-1,times(rho,p)),
      'p':scale(-1,times(rho,plus(u1range,scale(-4,times(a,p))))),
      'a':plus((1/rho[1],1/rho[0]),scale(F(1,3),times(rho,(p[1]**2,p[0]**2)))),
    }
    margins={}
    for name in outer:
        initial=tuple(F(v) for v in shell['boxes'][name+'_b']['exact'])
        tube=plus(initial,times((F(0),Z),rhs[name]))
        need(outer[name][0]<tube[0]<=tube[1]<outer[name][1],'FRESH_FIRST_EXIT:'+name)
        margins[name]={'flow_tube':tube,'lower':tube[0]-outer[name][0],'upper':outer[name][1]-tube[1]}
    need(rho[1]*Z<9,'POSITIVE_PROPER_DEPTH_MARGIN')
    hlow=u1range[1]/p[0];hhigh=u1range[0]/p[1]
    need(F(7,4)<hlow<hhigh<F(7,3),'H_ENDPOINT_RANGE')
    # Unlike the primary dependency-discarding bound, optimize the h,a terms jointly.
    need(4*F(7,4)-10*a[1]>0 and -10*F(7,4)+F(35,2)*a[1]<0,'Q_JOINT_MONOTONICITY')
    qlow=2*F(7,4)**2-10*a[1]*F(7,4)+F(35,4)*a[1]**2-F(25,8)-p[0]**2/6
    qhigh=2*F(7,3)**2-10*a[0]*F(7,3)+F(35,4)*a[0]**2-F(31,10)
    need(qlow>1 and qhigh<F(16,3),'INDEPENDENT_Q_RANGE')
    wlow=rho[0]**2-F(35,4)*ell**2
    whigh=rho[1]**2*F(16,3)
    need(wlow>16000 and whigh<92000<103000,'FULL_NEGATIVE_TIME_POTENTIAL')
    mu=F(103,10**9)
    amplitude=lambda tau:(1+ell*tau)/(1-mu*tau*tau/2)
    A=amplitude(F(1));delta=ell/2+mu*A/2
    share=(1-delta)**2/((1-delta)**2+(1+delta)**2)
    need(delta<F(1,10**6) and share>F(499999,10**6),'BOTH_CHARACTERISTIC_SHARES')
    # exp(4/3)<4 proves Q>1/4 throughout [5/4,7/4].
    x=F(4,3);term=F(1);partial=term
    for k in range(1,13):term*=x/k;partial+=term
    exponential_upper=partial+(term*x/13)/(1-x/14)
    need(exponential_upper<4,'CENTRAL_BUMP_LOWER_WITHOUT_BUMP_QUERY')
    corner_upper=-ell*(1-2*ell)/16+mu*amplitude(F(1,2))/8
    midpoint_lower=F(1,2)-ell/4-mu*amplitude(F(1,4))/32
    need(corner_upper<-F(4,10**8),'COMMON_CORNER_STRICT_NEGATIVE')
    need(midpoint_lower>F(499999,10**6),'BOTH_FACE_MIDPOINTS_STRICT_POSITIVE')
    goursat_norm=mu/16
    need(goursat_norm<1,'TRIANGULAR_GOURSAT_VOLERRA_CONTRACTION')
    for name,pin in PINS.items():need(hashlib.sha256((repo/name).read_bytes()).hexdigest()==pin,'SOURCE_PIN_AFTER')
    return encode({
      'schema':'HDBLAST_PAST_MASTER_INDEPENDENT_RATIONAL_AUDIT_V1','status':'PASS_CONDITIONAL_PAST_MASTER_MARGINS',
      'source_pins':PINS,'raw_shell_box_used_directly':True,'stored_wider_flow_intervals_used':False,
      'first_exit':margins,'static_radial_domain':[0,Z],'new_physical_time_domain':[-ell,0],
      'independent_coefficient_bounds':{'U1':u1range,'U2':u2range,'h':[hlow,hhigh],
        'Q_joint_lower':qlow,'Q_joint_upper':qhigh,'safe_Q':[1,F(16,3)],
        'strict_W':[16000,92000],'shared_W_absolute_cap':103000,'shared_scaled_cap':mu},
      'L2_route':{'Poincare_Q_over_Qprime':F(1,2),'U_norm_factor':A,'relative_J_norm_error':delta,
        'each_principal_weight_share_lower':share,'claimed_share_guard':F(499999,10**6)},
      'pointwise_route':{'exponential_4_over_3_upper':exponential_upper,
        'central_Q_lower':F(1,4),'corner_upper':corner_upper,'corner_guard':-F(4,10**8),
        'both_face_midpoint_lower':midpoint_lower,'midpoint_guard':F(499999,10**6)},
      'Goursat':{'past_corner':[F(-1,2),F(3,2)],'face_midpoints':[[F(-1,4),F(5,4)],[F(-1,4),F(7,4)]],
        'coordinates':'a=r+x-1,b=2-r+x; a,b>=0,a+b<=1',
        'equation':'U_ab=-(scaled_potential/4)U','Volterra_norm_upper':goursat_norm,
        'shared_corner_value_required':True,'each_face_has_an_interior_zero':True,
        'traces_are_restrictions_of_the_fixed_linear_Cauchy_solution':True},
      'scope':['Conditional actual linear master on the positive-depth collar before shell contact, -1<x<=0.',
        'Principal L2 weights are not physical energy/flux and are not invariant under spatial conjugation.',
        'A positive master conjugation preserves trace signs/zeros only, not all reconstructed component signs.',
        'No full Einstein-scalar characteristic constraints, nonlinear past preparation, source-selection law or uniform amplitude radius is proved.'],
      'activity':{'project_imports':0,'PDE_solves':0,'background_oracle_calls':0,'source_oracle_calls':0,'bump_queries':0},
    })

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    need(args.repo_root.is_absolute() and args.output.is_absolute(),'ABSOLUTE_ARGUMENTS')
    repo=args.repo_root.resolve(strict=True);out=args.output.absolute()
    need(all(not q.is_symlink() for q in [args.repo_root,*args.repo_root.parents,out,*out.parents]),'NO_SYMLINK_PATHS')
    need(not out.exists() and out.parent.is_dir(),'FRESH_OUTPUT')
    need(repo!=out and repo not in out.parents,'OUTPUT_OUTSIDE_REPOSITORY')
    result=calculate(repo);payload=(json.dumps(result,sort_keys=True,indent=2)+'\n').encode()
    with out.open('xb') as stream:stream.write(payload)
    print(json.dumps({'status':result['status'],'output_sha256':hashlib.sha256(payload).hexdigest()},sort_keys=True))

if __name__=='__main__':main()
