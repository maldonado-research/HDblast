#!/usr/bin/env python3
"""Exact independent arithmetic gates for the local characteristic precursor."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
PINS='59faf2ccbaa6f1ac3bbe876cf28b328edc565d540591b991fb93dc8ad47d85fc'
POTENTIAL='research/HDBLAST_INTERIOR_MASTER_JETS_20261015/potential/ACTUAL_POTENTIAL_JET_BOUNDS.json'

def need(ok,reason):
    if not ok:
        raise ValueError(reason)

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def exp_upper(x,n):
    partial=sum((x**j/F(math.factorial(j)) for j in range(n+1)),F(0))
    ratio=x/F(n+2)
    need(0<=ratio<1,'TAYLOR_TAIL_RATIO')
    return partial+x**(n+1)/F(math.factorial(n+1))/(1-ratio)

def regular(p):
    return p.is_file() and all(not q.is_symlink() for q in [p,*p.absolute().parents])

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo-root',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    repo,out=args.repo_root.absolute(),args.output.absolute()
    need(all(not q.is_symlink() for p in [repo,out] for q in [p,*p.parents]),'SYMLINK_PATH')
    repo=repo.resolve(strict=True)
    need(not out.exists() and out.parent.is_dir(),'FRESH_OUTPUT')
    need(not out.is_relative_to(repo) and not out.is_relative_to(HERE),'OUTPUT_OUTSIDE_SOURCE')
    need(regular(HERE/'SOURCE_PINS.json') and sha(HERE/'SOURCE_PINS.json')==PINS,'LOCAL_PINS')
    pins=json.loads((HERE/'SOURCE_PINS.json').read_bytes())['files']
    def authenticate():
        for rel,row in pins.items():
            p=repo/rel
            need(regular(p) and p.stat().st_size==row['bytes'] and sha(p)==row['sha256'],'SOURCE_BYTES')
    authenticate()
    checks=[]
    def check(ok,name):
        need(ok,name)
        checks.append(name)
    old=json.loads((repo/POTENTIAL).read_bytes())
    ell,mu,a=F(1,10**6),F(103,10**9),F(1,2)
    check(old['status']=='PASS_ACTUAL_RETAINED_POTENTIAL_THROUGH_WZZ','retained_actual_wide_collar_potential_status')
    check(F(old['domain']['ell'])==ell and F(old['domain']['Z'])==8*ell and F(old['domain']['T'])==3*ell,'wide_static_collar_not_narrow_legacy_interval')
    bounds=list(map(F,old['potential_intervals']['W']['exact']))
    check(max(map(abs,bounds))<F(old['simple_bounds']['absolute_W'])==102000<103000,'retained_potential_bound_implies_shared_conservative_M')
    check(mu==ell*ell*103000,'actual_master_scaled_potential_normalization')
    check(1-a==F(1,2) and 2+a==F(5,2) and 2+a<8,'half_time_support_avoids_shell_and_stays_in_actual_collar')
    beta=mu*a*a/2
    B=(1+a*ell)/(1-beta)
    error=beta*B
    check(beta==mu/8 and 0<beta<1,'backward_Volterra_contraction_factor')
    face=F(1,2)-ell/4-error
    face_guard=F(1,2)-F(1,10**6)
    check(face>face_guard>0,'both_actual_face_midpoints_strictly_positive')
    check(F(1,4)*F(3,4)==F(3,16) and 4-F(16,3)==F(-4,3),'central_half_bump_lower_exponent')
    e43=exp_upper(F(4,3),14)
    check(e43<4,'rational_exp_four_thirds_less_than_four')
    check(1-2*ell>0,'positive_exponential_linear_lower_bound')
    free_apex=-ell*(1-2*ell)/16
    apex=free_apex+error
    apex_guard=F(-1,25000000)
    check(apex<apex_guard<0,'actual_common_past_apex_strictly_negative')
    check(face_guard-apex_guard>F(499999,10**6),'opposite_signs_on_both_first_half_faces')
    slope=(face_guard-apex_guard)/F(1,2)
    check(slope>F(499999,500000) and 2*slope>F(499999,250000),'nonzero_principal_tangential_derivative_on_each_face')
    # In the diamond, a=u-1,b=vnull+2 are nonnegative with a+b<=1.
    # ab<=1/4 follows exactly from (a-b)^2>=0; 4 W_uv+vW=0.
    goursat_factor=mu/16
    check(goursat_factor<1 and F(1,4)*F(1,4)==F(1,16),'two_face_Goursat_contraction_factor')
    # Independent coarser energy route, for a quarter-time Cauchy slab.
    check(1+(1+ell*ell)/4<F(9,4),'initial_energy_bound_from_endpoint_Cauchy')
    e14=exp_upper(F(1,4),8)
    check((1+mu)/8<F(1,4) and e14<F(4,3),'quarter_time_independent_energy_growth_guard')
    lower=1-ell/2-mu/2
    check(lower>1-F(1,10**6)>0,'quarter_time_both_principal_L2_profiles_nonzero')
    # Exact null-coordinate affine maps and midpoint identification.
    check((F(-1,2)+F(3,2),F(-1,2)-F(3,2))==(1,-2),'common_past_corner_null_coordinates')
    check((F(-1,4)+F(5,4),F(-1,4)-F(5,4))==(1,F(-3,2)) and (F(-1,4)+F(7,4),F(-1,4)-F(7,4))==(F(3,2),-2),'two_face_midpoint_null_coordinates')
    authenticate()
    report={'schema':'HDBLAST_PAST_PRECURSOR_INDEPENDENT_ARITHMETIC_V1','status':'PASS_CONDITIONAL_LOCAL_CHARACTERISTIC_PRECURSOR_GATES','source_pins_sha256':PINS,'finite_checks':len(checks),'check_names':checks,'constants':{'ell':str(ell),'M':'103000','scaled_mu':str(mu),'half_time':'1/2','beta':str(beta),'Bhalf':str(B),'uniform_free_error':str(error),'face_midpoint_lower':str(face),'face_midpoint_guard':str(face_guard),'common_apex_upper':str(apex),'common_apex_guard':str(apex_guard),'face_tangential_slope_lower':str(slope),'Goursat_operator_guard':str(goursat_factor),'quarter_time_L2_lower_factor':str(lower),'exp_four_thirds_rational_upper':str(e43),'exp_quarter_rational_upper':str(e14)},'analytic_conditions':['Actual master potential is exactly even in time; only that master coefficient formula extends the retained time bound.','The regular actual static collar, first master reconstruction, linear smoothness and finite propagation remain retained analytic inputs.','Both characteristic traces are restrictions of the same unique local backward Cauchy solution and share the actual nonzero corner.'],'new_activity':{'project_scientific_imports':0,'background_field_bump_or_source_queries':0,'numerical_PDE_solves':0},'limits':['These are finite arithmetic/source-authentication gates, not a new full Einstein-scalar existence proof.','Principal derivatives and weights are fixed-master diagnostics, not physical Noether energy or stress flux.','Positive conjugation preserves the face signs/zeros, not the numerical principal derivative or weight bounds.','The local two-face precursor does not select a remote incoming source or complete nonlinear/second-order data.']}
    with out.open('xb') as stream:
        stream.write((json.dumps(report,sort_keys=True,indent=2)+'\n').encode())
    print(json.dumps({'status':report['status'],'finite_checks':len(checks),'output_sha256':sha(out)},sort_keys=True))

if __name__=='__main__':
    try:
        main()
    except (ValueError,ArithmeticError,OSError,KeyError,TypeError) as exc:
        print(json.dumps({'status':'FAIL','witness':str(exc)},sort_keys=True))
        sys.exit(2)
