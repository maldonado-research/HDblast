#!/usr/bin/env python3
"""Independent exact-rational join of the two new actual shell campaigns."""
import argparse,hashlib,json,sys
from fractions import Fraction as F
from pathlib import Path
HERE=Path(__file__).resolve().parent
PINS_SHA='b4a008faa3b8f65f481adbbcbb2c11c52ca375cb09a1a2ce3f345dd736677bc4'
class Failure(Exception):pass
def need(b,witness):
    if not b:raise Failure(witness)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def canonical(x):return (json.dumps(x,sort_keys=True,indent=2)+'\n').encode()
def regular(p):return p.is_file() and not p.is_symlink() and all(not x.is_symlink() for x in p.absolute().parents)
def pair(x):
    need(isinstance(x,list) and len(x)==2,'INTERVAL_PAIR');lo,hi=map(F,x);need(lo<=hi,'INTERVAL_ORDER');return lo,hi
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--gradient-report',type=Path,required=True);ap.add_argument('--kinematics-report',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
    pinpath=HERE/'HEADLINE_PINS.json';need(regular(pinpath) and sha(pinpath)==PINS_SHA,'HEADLINE_PINS');pins=json.loads(pinpath.read_bytes());own=sha(Path(__file__))
    output=a.output.absolute();need(not output.exists() and all(not p.is_symlink() for p in [output,*output.parents]),'FRESH_OUTPUT')
    repo=a.repo_root.resolve(strict=True)
    need(all(base!=output and base not in output.parents for base in [HERE.parent,repo]),'OUTPUT_OUTSIDE_CODE_AND_SOURCE')
    paths={'gradient':a.gradient_report,'kinematics':a.kinematics_report}
    for name,path in paths.items():need(regular(path) and sha(path)==pins['reports'][name]['sha256'],'INPUT_REPORT_'+name.upper())
    gradient=json.loads(a.gradient_report.read_bytes());kinematics=json.loads(a.kinematics_report.read_bytes())
    need(gradient['frozen_campaign']['continuous_postcontact_interval']==['31/25','63/50']
         and kinematics['frozen_campaign']['continuous_target_interval']==['31/25','63/50'],'SAME_CONTINUOUS_INTERVAL')
    need(gradient['actual_oracle_calls']==9 and kinematics['actual_oracle_calls']==73,'FIXED_ACTUAL_CAMPAIGN_COUNTS')
    c=pair(kinematics['continuous_interval_hulls']['curvature_defect_second'])
    qb=pair(kinematics['continuous_interval_hulls']['buchert_Q_second'])
    q=pair(gradient['continuous_interval_hulls']['mean_Q2'])
    ratio=pair(gradient['continuous_interval_hulls']['spatial_to_temporal_ratio'])
    need(F(3,10**24)<=c[0]<=c[1]<=F(2,10**22),'CONTINUOUS_CURVATURE_DEFECT_HEADLINE')
    need(F(1,10**9)<=qb[0]<=qb[1]<=F(2,10**8),'CONTINUOUS_BUCHERT_HEADLINE')
    need(q[1]<-4*10**17,'CONTINUOUS_NEGATIVE_GRADIENT_HEADLINE')
    need(0<ratio[0]<=ratio[1]<F(1,10**12),'CONTINUOUS_TEMPORAL_DOMINATION')
    R=ratio[1];r=R/(R+2);rcoarse=F(1,2*10**12);sqrt_compare=F(3,4*10**6)
    need(0<r<R/2<rcoarse<sqrt_compare**2,'ANGULAR_RATIONAL_SQUARE_CHAIN')
    fraction=4*sqrt_compare/3;need(fraction==F(1,10**6),'ANGULAR_FRACTION_HEADLINE')
    result={'schema':'HDBLAST_INTRINSIC_SHELL_BUDGET_HEADLINE_JOIN_V1','status':'PASS_EXACT_CONTINUOUS_SHELL_BUDGET_HEADLINES',
        'checker_sha256':own,'pins_sha256':PINS_SHA,'input_report_sha256':{name:row['sha256'] for name,row in pins['reports'].items()},
        'continuous_scaled_x':['31/25','63/50'],'ell':'1/1000000','all_comparisons_exact_Fraction':True,'new_actual_oracle_calls':0,
        'raw_continuous_hulls':{'curvature_defect_second':list(map(str,c)),'buchert_Q_second':list(map(str,qb)),'mean_Q2':list(map(str,q)),'spatial_to_temporal_ratio':list(map(str,ratio))},
        'verified_loose_bounds':{'curvature_defect_second':['3/1000000000000000000000000','1/5000000000000000000000'],
                                 'buchert_Q_second':['1/1000000000','1/50000000'],'mean_Q2_strict_upper':'-400000000000000000','spatial_to_temporal_ratio_strict_upper':'1/1000000000000'},
        'angular_fraction':{'R_upper':str(R),'r_upper':str(r),'r_strict_coarse_upper':str(rcoarse),
                            'sqrt_r_strict_upper_comparison_square':str(sqrt_compare**2),'pi_lower_premise':'3',
                            'non_timelike_coefficient_fraction_strict_upper':str(fraction),'timelike_coefficient_fraction_strict_lower':str(1-fraction),
                            'uses_independent_geometry_angular_measure_proof':True},
        'scope':{'conditional_exact_family_and_actual_first_field_premises':True,'uniform_continuous_background_time_interval':True,
                 'coefficient_angular_fraction_is_probability_normalized_round_S3':True,
                 'proper_volume_finite_parameter_consequence':'For each fixed admissible C2 completion, uniform convergence gives the stated loose angular fraction for sufficiently small nonzero parameter.',
                 'effective_parameter_radius':False,'completion_uniform_parameter_radius':False,'full_scalar_monopole_or_acceleration_response':False,
                 'SI_units_energy_force_hot_matter_or_cosmological_origin_claim':False}}
    need(sha(Path(__file__))==own and sha(pinpath)==PINS_SHA and all(sha(paths[name])==row['sha256'] for name,row in pins['reports'].items()),'HEADLINE_POSTCHECK')
    with output.open('xb') as handle:handle.write(canonical(result))
    print(json.dumps({'status':'PASS','witness':'EXACT_CONTINUOUS_HEADLINE_JOIN','report_sha256':sha(output)},sort_keys=True))
if __name__=='__main__':
    try:main()
    except Failure as exc:
        print(json.dumps({'status':'FAIL','witness':str(exc)},sort_keys=True));sys.exit(2)
