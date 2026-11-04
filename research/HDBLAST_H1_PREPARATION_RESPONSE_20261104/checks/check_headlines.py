#!/usr/bin/env python3
"""Join the conditional H1 energy estimates, complete observables and retained witness."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

CAP=Path(__file__).resolve().parents[1]
BASE='research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/response/references/FINAL_NORMAL.json'
BASE_PIN='1a9ebffb6bee9332bc0fd0df0f5d1d9f7f2a796825c2a563e4c644974ba4d8a3'
PREMISES='research/HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/EXPLICIT_PREMISES.json'
PREMISES_PIN='436305e19f375ddead2958dd44ea48dd080302852ef056631efd3ae2e6b433c1'
WITNESS='research/HDBLAST_INTRINSIC_PREPARATION_ROBUSTNESS_20261104/independent/references/CONTROLS.artifacts/WITNESS_NORMAL.json'
WITNESS_PIN='ba3973210fb364f40cafaada84d3dd04abb29b385c594d70603903817ad31e5c'
IMAGE_PIN='af2e1c9e54a04e603a6b721a71c773200bbbf67f6dd0f32fedd924613672dd75'

class Rejected(Exception):pass
def need(ok,witness):
 if not ok:raise Rejected(witness)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def regular(p):return p.is_file() and all(not q.is_symlink() for q in [p,*p.absolute().parents])
def read(p):
 need(regular(p),'REGULAR_INPUT');return json.loads(p.read_bytes())
def box(value):
 lo,hi=map(F,value['exact'] if isinstance(value,dict) else value);need(lo<=hi,'ORDERED_INTERVAL');return lo,hi
def record(value):return [str(v) for v in value] if isinstance(value,tuple) else str(value)

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 for name in ['repo-root','energy-report','independent-report','image-report','output']:ap.add_argument('--'+name,required=True,type=Path)
 a=ap.parse_args();repo=a.repo_root.absolute();out=a.output.absolute()
 need(all(not p.is_symlink() for p in [repo,*repo.parents,out,*out.parents]),'PATH_SYMLINK');repo=repo.resolve(strict=True);out=out.resolve()
 need(not out.exists(),'FRESH_OUTPUT');need(out.parent.is_dir(),'OUTPUT_PARENT');need(all(out!=p and p not in out.parents for p in [repo,CAP]),'OUTPUT_OUTSIDE_CODE_AND_SOURCE')
 for rel,pin in [(BASE,BASE_PIN),(PREMISES,PREMISES_PIN),(WITNESS,WITNESS_PIN)]:need(regular(repo/rel) and sha(repo/rel)==pin,'RETAINED_SOURCE_PIN:'+rel)
 base=read(repo/BASE);premises=read(repo/PREMISES);witness=read(repo/WITNESS);p=read(a.energy_report);i=read(a.independent_report);image=read(a.image_report)
 need(p['status']=='PASS_CONDITIONAL_H1_PREPARATION_RESPONSE_BOUNDS' and i['status']=='PASS_INDEPENDENT_H1_PREPARATION_RESPONSE','PRIMARY_AND_INDEPENDENT_STATUS')
 need(sha(a.image_report)==IMAGE_PIN and image['status']=='PASS_EXACT_FORMAL_ALGEBRA_AND_CONTINUOUS_BOUND_CONSTANTS','FORMAL_ACTUAL_IMAGE_REFERENCE')
 need(premises['status']=='EXPLICIT_CONDITIONAL_PREMISES_NOT_REPROVED_BY_PORTABLE_CHECKER' and p['continuous_operator_caps_are_retained_conditional_assertions'] is True,'RETAINED_CONDITIONAL_OPERATOR_ASSERTIONS')
 s=i['shared_projection'];need(s['conditional_operator_assertions'] is True and i['primary_comparison']['shared_projection_matches'] is True and i['primary_comparison']['unknown_top_level_fields_rejected'] is True,'INDEPENDENT_SCOPE_AND_COMPARISON')
 for key in ['mean_Gram','new_domain','retained_domain','observable_difference_caps_per_R','whole_class_absolute_difference_caps','whole_class_exact_intervals','whole_class_outward_intervals']:need(p[key]==s[key],'SHARED_MATHEMATICAL_JOIN:'+key)
 need(p['class']['norm']==s['class_norm'] and p['class']['parameter_normalization_fixed'] is True,'COMPLETE_INCREMENT_NORM')
 radius=F(p['class']['radius_closed']);need(radius==F(s['class_radius'])==2*10**31,'PREDECLARED_CLOSED_H1_RADIUS')
 shared={'tilt_k':F(1,2),'boundary_alpha':F(1,2),'growth_guard':F(201,100),'exponential_guard':F(13),'initial_energy_guard':F(1,5000000),'boundary_scalar_factor':F(1,300),'bulk_growth_exact_upper':F(1000351,500000),'boundary_growth_exact_upper':F(2003001,1000000),'exponential_argument':F(201,80),'boundary_value_square_guard':F(13,1250000)}
 for key,value in shared.items():need(F(p['energy'][key])==F(s['energy'][key])==value,'SHARED_ENERGY_GUARD:'+key)
 need(p['mean_Gram']==[1,2,1,3,1],'COMPLETE_FIVE_MEAN_GRAM')
 metric_factor=F(3,2);scalar_pullback_factor=F(-3)
 need(2*metric_factor+scalar_pullback_factor==0,'ACTUAL_SECOND_DERIVATIVE_GEOMETRY_CANCELLATION')
 c={k:F(v) for k,v in p['initial_L2_caps_per_R'].items()};d={k:F(v) for k,v in i['independent_initial_L2_caps_per_R'].items()}
 need(c['Gram_U0']==2*c['h']+c['w'] and c['Gram_U0_r']==2*c['h_r']+c['w_r'] and c['Gram_U0_x']==2*c['V_L'],'PRIMARY_COMPLETE_INITIAL_GRAM_CHANNELS')
 need(c['initial_boundary_value']==c['initial_tilt_cross_term']==0,'INITIAL_BOUNDARY_AND_ORTHOGONAL_TILT')
 E0=F(p['energy']['initial_energy_per_R_squared']);Ei=F(i['independent_initial_energy_per_R_squared'])
 need(E0==(c['Gram_U0']**2+c['Gram_U0_r']**2+c['Gram_U0_x']**2)/2,'PRIMARY_INITIAL_ENERGY_RECONSTRUCTION')
 need(Ei==(4*d['h']**2+d['w']**2+4*d['h_r']**2+d['w_r']**2+2*d['V_L']**2)/2,'INDEPENDENT_INITIAL_ENERGY_RECONSTRUCTION')
 need(0<E0<shared['initial_energy_guard'] and 0<Ei<shared['initial_energy_guard'],'DISTINCT_VALID_INITIAL_ENERGY_MAJORANTS')
 need(F(p['energy']['exponential_rational_upper'])<13 and F(i['independent_exponential']['rational_upper'])<13,'BOTH_RATIONAL_EXPONENTIAL_GUARDS')
 need(max(shared['bulk_growth_exact_upper'],shared['boundary_growth_exact_upper'])<shared['growth_guard'],'BULK_AND_NONSYMMETRIC_ROBIN_GROWTH')
 need(2/shared['boundary_alpha']*shared['exponential_guard']*shared['initial_energy_guard']==shared['boundary_value_square_guard']<shared['boundary_scalar_factor']**2,'BOUNDARY_SCALAR_TRACE_GUARD')
 ell=F(1,10**6);delta=F(1,4);rhoM=F(131);force=F(1,1000)
 Cv=(rhoM*ell)**2*force*delta**2/(4*300);Ct=rhoM*ell*force*delta/(2*300)
 need((Cv,Ct)==(F(p['observable_difference_caps_per_R']['v2']),F(p['observable_difference_caps_per_R']['Theta2'])),'PROPER_CLOCK_GREEN_TRANSFER_FACTORS')
 bounds={};displays={'v2':('-7.2400e16','-8.6400e15'),'Theta2':('-1.0113e22','-3.3540e21')}
 for key,oldkey,C in [('v2','proper_v2',Cv),('Theta2','proper_Theta2',Ct)]:
  old=box(base[oldkey]);new=(old[0]-C*radius,old[1]+C*radius);need(new==box(p['whole_class_exact_intervals'][key]),'EXACT_ABSOLUTE_CLASS_INTERVAL:'+key)
  need(F(p['whole_class_absolute_difference_caps'][key])==C*radius,'WHOLE_CLASS_DIFFERENCE_CAP:'+key)
  lo,hi=map(F,displays[key]);need((lo,hi)==box(p['whole_class_outward_intervals'][key]) and lo<=new[0]<=new[1]<=hi<0,'STRICT_OUTWARD_PRINTED_INTERVAL:'+key)
  bounds[key]={'retained_absolute_base':record(old),'whole_closed_H1_class':record(new),'printed_outward':list(displays[key]),'difference_cap':str(C*radius)}
 trace=F(21,10);need(1/delta+delta==F(17,4)<trace**2,'ENDPOINT_H1_TRACE_FACTOR')
 need(witness['status']=='PASS_INDEPENDENT_WITNESS_EXCLUSION' and F(witness['lambda'])==-10**21 and F(witness['target'])==F(5,4),'RETAINED_SCALED_WITNESS_IDENTITY')
 oldD=F(witness['witness_absolute_D_lower']);lower=oldD/trace
 need(lower==F(i['old_witness_H1_radius_lower']) and oldD>F('3.7621e33') and lower>F('3.7621e33')/trace>F('1.7914e33')>radius,'OLD_WITNESS_OUTSIDE_CLOSED_H1_CLASS')
 need(p['verification']['finite_checks']==len(p['verification']['check_names'])==24 and i['verification']['finite_checks']==len(i['verification']['check_names'])==25,'OBSERVED_FINITE_GATE_ROSTERS')
 need(p['verification']['new_PDE_solves']==p['verification']['new_background_or_evolved_field_queries']==p['verification']['new_bump_queries']==0 and all(type(v)is int and v==0 for v in i['activity'].values()) and all(type(v)is int and v==0 for v in image['new_activity'].values()),'ZERO_NEW_SCIENTIFIC_QUERIES')
 need(image['exact_bump_jets_at_zero']==['0','0','0','6'] and image['m_subsequence_L2_exponents_j0_to_3']==['-6','-4','-2','0'],'ACTUAL_IMAGE_CONCENTRATION_IDENTITIES')
 result={'schema':'HDBLAST_H1_PREPARATION_RESPONSE_HEADLINES_V1','status':'PASS_EXACT_CONDITIONAL_H1_PREPARATION_RESPONSE_HEADLINES',
  'conditional_operator_norm_assertions':{'source':PREMISES,'source_sha256':PREMISES_PIN,'status':premises['status'],'new_matrix_operator_bound_proof':False},
  'geometry':{'coefficient_convention':'actual second parameter derivative; no extra Taylor half','profile_formula':'D=3*J/2-3*a*w/p','metric_factor':str(metric_factor),'scalar_pullback_factor':str(scalar_pullback_factor),'background_scalar_label_derivative':'D_r=(-ell*rho*p)*D_v','specified_physical_Cauchy_slice_required':True,'selected_scalar_level_component_required':True,'class_restricted_to_retained_admissible_image':True,'profile_map_asserted_onto_arbitrary_H1_data':False,'causal_segment_determines_exterior_preparations':False},
  'complete_five_mean_initial_energy':{'order':['K','L','N','J','w'],'Gram':[1,2,1,3,1],'initial_slots':['h','0','0','h','w'],'velocity_slots':['0','V_L','0','0','0'],'radial_slots':['h_r','0','0','h_r','w_r'],'tilted_initial_cross_term':'0','primary_majorant_per_R_squared':str(E0),'independent_majorant_per_R_squared':str(Ei),'majorants_required_equal':False,'shared_guard':str(shared['initial_energy_guard'])},
  'shared_energy_guards':{k:str(v) for k,v in shared.items()},
  'proper_clock_transfer_per_H1_radius':{'volume':str(Cv),'expansion':str(Ct),'volume_denominator_factor':'4*300','expansion_denominator_factor':'2*300','scalar_boundary_factor':'1/300','Green_kernel_upper':'T-x','Green_kernel_target_derivative_upper':'1'},
  'whole_closed_class':{'norm':p['class']['norm'],'radius':str(radius),'initial_scalar_label_interval':['1','5/4'],'bounds':bounds},
  'retained_scaled_witness':{'source':WITNESS,'source_sha256':WITNESS_PIN,'lambda':str(F(witness['lambda'])),'absolute_D_lower':str(oldD),'H1_norm_lower':str(lower),'H1_norm_lower_relation':'>=','strict_display_lower':'1.7914e33','strictly_outside_new_H1_class':True,'sample_is_initial_scalar_label_not_future_time':True,'old_analytic_construction_replayed':False},
  'actual_image_support':{'formal_certificate_sha256':IMAGE_PIN,'bounded_H1_does_not_bound_old_C1_or_initial_sup_norm_route':True,'response_instability_claimed':False,'set_containment_of_fixed_radius_balls_claimed':False},
  'finite_gates':{'primary':24,'independent':25,'are_subprocess_negative_controls':False},
  'scope':{'conditional_coefficient_class_result':True,'tilted_energy_identified_as_positive_physical_energy':False,'physical_class_selection_established':False,'uniform_effective_nonlinear_epsilon_window_established':False,'new_weak_solution_existence_established':False,'absolute_shared_first_quadratics_discarded':False,'first_only_bulk_shell_constraint_clock_measure_Gauss_terms_cancel_only_in_differences':True},
  'new_activity':{'oracle_calls':0,'background_queries':0,'PDE_solves':0,'bump_queries':0,'analytic_target_bump_replays':0}}
 need(regular(CAP/'DERIVATION.md'),'ROOT_DERIVATION_REGULAR');result['geometry']['root_derivation_sha256']=sha(CAP/'DERIVATION.md')
 payload=(json.dumps(result,sort_keys=True,indent=2)+'\n').encode()
 with out.open('xb') as stream:stream.write(payload)
 print(json.dumps({'status':'PASS','witness':'EXACT_CONDITIONAL_H1_PREPARATION_RESPONSE_HEADLINES','output_sha256':hashlib.sha256(payload).hexdigest()},sort_keys=True))

if __name__=='__main__':
 try:main()
 except (Rejected,OSError,ValueError,KeyError,TypeError,ZeroDivisionError) as error:
  print(json.dumps({'status':'FAIL','witness':str(error) if isinstance(error,Rejected) else 'HEADLINE_INPUT_FAILURE'},sort_keys=True));sys.exit(2)
