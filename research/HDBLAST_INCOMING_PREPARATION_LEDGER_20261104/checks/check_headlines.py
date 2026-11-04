#!/usr/bin/env python3
"""Exact joins for the future source ledger, partial null constraint and image-only inverse."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

CAP=Path(__file__).resolve().parents[1]
PREMISES='research/HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/EXPLICIT_PREMISES.json'
PREMISES_PIN='436305e19f375ddead2958dd44ea48dd080302852ef056631efd3ae2e6b433c1'
RAW='research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/ACTUAL_SOURCES_NORMAL.json'
RAW_PIN='45a5550d2e4647bdfb84f9c4a36bf9861ef500a1de89a76f978c5afe9c51a60f'
class Rejected(Exception):pass
def need(ok,witness):
 if not ok:raise Rejected(witness)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def regular(p):return p.is_file() and all(not q.is_symlink() for q in [p,*p.absolute().parents])
def read(p):
 need(regular(p),'REGULAR_INPUT');return json.loads(p.read_bytes())
def fractions(values):return {k:F(v) for k,v in values.items()}
def encode(value):
 if isinstance(value,F):return str(value)
 if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
 if isinstance(value,(tuple,list)):return [encode(v) for v in value]
 return value

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 for name in ['repo-root','source-report','theory-report','independent-report','output']:ap.add_argument('--'+name,required=True,type=Path)
 a=ap.parse_args();repo=a.repo_root.absolute();out=a.output.absolute()
 need(all(not p.is_symlink() for p in [repo,*repo.parents,out,*out.parents]),'PATH_SYMLINK');repo=repo.resolve(strict=True);out=out.resolve()
 need(not out.exists(),'FRESH_OUTPUT');need(out.parent.is_dir(),'OUTPUT_PARENT');need(all(out!=p and p not in out.parents for p in [repo,CAP]),'OUTPUT_OUTSIDE_CODE_AND_SOURCE')
 for rel,pin in [(PREMISES,PREMISES_PIN),(RAW,RAW_PIN)]:need(regular(repo/rel) and sha(repo/rel)==pin,'RETAINED_SOURCE_PIN:'+rel)
 premises=read(repo/PREMISES);raw=read(repo/RAW);s=read(a.source_report);t=read(a.theory_report);i=read(a.independent_report)
 need(s['status']=='PASS_EXACT_ZERO_ORACLE_SOURCE_LEDGER' and t['status']=='PASS_CONDITIONAL_INCOMING_PREPARATION_IDENTITIES' and i['status']=='PASS_INDEPENDENT_INCOMING_PREPARATION_RECONSTRUCTION','THREE_COMPONENT_STATUSES')
 need(i['candidate_sha256']=={'source':sha(a.source_report),'theory':sha(a.theory_report)},'FRESH_PRIMARY_CANDIDATE_JOIN')
 need(premises['status']=='EXPLICIT_CONDITIONAL_PREMISES_NOT_REPROVED_BY_PORTABLE_CHECKER' and t['operator_caps_remain_conditional_future_assertions'] is True and s['analytic_identity_premises_reproved'] is False,'RETAINED_CONDITIONAL_PREMISES')
 gram=[1,2,1,3,1];need(premises['mean_gram']==s['mean_gram']==gram and s['mean_order']==['K0','L0','N0','J0','w0'],'COMPLETE_FIVE_MEAN_GRAM')
 ell,T,delta,k,alpha,m=F(1,10**6),F(5,4),F(1,4),F(1,2),F(1,2),F(1,2000)
 C,Cf,Cg=F(301,100),F(1,2),1+(F(1,2)+m)**2;exp_guard=F(44)
 expected={'tilt_k':k,'boundary_alpha':alpha,'homogeneous_bulk_rate':F(1000351,500000),'homogeneous_boundary_rate':F(2003001,1000000),'forced_growth_guard':C,'bulk_square_factor':Cf,'boundary_square_factor':Cg,'exponential_guard':exp_guard,'scalar_trace_square_factor':F(4)}
 need(fractions(i['shared_projection']['forced_energy'])==expected,'INDEPENDENT_FORCED_GUARDS')
 need(all(F(t['forced_energy'][key])==value for key,value in expected.items()),'PRIMARY_FORCED_GUARDS')
 need(max(expected['homogeneous_bulk_rate']+1,expected['homogeneous_boundary_rate']+1)<C and Cg==F(5002001,4000000),'NONSYMMETRIC_FORCING_GROWTH_AND_BOUNDARY_SQUARE')
 need(F(t['forced_energy']['exponential_argument'])==C*T==F(301,80) and F(t['forced_energy']['exponential_rational_upper'])<44 and F(i['independent_route']['order_24_exponential_upper'])<44,'DISTINCT_RATIONAL_EXPONENTIAL_GUARDS')
 need(F(t['forced_energy']['scalar_trace_square_factor'])==2/alpha==4,'BOUNDARY_SCALAR_FROM_ENERGY')
 need(fractions(s['energy_coefficients_conditional_on_accompanying_forced_proof'])=={'C':C,'Cf':Cf,'Cg':Cg,'exp_CT_upper':exp_guard},'SOURCE_LEDGER_FORCED_PROOF_JOIN')
 geometry={'T':T,'cone_area':F(25,32),'active_area':F(9,32),'boundary_length':T,'active_boundary_length':delta}
 need(fractions(s['exact_geometry'])==geometry,'EXACT_FUTURE_CONE_AND_ZERO_WEDGE')
 FH,GH=F(s['registered_caps']['F_H']),F(s['registered_caps']['G_H']);need((FH,GH)==(F(17387513156573),F(3851400)),'REGISTERED_SOURCE_CAPS')
 need(s['retained_actual_source_report_sha256']==RAW_PIN and raw['global_cover']['cells']==len(raw['global_cover']['cell_enclosures'])==81,'RETAINED_COMPLETE_81_CELL_REPORT')
 bulk=s['retained_bulk_cells_intersecting_cone'];shell=s['retained_shell_cells_intersecting_boundary'];need(len(bulk)==6,'SIX_CONE_INTERSECTIONS')
 active={tuple(row['old_cell_index']):F(row['active_area']) for row in bulk if F(row['active_area'])}
 need(active=={(0,1):F(3,32),(1,0):F(3,32),(0,2):F(1,32),(1,1):F(1,32),(2,0):F(1,32)},'FIVE_EXACT_ACTIVE_CLIPPED_AREAS')
 need(sum((F(row['cone_area']) for row in bulk),F(0))==geometry['cone_area'] and sum(active.values(),F(0))==geometry['active_area'],'BULK_AREA_SUMS')
 need(sum((F(row['length']) for row in shell),F(0))==T and sum((F(row['active_length']) for row in shell),F(0))==delta,'BOUNDARY_LENGTH_SUMS')
 for row in bulk:need(F(row['chosen_square_cap'])==min(F(row['direct_mean_square_cap']),F(row['full_angular_H_square_cap']),FH**2),'REPORTED_SOURCE_COMPONENT_MINIMUM')
 for row in shell:need(F(row['chosen_square_cap'])==min(F(row['direct_mean_square_cap']),GH**2),'REPORTED_BOUNDARY_COMPONENT_MINIMUM')
 integral_bulk=sum((F(row['active_area'])*F(row['chosen_square_cap']) for row in bulk),F(0));integral_boundary=sum((F(row['active_length'])*F(row['chosen_square_cap']) for row in shell),F(0))
 budgets={name:fractions(values) for name,values in s['budgets'].items()};indbudgets={name:fractions(values) for name,values in i['shared_projection']['source_budgets'].items()};need(budgets==indbudgets,'INDEPENDENT_EXACT_SOURCE_LEDGER_JOIN')
 active_budget=budgets['wedge_zero_retained_cell_caps'];Bsrc=Cf*integral_bulk+Cg*integral_boundary
 need(active_budget=={'bulk_L2_squared_upper':integral_bulk,'boundary_L2_squared_upper':integral_boundary,'energy_source_ledger_upper':Bsrc,'conditional_zero_initial_energy_upper':44*Bsrc,'conditional_zero_initial_mean_trace_squared_upper':176*Bsrc},'RECONSTRUCTED_CLIPPED_SOURCE_LEDGER')
 displayed_bulk=F(8526418129275428520152526);displayed_boundary=F(3708320490000);displayed_B=F(4263209064642351515776089);displayed_trace=F(27392057158547)
 need(integral_bulk<=displayed_bulk and integral_boundary==displayed_boundary and Bsrc<=displayed_B and 176*Bsrc<=displayed_trace**2,'PRINTED_OUTWARD_FUTURE_LEDGER_AND_AUXILIARY_TRACE')
 need(displayed_trace>F('5.4e12'),'NO_NUMERICAL_IMPROVEMENT_CLAIM_FROM_AUXILIARY_BOUND')
 common_imbalance=ell**2/4;need(F(i['shared_projection']['packet_common_principal_imbalance_upper'])==common_imbalance==F(1,4000000000000),'PRINCIPAL_PACKET_COMMON_GUARD')
 need(F(i['independent_route']['stronger_principal_imbalance_upper'])<=common_imbalance and ell*(1-F(3,2)*ell)>0,'INDEPENDENT_PACKET_GUARD_AND_NONZERO_BOTH_DIRECTIONS')
 packet=t['packet_diagnostic'];need(packet['not_current_time_symmetric_or_pure_oneway'] is True and packet['is_physical_energy_flux_measurement'] is False and packet['past_incoming_only_state_excluded'] is False,'PRESENT_PACKET_SCOPE')
 need(F(packet['Q_center'])==1 and F(packet['Q_r_center'])==0 and F(packet['central_scaled_radius'])==F(3,2),'ACTUAL_PACKET_CENTER')
 inverse=t['inverse_history'];cM=F(131,10500000000);bM=F(917,62500000);Hw=cM/3;Hrw=(3*ell**2+bM)*Hw+cM;Vw=3*ell**2*Hw+cM;N=max((1+2*Hw+2*T*Vw)/3,1+2*Hrw+2*Vw);Nguard=F(20000001,20000000);q=F(1201,3200000);kappa=q*Nguard/(1-q);kappaguard=F(47,125000)
 expected_inverse={'cM':cM,'Hw':Hw,'Hrw':Hrw,'Vw':Vw,'N_exact_majorant':N,'N_guard':Nguard,'q_retained':q,'kappa_exact':kappa,'kappa_guard':kappaguard}
 need(all(F(inverse[key])==value for key,value in expected_inverse.items()) and fractions(i['shared_projection']['inverse_history'])==expected_inverse,'COMPLETE_FIVE_MEAN_HISTORY_INVERSE_JOIN')
 need(N<Nguard and kappa<kappaguard<1,'STRICT_IMAGE_ONLY_INVERSE_GUARDS')
 need(inverse['on_admissible_h_only_image'] is True and inverse['onto_arbitrary_histories'] is False and inverse['past_state_constructed'] is False,'IMAGE_NOT_ONTO_OR_PAST_CONSTRUCTION')
 null=t['null_constraint'];need(t['background_null_identity']['screen_dimension']==i['shared_projection']['null_screen_dimension']==3 and null['integrating_factor']==i['shared_projection']['null_integrating_factor']=='R0^2','AFFINE_NULL_DIMENSION_AND_INTEGRATING_FACTOR')
 need(null['general_second_scalar_cross_term']=='-2 Phi0prime Phi2prime retained before setting pulled-back Phi2=0' and null['actual_area_second_formula']==i['shared_projection']['area_second_formula']=='average(L2+L1^2); no sign from L2 alone','GENERAL_NULL_SECOND_VARIATION_BEFORE_NEW_STATE_CHOICE')
 need(null['full_characteristic_solution_constructed'] is False and null['current_shell_expansion_or_H1_bound_established'] is False and t['forced_energy']['past_source_history_computed'] is False and t['background_null_identity']['past_coefficient_enclosures_computed'] is False,'NO_PAST_TRANSFER_OR_FULL_CHARACTERISTIC_CONSTRUCTION')
 need(F(-1)<0 and F(-1)+F(2)**2>0,'LOG_DENSITY_SIGN_DOES_NOT_FIX_AREA_SIGN')
 need(t['verification']['finite_checks']==len(t['verification']['check_names'])==32 and i['verification']['finite_checks']==len(i['verification']['check_names']),'OBSERVED_FINITE_GATE_ROSTERS')
 need(all(type(v)is int and v==0 for v in s['scientific_activity'].values()) and t['verification']['new_PDE_solves']==t['verification']['new_background_or_evolved_field_queries']==t['verification']['new_bump_queries']==0 and i['verification']['new_PDE_solves']==i['verification']['new_field_background_bump_or_source_queries']==i['verification']['project_imports']==0,'ZERO_NEW_SCIENTIFIC_ACTIVITY')
 need(regular(CAP/'DERIVATION.md'),'ROOT_DERIVATION_REGULAR')
 result={'schema':'HDBLAST_INCOMING_PREPARATION_LEDGER_HEADLINES_V1','status':'PASS_EXACT_CONDITIONAL_INCOMING_PREPARATION_LEDGER_HEADLINES','root_derivation_sha256':sha(CAP/'DERIVATION.md'),
  'retained_inputs':{'conditional_continuous_operator_premises':{'path':PREMISES,'sha256':PREMISES_PIN,'status':premises['status']},'continuous_source_report':{'path':RAW,'sha256':RAW_PIN,'old_cells':81,'proof_reproved':False},'new_matrix_operator_bound_proof':False},
  'forced_future_energy':{'mean_order':['K','L','N','J','w'],'Gram':gram,'constants':expected,'forced_growth':C,'auxiliary_zero_Cauchy_energy_cap':44*Bsrc,'auxiliary_zero_Cauchy_scalar_trace_squared_cap':176*Bsrc,'PDE_source_metric_multiplier':-2,'PDE_source_scalar_multiplier':1,'bulk_one_half_is_Young_factor_not_PDE_Hessian_half':True,'increasing_radial_boundary_convention':True,'nonsymmetric_Robin_feedback_retained':True},
  'future_source_ledger':{'geometry':geometry,'cone_intersecting_cells':6,'active_clipped_cells':5,'bulk_integral_cap':integral_bulk,'boundary_integral_cap':integral_boundary,'Bsrc':Bsrc,'outward_printed_caps':{'bulk':displayed_bulk,'boundary':displayed_boundary,'Bsrc':displayed_B,'auxiliary_zero_Cauchy_trace':displayed_trace},'numerical_improvement_over_stated_old_5_4e12_claimed':False,'physical_solution_alone_identified':False},
  'actual_first_packet':{'principal_master_only':True,'normalized_imbalance_interval':['0',str(common_imbalance)],'imbalance_lower_endpoint_open':True,'toward_shell_plus_fraction_bounds':[str(F(1,2)),str(F(1,2)+common_imbalance/2)],'plus_fraction_lower_open':True,'away_from_shell_minus_fraction_bounds':[str(F(1,2)-common_imbalance/2),str(F(1,2))],'minus_fraction_upper_open':True,'principal_transport_speeds':{'plus':-1,'minus':1},'both_present_directions_nonzero':True,'physical_energy_or_radiation_flux_measurement':False,'past_incoming_history_excluded':False,'background_time_symmetry_implies_first_tangent_symmetry':False},
  'partial_null_constraint':{'screen_dimension':3,'integrating_factor':'R0^2','new_state_conditions':t['candidate_null_state_assumptions'],'background_scalar_acceleration_cross_term_retained_before_Phi2_cut_zero':True,'actual_second_area_formula':'average(L2+L1^2)','area_sign_from_second_log_density_alone':False,'full_Goursat_or_shell_H1_bound_established':False},
  'image_only_scalar_history_inverse':{'constants':expected_inverse,'inverse_C1_factor':1/(1-kappa),'norm':inverse['norm'],'observable':inverse['observable'],'homogeneous_difference_of_admitted_families_sharing_first':True,'normalized_scalar_mean_not_unnormalized_integral':True,'same_first_quadratic_clock_measure_terms_cancel_only_in_difference':True,'causal_segment':['1','5/4'],'other_canonical_sectors_fixed':True,'injective_on_retained_admissible_image':True,'onto_arbitrary_histories':False,'exterior_preparation_determined':False,'past_state_constructed':False},
  'finite_gates':{'theory':32,'independent':i['verification']['finite_checks'],'are_negative_subprocess_controls':False},
  'scope':{'forced_future_ledger_not_past_selection_law':True,'physical_energy_identified':False,'all_zero_physical_second_Cauchy_data_admissibility_established':False,'future_to_past_norm_transfer_proved':False,'source_cell_proof_reproved':False,'absolute_proper_observable_sign_or_first_quadratic_discarding_claimed':False,'physical_class_or_radius_selection_established':False,'uniform_effective_nonlinear_parameter_window_established':False},
  'new_activity':{'oracle_calls':0,'background_queries':0,'PDE_solves':0,'bump_queries':0,'project_imports':0,'first_jet_calls':0}}
 payload=(json.dumps(encode(result),sort_keys=True,indent=2)+'\n').encode()
 with out.open('xb') as stream:stream.write(payload)
 print(json.dumps({'status':'PASS','witness':'EXACT_CONDITIONAL_INCOMING_PREPARATION_LEDGER_HEADLINES','output_sha256':hashlib.sha256(payload).hexdigest()},sort_keys=True))

if __name__=='__main__':
 try:main()
 except (Rejected,OSError,ValueError,KeyError,TypeError,ArithmeticError) as error:
  print(json.dumps({'status':'FAIL','witness':str(error) if isinstance(error,Rejected) else 'HEADLINE_INPUT_FAILURE'},sort_keys=True));sys.exit(2)
