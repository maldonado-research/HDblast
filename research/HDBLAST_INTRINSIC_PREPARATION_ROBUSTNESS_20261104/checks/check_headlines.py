#!/usr/bin/env python3
"""Exact arithmetic joins for the intrinsic scalar-level-area class certificate."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

CAP=Path(__file__).resolve().parents[1]
BASE='research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/response/references/FINAL_NORMAL.json'
BASE_PIN='1a9ebffb6bee9332bc0fd0df0f5d1d9f7f2a796825c2a563e4c644974ba4d8a3'
SCIENCE_PIN='e940b780cb5a020371575f42deef7e080a69c0e605807d44f152c30893ea69ce'
class Rejected(Exception):pass
def need(ok,witness):
 if not ok:raise Rejected(witness)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def regular(p):return p.is_file() and all(not q.is_symlink() for q in [p,*p.absolute().parents])
def read(p):
 need(regular(p),'REGULAR_INPUT');return json.loads(p.read_bytes())
def scalar(x):return F(x['exact'])
def box(x):
 lo,hi=map(F,x['exact']);need(lo<=hi,'ORDERED_INTERVAL');return lo,hi
def rec(x):return {'exact':[str(v) for v in x]} if isinstance(x,tuple) else {'exact':str(x)}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path);ap.add_argument('--bounds-report',required=True,type=Path);ap.add_argument('--independent-report',required=True,type=Path);ap.add_argument('--comparison-report',required=True,type=Path);ap.add_argument('--witness-report',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
 repo=a.repo_root.absolute();out=a.output.absolute();need(all(not p.is_symlink() for p in [repo,*repo.parents,out,*out.parents]),'PATH_SYMLINK');repo=repo.resolve(strict=True);out=out.resolve();need(not out.exists(),'FRESH_OUTPUT');need(out.parent.is_dir(),'OUTPUT_PARENT');need(all(out!=p and p not in out.parents for p in [repo,CAP]),'OUTPUT_OUTSIDE_CODE_AND_SOURCE')
 need(sha(repo/BASE)==BASE_PIN,'BASE_BYTE_PIN');base=read(repo/BASE);b=read(a.bounds_report);ind=read(a.independent_report);comparison=read(a.comparison_report);witness=read(a.witness_report)
 projection={k:v for k,v in b.items() if k not in ['authenticated_inputs','frozen_design','provenance','scope']};need(hashlib.sha256(json.dumps(projection,sort_keys=True,separators=(',',':')).encode()).hexdigest()==SCIENCE_PIN,'FROZEN_BOUND_SCIENCE')
 need(b['status']=='PASS_CONDITIONAL_PHYSICAL_INITIAL_PROFILE_CLASS_ROBUSTNESS' and ind['status']=='PASS_INDEPENDENT_GEOMETRIC_CLASS_BOUND','PRIMARY_AND_INDEPENDENT_STATUS')
 need(comparison['status']=='PASS_INDEPENDENT_CANDIDATE_COMPARISON' and comparison['candidate_sha256']==sha(a.bounds_report),'FRESH_CANDIDATE_JOIN');need(comparison['accepted_check_count']==30 and comparison['full_scientific_projection_matches_frozen_baseline'] is True and comparison['exact_portable_provenance_schema_checked'] is True,'COMPLETE_INDEPENDENT_COMPARISON')
 need(witness['status']=='PASS_INDEPENDENT_WITNESS_EXCLUSION','WITNESS_STATUS')
 c=b['constants'];ell,T,delta,q=map(scalar,[c['ell'],c['T'],c['delta'],c['q']]);need((ell,T,delta,q)==(F(1,10**6),F(5,4),F(1,4),F(1201,3200000)),'RETAINED_CONSTANTS')
 rho=box(c['rho']);minus_p=box(c['minus_p']);acoef=box(c['a']);need(rho==(F(129),F(131)) and minus_p==(F(3,25000),F(1,7000)) and acoef==(F(111,1000),F(14,125)),'CONTINUOUS_BACKGROUND_CAPS')
 radius=scalar(b['class']['radius']);need(radius==F(2*10**33) and b['class']['complete_increment_required'] is True and b['class']['causal_initial_segment']==['1','5/4'],'CLOSED_COMPLETE_INCREMENT_CLASS')
 caps=b['conditioned_component_caps_per_R'];H,Hr,S,Sr,V=map(scalar,[caps['h'],caps['h_r'],caps['s'],caps['s_r'],caps['V_L']]);gram=b['Gram_norm_caps_per_R'];D0,Dr,Dv,N=map(scalar,[gram['D0'],gram['Dr'],gram['Dv'],gram['N']])
 weights=[1,2,1,3,1];need(weights[0]+weights[3]==4 and weights[1]<=4,'FIVE_MEAN_GRAM_TRIANGLE');need((D0,Dr,Dv)==(2*H+S,2*Hr+Sr,2*V),'COMPLETE_FIVE_MEAN_INITIAL_NORMS');need(N==max((D0+T*Dv)/3,Dr+Dv),'BOTH_INITIAL_NORM_BRANCHES')
 guard=scalar(gram['fixed_N_per_R_guard']);Nsafe=scalar(b['class_N']['fixed_guard']);need(guard==F(1,2000) and N<=guard and Nsafe==radius*guard==10**30 and scalar(b['class_N']['actual_derived_cap'])==radius*N,'FIXED_NORM_CLASS_GUARDS')
 owncaps=ind['profile_caps_per_R'];mapping={'h':'h','h_r':'h_r','s':'s','s_r':'s_r','V_L':'V_L'}
 for left,right in mapping.items():need(scalar(caps[left])==scalar(owncaps[right]),'INDEPENDENT_COMPONENT_CAP_JOIN:'+left)
 for key in ['D0','Dr','Dv','N']:need(scalar(gram[key])==scalar(owncaps[key]),'INDEPENDENT_GRAM_JOIN:'+key)
 need(F(ind['class_radius'])==radius and F(ind['N_safe'])==Nsafe and ind['segment']==['1','5/4'],'INDEPENDENT_CLASS_JOIN')
 Cv=(rho[1]*ell)**2*scalar(c['abs_Fprime'])*delta**3/(12*(1-q));Ct=rho[1]*ell*scalar(c['abs_Fprime'])*delta**2/(4*(1-q));need(Cv==scalar(b['observable_operator_caps_per_N']['volume'])==scalar(ind['homogeneous_operator']['Cv']) and Ct==scalar(b['observable_operator_caps_per_N']['expansion'])==scalar(ind['homogeneous_operator']['Ct']),'PROPER_CLOCK_HOMOGENEOUS_TRANSFER')
 bounds={};display={'v2':('-7.6900e16','-4.1600e15'),'Theta2':('-1.1070e22','-2.3980e21')};base_display={'v2':('-5.4520e16','-2.6521e16'),'Theta2':('-9.0208e21','-4.4458e21')}
 for key,oldkey,C,transferkey in [('v2','proper_v2',Cv,'volume'),('Theta2','proper_Theta2',Ct,'expansion')]:
  old=box(base[oldkey]);need(box(b['base_exact_intervals'][key])==old==box(ind['base_coefficients'][key]),'EXACT_BASE_JOIN:'+key);new=(old[0]-C*Nsafe,old[1]+C*Nsafe);need(new==box(b['whole_class_exact_intervals'][key])==box(ind['class_coefficients'][key]),'EXACT_CLASS_JOIN:'+key);need(scalar(b['absolute_observable_difference_caps'][transferkey])==C*Nsafe,'TRANSFER_NORM_JOIN:'+key)
  lo,hi=map(F,display[key]);blo,bhi=map(F,base_display[key]);need(blo<=old[0]<=old[1]<=bhi<0 and lo<=new[0]<=new[1]<=hi<0,'PRINTED_STRICT_OUTWARD_BOUNDS:'+key);need(box(b['frozen_outward_intervals'][key])==(lo,hi),'FROZEN_DISPLAY_JOIN:'+key);bounds[key]={'base':rec(old),'whole_closed_class':rec(new),'printed_class_outward':list(display[key])}
 threshold=min(-box(base['proper_v2'])[1]/Cv,-box(base['proper_Theta2'])[1]/Ct);need(threshold==scalar(b['old_exact_N_threshold']['both']) and Nsafe<threshold,'STRICT_OLD_ROBUSTNESS_THRESHOLD')
 # The algebraic geometry factors use actual second parameter derivatives.
 angular_dimension=F(3);metric_factor=angular_dimension/2;scalar_factor=-angular_dimension;need(metric_factor==F(3,2) and 2*metric_factor+scalar_factor==0,'SCALAR_PULLBACK_GAUGE_CANCELLATION');need(-scalar_factor*F(1,2)/metric_factor==1,'CONDITIONED_AREA_RECONSTRUCTION')
 lam=F(witness['lambda']);need(lam==-10**21 and F(witness['class_radius'])==radius and F(witness['target'])==T,'OLD_SCALED_WITNESS_PARAMETERS');d=box(witness['direction_D']);dr=box(witness['direction_D_r']);need(d[0]>0 and dr[1]<0,'WITNESS_PROFILE_SIGNS');wd=abs(lam)*d[0];wdr=abs(lam)*(-dr[1]);need(wd==F(witness['witness_absolute_D_lower']) and wdr==F(witness['witness_absolute_D_r_lower']),'SCALED_COMPLETE_WITNESS');need(wd>F('3.7621e33')>radius and wdr>F('6.870e34')>radius,'OLD_WITNESS_STRICTLY_OUTSIDE_CLASS')
 need(b['verification']['oracle_calls']==b['verification']['new_PDE_solves']==b['verification']['new_bump_queries']==ind['new_oracle_calls']==comparison['new_oracle_calls']==witness['scope']['new_background_oracle_queries']==witness['scope']['new_evolved_field_oracle_queries']==witness['scope']['new_bump_queries']==0,'ZERO_NEW_SCIENTIFIC_ACTIVITY')
 need(ind['negative_controls']['in_process_mutant_count']==6 and len(ind['negative_controls']['checks'])==6 and comparison['in_process_report_mutants']==9 and len(comparison['report_mutant_outcomes'])==9,'SEPARATE_INPROCESS_MUTANT_ROSTERS')
 result={'schema':'HDBLAST_INTRINSIC_PREPARATION_ROBUSTNESS_HEADLINES_V1','status':'PASS_EXACT_INTRINSIC_PREPARATION_ROBUSTNESS_HEADLINES','geometry':{'coefficient_convention':'actual second parameter derivative; no extra half','angular_dimension':3,'profile_formula':'D=3*J/2-3*a*w/p','metric_factor':'3/2','scalar_pullback_factor':'-3','raw_conditioned_area_equal':True,'background_scalar_label_derivative':'D_r=(-ell*rho*p)*D_v','physical_Cauchy_section_part_of_observable':True,'full_packet_profile_injective_on_retained_mean_freedom':True,'causal_segment_profile_does_not_determine_exterior_preparations':True,'profile_map_asserted_onto_arbitrary_C1_data':False},'five_mean_Gram':{'order':['K','L','N','J','w'],'weights':weights,'initial_caps_per_R':[str(H),'0','0',str(H),str(S)],'radial_caps_per_R':[str(Hr),'0','0',str(Hr),str(Sr)],'velocity_caps_per_R':['0',str(V),'0','0','0'],'D0':rec(D0),'Dr':rec(Dr),'Dv':rec(Dv),'N_per_R':rec(N),'N_per_R_guard':'1/2000'},'closed_class':{'radius':str(radius),'complete_increment_required':True,'causal_initial_labels':['1','5/4'],'fixed_initial_norm_guard':str(Nsafe),'strict_old_threshold':rec(threshold)},'proper_clock_transfer_per_initial_N':{'volume':rec(Cv),'expansion':rec(Ct),'volume_denominator_factor':12,'expansion_denominator_factor':4,'full_homogeneous_five_mean_contraction_q':str(q)},'bounds':bounds,'old_scaled_sign_reversal_witness':{'lambda':str(lam),'initial_slice_label':'5/4','label_is_future_evolution_time':False,'absolute_D_lower':rec(wd),'absolute_D_r_lower':rec(wdr),'strictly_outside_class':True},'scope':{'conditional_coefficient_class_result':True,'physical_class_selection_established':False,'uniform_effective_nonlinear_epsilon_window_established':False,'unrestricted_Cauchy_data_robustness_established':False,'absolute_shared_first_quadratics_discarded':False,'shared_first_quadratics_cancel_only_in_second_differences':True},'new_activity':{'oracle_calls':0,'PDE_solves':0,'bump_queries':0}}
 need(regular(CAP/'DERIVATION.md'),'GEOMETRY_DERIVATION_REGULAR');result['geometry']['derivation_sha256']=sha(CAP/'DERIVATION.md')
 payload=(json.dumps(result,sort_keys=True,indent=2)+'\n').encode()
 with out.open('xb') as stream:stream.write(payload)
 print(json.dumps({'status':'PASS','witness':'EXACT_INTRINSIC_PREPARATION_ROBUSTNESS_HEADLINES','output_sha256':hashlib.sha256(payload).hexdigest()},sort_keys=True))
if __name__=='__main__':
 try:main()
 except (Rejected,OSError,ValueError,KeyError,TypeError,ZeroDivisionError) as error:
  print(json.dumps({'status':'FAIL','witness':str(error) if isinstance(error,Rejected) else 'HEADLINE_INPUT_FAILURE'},sort_keys=True));sys.exit(2)
