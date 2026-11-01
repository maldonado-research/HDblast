#!/usr/bin/env python3
"""Exact end-to-end joins of the new canonical scientific artifacts."""
import argparse,importlib.util,json,sys
from fractions import Fraction as F
from pathlib import Path
HERE=Path(__file__).resolve().parent;CAP=HERE.parent
def main():
 p=argparse.ArgumentParser();p.add_argument('--repo-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 spec=importlib.util.spec_from_file_location('root_checkpoint_guards',HERE/'run_checkpoint.py');g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
 repo=a.repo_root.resolve(strict=True);out=a.output.absolute()
 g.need(not out.exists()and all(not x.is_symlink()for x in [out,*out.parents]),'FRESH_OUTPUT')
 out=out.resolve()
 g.need(all(out!=x and x not in out.parents for x in [repo,CAP]),'OUTPUT_OUTSIDE_CODE_AND_SOURCE')
 own=g.authenticate(CAP);plan=g.read(HERE/'REPLAY_PLAN.json');before=g.sources(repo,plan)
 initial=g.read(CAP/'initial/references/RAW_BOXES.json');conditioned=g.read(CAP/'initial/references/CONDITIONED_RAW_BOXES.json')
 r=g.read(CAP/'response/references/FINAL_NORMAL.json');science=g.read(CAP/'response/FROZEN_SCIENTIFIC_REFERENCE.json')
 g.need({k:v for k,v in r.items()if k not in ('source_authentication','actual_first_oracle_calls','unchanged_frozen_scientific_values')}==science,'ORIGINAL_SCIENTIFIC_FIELDS')
 g.need(initial['partition']=={'N':1024,'cells':994,'tail_width':'1/64'},'FROZEN_N1024_PARTITION')
 g.need(len(initial['profile_cells'])==len(conditioned['profile_cells'])==994,'CONTINUOUS_CELL_COVERAGE')
 g.need([x['r']for x in initial['targets']]==['3/2','5/4'],'FROZEN_INITIAL_TARGETS')
 g.need(r['conditioned_adapter_independent_inclusion_checks']==14910,'ALL_CONDITIONED_COMPONENTS')
 q=F(r['q']['exact']);B=F(r['free_X_cap']['exact']);S=F(r['source_history_X_cap']['exact']);E=F(r['full_resolvent_minus_free_X_cap']['exact'])
 g.need(q==F(1201,3200000)and E==(q*B+S)/(1-q),'NEW_FAMILY_RESOLVENT')
 g.need(S==F(5,4)*17387513156573+3851400,'COMPLETE_SOURCE_HISTORY_ONCE')
 intervals={k:list(map(F,r[k]['exact']))for k in ('proper_M2_target','proper_v2','proper_Theta2')}
 g.need(intervals['proper_M2_target'][0]>0 and intervals['proper_v2'][1]<0 and intervals['proper_Theta2'][1]<0,'STRICT_REGISTERED_SIGNS')
 outward={'proper_M2_target':['1.1296e30','1.9779e30'],'proper_v2':['-5.4520e16','-2.6521e16'],'proper_Theta2':['-9.0208e21','-4.4458e21']}
 loose={'proper_M2_target':['1.1e30','2e30'],'proper_v2':['-6e16','-2e16'],'proper_Theta2':['-1e22','-4e21']}
 for bounds in (outward,loose):
  for key,ends in bounds.items():g.need(F(ends[0])<=intervals[key][0]<=intervals[key][1]<=F(ends[1]),'EXACT_OUTWARD_HEADLINE_ENCLOSURE')
 g.need(F(r['cone']['T']['exact'])==F(5,4)and F(r['first']['rho']['exact'][0])>0,'POSITIVE_BACKGROUND_PROPER_CLOCK')
 g.need(r['first']['M2_correction']['exact']!='0'and r['first']['Q1']['exact']!='0','FIRST_QUADRATIC_TERMS_RETAINED')
 g.need(r['actual_first_oracle_calls']==1 and bool(r['first']['raw_actual_shell_2jet']['m']),'ACTUAL_FIRST_JETS_RETAINED')
 g.need(g.authenticate(CAP,own['manifest_sha256'])==own and g.sources(repo,plan)==before,'SOURCES_UNCHANGED')
 result={'status':'PASS_NEW_PREPARED_MEAN_HEADLINES','T':'5/4','N':1024,'continuous_initial_cells':994,
  'conditioned_component_inclusions':14910,'full_five_mean_block':r['mean_operator'],
  'q':str(q),'new_family_free_X_cap':str(B),'complete_source_history_X_cap':str(S),'full_resolvent_error_X_cap':str(E),
  'proper_M2_target':r['proper_M2_target']['exact'],'proper_v2':r['proper_v2']['exact'],'proper_Theta2':r['proper_Theta2']['exact'],
  'original_scientific_fields_unchanged':True,'all_three_strict_registered_signs':True,
  'safe_outward_headline_enclosures':outward,'safe_looser_exact_enclosures':loose,
  'point_target_only':True,'continuous_signed_history_claimed':False,'baseline_Theta0_positive':True,
  'different_explicit_preparation_from_old_weighted_only_family':True,
  'same_physical_second_jet_under_conditioning_is_retained_analytic_input':True,
  'effective_nonzero_parameter_radius':None,'SI_units_conversion_supplied':False,
  'physical_origin_or_radiation_established':False,'total_contraction_established':False,
  'source_authentication_before_and_after':True,'inherited_source_capsules_authenticated':len(before),
  'retained_nonlinear_premises_proved_by_byte_authentication':False,'old_family_total_cap_reused':False,
  'scope':'Formal second changes of the explicitly prepared family at the same anchored proper clock; no total contraction, negative total volume, radiation, finite parameter radius or cosmological-origin conclusion.'}
 out.write_bytes(g.canonical(result));print(json.dumps({'status':'PASS','witness':'NEW_PREPARED_MEAN_HEADLINES'},sort_keys=True))
if __name__=='__main__':
 try:main()
 except Exception as exc:
  print(json.dumps({'status':'FAIL','witness':str(exc)},sort_keys=True));sys.exit(2)
