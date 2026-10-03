"""Invented decimal-output accounting fixtures only; no model/source/array calls."""
import copy,hashlib,json,sys,tempfile
from decimal import Decimal,localcontext
from fractions import Fraction as Q
from pathlib import Path
from review_outputs import CASES,audit_values,witness,margins,authenticate_fresh_package_manifest,sha


def need(ok,message):
 if not ok:raise RuntimeError(message)
def text(x):
 with localcontext() as ctx:
  ctx.prec=100;return str(Decimal(x.numerator)/Decimal(x.denominator))


def fixtures():
 primary={'records':[]};independent={'rows':[]}
 for s in ('positive_B','signed_uB'):
  for t in ('coarse','fine'):
   record={'source':s,'setting':t,'controls':[]}
   n=129 if t=='coarse' else 257
   ds=Q('3e-6')*(1 if s=='positive_B' else -1);mom=Q('2e-10')
   for order in (24,32):
    control={'quadrature_order':order,'levels':{}}
    for dps in ('80','100'):
     control['levels'][dps]=[]
     for K in (64,128,256):
      values={'DeltaR':ds,'S_ab':Q(0),'I_ab':ds,'I_A':ds-mom,'D_S':ds,'D_cont':Q(0),'D_cont_A':mom,'E_Q':ds,'E_Q_A':ds-mom,
              'E_momentum':mom,'E_flow':Q(0),'E_operator':Q(0),'E_operator_initial':Q(0),'E_operator_end':Q(0),'E_reconstruction':Q(0),
              'canonical_constant_drift':Q(0),'source_primitive_residual':Q(0),'E_flow_triangle_bound':Q('1e-12'),'midpoint_E_flow':Q(0),'midpoint_flow_triangle_bound':Q('1e-12')}
      values.update({k:Q(0) for k in ('midpoint_contact_plus_flow_density_gap','midpoint_E_momentum','midpoint_E_operator')})
      row={'K':K,**{k:text(v) for k,v in values.items()}}
      row.update({f:['0']*n for f in ('R','P','F','contact_R','contact_P','contact_F')})
      row.update({f:['0']*3 for f in ('discrete_contact_R_anchors','discrete_contact_P_anchors','discrete_rho0_anchors','discrete_p0_anchors')})
      control['levels'][dps].append(row)
    record['controls'].append(control)
   primary['records'].append(record)
   for K in (64,128,256):
    fields={'DeltaR':ds,'S_ab':Q(0),'I_ab':ds,'D_S':ds,'D_cont':Q(0),'E_Q':ds,'E_flow':Q(0),'triangle_bound':Q('1e-12'),
            'E_operator':Q(0),'E_evolution_ledger':Q(0),'midpoint_E_flow':Q(0),'midpoint_triangle_bound':Q('1e-12'),
            'forcing_control_gap':Q('1e-10'),'ledger_control_gap':Q('1e-10'),'joint_control_gap':Q('2e-10')}
    fields.update({k:Q(0) for k in ('midpoint_density_gap','midpoint_E_operator','midpoint_projection_closure','initial_R_reconstruction_error','endpoint_R_reconstruction_error')})
    if t=='fine':fields.update(S_ab_double_global=Q(0),S_native_minus_double=Q(0))
    control_values={name:{'I_ab':text(ds-Q(delta))} for name,delta in [('16_16','2e-10'),('16_24','1e-10'),('24_24','0')]}
    independent['rows'].append({'source':s,'setting':t,'K':K,'times':[text(Q(-9,2)+Q(j,n-1)) for j in range(n)],
      'precisions':{dps:{'fields':{k:text(v) for k,v in fields.items()},'controls':copy.deepcopy(control_values)} for dps in ('80','100')},
      'control_profiles':{name:{f:['0']*n for f in ('R','P','F')} for name in ('16_16','16_24','24_24')},
      'contact_profiles':{f:['0']*n for f in ('R','P','F')},'baseline_profiles':{f:['0']*n for f in ('r0','p0')}})
 validation={'failures':[],'classification':'LEDGER_ERROR_DEMONSTRATED','attribution_cases':[{'source':s,'setting':t,'K':k} for s,t,k in CASES]}
 return primary,independent,validation


def main():
 checks=[]
 need(not witness(Q('2e-6'),Q(0),Q('2e-6')),'Strict magnitude boundary accepted');checks.append('strict_threshold_equality_rejected')
 need(witness(Q('-3e-6'),Q('3e-7'),Q('-2.7e-6')),'Signed inclusive fractional witness rejected');checks.append('negative_signed_fraction_boundaries_accepted')
 p,i,v=fixtures();a=audit_values(p,i,v)
 need(not a['accounting_audit_errors'] and len(a['case_reviews'])==12 and len(a['cross_control_matrix'])==144,'Invented accounting contract failed');checks.append('all12_cases_all144_cross_pairs')
 bad=copy.deepcopy(p);bad['records'][0]['controls'][1]['levels']['100'][0]['E_momentum']='-0.0000000002'
 b=audit_values(bad,i,v);need(any('correction' in e['check'] for e in b['accounting_audit_errors']),'Momentum sign mutation hidden');checks.append('momentum_sign_mutation_detected')
 bad=copy.deepcopy(i);bad['rows'][0]['precisions']['100']['fields']['midpoint_E_operator']='0.0000001'
 b=audit_values(p,bad,v);need(any('midpoint' in e['check'] for e in b['accounting_audit_errors']),'Midpoint operator mutation hidden');checks.append('midpoint_operator_mutation_detected')
 badv=copy.deepcopy(v);badv.update(classification='CONSISTENCY_FAILURE',attribution_cases=[],failures=[{'route':'fabricated','gate':'invented_consistency_failure','exact_rational_value':'3/10000000','threshold':'1/5000000'}])
 c=audit_values(p,i,badv);need(c['classification']=='CONSISTENCY_FAILURE' and len(c['canonical_same_case_candidates'])==12,'Failures did not block accepted attribution');checks.append('global_failure_blocks_candidate_witnesses')
 # Authenticate package transition against actual invented manifest bytes/payload.
 with tempfile.TemporaryDirectory(dir='/tmp') as td:
  root=Path(td)/'original';pkg=Path(td)/'package';root.mkdir();pkg.mkdir()
  (root/'source.txt').write_text('invented registered source; no physics')
  reg={'files':{'source.txt':sha(root/'source.txt')}}
  (root/'FULL_REGISTRATION.json').write_text(json.dumps(reg))
  (root/'FREEZE_RECEIPT.json').write_text(json.dumps({'public_freeze_commit':'a'*40,'registration_sha256':sha(root/'FULL_REGISTRATION.json')}))
  for name in ('source.txt','FULL_REGISTRATION.json','FREEZE_RECEIPT.json'):(pkg/name).write_bytes((root/name).read_bytes())
  manifest=pkg/'MANIFEST.json';manifest.write_text(json.dumps({'files':{name:sha(pkg/name) for name in ('source.txt','FULL_REGISTRATION.json','FREEZE_RECEIPT.json')}}))
  proof=authenticate_fresh_package_manifest(manifest,root,sha(root/'FULL_REGISTRATION.json'),sha(manifest))
  need(proof['payload_files_verified']==3,'Package payload universe differs');checks.append('actual_manifest_and_payload_byte_authentication')
  (pkg/'source.txt').write_text('altered payload')
  try:authenticate_fresh_package_manifest(manifest,root,sha(root/'FULL_REGISTRATION.json'),sha(manifest))
  except ValueError:checks.append('package_payload_corruption_rejected')
  else:raise RuntimeError('Changed payload accepted')
  try:authenticate_fresh_package_manifest(manifest,root,sha(root/'FULL_REGISTRATION.json'),'b'*64)
  except ValueError:checks.append('reported_hex_without_matching_manifest_rejected')
  else:raise RuntimeError('Unbound manifest hash accepted')
 result={'status':'PASS_INVENTED_ACCOUNTING_AND_PACKAGE_PROVENANCE_REVIEW','physical_evaluations':0,'actual_result_paths_read':False,
         'fixture_scope':'Invented decimal dictionaries test advisory accounting, not a claimed physically valid experiment.','checks':checks,'python_optimization':sys.flags.optimize,
         'review_tool_sha256':sha(Path(__file__).with_name('review_outputs.py'))}
 name='ACCOUNTING_TESTS_OPTIMIZED.json' if sys.flags.optimize else 'ACCOUNTING_TESTS_NORMAL.json'
 Path(__file__).with_name(name).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
