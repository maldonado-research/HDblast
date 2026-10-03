#!/usr/bin/env python3
"""Advisory post-public-freeze audit of JSON outputs; never imports physics or reads arrays."""
from __future__ import annotations
import argparse,csv,hashlib,json,re,shutil
from collections import Counter
from decimal import Decimal
from fractions import Fraction as Q
from pathlib import Path
from scientific_projection import compare,scientific_digest

SCIENCE=Q('2e-7');ARITHMETIC=Q('1e-12');ATTRIBUTION=Q('2e-6')
CASES=tuple((s,t,k) for s in ('positive_B','signed_uB') for t in ('coarse','fine') for k in (64,128,256))


def need(ok,message):
 if not ok:raise ValueError(message)
def number(value):
 need(isinstance(value,str) and 0<len(value)<=4096,'Scientific audit scalar must remain a decimal string')
 d=Decimal(value);need(d.is_finite() and abs(d.as_tuple().exponent)<=10000,'Invalid scientific scalar')
 return Q(d)
def rational(value):return Q(value)
def read(path):
 def unique(pairs):
  out={}
  for k,v in pairs:need(k not in out,'Duplicate JSON key');out[k]=v
  return out
 def invalid(value):raise ValueError('Nonfinite JSON literal '+value)
 return json.loads(Path(path).read_text(),object_pairs_hook=unique,parse_constant=invalid)
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def gap(values):return max((abs(x) for x in values),default=Q(0))
def margins(ds,dc,eq):
 return {'strict_magnitude_margin':abs(ds)-ATTRIBUTION,'continuum_fraction_margin':abs(ds)/10-abs(dc),'sampling_fraction_margin':abs(eq)-9*abs(ds)/10}
def witness(ds,dc,eq):
 m=margins(ds,dc,eq)
 return m['strict_magnitude_margin']>0 and m['continuum_fraction_margin']>=0 and m['sampling_fraction_margin']>=0
def fraction_map(values):return {k:str(v) for k,v in values.items()}
def exact_close(value,label,errors):
 if abs(value)>ARITHMETIC:errors.append({'check':label,'exact_signed_error':str(value),'threshold':str(ARITHMETIC)})


def lookups(primary,independent):
 p={};i={}
 for record in primary['records']:
  for control in record['controls']:
   for dps,rows in control['levels'].items():
    for row in rows:
     key=(record['source'],record['setting'],row['K'],control['quadrature_order'],dps)
     need(key not in p,'Duplicate primary case/control/context');p[key]=row
 for row in independent['rows']:
  key=(row['source'],row['setting'],row['K']);need(key not in i,'Duplicate independent case');i[key]=row
 expected={(s,t,k,q,dps) for s,t,k in CASES for q in (24,32) for dps in ('80','100')}
 need(set(p)==expected and set(i)==set(CASES),'Declared twelve-case/control/context universe is incomplete')
 return p,i


def audit_values(primary,independent,validation):
 p,i=lookups(primary,independent);errors=[];rows=[];cross_matrix=[];witnesses=[]
 for case in CASES:
  source,setting,K=case;pr=p[(*case,32,'100')];ir=i[case];iv=ir['precisions']['100']['fields']
  pv={k:number(pr[k]) for k in ('DeltaR','S_ab','I_ab','I_A','D_S','D_cont','D_cont_A','E_Q','E_Q_A','E_momentum','E_flow','E_operator','E_operator_initial','E_operator_end','E_reconstruction','canonical_constant_drift','source_primitive_residual','E_flow_triangle_bound','midpoint_E_flow','midpoint_flow_triangle_bound')}
  v={k:number(iv[k]) for k in ('DeltaR','S_ab','I_ab','D_S','D_cont','E_Q','E_flow','triangle_bound','E_operator','E_evolution_ledger','midpoint_E_flow','midpoint_triangle_bound','forcing_control_gap','ledger_control_gap','joint_control_gap')}
  pv.update({k:number(pr[k]) for k in ('midpoint_contact_plus_flow_density_gap','midpoint_E_momentum','midpoint_E_operator')})
  v.update({k:number(iv[k]) for k in ('midpoint_density_gap','midpoint_E_operator','midpoint_projection_closure','initial_R_reconstruction_error','endpoint_R_reconstruction_error')})
  for name,value in {'P signedledger':pv['D_S']-pv['D_cont']-pv['E_Q'],'I signedledger':v['D_S']-v['D_cont']-v['E_Q'],
                     'Primary discrete correction':pv['I_ab']-pv['I_A']-pv['E_momentum'],
                     'Primary raw-to-matched residual':pv['D_cont_A']-pv['D_cont']-pv['E_momentum'],
                     'Primary raw-to-matched sampling':pv['E_Q']-pv['E_Q_A']-pv['E_momentum'],
                     'Primary signedoperator':pv['E_operator']-pv['E_operator_end']+pv['E_operator_initial'],
                     'Primary reconstruction components':pv['E_reconstruction']-pv['canonical_constant_drift']-pv['source_primitive_residual'],
                     'Primary continuous components':pv['D_cont']-pv['E_flow']-pv['E_operator']-pv['E_reconstruction'],
                     'Independent continuous components':v['D_cont']-v['E_flow']-v['E_operator']-v['E_evolution_ledger'],
                     'Independent operator brackets':v['E_operator']-v['endpoint_R_reconstruction_error']+v['initial_R_reconstruction_error'],
                     'Primary midpoint components':pv['midpoint_contact_plus_flow_density_gap']-pv['midpoint_E_flow']-pv['midpoint_E_momentum']-pv['midpoint_E_operator'],
                     'Independent midpoint components':v['midpoint_density_gap']-v['midpoint_E_flow']-v['midpoint_E_operator'],
                     'Independent midpoint reported closure':v['midpoint_projection_closure']-v['midpoint_density_gap']+v['midpoint_E_flow']+v['midpoint_E_operator'],
                     'Identical saved endpoint change':pv['DeltaR']-v['DeltaR'],'Identical sampled Simpson':pv['S_ab']-v['S_ab']}.items():exact_close(value,f'{case}: {name}',errors)
  for route,signed,bound in [('primary',pv['E_flow'],pv['E_flow_triangle_bound']),('independent',v['E_flow'],v['triangle_bound']),('primary midpoint',pv['midpoint_E_flow'],pv['midpoint_flow_triangle_bound']),('independent midpoint',v['midpoint_E_flow'],v['midpoint_triangle_bound'])]:
   need(bound>=0,'Negative reported flow envelope')
   exact_close(max(Q(0),abs(signed)-bound),f'{case}: {route} flow exceeds envelope',errors)
  controls=ir['precisions']['100']['controls']
  q1616=number(controls['16_16']['I_ab']);q1624=number(controls['16_24']['I_ab']);q2424=number(controls['24_24']['I_ab'])
  exact_close(v['forcing_control_gap']-q2424+q1624,f'{case}: forcing-control sign',errors)
  exact_close(v['ledger_control_gap']-q1624+q1616,f'{case}: ledger-control sign',errors)
  exact_close(v['joint_control_gap']-q2424+q1616,f'{case}: joint-control sign',errors)
  both_witness=witness(pv['D_S'],pv['D_cont'],pv['E_Q']) and witness(v['D_S'],v['D_cont'],v['E_Q'])
  if both_witness:witnesses.append({'source':source,'setting':setting,'K':K})
  row={'source':source,'setting':setting,'K':K,'canonical_same_case_witness_conditions':both_witness,
       'primary_attribution_margins':fraction_map(margins(pv['D_S'],pv['D_cont'],pv['E_Q'])),
       'independent_attribution_margins':fraction_map(margins(v['D_S'],v['D_cont'],v['E_Q'])),
       'canonical_integral_gap_exact_rational':str(pv['I_ab']-v['I_ab']),
       'primary_terms':fraction_map(pv),'independent_terms':fraction_map(v),
       'flow_cancellation_fraction':{'primary':str(abs(pv['E_flow'])/pv['E_flow_triangle_bound']) if pv['E_flow_triangle_bound'] else None,
                                     'independent':str(abs(v['E_flow'])/v['triangle_bound']) if v['triangle_bound'] else None},
       'control_terms':{'forcing_gap':str(q2424-q1624),'ledger_gap':str(q1624-q1616),'joint_gap':str(q2424-q1616)}}
  if setting=='fine':
   row['same_trajectory_Simpson_native_minus_double']=str(number(iv['S_ab'])-number(iv['S_ab_double_global']))
   exact_close(number(iv['S_native_minus_double'])-number(iv['S_ab'])+number(iv['S_ab_double_global']),f'{case}: fine-sampling sign',errors)
  rows.append(row)
  # Every predefined cross-control pair at both contexts; contact-removed full mode profiles.
  for dps in ('80','100'):
   for pq in (24,32):
    prow=p[(*case,pq,dps)]
    for iq in ('16_16','16_24','24_24'):
     ic=ir['precisions'][dps]['controls'][iq]
     profiles={f:gap((number(a)-number(c))-(number(b)-number(d)) for a,c,b,d in zip(prow[f],prow['contact_'+f],ir['control_profiles'][iq][f],ir['contact_profiles'][f])) for f in ('R','P','F')}
     integral_gap=number(prow['I_ab'])-number(ic['I_ab'])
     cross_matrix.append({'source':source,'setting':setting,'K':K,'dps':dps,'primary_order':pq,'independent_control':iq,
                          'matched_integral_signed_gap':str(integral_gap),'mode_profile_max_gaps':fraction_map(profiles),
                          'largest_gap_over_registered_science_gate':str(max(abs(integral_gap),*profiles.values())/SCIENCE)})
  # Only these three direct discrete contact/baseline knots are cross-target audited.
  n=len(ir['times']);mid=n//2
  row['discrete_knot_gaps']={}
  for knot,j in zip(('initial','midpoint','endpoint'),(0,mid,n-1)):
   knot_index={'initial':0,'midpoint':1,'endpoint':2}[knot];L=-1/number(ir['times'][j])
   row['discrete_knot_gaps'][knot]=fraction_map({
    'density_contact':number(pr['discrete_contact_R_anchors'][knot_index])-number(ir['contact_profiles']['R'][j]),
    'pressure_contact':number(pr['discrete_contact_P_anchors'][knot_index])-number(ir['contact_profiles']['P'][j]),
    'conformal_baseline_density':L**4*number(pr['discrete_rho0_anchors'][knot_index])-number(ir['baseline_profiles']['r0'][j]),
    'conformal_baseline_pressure':L**4*number(pr['discrete_p0_anchors'][knot_index])-number(ir['baseline_profiles']['p0'][j])})
 failures=validation['failures'];need(type(failures) is list,'Missing authoritative failure universe')
 expected='CONSISTENCY_FAILURE' if failures else 'LEDGER_ERROR_DEMONSTRATED' if witnesses else 'NO_GATE_SCALE_ATTRIBUTION'
 need(validation['classification']==expected,'Canonical witness/failure interpretation differs from authoritative classification')
 need(validation['attribution_cases']==([] if failures else witnesses),'Authoritative same-case attribution list differs')
 severities=[]
 for failure in failures:
  value=abs(rational(failure['exact_rational_value']));threshold=rational(failure['threshold'])
  need(threshold>0,'Invalid failure threshold')
  severities.append({'identity':{k:v for k,v in failure.items() if k not in ('exact_rational_value','threshold')},'absolute_value':str(value),'threshold':str(threshold),'multiple_of_threshold':str(value/threshold)})
 severities.sort(key=lambda r:rational(r['multiple_of_threshold']),reverse=True)
 return {'classification':validation['classification'],'accounting_audit_errors':errors,'case_reviews':rows,'cross_control_matrix':cross_matrix,
         'failure_counts_by_gate':dict(Counter(f['gate'] for f in failures)),'failure_counts_by_route':dict(Counter(f.get('route','unspecified') for f in failures)),
         'failure_severity':severities,'canonical_same_case_candidates':witnesses,
         'scope':'Independent advisory recomputation from serialized outputs only; the frozen validator remains authoritative.'}


def verify_bundle(root,freeze,registration_pin,paths):
 need(re.fullmatch('[0-9a-f]{40}',freeze) and re.fullmatch('[0-9a-f]{64}',registration_pin),'Explicit public freeze/registration pins required')
 need(sha(root/'FULL_REGISTRATION.json')==registration_pin,'Local registration hash differs')
 reg=read(root/'FULL_REGISTRATION.json');receipt=read(root/'FREEZE_RECEIPT.json')
 need(receipt['public_freeze_commit']==freeze and receipt['registration_sha256']==registration_pin,'Freeze receipt differs')
 for name,pin in reg['files'].items():
  path=root/name;need(path.is_file() and not path.is_symlink() and path.resolve().is_relative_to(root.resolve()),'Invalid frozen path')
  need(sha(path)==pin,'Registered source/input bytes changed')
 p,i,v,vo,execution=[read(paths[x]) for x in ('primary','independent','validation','validation_optimized','execution')]
 need(p['producer_sha256']==reg['files']['code/diagnostic_primary.py'] and i['producer_sha256']==reg['files']['independent/diagnostic_independent.py'],'Producer pin differs')
 need(p['input_manifest_sha256']==i['input_manifest_sha256']==reg['files']['inputs/INPUT_MANIFEST.json'],'Input lineage differs')
 need(p['freeze_commit']==freeze and p['registration_sha256']==registration_pin,'Primary freeze differs')
 need(i['provenance']['freeze_commit']==freeze and i['provenance']['registration_sha256']==registration_pin,'Independent freeze differs')
 need(v['status']==vo['status']=='VERIFIED_COMPLETE_DIAGNOSTIC','No complete accepted validation')
 need({k:x for k,x in v.items() if k!='python_optimization'}=={k:x for k,x in vo.items() if k!='python_optimization'},'Normal/optimized validation differs')
 need(v['result_sha256']=={'primary':sha(paths['primary']),'independent':sha(paths['independent'])},'Validator result hashes differ')
 need(v['registration_sha256']==registration_pin and v['public_freeze_commit']==freeze,'Validation freeze differs')
 need(v['old_metric_status']==p['old_metric_status']==i['old_metric_status']=='FAIL','Historical FAIL was altered')
 physical=[r for r in execution if r.get('physical_route') is True]
 need([r['name'] for r in physical]==['primary','independent'],'Expected two complete physical executions')
 for r in physical:
  need(r['status']=='PASS_EXECUTION' and r['exit_code']==0 and not r['timed_out'] and not r['memory_exceeded'],'Execution did not pass')
  need(0<=r['elapsed_seconds']<=900 and 0<r['peak_rss_kib']<=262144,'Complete-route resource limit exceeded')
 return p,i,v,physical


def authenticate_fresh_package_manifest(path,root,registration_pin,reported_pin):
 need(path is not None and path.is_file() and not path.is_symlink(),'A fresh package manifest file is required for a package provenance transition')
 need(sha(path)==reported_pin,'Reported fresh package manifest pin does not match actual manifest bytes')
 package_root=path.parent.absolute()
 for parent in (package_root,*package_root.parents):need(not parent.is_symlink(),'Package root symlink ancestry forbidden')
 manifest=read(path);files=manifest['files'];reg=read(root/'FULL_REGISTRATION.json')
 need(files.get('FULL_REGISTRATION.json')==registration_pin,'Packaged registration differs')
 need(files.get('FREEZE_RECEIPT.json')==sha(root/'FREEZE_RECEIPT.json'),'Packaged freeze receipt differs')
 need(all(files.get(name)==pin for name,pin in reg['files'].items()),'Packaged registered source/input pins differ')
 for name,pin in files.items():
  need(isinstance(name,str) and not Path(name).is_absolute() and '..' not in Path(name).parts and '\\' not in name,'Unsafe manifest member')
  target=package_root/name
  need(target.is_file() and not target.is_symlink() and target.resolve().is_relative_to(package_root.resolve()),'Unsafe package payload path')
  need(all(not parent.is_symlink() for parent in target.parents if parent.is_relative_to(package_root)),'Package payload symlink ancestry forbidden')
  need(sha(target)==pin,'Package manifest payload bytes differ: '+name)
 return {'manifest_sha256':reported_pin,'payload_files_verified':len(files),'registration_sha256':registration_pin,'freeze_receipt_sha256':files['FREEZE_RECEIPT.json'],
         'scope':'Actual manifest bytes and every listed payload byte verified; frozen registration/input/source family preserved.'}


def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--review-after-public-freeze',action='store_true',required=True)
 parser.add_argument('--checkpoint-root',type=Path,required=True);parser.add_argument('--freeze-commit',required=True);parser.add_argument('--registration-sha256',required=True)
 for field in ('primary','independent','validation','validation-optimized','execution','output-dir'):parser.add_argument('--'+field,type=Path,required=True)
 for field in ('fresh-primary','fresh-independent','fresh-validation','fresh-package-manifest'):parser.add_argument('--'+field,type=Path)
 args=parser.parse_args();root=args.checkpoint_root.absolute();output=args.output_dir.absolute()
 need(not output.exists() and not output.resolve().is_relative_to(root.resolve()),'Choose a fresh external review directory')
 paths={key:getattr(args,key) for key in ('primary','independent','validation','validation_optimized','execution')}
 for path in paths.values():need(path.is_file() and path.suffix=='.json' and not path.is_symlink(),'Only explicit nonsymlink result JSON paths may be reviewed')
 p,i,v,executions=verify_bundle(root,args.freeze_commit,args.registration_sha256,paths)
 audit=audit_values(p,i,v)
 audit.update(schema_version=1,status='ADVISORY_AUDIT_CONSISTENT' if not audit['accounting_audit_errors'] else 'ADVISORY_ACCOUNTING_DIFFERENCE',old_metric_status='FAIL',
              freeze_commit=args.freeze_commit,registration_sha256=args.registration_sha256,
              actual_source_evaluations_by_review=0,input_array_values_read_by_review=False,
              raw_result_sha256={key:sha(path) for key,path in paths.items()},
              scientific_projection_sha256={kind:scientific_digest(value,kind) for kind,value in [('primary',p),('independent',i),('validation',v)]},
              execution_resources=[{'route':r['name'],'seconds':r['elapsed_seconds'],'peak_rss_kib':r['peak_rss_kib'],'wall_headroom_seconds':900-r['elapsed_seconds'],'memory_headroom_kib':262144-r['peak_rss_kib']} for r in executions])
 fresh=[args.fresh_primary,args.fresh_independent,args.fresh_validation]
 need(not any(fresh) or all(fresh),'Supply all three fresh result paths together')
 if all(fresh):
  fresh_objects=[read(path) for path in fresh]
  need(fresh_objects[2]['status']=='VERIFIED_COMPLETE_DIAGNOSTIC' and fresh_objects[2]['old_metric_status']=='FAIL','Fresh authoritative validation absent')
  need(fresh_objects[2]['result_sha256']=={'primary':sha(fresh[0]),'independent':sha(fresh[1])},'Fresh validator result hashes differ')
  audit['fresh_replay_comparison']=[{**compare(original,new,kind),'original_raw_sha256':sha(paths[kind]),'fresh_raw_sha256':sha(path)} for original,new,kind,path in zip((p,i,v),fresh_objects,('primary','independent','validation'),fresh)]
  if any(not row['science_equal'] or not row['stable_frame_equal'] for row in audit['fresh_replay_comparison']):audit['status']='ADVISORY_REPLAY_DIFFERENCE'
  transition=next(r for r in audit['fresh_replay_comparison'] if r['kind']=='independent')['allowed_package_manifest_transition']
  if transition is not None:audit['fresh_package_manifest_authentication']=authenticate_fresh_package_manifest(args.fresh_package_manifest,root,args.registration_sha256,transition['fresh'])
  audit['fresh_replay_scope']='Frozen validator outputs and exact scientific/stable-frame projections audited. Fresh outer receipt origin/remote public verification remain separately supplied by root; this tool executes no producer.'
 output.mkdir()
 for key,path in paths.items():shutil.copyfile(path,output/(key.upper()+'_RAW.json'))
 if all(fresh):
  for key,path in zip(('PRIMARY','INDEPENDENT','VALIDATION'),fresh):shutil.copyfile(path,output/('FRESH_'+key+'_RAW.json'))
  if args.fresh_package_manifest:shutil.copyfile(args.fresh_package_manifest,output/'FRESH_PACKAGE_MANIFEST_RAW.json')
 (output/'RESULT_REVIEW.json').write_text(json.dumps(audit,indent=2,allow_nan=False)+'\n')
 with (output/'CASE_ATTRIBUTION_MARGINS.csv').open('w',newline='') as f:
  writer=csv.writer(f);writer.writerow(['source','setting','K','route','strict_magnitude_margin','continuum_fraction_margin','sampling_fraction_margin'])
  for row in audit['case_reviews']:
   for route in ('primary','independent'):
    m=row[route+'_attribution_margins'];writer.writerow([row['source'],row['setting'],row['K'],route,*[m[k] for k in ('strict_magnitude_margin','continuum_fraction_margin','sampling_fraction_margin')]])
 print(json.dumps({'status':audit['status'],'classification':audit['classification'],'cases':len(audit['case_reviews']),'accounting_errors':len(audit['accounting_audit_errors'])}))
 return 0 if audit['status']=='ADVISORY_AUDIT_CONSISTENT' else 1
if __name__=='__main__':raise SystemExit(main())
