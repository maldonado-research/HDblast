#!/usr/bin/env python3
"""Portable real-subprocess checks for the independent Cauchy response audit."""
import argparse,copy,hashlib,json,shutil,subprocess,sys
from fractions import Fraction as Q
from pathlib import Path

def need(ok,why):
 if not ok:raise RuntimeError(why)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,d):p.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--source-root',required=True,type=Path)
 for role in ('same-response','absolute','velocity','independent-velocity'):ap.add_argument('--'+role+'-root',required=True,type=Path)
 ap.add_argument('--producer-source',required=True,type=Path);ap.add_argument('--producer-freeze',required=True,type=Path)
 ap.add_argument('--producer-report',required=True,type=Path);ap.add_argument('--producer-report-sha256',required=True)
 ap.add_argument('--output-dir',required=True,type=Path);a=ap.parse_args()
 need(not a.output_dir.exists(),'fresh Cauchy controls directory');a.output_dir.mkdir(parents=True)
 checker=Path(__file__).with_name('check_cauchy_response.py');checker_before=sha(checker)
 roots={role:getattr(a,role+'_root') for role in ('same_response','absolute','velocity','independent_velocity')}
 original_paths=dict(checker=checker,source=a.source_root,roots=roots,producer_source=a.producer_source,producer_freeze=a.producer_freeze,report=a.producer_report)
 def run(paths,out,opt=False,reportsha=None):
  cmd=[sys.executable]+(['-O'] if opt else [])+[str(paths['checker']),'--source-root',str(paths['source'])]
  for role,root in paths['roots'].items():cmd+=['--'+role.replace('_','-')+'-root',str(root)]
  cmd+=['--producer-source',str(paths['producer_source']),'--producer-freeze',str(paths['producer_freeze']),
   '--producer-report',str(paths['report']),'--producer-report-sha256',reportsha or sha(paths['report']),'--output',str(out)]
  return subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
 outputs=[]
 for opt in (False,True):
  out=a.output_dir/('INDEPENDENT_RESPONSE_OPTIMIZED.json' if opt else 'INDEPENDENT_RESPONSE_NORMAL.json')
  p=run(original_paths,out,opt,a.producer_report_sha256);need(p.returncode==0,p.stderr);outputs.append(out)
 need(outputs[0].read_bytes()==outputs[1].read_bytes(),'normal optimized canonical equality')
 result=json.loads(outputs[0].read_text());freeze=json.loads(a.producer_freeze.read_text())
 stage=a.output_dir/'relocated';source=stage/'source';proofs=stage/'proofs';producer=stage/'producer';checks=stage/'checker'
 for d in (source,proofs,producer,checks):d.mkdir(parents=True)
 for rel in result['selected_source_sha256']:
  if rel.startswith('research/'):
   p=source/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(a.source_root/rel,p)
 staged_roots={role:proofs/role for role in roots}
 for role,pins in freeze['input_pins'].items():
  for rel in pins:
   p=staged_roots[role]/rel;p.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(roots[role]/rel,p)
 new_checker=checks/checker.name;shutil.copyfile(checker,new_checker)
 new_producer=producer/a.producer_source.name;shutil.copyfile(a.producer_source,new_producer)
 new_freeze=producer/'FREEZE.json';shutil.copyfile(a.producer_freeze,new_freeze)
 report=producer/'RESPONSE.json';shutil.copyfile(a.producer_report,report)
 paths=dict(checker=new_checker,source=source,roots=staged_roots,producer_source=new_producer,producer_freeze=new_freeze,report=report)
 for opt in (False,True):
  out=a.output_dir/('PORTABLE_OPTIMIZED.json' if opt else 'PORTABLE_NORMAL.json');p=run(paths,out,opt,a.producer_report_sha256)
  need(p.returncode==0,p.stderr);need(out.read_bytes()==outputs[0].read_bytes(),'relocated canonical equality')
 baseline=json.loads((a.source_root/'research/HDBLAST_EXACT_PACKET_INITIAL_20261031/reference/INITIAL.json').read_text())['campaigns'][0]
 original_report=report.read_bytes()
 def old_abs(d):
  for n in ('v','v_r'):d['candidate_initial_absolute_H_caps'][n]=baseline['candidate_initial_absolute_H_caps'][n]
 def old_errors(d):
  for n in ('v0','v0_r'):d['initial_error_H_caps'][n]=baseline['initial_error_H_caps'][n]
 def altered_witness(d):
  row=d['formula_and_stale_velocity_corruption_witnesses']['omit_original_Fx'];row['wrong_minus_correct_rr']=str(Q(row['wrong_minus_correct_rr'])+1)
 cases=[('stale_velocity_absolute',old_abs,'fresh analytic absolute velocity inputs required'),
 ('stale_velocity_errors',old_errors,'fresh paired analytic velocity errors required'),
 ('false_absolute_cap',lambda d:d['candidate_initial_absolute_H_caps'].__setitem__('u_rr','0'),'fresh analytic absolute velocity inputs required'),
 ('false_error_cap',lambda d:d['initial_error_H_caps'].__setitem__('v0_r','0'),'fresh paired analytic velocity errors required'),
 ('false_original_Fx',lambda d:d['constants'].__setitem__('Fx','0'),'every Cauchy operator/source/commutator constant'),
 ('false_ledger',lambda d:d['ledger']['second_coordinate_error_H_caps'].__setitem__('rr',str(Q(d['ledger']['second_coordinate_error_H_caps']['rr'])+1)),'every Cauchy ledger node independently reproduced'),
 ('false_scope_flag',lambda d:d.__setitem__('no_global_true_error_factor_claim',False),'advertised Cauchy boolean flag: no_global_true_error_factor_claim'),
 ('false_Gram',lambda d:d['Gram'].__setitem__(0,'2'),'advertised response domain/full23 norm/baseline count'),
 ('false_selected_map',lambda d:d['selected_input_sha256_after'].pop('FREEZE.json'),'every advertised selected source pin'),
 ('false_comparison',lambda d:d['node_comparisons']['complete_C1_error_X_cap'].__setitem__('old_over_new','0'),'every advertised baseline-to-new node comparison'),
 ('false_witness',altered_witness,'every advertised formula/stale-velocity corruption witness'),
 ('false_field_role',lambda d:d['field_roles'].__setitem__('v/v_r','old P1 velocity'),'every advertised initial field role'),
 ('false_status',lambda d:d.__setitem__('status','PASS_OTHER'),'advertised response schema and status')]
 records=[]
 def reject(name,marker,opt,altered_paths=paths):
  out=a.output_dir/(name+('_O' if opt else '')+'_unexpected.json');p=run(altered_paths,out,opt)
  need(p.returncode!=0 and marker in p.stderr and not out.exists(),'effective Cauchy rejection:'+name+':'+p.stderr)
  records.append(dict(case=name,optimized=opt,expected_rejection=marker,nonzero_exit=True,no_success_artifact=True))
 for name,change,marker in cases:
  d=json.loads(original_report);change(d);write(report,d)
  for opt in (False,True):reject(name,marker,opt)
  report.write_bytes(original_report)
 for name,path,marker in [
 ('native_velocity_table_disk',staged_roots['velocity']/'INITIAL_NORMAL.json','selected source hash: velocity/INITIAL_NORMAL.json'),
 ('Cauchy_recipe_disk',new_producer,'selected source hash: producer/response_cauchy.py'),
 ('original_Fx_disk',source/'research/HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/FX.json','selected source hash: research/HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/reference/FX.json')]:
  raw=path.read_bytes();path.write_bytes(raw+b'\n')
  for opt in (False,True):reject(name,marker,opt)
  path.write_bytes(raw)
 alias=producer/'alias.json';alias.symlink_to(report.name);altered=dict(paths,report=alias)
 for opt in (False,True):reject('report_symlink','regular symlink-free input: alias.json',opt,altered)
 need(sha(checker)==checker_before,'independent checker unchanged')
 receipt=dict(status='PASS_INDEPENDENT_CAUCHY_RESPONSE_SUBPROCESS_CONTROLS',checker_sha256=checker_before,
  producer_report_sha256=a.producer_report_sha256,canonical_report_sha256=sha(outputs[0]),
  normal_optimized_byte_identical=True,relocated_normal_optimized_byte_identical=True,
  rehashed_report_corruption_cases=len(cases),rejected_subprocesses=len(records),controls=records,
  every_rejection_nonzero_without_success_artifact=True)
 write(a.output_dir/'CONTROLS.json',receipt)
 print(json.dumps({'status':receipt['status'],'rejected_subprocesses':len(records),'canonical_report_sha256':receipt['canonical_report_sha256']}))

if __name__=='__main__':main()
