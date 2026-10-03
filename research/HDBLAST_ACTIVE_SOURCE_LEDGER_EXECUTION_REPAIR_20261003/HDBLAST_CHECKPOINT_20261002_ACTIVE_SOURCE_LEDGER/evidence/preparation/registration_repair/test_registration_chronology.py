"""Fabricated metadata and opaque bytes only; no scientific result or array reads."""
import ast,copy,hashlib,json,sys,tempfile
from contextlib import contextmanager
from pathlib import Path
import registration_chronology as guard

def need(condition,message):
 if not condition:raise RuntimeError(message)
def write(path,value):
 path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(value,sort_keys=True,indent=2)+'\n')
def put(root,name,raw):
 p=root/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(raw);return guard.sha(p)

@contextmanager
def fixture(historical_marker=0,remote_verified=True):
 saved=(guard.ORIGINAL_REGISTRATION_SHA256,guard.ORIGINAL_FREEZE_COMMIT,guard.ORIGINAL_RECEIPT_SHA256)
 with tempfile.TemporaryDirectory(dir='/tmp') as temporary:
  root=Path(temporary)
  fake_core=guard.IMPORT_ANCHOR+b'\ndef verify(reg):\n'+guard.OLD_MARKER+b' return "invented numerical tail, never executed"\n'
  original_files={guard.CORE_PATH:guard.sha_bytes(fake_core),'EXPERIMENT.json':put(root,'EXPERIMENT.json',b'{"invented_choices":true}'),'code/diagnostic_primary.py':put(root,'code/diagnostic_primary.py',b'# invented retained method bytes\n')}
  original={'schema_version':1,'new_followup_physical_diagnostic_evaluations_before_freeze':historical_marker,'scientific_outcome_before_freeze':'NEW_FOLLOWUP_DIAGNOSTIC_UNCOMPUTED','files':original_files,'experiment_sha256':original_files['EXPERIMENT.json'],'input_manifest_sha256':'a'*64,'frozen_configuration':{'invented_choices':True}}
  write(root/guard.EVIDENCE_PATHS['original_registration'],original)
  guard.ORIGINAL_REGISTRATION_SHA256=guard.sha(root/guard.EVIDENCE_PATHS['original_registration']);guard.ORIGINAL_FREEZE_COMMIT='b'*40
  remote={name:{'sha256':pin,'independent_public_blob_download_sha256_verified':remote_verified} for name,pin in {**original_files,'FULL_REGISTRATION.json':guard.ORIGINAL_REGISTRATION_SHA256}.items()}
  receipt={'public_freeze_commit':guard.ORIGINAL_FREEZE_COMMIT,'registration_sha256':guard.ORIGINAL_REGISTRATION_SHA256,'remote_complete_tree_sha':guard.ORIGINAL_FREEZE_COMMIT,'all_registered_public_blobs_independently_downloaded':True,'new_followup_physical_evaluations_before_verification':0,'remote_files_verified':len(remote),'remote_registered_files':remote}
  write(root/guard.EVIDENCE_PATHS['original_freeze_receipt'],receipt);guard.ORIGINAL_RECEIPT_SHA256=guard.sha(root/guard.EVIDENCE_PATHS['original_freeze_receipt'])
  write(root/guard.EVIDENCE_PATHS['executions'],[{'name':'primary','physical_route':True,'status':'PASS_EXECUTION','exit_code':0},{'name':'independent','physical_route':True,'status':'FAIL_EXECUTION','exit_code':1}])
  write(root/guard.EVIDENCE_PATHS['failed_replay'],{'status':'FAIL_REPLAY','classification':'UNCOMPUTED','plan_only':False,'physical_routes_completed':1})
  put(root,guard.EVIDENCE_PATHS['primary_output'],b'{"invented_result":"opaque bytes, not interpreted"}\n')
  counts={'schema_version':1,'registration_kind':'EXECUTION_SCHEMA_REPAIR','combined_classification':'UNCOMPUTED','original_frozen_scientific_choices_unchanged':True,'prior_physical_execution_counts_before_repair':{'primary':1,'independent':0},'prior_primary_completed_cases':12,'independent_initial_guard_exit':1,'original_public_freeze_commit':guard.ORIGINAL_FREEZE_COMMIT,'original_registration_sha256':guard.ORIGINAL_REGISTRATION_SHA256,'evidence_paths':guard.EVIDENCE_PATHS,'evidence_sha256':{label:guard.sha(root/name) for label,name in guard.EVIDENCE_PATHS.items()}}
  write(root/guard.COUNTS_PATH,counts)
  files={**original_files,guard.CORE_PATH:put(root,guard.CORE_PATH,guard.patch_core(fake_core)),guard.HELPER_PATH:put(root,guard.HELPER_PATH,Path(guard.__file__).read_bytes()),guard.COUNTS_PATH:guard.sha(root/guard.COUNTS_PATH),**{name:guard.sha(root/name) for name in guard.EVIDENCE_PATHS.values()}}
  registration={'schema_version':1,'registration_kind':'EXECUTION_SCHEMA_REPAIR','pre_repair_counts_sha256':files[guard.COUNTS_PATH],'files':files,**{k:original[k] for k in ('experiment_sha256','input_manifest_sha256','frozen_configuration')}}
  try:yield root,registration
  finally:guard.ORIGINAL_REGISTRATION_SHA256,guard.ORIGINAL_FREEZE_COMMIT,guard.ORIGINAL_RECEIPT_SHA256=saved

def update_counts(root,registration,mutate):
 p=root/guard.COUNTS_PATH;a=guard.read(p);mutate(a);write(p,a)
 registration['pre_repair_counts_sha256']=registration['files'][guard.COUNTS_PATH]=guard.sha(p)
def rejected(root,registration):
 try:guard.verify_registration_chronology(root,registration)
 except (RuntimeError,KeyError,TypeError):return True
 return False

def main():
 checks=[]
 guard.verify_registration_chronology(Path('/tmp'),{'physical_diagnostic_evaluations_before_freeze':0});checks.append('legacy_exact_integer_zero')
 for name,value in [('boolean_false',False),('float_zero',0.0),('missing',None),('nonzero',1)]:
  need(rejected(Path('/tmp'),{'physical_diagnostic_evaluations_before_freeze':value}),'Bad legacy marker accepted');checks.append('legacy_'+name+'_rejected')
 need(rejected(Path('/tmp'),{'new_followup_physical_diagnostic_evaluations_before_freeze':0}),'Missing marker fallback accepted');checks.append('missing_key_fallback_rejected')
 need(rejected(Path('/tmp'),{'physical_diagnostic_evaluations_before_freeze':0,'pre_repair_counts_sha256':'a'*64}),'Repair masquerade accepted');checks.append('legacy_repair_masquerade_rejected')
 with fixture() as (root,registration):
  guard.verify_registration_chronology(root,registration);checks.append('authenticated_repair_chronology')
 mutations=[('primary_zero',lambda a:a['prior_physical_execution_counts_before_repair'].update(primary=0)),('primary_boolean',lambda a:a['prior_physical_execution_counts_before_repair'].update(primary=True)),('independent_boolean',lambda a:a['prior_physical_execution_counts_before_repair'].update(independent=False)),('independent_one',lambda a:a['prior_physical_execution_counts_before_repair'].update(independent=1)),('changed_cases',lambda a:a.update(prior_primary_completed_cases=11)),('combined_claim',lambda a:a.update(combined_classification='LEDGER_ERROR_DEMONSTRATED')),('choices_flag_integer',lambda a:a.update(original_frozen_scientific_choices_unchanged=1)),('original_freeze_change',lambda a:a.update(original_public_freeze_commit='c'*40)),('extra_counts_key',lambda a:a.update(unregistered_fact=True))]
 for name,mutation in mutations:
  with fixture() as (root,registration):
   update_counts(root,registration,mutation);need(rejected(root,registration),'Bad repair metadata accepted: '+name);checks.append(name+'_rejected')
 for name in ('bad_count_hash','helper_change','retained_method_change','numeric_core_change','frozen_choices_change','evidence_corruption','false_zero_marker'):
  with fixture() as (root,registration):
   if name=='bad_count_hash':registration['pre_repair_counts_sha256']='c'*64
   elif name=='helper_change':put(root,guard.HELPER_PATH,b'# unregistered helper mutation\n')
   elif name=='retained_method_change':registration['files']['code/diagnostic_primary.py']=put(root,'code/diagnostic_primary.py',b'# altered retained method\n')
   elif name=='numeric_core_change':
    p=root/guard.CORE_PATH;p.write_bytes(p.read_bytes()+b'# non-guard numerical core alteration\n');registration['files'][guard.CORE_PATH]=guard.sha(p)
   elif name=='frozen_choices_change':registration['frozen_configuration']={'invented_choices':False}
   elif name=='evidence_corruption':put(root,guard.EVIDENCE_PATHS['primary_output'],b'altered opaque output bytes')
   else:registration['new_followup_physical_diagnostic_evaluations_before_freeze']=0
   need(rejected(root,registration),'Bad repair accepted: '+name);checks.append(name+'_rejected')
 with fixture(historical_marker=False) as (root,registration):need(rejected(root,registration),'Boolean historical zero accepted');checks.append('historical_boolean_zero_rejected')
 with fixture(remote_verified=1) as (root,registration):need(rejected(root,registration),'Integer public proof flag accepted');checks.append('public_blob_boolean_type_rejected')
 original_path=Path(__file__).with_name('ORIGINAL_CORE.py.txt')
 original=original_path.read_bytes();patched=Path(__file__).with_name('run_independent.py').read_bytes()
 need(guard.patch_core(original)==patched and guard.original_core_from_patch(patched)==original,'Staged patch altered core methods')
 need(ast.dump(ast.parse(original))==ast.dump(ast.parse(guard.original_core_from_patch(patched))),'Core AST changed outside exact guard patch');checks.append('actual_core_exact_byte_transform_and_ast_invariance')
 receipt={'schema_version':1,'status':'PASS_FABRICATED_REGISTRATION_CHRONOLOGY_GUARDS','physical_evaluations':0,'actual_input_arrays_read':False,'actual_scientific_result_values_interpreted':False,'fixture_scope':'Invented metadata and opaque temporary payloads; actual source bytes examined only to prove deterministic guard-only transform.','python_optimization':sys.flags.optimize,'checks':checks,'helper_sha256':guard.sha(Path(guard.__file__)),'staged_core_sha256':guard.sha(Path(__file__).with_name('run_independent.py')),'original_core_sha256':guard.sha(original_path)}
 path=Path(__file__).with_name('GUARDS_OPTIMIZED.json' if sys.flags.optimize else 'GUARDS_NORMAL.json');write(path,receipt);print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
