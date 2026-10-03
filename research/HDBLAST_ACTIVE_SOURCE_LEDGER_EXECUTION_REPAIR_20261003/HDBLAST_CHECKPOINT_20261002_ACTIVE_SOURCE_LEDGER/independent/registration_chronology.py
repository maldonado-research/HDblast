"""Authenticated execution-schema repair; stdlib only, no numerical evaluation.

The historical zero applies only to the original public scientific freeze.
The repair follows one completed primary run and zero independent evaluations.
"""
from __future__ import annotations
import hashlib,json,re
from pathlib import Path

ORIGINAL_REGISTRATION_SHA256='8c772f5ef14970cee8e941bdbdb49a6b709e6aa4364f30afc9211f236e6f8791'
ORIGINAL_FREEZE_COMMIT='2e7dfddaa6b7aa2b43bcb12c7ce60a17eed39d3c'
ORIGINAL_RECEIPT_SHA256='b912e1918b5b77bca98b6b7217cc40800eb51a3a071913a1d8e80be3e5f3f2c0'
COUNTS_PATH='evidence/original_registration/PRE_REPAIR_COUNTS.json'
HELPER_PATH='independent/registration_chronology.py'
CORE_PATH='independent/run_independent.py'
OLD_MARKER=b" require(reg.get('physical_diagnostic_evaluations_before_freeze')==0,'Prospective evaluation marker mismatch')\n"
NEW_MARKER=b' verify_registration_chronology(root,reg)\n'
NEW_IMPORT=b'from registration_chronology import verify_registration_chronology\n'
IMPORT_ANCHOR=b'from engine import prepare,physical_sources,trajectory,cut_prefix\n'
COUNT_FIELDS={'combined_classification','evidence_paths','evidence_sha256','independent_initial_guard_exit','original_frozen_scientific_choices_unchanged','original_public_freeze_commit','original_registration_sha256','prior_physical_execution_counts_before_repair','prior_primary_completed_cases','registration_kind','schema_version'}
EVIDENCE_PATHS={'executions':'evidence/initial_attempt/EXECUTION.json','failed_replay':'evidence/initial_attempt/REPLAY.json','original_freeze_receipt':'evidence/original_registration/FREEZE_RECEIPT.json','original_registration':'evidence/original_registration/FULL_REGISTRATION.json','primary_output':'outputs/initial_attempt/primary/diagnostic.json'}

def require(condition,message):
 if not condition:raise RuntimeError(message)
def sha_bytes(raw):return hashlib.sha256(raw).hexdigest()
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def read(path):
 def unique(pairs):
  result={}
  for key,value in pairs:
   require(key not in result,'Duplicate repair JSON key');result[key]=value
  return result
 return json.loads(path.read_text(),object_pairs_hook=unique)
def hexadecimal(value,n=64):return type(value) is str and re.fullmatch('[0-9a-f]{%d}'%n,value) is not None
def safe_file(root,name):
 require(type(name) is str and name and not Path(name).is_absolute() and '..' not in Path(name).parts and '\\' not in name,'Unsafe repair evidence path')
 path=root/name
 require(path.is_file() and not path.is_symlink(),'Missing or symlink repair evidence')
 require(all(not p.is_symlink() for p in (path,*path.parents)),'Repair path symlink ancestry forbidden')
 require(path.resolve().is_relative_to(root.resolve()),'Repair evidence escaped checkpoint')
 return path
def exact_int(value,expected):return type(value) is int and value==expected

def patch_core(original):
 """The sole permissible core change is an import and the registration guard."""
 require(original.count(OLD_MARKER)==1 and original.count(IMPORT_ANCHOR)==1,'Original core patch anchors differ')
 require(NEW_IMPORT not in original and NEW_MARKER not in original,'Original core already repaired')
 return original.replace(IMPORT_ANCHOR,IMPORT_ANCHOR+NEW_IMPORT).replace(OLD_MARKER,NEW_MARKER)

def original_core_from_patch(patched):
 require(patched.count(NEW_IMPORT)==1 and patched.count(NEW_MARKER)==1 and OLD_MARKER not in patched,'Repair core patch anchors differ')
 original=patched.replace(NEW_IMPORT,b'').replace(NEW_MARKER,OLD_MARKER)
 require(patch_core(original)==patched,'Noncanonical repair core transform')
 return original

def verify_registration_chronology(root,registration):
 root=Path(root).absolute()
 require(type(registration) is dict,'Registration object required')
 kind=registration.get('registration_kind')
 if kind is None:
  require(exact_int(registration.get('physical_diagnostic_evaluations_before_freeze'),0),'Legacy prospective marker must be exact integer zero')
  require('pre_repair_counts_sha256' not in registration,'Repair proof cannot masquerade as legacy registration')
  return
 require(kind=='EXECUTION_SCHEMA_REPAIR','Unsupported registration kind')
 require('physical_diagnostic_evaluations_before_freeze' not in registration and 'new_followup_physical_diagnostic_evaluations_before_freeze' not in registration,'Repair must not claim a new before-freeze zero')
 files=registration.get('files');require(type(files) is dict and files,'Repair registered inventory missing')
 require(files.get(HELPER_PATH)==sha(safe_file(root,HELPER_PATH)),'Repair helper source not authenticated')
 countpath=safe_file(root,COUNTS_PATH)
 require(hexadecimal(registration.get('pre_repair_counts_sha256')) and sha(countpath)==registration['pre_repair_counts_sha256']==files.get(COUNTS_PATH),'Repair count receipt hash differs')
 counts=read(countpath)
 require(type(counts) is dict and set(counts)==COUNT_FIELDS,'Repair count receipt schema differs')
 require(exact_int(counts['schema_version'],1) and counts['registration_kind']=='EXECUTION_SCHEMA_REPAIR','Unsupported repair count schema')
 require(counts['combined_classification']=='UNCOMPUTED','Initial combined outcome must remain uncomputed')
 require(counts['original_frozen_scientific_choices_unchanged'] is True,'Repair choices changed')
 prior=counts['prior_physical_execution_counts_before_repair']
 require(type(prior) is dict and set(prior)=={'primary','independent'} and exact_int(prior['primary'],1) and exact_int(prior['independent'],0),'Repair chronology must record primary1/independent0')
 require(exact_int(counts['prior_primary_completed_cases'],12) and exact_int(counts['independent_initial_guard_exit'],1),'Original attempted-route facts differ')
 require(counts['original_public_freeze_commit']==ORIGINAL_FREEZE_COMMIT and counts['original_registration_sha256']==ORIGINAL_REGISTRATION_SHA256,'Original scientific freeze identity differs')
 require(counts['evidence_paths']==EVIDENCE_PATHS and set(counts['evidence_sha256'])==set(EVIDENCE_PATHS),'Original evidence inventory differs')
 for label,name in EVIDENCE_PATHS.items():
  pin=counts['evidence_sha256'][label]
  require(hexadecimal(pin) and files.get(name)==pin and sha(safe_file(root,name))==pin,'Original evidence bytes differ: '+label)
 require(counts['evidence_sha256']['original_registration']==ORIGINAL_REGISTRATION_SHA256 and counts['evidence_sha256']['original_freeze_receipt']==ORIGINAL_RECEIPT_SHA256,'Original registration/receipt anchors differ')
 original=read(safe_file(root,EVIDENCE_PATHS['original_registration']))
 receipt=read(safe_file(root,EVIDENCE_PATHS['original_freeze_receipt']))
 require(exact_int(original.get('schema_version'),1) and exact_int(original.get('new_followup_physical_diagnostic_evaluations_before_freeze'),0),'Historical first-registration zero not authenticated')
 require(original.get('scientific_outcome_before_freeze')=='NEW_FOLLOWUP_DIAGNOSTIC_UNCOMPUTED','Historical prefreeze outcome differs')
 original_files=original.get('files');require(type(original_files) is dict and original_files,'Original file inventory absent')
 require(receipt.get('public_freeze_commit')==ORIGINAL_FREEZE_COMMIT and receipt.get('registration_sha256')==ORIGINAL_REGISTRATION_SHA256 and receipt.get('remote_complete_tree_sha')==ORIGINAL_FREEZE_COMMIT,'Original public receipt identity differs')
 require(receipt.get('all_registered_public_blobs_independently_downloaded') is True and exact_int(receipt.get('new_followup_physical_evaluations_before_verification'),0),'Original public verification did not precede source evaluation')
 remote=receipt.get('remote_registered_files');expected={**original_files,'FULL_REGISTRATION.json':ORIGINAL_REGISTRATION_SHA256}
 require(type(remote) is dict and set(remote)==set(expected) and exact_int(receipt.get('remote_files_verified'),len(expected)),'Original complete public proof inventory differs')
 for name,pin in expected.items():
  proof=remote[name]
  require(type(proof) is dict and proof.get('sha256')==pin and proof.get('independent_public_blob_download_sha256_verified') is True,'Original public blob proof differs: '+name)
 for name,pin in original_files.items():
  if name==CORE_PATH:continue
  require(files.get(name)==pin,'Previously frozen non-guard bytes changed: '+name)
 current_core=safe_file(root,CORE_PATH).read_bytes()
 require(files.get(CORE_PATH)==sha_bytes(current_core) and sha_bytes(original_core_from_patch(current_core))==original_files.get(CORE_PATH),'Numerical core changed beyond the exact guard patch')
 for name in ('experiment_sha256','input_manifest_sha256','frozen_configuration'):
  require(registration.get(name)==original.get(name),'Previously frozen scientific registration differs: '+name)
 replay=read(safe_file(root,EVIDENCE_PATHS['failed_replay']))
 require(replay.get('status')=='FAIL_REPLAY' and replay.get('classification')=='UNCOMPUTED' and replay.get('plan_only') is False and exact_int(replay.get('physical_routes_completed'),1),'Initial failed combined replay facts differ')
 executions=read(safe_file(root,EVIDENCE_PATHS['executions']))
 require(type(executions) is list and len(executions)==2,'Initial two route attempts not preserved')
 for row,name,status,code in zip(executions,('primary','independent'),('PASS_EXECUTION','FAIL_EXECUTION'),(0,1)):
  require(type(row) is dict and row.get('name')==name and row.get('physical_route') is True and row.get('status')==status and exact_int(row.get('exit_code'),code),'Initial execution chronology differs')
