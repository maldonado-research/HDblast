"""Seal an independent prospective-only review after all final evidence exists."""
from pathlib import Path
import hashlib
import json
import os
import stat

HERE=Path(__file__).resolve().parent
BASE=HERE.parent
REGISTRATION_SHA='33eded898f47a96142f85002b896edd4c294180653f50f2b38e83f22d5a75796'
PUBLIC_HELPER_SHA='f5269a494287a78defa6cf886c533d7cc081e5a9bbe3c7f4d8aa1537cb0c66b5'
def require(ok,message):
 if not ok:raise ValueError(message)
def pin(path):
 raw=path.read_bytes();return {'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}
def read(name):return json.loads((HERE/name).read_bytes())
def inventory(root):
 files={};directories=set()
 def error(exc):raise exc
 for parent,names,leaves in os.walk(root,followlinks=False,onerror=error):
  for name in names+leaves:
   path=Path(parent)/name;metadata=path.lstat();relative=str(path.relative_to(root))
   require(not stat.S_ISLNK(metadata.st_mode),'source link')
   if stat.S_ISDIR(metadata.st_mode):directories.add(relative)
   else:
    require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink==1,'unaliased regular source required')
    files[relative]=pin(path)
 return files,directories

def main():
 receipt_path=HERE/'FINAL_REVIEW_RECEIPT_005.json';manifest_path=HERE/'FINAL_REVIEW_MANIFEST_005.json'
 require(not receipt_path.exists() and not manifest_path.exists(),'review already sealed')
 closure=read('SOURCE_CLOSURE_005.json');pair=read('BOUNDED_PAIR_005.json');allnode=read('FINAL_005_INDEPENDENT_ALLNODE_READBACK.json')
 require(closure['registration_sha256']==REGISTRATION_SHA and pair['registration']['sha256']==REGISTRATION_SHA,'registration bindings')
 require(pair['status']=='PASS_INDEPENDENT_BOUNDED_PAIR_AND_SOURCE_BYTES' and not pair['extra_empty_source_directories'],'final bounded pair/closure')
 require(allnode['status']=='PASS_INDEPENDENT_MANUFACTURED_ALLNODE_READBACK' and allnode['nodes']==49152 and allnode['comparisons']==98304 and allnode['prefixes']==24 and allnode['decoded_real_slots']==688152,'all-node coverage')
 for name,digest in allnode['input_output_files'].items():require(pair['stable_scientific_files'][name]['sha256']==digest,'all-node input pin')
 for relative in ['implementation','implementation-frozen-005','endpoint-independent-review/frozen-source-005']:
  files,dirs=inventory(BASE/relative)
  require(files==closure['source_files']==pair['source_files'] and dirs==set(closure['source_directories']),'source changed after review:'+relative)
 require(pin(BASE/'manufactured-registration-005.json')['sha256']==REGISTRATION_SHA,'registration changed')
 require(pin(BASE/'root/public_endpoint_readback_go.py')['sha256']==PUBLIC_HELPER_SHA,'public helper changed')
 for name in ['directory-scanner-frozen-005-normal/RESULT.json','directory-scanner-frozen-005-optimized/RESULT.json',
              'authenticated-bootstrap-frozen-005/RESULT.json','authenticated-bootstrap-frozen-005-optimized/RESULT.json',
              'unreadable-directory-frozen-005/RESULT.json','public-readback-controls-005/RESULT.json',
              'EXACT_CONTROLS_NORMAL.json','EXACT_CONTROLS_OPTIMIZED.json','SOURCE_ERROR_ROUNDING_CONTROL.json',
              'CUSTODIAN_005_EVIDENCE_BYTE_CHECK.json','producer-controls-frozen-005/BYTE_VERIFICATION.json']:
  require(read(name)['status'].startswith('PASS'),'nonpassing relevant control:'+name)
 require(read('public-readback-controls-005/RESULT.json')['helper_sha256']==PUBLIC_HELPER_SHA,'helper control binding')
 # Older exact numerical controls remain applicable only by identical source bytes.
 for mode in ['NORMAL','OPTIMIZED']:
  for name,digest in read('EXACT_CONTROLS_'+mode+'.json')['source_pins'].items():
   require(closure['source_files'][name]['sha256']==digest,'numeric source changed since oracle controls')
 notes=HERE/'REVIEW_NOTES.md';text=notes.read_text()
 text=text.replace('These notes are an intermediate review, not a final source acceptance or GO.\nThe root must receive a final receipt binding the complete frozen source tree,\nthe registration bytes, current mathematical and kernel control receipts,\nsuccessful complete bounded manufactured execution, and the independent final\nall-node readback. Any source edit after that snapshot requires renewed review.\nHistorical FAIL and UNRESOLVED outcomes remain unchanged.',
 '''Frozen005 is accepted for prospective public source publication only. Its
41 files and11 exact parent directories match registration
33eded898f47a96142f85002b896edd4c294180653f50f2b38e83f22d5a75796.
The directory scanner propagates errors and holds no-follow descriptors while
checking directory identities. The worker captures and authenticates all nine
consumed local modules before executing any, retains their dependency identities,
and keeps candidate directories outside sys.path throughout execution.

Final controls include11 independent scanner/error/race controls per mode and
12 independent entrypoint alias/cache/native controls per mode. The preserved
004 execute-only package reproduction now rejects before import. The producer's
23 semantic/GO checks and14 worker/captured-source checks pass per mode; the
launcher has23 methods per mode with54 final authoritative child outcomes.
All298 compact launcher evidence artifacts were independently rehashed.

Both full bounded manufactured005 runs pass. Their nine scientific files are
byte identical. An independent reader importing no production modules checks
49,152 nodes,98,304 endpoint comparisons, all24 finite aggregates,589,824 stored
normalized scalar components and688,128 sign flags (688,152 total scalar slots).
The fixed export radius is2^-96 for U and W. The public-readback helper passes
24 entirely offline protocol controls and is separately byte pinned.

This acceptance is not an actual source/decode GO. The root owns immutable
public byte readback and any later authorization. Any source edit requires a
new freeze and review. Frozen003 and004 remain rejected for prospective use;
their successful numerical runs are development evidence. Historical physical
FAIL and UNRESOLVED outcomes remain unchanged.''')
 notes.write_text(text)
 # Compact artifacts only: never copied candidate trees, caches or fixture GO files.
 selected={}
 for path in HERE.iterdir():
  if path.is_file() and path.name not in {receipt_path.name,manifest_path.name}:
   selected[path.name]=pin(path)
 allowed_prefixes=('bootstrap-cache-','bootstrap-python-alias-','authenticated-bootstrap-',
 'final-bootstrap-','unreadable-directory-','directory-scanner-')
 for directory in HERE.iterdir():
  if not directory.is_dir() or directory.is_symlink():continue
  if directory.name.startswith(allowed_prefixes) or directory.name in ('final-boundary-audit','producer-controls-frozen-005'):
   for path in directory.iterdir():
    if path.is_file() and not path.is_symlink():selected[str(path.relative_to(HERE))]=pin(path)
  if directory.name in ('public-readback-controls-final','public-readback-controls-005'):
   for path in directory.iterdir():
    if path.is_file():selected[str(path.relative_to(HERE))]=pin(path)
    elif path.is_dir() and (path/'CONTROL.log').is_file():selected[str((path/'CONTROL.log').relative_to(HERE))]=pin(path/'CONTROL.log')
 require(not any('__pycache__' in name or name.endswith(('.pyc','.pdf')) or 'MOCK_GO' in name for name in selected),'private/fixture artifact selected')
 external_names=['bounded-review/FINAL_CONTROL_MANIFEST_005.json','bounded-review/FINAL_CONTROL_REVIEW_RECEIPT_005.json']
 external={name:pin(BASE/name) for name in external_names}
 report={'schema_version':1,'status':'GO_FOR_PROSPECTIVE_PUBLIC_FREEZE','review_scope':'INDEPENDENT_FROZEN005_ENDPOINT_CANDIDATE_PREFLIGHT',
 'registration_sha256':REGISTRATION_SHA,'registered_file_pins':closure['source_files'],'registered_source_directories':closure['source_directories'],
 'independent_review_pass':True,'physical_execution_authorized_by_this_receipt':False,'source_go':False,'decode_go':False,'actual_public_GO_issued':False,
 'new_target_evaluations_before_freeze':{name:0 for name in ['physical_source_callbacks','retained_array_decodes','stored_endpoint_comparisons','observational_likelihood_evaluations']},
 'manufactured_full_pipeline':{'nodes':49152,'endpoint_comparisons':98304,'finite_aggregates':24,'scalar_slots':688152,'nine_scientific_files_byte_identical':True,
 'scientific_file_pins':pair['stable_scientific_files'],'bounded_runs':{mode:{key:value for key,value in run.items() if key!='file_pins'} for mode,run in pair['runs'].items()}},
 'controls':{'independent_exact_polynomial_and_stress_checks_each_mode':18,'independently_rerun_binary80_tests_each_mode':36,
 'independent_directory_scan_controls_each_mode':11,'independent_entrypoint_boundary_controls_each_mode':12,
 'independent_realistic_odd_denominator_source_error_rounding':True,'producer_semantic_GO_controls_each_mode':23,
 'producer_worker_bootstrap_captured_source_checks_each_mode':14,'launcher_methods_each_mode':23,
 'launcher_final_authoritative_child_outcomes':54,'launcher_preserved_initial_child_outcomes':54,'launcher_artifacts_rehashed':298,
 'offline_public_readback_protocol_checks':24},
 'public_readback_helper_sha256':PUBLIC_HELPER_SHA,'external_evidence_pins':external,'compact_review_evidence_pins':selected,
 'preserved_rejections':{'source003':'preauthentication Python package/stdlib alias execution','source004':'unreadable-directory omission from source inventory'},
 'scientific_limits':{'target':'Exact forced flow starting from exact represented saved t1=-9/2; endpoint comparisons at -4 and -7/2',
 'input_selection':'36 selected arrays; all76 capsule archive members authenticated before selected decoding; exact epsilon normalization once',
 'source_calls_for_future_authorized_run':128,'incoming_error_included':False,'historical_unsaved_trajectory':'UNRECOVERABLE_FROM_RETAINED_MODE_SNAPSHOTS',
 'prior_metric_calibration':'FAIL_UNCHANGED','continuous_momentum_pressure_contact_time_UV_certificate':'UNRESOLVED',
 'aggregate_output_bound':'POLL_AND_FINAL_ACCEPTANCE_CHECK_WITH_RECEIPT_RESERVE; not a filesystem quota'},
 'required_next_gate':'Root-owned immutable public byte readback and separate actual GO; this receipt authorizes prospective source publication only.'}
 receipt_path.write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
 selected[receipt_path.name]=pin(receipt_path)
 manifest={'schema_version':1,'status':'COMPLETE_COMPACT_INDEPENDENT_REVIEW_EVIDENCE_005','files':selected,
 'file_count':len(selected),'total_bytes':sum(p['bytes'] for p in selected.values()),'excluded':'Copied candidate/source trees, cache/native fixtures, large manufactured node streams, mock GO fixtures, private literature; retained locally and/or bound by hashes.'}
 manifest_path.write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
 print(json.dumps({'receipt':pin(receipt_path),'manifest':pin(manifest_path),'compact_file_count':len(selected),'total_bytes':manifest['total_bytes']},sort_keys=True))
if __name__=='__main__':main()
