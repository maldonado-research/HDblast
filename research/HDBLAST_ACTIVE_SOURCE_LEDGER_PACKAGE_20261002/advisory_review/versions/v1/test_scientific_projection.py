"""Invented output fixtures only: ensure omitted metadata cannot hide science/pin changes."""
import copy,hashlib,importlib.util,json,sys
from pathlib import Path
from scientific_projection import policy,compare,stable_diagnostic_digest

AUTOMATION_HELPER=Path(__file__).parent.parent/'automation'/'scientific_projection.py'
AUTOMATION_PIN='2fe526a120544766e52a4dbad0886c32aafcc2b8e72bc48e006195828d070392'


def need(ok,message):
 if not ok:raise RuntimeError(message)


def fixture(kind):
 p=policy();r={k:None for k in p[kind+'_scientific_allowlist']+p[kind+'_excluded_top_keys']}
 r['resources']={('seconds' if kind=='primary' else 'elapsed_seconds'):123.4,'peak_rss_kib':123456,'scope':'fabricated full route','limit_seconds':900,'limit_peak_rss_kib':262144} if kind!='validation' else None
 if kind=='validation':r.pop('resources')
 if kind=='primary':
  r.update(producer_sha256='1'*64,registration_sha256='2'*64,freeze_commit='3'*40,input_manifest_sha256='4'*64,source_unchanged_during_execution=True,
           schema_version=1,fixed_cases=12,old_metric_status='FAIL',status='COMPLETED_SAVED_DATA_DIAGNOSTIC',records=[{'invented_numeric_value':'0.00000017','unicode_scope':'invented μ'}])
 elif kind=='independent':
  r.update(producer_sha256='1'*64,input_manifest_sha256='4'*64,schema_version=1,fixed_cases=12,old_metric_status='FAIL',status='COMPLETED_DIAGNOSTIC_WITHOUT_RECLASSIFYING_OLD_FAIL',resource_limits={'wall_seconds':900,'peak_rss_kib':262144},rows=[{'invented_numeric_value':'0.00000017','unicode_scope':'invented μ'}],
           provenance={'freeze_commit':'3'*40,'registration_sha256':'2'*64,'registered_files_checked':1,'input_manifest_sha256':'4'*64,'post_run_frozen_verification':'PASS',
                       'full_frozen_verification':{'public_freeze_commit':'3'*40,'registration_sha256':'2'*64,'freeze_receipt_sha256':'a'*64,'package_manifest_sha256':None,'frozen_files':{'inputs/INPUT_MANIFEST.json':'4'*64},
                                                  'input_capsule_verification':{'capsules_verified':4,'members_verified':76,'array_payload_values_interpreted':False,'input_manifest_sha256':'4'*64}}})
 else:r.update(registration_sha256='2'*64,public_freeze_commit='3'*40,result_sha256={'primary':'5'*64,'independent':'6'*64},executions=[{'elapsed_seconds':99,'peak_rss_kib':111111}],python_optimization=0)
 return r


def rejected(original,fresh,kind):
 try:compare(original,fresh,kind)
 except ValueError:return True
 return False


def main():
 checks=[]
 for kind in ('primary','independent','validation'):
  a=fixture(kind);b=copy.deepcopy(a)
  need(compare(a,b,kind)['science_equal'],'Identical fabricated science differs');checks.append(kind+'_identical')
  if kind!='validation':b['resources']['seconds' if kind=='primary' else 'elapsed_seconds']+=100;b['resources']['peak_rss_kib']+=111
  else:b['executions'][0]['elapsed_seconds']+=100;b['python_optimization']=1;b['result_sha256']['primary']='7'*64
  need(compare(a,b,kind)['science_equal'] and compare(a,b,kind)['stable_frame_equal'],'Runtime/raw hash metadata was treated as scientific difference');checks.append(kind+'_metadata_ignored')
  if kind!='validation':
   c=copy.deepcopy(a);c['resources']['scope']='changed scope must remain visible'
   need(compare(a,c,kind)['stable_frame_equal'] is False,'Changed resource scope hidden');checks.append(kind+'_resource_scope_mutation_detected')
  b=copy.deepcopy(a)
  b['records' if kind=='primary' else 'rows' if kind=='independent' else 'interpretation']=[{'invented_numeric_value':'0.000000170000000000000001'}]
  need(compare(a,b,kind)['science_equal'] is False,'Small serialized science mutation hidden');checks.append(kind+'_science_mutation_detected')
  b=copy.deepcopy(a);b['unexpected_new_key']='must never disappear'
  need(rejected(a,b,kind),'Unexpected top-level key silently discarded');checks.append(kind+'_unknown_key_rejected')
  b=copy.deepcopy(a)
  if kind=='primary':b['input_manifest_sha256']='8'*64
  elif kind=='independent':b['provenance']['full_frozen_verification']['frozen_files']['inputs/INPUT_MANIFEST.json']='8'*64
  else:b['registration_sha256']='8'*64
  need(rejected(a,b,kind),'Changed frozen pin hidden by projection');checks.append(kind+'_pin_mutation_rejected')
 a=fixture('independent');b=copy.deepcopy(a);b['provenance']['full_frozen_verification']['package_manifest_sha256']='9'*64
 c=compare(a,b,'independent');need(c['science_equal'] and c['allowed_package_manifest_transition'] is not None,'Declared None-to-package transition rejected');checks.append('explicit_package_transition')
 a['provenance']['full_frozen_verification']['package_manifest_sha256']='a'*64
 need(rejected(a,b,'independent'),'Two different nonempty manifest pins accepted');checks.append('different_nonempty_package_pins_rejected')
 a=fixture('independent');b=copy.deepcopy(a);b['producer_sha256']='a'*64
 need(rejected(a,b,'independent'),'Changed producer pin hidden');checks.append('producer_mutation_rejected')
 need(hashlib.sha256(AUTOMATION_HELPER.read_bytes()).hexdigest()==AUTOMATION_PIN,'Automation projector changed; review compatibility explicitly')
 spec=importlib.util.spec_from_file_location('automation_projection_review_fixture',AUTOMATION_HELPER)
 module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
 for kind in ('primary','independent'):
  a=fixture(kind);b=copy.deepcopy(a)
  b['resources']['seconds' if kind=='primary' else 'elapsed_seconds']+=100
  b['resources']['peak_rss_kib']+=11
  if kind=='independent':b['provenance']['full_frozen_verification']['package_manifest_sha256']='9'*64
  for name,value in (('original',a),('fresh',b)):
   need(stable_diagnostic_digest(value,kind)==module.project(value,kind)[0],'Full-frame digest disagrees with independent automation policy')
   checks.append(kind+'_'+name+'_automation_digest_compatible')
 result={'status':'PASS_SCIENTIFIC_PROJECTION_FABRICATED_MUTATIONS','physical_evaluations':0,'actual_result_paths_read':False,'checks':checks,
         'python_optimization':sys.flags.optimize,'automation_projector_sha256':AUTOMATION_PIN,'allowlist_sha256':hashlib.sha256(Path(__file__).with_name('SCIENTIFIC_PROJECTION_ALLOWLIST.json').read_bytes()).hexdigest()}
 name='PROJECTION_TESTS_OPTIMIZED.json' if sys.flags.optimize else 'PROJECTION_TESTS_NORMAL.json'
 Path(__file__).with_name(name).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
