"""Manufactured byte-closure controls; no fake public GO or scientific calls."""
from pathlib import Path
import hashlib,importlib.util,json,os,shutil,socket,sys,tempfile,types,unittest
from unittest import mock
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
SRC=HERE/'replay_checkpoint.py'
def pin(raw):return {'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}
def encode(obj):return (json.dumps(obj,sort_keys=True,indent=2)+'\n').encode()
class Controls(unittest.TestCase):
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory(prefix='portable-fabricated-');self.root=Path(self.tmp.name)/'package';self.root.mkdir()
  shutil.copyfile(SRC,self.root/'replay_checkpoint.py')
  spec=importlib.util.spec_from_file_location('manufactured_portable_helper',self.root/'replay_checkpoint.py');self.m=importlib.util.module_from_spec(spec);sys.modules[spec.name]=self.m;spec.loader.exec_module(self.m)
  self.payload={}
  def add(name,raw):self.payload[name]=raw
  add('PUBLIC_GO.json',b'{"status":"NOT_AN_AUTHORIZATION_RECEIPT"}\n')
  self.contract={'resources':self.m.LIMITS};add(self.m.PREFIX+'/REGISTRATION_CONTRACT.json',encode(self.contract))
  add(self.m.PREFIX+'/FULL_REGISTRATION.json',b'{"status":"MANUFACTURED_NONEXECUTABLE_FIXTURE"}\n')
  self.spec={'manifest_repository_path':'research/prior/INPUT_MANIFEST.json','producer_repository_path':'research/prior/producer.py','prior_static_audit':{'repository_path':'research/prior/audit.json'},'capsules':[]}
  for i in range(4):
   compact=f'research/inputs/compact{i}.npz';original=f'research/inputs/original{i}.npz'
   add(compact,('MANUFACTURED-COMPACT-'+str(i)).encode());add(original,('MANUFACTURED-ORIGINAL-'+str(i)).encode())
   self.spec['capsules'].append({'repository_path':compact,'reader_inputspec':{'file_sha256':pin(self.payload[compact])['sha256']},'original':{'repository_path':original,'sha256':pin(self.payload[original])['sha256']}})
  for path in [self.spec['manifest_repository_path'],self.spec['producer_repository_path'],self.spec['prior_static_audit']['repository_path']]:add(path,b'ORIGINAL MANUFACTURED STATIC BYTES ONLY\n')
  self.spec['manifest_sha256']=pin(self.payload[self.spec['manifest_repository_path']])['sha256'];self.spec['producer_sha256']=pin(self.payload[self.spec['producer_repository_path']])['sha256'];self.spec['prior_static_audit']['sha256']=pin(self.payload[self.spec['prior_static_audit']['repository_path']])['sha256']
  add(self.m.PREFIX+'/provenance/INPUT_SPEC.json',encode(self.spec))
  # A source sentinel proves tests never import even externally pinned fixture code.
  add(self.m.PREFIX+'/execution/registration_guard.py',b'raise RuntimeError("SCIENTIFIC_IMPORT_FORBIDDEN_IN_BYTE_CLOSURE_TEST")\n')
  for mode in ('normal','optimized'):
   for name in self.m.STABLE: add('evidence/'+mode+'/'+name,('MANUFACTURED SAVED EXPORT '+name+'\n').encode())
   for name in ('RESOURCE.json','ENTRY_RECEIPT.json','EXECUTION.json','child.log'):add('evidence/'+mode+'/'+name,('MANUFACTURED '+mode+' '+name+'\n').encode())
  for name,raw in self.payload.items():p=self.root/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(raw)
  self.manifest={'schema_version':1,'files':{name:pin(raw) for name,raw in self.payload.items()}}
  self.pins={'schema_version':1,'sealed':True,'checkpoint_path':self.m.PREFIX,'public_go_path':'PUBLIC_GO.json','registration_sha256':pin(self.payload[self.m.PREFIX+'/FULL_REGISTRATION.json'])['sha256'],'public_go_sha256':pin(self.payload['PUBLIC_GO.json'])['sha256'],'payload_manifest_sha256':'','helper_sha256':pin(SRC.read_bytes())['sha256'],'runtime':self.m.RUNTIME,'saved_runs':{'normal':'evidence/normal','optimized':'evidence/optimized'},'stable_files':list(self.m.STABLE)}
  self.seal_fixture()
  self.socket_patch=mock.patch.object(socket,'socket',side_effect=AssertionError('REAL_NETWORK_FORBIDDEN'));self.socket_patch.start()
 def tearDown(self):self.socket_patch.stop();self.tmp.cleanup()
 def seal_fixture(self):
  raw=encode(self.manifest);(self.root/'PAYLOAD_MANIFEST.json').write_bytes(raw);self.pins['payload_manifest_sha256']=pin(raw)['sha256'];(self.root/'REPLAY_PINS.json').write_bytes(encode(self.pins))
  self.external=(self.pins['payload_manifest_sha256'],self.pins['helper_sha256'],pin((self.root/'REPLAY_PINS.json').read_bytes())['sha256'])
 def authenticate(self):return self.m.authenticate_package(self.root,*self.external)
 def test_complete_byte_closure_and_all8_inputs_no_scientific_import(self):
  with mock.patch.object(self.m,'study_module',side_effect=AssertionError('STUDY_IMPORT_BEFORE_AUTH')):
   v=self.authenticate();self.assertEqual(set(v['manifest']['files']),set(self.payload))
  self.assertNotIn('registration_guard',sys.modules)
 def test_changed_compact_or_original_input_fails(self):
  for name in ['research/inputs/compact0.npz','research/inputs/original3.npz']:
   original=(self.root/name).read_bytes();(self.root/name).write_bytes(original+b'changed')
   with self.assertRaisesRegex(self.m.Stop,'PACKAGE_MEMBER_PIN_DIFFERS'):self.authenticate()
   (self.root/name).write_bytes(original)
 def test_extra_python_before_any_import_fails(self):
  (self.root/'rogue.py').write_text('raise AssertionError("must not import")')
  with self.assertRaisesRegex(self.m.Stop,'COMPLETE_PACKAGE_MEMBERSHIP_DIFFERS'):self.authenticate()
 def test_extra_empty_directory_fails(self):
  (self.root/'unregistered').mkdir()
  with self.assertRaisesRegex(self.m.Stop,'UNREGISTERED_EMPTY_PACKAGE_DIRECTORY'):self.authenticate()
 def test_symlink_hardlink_fifo_rejected_without_blocking(self):
  extra=self.root/'extra'
  for create in [lambda:extra.symlink_to(self.root/'PUBLIC_GO.json'),lambda:os.link(self.root/'PUBLIC_GO.json',extra),lambda:os.mkfifo(extra)]:
   create()
   with self.assertRaises(self.m.Stop):self.authenticate()
   extra.unlink()
 def test_symlink_ancestor_rejected(self):
  alias=Path(self.tmp.name)/'alias';alias.symlink_to(self.root,target_is_directory=True)
  with self.assertRaisesRegex(self.m.Stop,'REAL_DIRECTORY'):self.m.authenticate_package(alias,*self.external)
 def test_pending_replay_seal_rejected(self):
  self.pins['sealed']=False;self.seal_fixture()
  with self.assertRaisesRegex(self.m.Stop,'FINAL_SEALED_REPLAY_PINS_REQUIRED'):self.authenticate()
 def test_external_helper_manifest_or_replay_pin_mismatch_fails(self):
  for i in range(3):
   pins=list(self.external);pins[i]='0'*64
   with self.assertRaisesRegex(self.m.Stop,'EXTERNAL_PIN_DIFFERS'):self.m.authenticate_package(self.root,*pins)
 def test_duplicate_nonfinite_and_float_manifest_rejected(self):
  for raw in [b'{"schema_version":1,"schema_version":1,"files":{}}',b'{"x":NaN}',b'{"x":1e309}']:
   with self.assertRaises((ValueError,self.m.Stop)):self.m.load(raw)
  self.manifest['files']['PUBLIC_GO.json']['bytes']=float(self.manifest['files']['PUBLIC_GO.json']['bytes']);self.seal_fixture()
  with self.assertRaisesRegex(self.m.Stop,'PACKAGE_MEMBER_SCHEMA_DIFFERS'):self.authenticate()
 def test_unsafe_and_circular_manifest_names_rejected(self):
  for name in ['../rogue','/rogue','a//b','a\\b','replay_checkpoint.py']:
   self.manifest['files'][name]=pin(b'');self.seal_fixture()
   with self.assertRaises(self.m.Stop):self.authenticate()
   del self.manifest['files'][name]
 def test_wrong_runtime_rejected(self):
  self.pins['runtime']={'python_version':'3.12.14','platform':'linux','packages':{}};self.seal_fixture()
  with self.assertRaisesRegex(self.m.Stop,'REPLAY_PIN_BINDING_OR_POLICY'):self.authenticate()
 def test_saved_scientific_pin_disagreement_rejected(self):
  name='evidence/optimized/DATA.json';(self.root/name).write_bytes(b'changed fabricated data');self.manifest['files'][name]=pin((self.root/name).read_bytes());self.seal_fixture()
  with self.assertRaisesRegex(self.m.Stop,'SAVED_NORMAL_OPTIMIZED_SCIENCE_PINS_DIFFER'):self.authenticate()
 def test_resource_time_and_RSS_not_byte_identity(self):
  self.authenticate();self.assertNotEqual(self.manifest['files']['evidence/normal/RESOURCE.json'],self.manifest['files']['evidence/optimized/RESOURCE.json'])
 def test_missing_original_input_binding_rejected(self):
  self.spec['capsules'][0]['original']['repository_path']='research/inputs/compact0.npz';self.spec['capsules'][0]['original']['sha256']=self.spec['capsules'][0]['reader_inputspec']['file_sha256'];name=self.m.PREFIX+'/provenance/INPUT_SPEC.json';(self.root/name).write_bytes(encode(self.spec));self.manifest['files'][name]=pin((self.root/name).read_bytes());self.seal_fixture()
  with self.assertRaisesRegex(self.m.Stop,'ALL_EIGHT_INPUTS_AND_THREE'):self.authenticate()
 def test_actual_authorization_failure_never_reads_saved_science_or_launches(self):
  v=self.authenticate();guard=types.SimpleNamespace(authenticate=mock.Mock(side_effect=ValueError('NO_ACTUAL_PUBLIC_GO')))
  with mock.patch.object(self.m.sys,'flags',types.SimpleNamespace(isolated=True)),mock.patch.object(self.m,'study_module',return_value=guard) as loader:
   with self.assertRaisesRegex(ValueError,'NO_ACTUAL_PUBLIC_GO'):self.m.replay(v,True)
   self.assertEqual(loader.call_count,1);self.assertEqual(guard.authenticate.call_count,1)
 def test_changed_registered_resource_gate_rejected(self):
  contract={'resources':dict(self.m.LIMITS,wall_seconds=901)};name=self.m.PREFIX+'/REGISTRATION_CONTRACT.json';(self.root/name).write_bytes(encode(contract));self.manifest['files'][name]=pin((self.root/name).read_bytes());self.seal_fixture()
  with self.assertRaisesRegex(self.m.Stop,'REGISTERED_FIXED_RESOURCE_POLICY_DIFFERS'):self.authenticate()
 def test_registered_commands_have_fixed_auth_paths_mode_and_isolation(self):
  v=self.authenticate();output=Path(self.tmp.name)/'fresh'/'normal'
  for mode in ('normal','optimized'):
   command=self.m.registered_command(v,output,mode)
   self.assertEqual(command[1:3],['-I','-B'])
   self.assertEqual(command.count('--registration-sha256'),1)
   self.assertEqual(command.count('--public-go-sha256'),1)
   self.assertNotIn('--fabricated-only',command)
   self.assertEqual('--optimized' in command,mode=='optimized')
   self.assertEqual(command[command.index('--repository-root')+1],str(self.root))
 def test_foreign_imported_study_module_rejected(self):
  with mock.patch.dict(sys.modules,{'registration_guard':types.SimpleNamespace()}):
   with self.assertRaisesRegex(self.m.Stop,'FOREIGN_STUDY_MODULE_ALREADY_IMPORTED'):
    self.m.study_module('registration_guard',self.root/self.m.PREFIX/'execution/registration_guard.py')
 def custody_fixture(self,mode='normal'):
  out=Path(self.tmp.name)/('manufactured-custody-'+mode);out.mkdir();(out/'ENTRY_RECEIPT.json').write_bytes(b'NONAUTHORIZING MANUFACTURED ENTRY BYTES');(out/'child.log').write_bytes(b'MANUFACTURED LOG')
  self.custody={'status':'PASS_BOUNDED_REGISTERED_EXECUTION','fabricated_only':False,'optimized':mode=='optimized','exit_code':0,'wrapper_exit_code':0,'wait_status':0,'stop_reason':None,'monitor_error':None,'limits':{'wall_seconds':900,'cpu_seconds':900,'address_space_bytes':536870912,'aggregate_output_bytes':134217728,'file_size_bytes':134217728,'process_limit':0,'core_bytes':0,'custodian_receipt_reserve_bytes':65536},'entry_receipt':pin((out/'ENTRY_RECEIPT.json').read_bytes()),'captured_child_log_sha256':pin((out/'child.log').read_bytes())['sha256'],'observed_child_log_sha256':pin((out/'child.log').read_bytes())['sha256'],'preserved_child_log':None,'uid':1000,'peak_rss_kib_wait4':1000,'wall_seconds':1.5,'cpu_seconds_wait4':1.0,'worker_output':{'bytes':1000,'file_count':12}}
  return out
 def write_custody(self,out): (out/'EXECUTION.json').write_bytes(encode(self.custody))
 def test_custody_validates_each_mode_without_public_GO(self):
  for mode in ('normal','optimized'):
   out=self.custody_fixture(mode);self.write_custody(out);self.assertEqual(self.m.validate_custodian(out,mode)['mode'],mode)
 def test_custody_failed_wrong_mode_or_fabricated_rejected(self):
  out=self.custody_fixture();original=dict(self.custody)
  for key,value in [('status','EXECUTION_FAILED'),('optimized',True),('fabricated_only',True)]:
   self.custody=dict(original,**{key:value});self.write_custody(out)
   with self.assertRaisesRegex(self.m.Stop,'CUSTODY_SCOPE_OR_MODE'):self.m.validate_custodian(out,'normal')
 def test_custody_wait_outcome_types_and_limits_rejected(self):
  out=self.custody_fixture();original=dict(self.custody)
  for key,value in [('exit_code',False),('wrapper_exit_code',0.0),('wait_status',256),('stop_reason','TIMEOUT'),('monitor_error','unsafe')]:
   self.custody=dict(original,**{key:value});self.write_custody(out)
   with self.assertRaisesRegex(self.m.Stop,'WAIT4_OUTCOME'):self.m.validate_custodian(out,'normal')
  self.custody=dict(original,limits=dict(original['limits'],process_limit=False));self.write_custody(out)
  with self.assertRaisesRegex(self.m.Stop,'CUSTODIAN_LIMITS'):self.m.validate_custodian(out,'normal')
 def test_custody_entry_and_log_linkage_rejected(self):
  out=self.custody_fixture();self.write_custody(out);(out/'ENTRY_RECEIPT.json').write_bytes(b'changed entry')
  with self.assertRaisesRegex(self.m.Stop,'ENTRY_OR_LOG_LINKAGE'):self.m.validate_custodian(out,'normal')
 def test_custody_measured_RSS_and_CPU_bound_rejected(self):
  out=self.custody_fixture();original=dict(self.custody)
  for key,value in [('uid',0),('peak_rss_kib_wait4',524289),('wall_seconds',901),('cpu_seconds_wait4',True),('worker_output',{'bytes':134217728,'file_count':12})]:
   self.custody=dict(original,**{key:value});self.write_custody(out)
   with self.assertRaises(self.m.Stop):self.m.validate_custodian(out,'normal')
 def test_helper_path_cannot_be_substituted(self):
  with mock.patch.object(self.m,'__file__',str(self.root/'other.py')):
   with self.assertRaisesRegex(self.m.Stop,'EXECUTED_HELPER_MUST_BE_PACKAGE_HELPER'):self.authenticate()

if __name__=='__main__':
 result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Controls))
 report={'status':'PASS' if result.wasSuccessful() else 'FAIL','tests_run':result.testsRun,'failures':[str(t) for t,_ in result.failures],'errors':[str(t) for t,_ in result.errors],'skipped':len(result.skipped),'helper_sha256':pin(SRC.read_bytes())['sha256'],'suite_sha256':pin(Path(__file__).read_bytes())['sha256'],'real_study_callbacks':0,'real_array_decodes':0,'fake_public_GO_created':False,'real_network':'FORBIDDEN','fixture_basis':'NONEXECUTABLE MANUFACTURED BYTE-CLOSURE PACKAGE ONLY'}
 (HERE/('CONTROLS_OPTIMIZED.json' if sys.flags.optimize else 'CONTROLS_NORMAL.json')).write_text(json.dumps(report,indent=2)+'\n')
 raise SystemExit(0 if result.wasSuccessful() else 1)
