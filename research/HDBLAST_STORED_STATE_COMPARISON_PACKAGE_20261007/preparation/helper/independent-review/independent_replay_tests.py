#!/usr/bin/env python3
"""Independent manufactured closure/orchestration controls; no real study APIs."""
import ast
from contextlib import ExitStack
import hashlib
import importlib.util
import json
from pathlib import Path
import signal
import socket
import ssl
import sys
import tempfile
import types
import unittest
from unittest import mock
from urllib import request

HERE = Path(__file__).resolve().parent
DELIVERY = HERE.parent
W2 = DELIVERY.parent
SNAPSHOT = HERE/'replay_checkpoint_5d50ce6c_snapshot.py'
EXPECTED_SHA = '5d50ce6c24ebf7f9faeb0704201fee13bf3ce2dc94f5d56eaa7325c010ece0f0'
if hashlib.sha256(SNAPSHOT.read_bytes()).hexdigest() != EXPECTED_SHA:
    raise RuntimeError('REPLAY_HELPER_SNAPSHOT_PIN_DIFFERS')
sys.dont_write_bytecode = True

def forbidden(*args,**kwargs):
    raise AssertionError('REAL_NETWORK_PROCESS_OR_STUDY_IMPORT_FORBIDDEN')
socket.socket = forbidden
socket.create_connection = forbidden
request.urlopen = forbidden
request.OpenerDirector.open = forbidden

def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec)
    sys.modules[name]=result
    spec.loader.exec_module(result)
    return result
F=module(DELIVERY/'test_replay_checkpoint.py','independent_nonexecutable_package_fixture')
F.SRC=SNAPSHOT
work=HERE/'temporary-nonauthorization-fixtures'
work.mkdir(exist_ok=True)
tempfile.tempdir=str(work)

class Independent(F.Controls):
    def setUp(self):
        super().setUp()
        self.calls=[]
        self.context={'scope':'MANUFACTURED_STUB_NOT_ACTUAL_PUBLIC_GO_AUTHORIZATION'}
        for mode in ('normal','optimized'):
            name='evidence/'+mode+'/EXECUTION.json'
            raw=F.encode(self.execution_receipt(mode,self.root/'evidence'/mode))
            (self.root/name).write_bytes(raw)
            self.manifest['files'][name]=F.pin(raw)
        self.seal_fixture()

    def execution_receipt(self,mode,run):
        entry=run/'ENTRY_RECEIPT.json'
        log=run/'child.log'
        return {'status':'PASS_BOUNDED_REGISTERED_EXECUTION','fabricated_only':False,
                'optimized':mode=='optimized','uid':1234,'exit_code':0,'wait_status':0,
                'wrapper_exit_code':0,'stop_reason':None,'monitor_error':None,
                'entry_receipt':F.pin(entry.read_bytes()) if entry.exists() else F.pin(b'fake-entry'),
                'captured_child_log_sha256':F.pin(log.read_bytes())['sha256'] if log.exists() else F.pin(b'fake-log')['sha256'],
                'observed_child_log_sha256':F.pin(log.read_bytes())['sha256'] if log.exists() else F.pin(b'fake-log')['sha256'],
                'wall_seconds':0.1,'cpu_seconds_wait4':0.1,'peak_rss_kib_wait4':1024,
                'limits':{'wall_seconds':900,'cpu_seconds':900,'address_space_bytes':536870912,
                          'aggregate_output_bytes':134217728,'file_size_bytes':134217728,
                          'process_limit':0,'core_bytes':0,'custodian_receipt_reserve_bytes':65536}}

    def change_saved_receipt(self,mode,change):
        name='evidence/'+mode+'/EXECUTION.json'
        value=json.loads((self.root/name).read_bytes())
        change(value)
        (self.root/name).write_bytes(F.encode(value))
        self.manifest['files'][name]=F.pin((self.root/name).read_bytes())
        self.seal_fixture()

    def fake_flow(self,validate_hook=None,compare_hook=None,launch_hook=None):
        def authenticate(*args):
            self.calls.append(('authenticate',args))
            return self.context
        def reauthenticate(verified):
            self.assertIs(verified,self.context)
            self.calls.append(('guard_reauthenticate',))
        def validate(path,candidate,fabricated,verified):
            self.calls.append(('validate',path,candidate,fabricated,verified))
            self.assertFalse(fabricated)
            self.assertIs(verified,self.context)
            if validate_hook:validate_hook(path)
            return {'status':'MANUFACTURED_SAVED_READBACK_STUB_ONLY'}
        def compare(first,second):
            self.calls.append(('compare',first,second))
            if compare_hook:compare_hook(first,second)
            return {'status':'MANUFACTURED_COMPARE_STUB_ONLY'}
        def launch(v,run,mode,logs):
            self.calls.append(('launch',mode,run))
            if launch_hook:return launch_hook(run,mode)
            run.mkdir()
            receipt=self.execution_receipt(mode,run)
            (run/'EXECUTION.json').write_bytes(F.encode(receipt))
            return receipt
        guard=types.SimpleNamespace(authenticate=authenticate,reauthenticate=reauthenticate)
        reader=types.SimpleNamespace(validate_output=validate,compare_replays=compare)
        def loader(name,path):
            self.calls.append(('import',name,path))
            if name=='registration_guard':return guard
            if name=='readback':return reader
            raise AssertionError('UNEXPECTED_MANUFACTURED_STUDY_MODULE')
        stack=ExitStack()
        stack.enter_context(mock.patch.object(self.m.sys,'flags',types.SimpleNamespace(isolated=True)))
        stack.enter_context(mock.patch.object(self.m,'study_module',side_effect=loader))
        stack.enter_context(mock.patch.object(self.m,'launch_registered',side_effect=launch))
        stack.enter_context(mock.patch.object(self.m.os,'getuid',return_value=1234))
        stack.enter_context(mock.patch.object(self.m.os,'geteuid',return_value=1234))
        return stack

    def test_independent_changed_helper_bytes_reject_before_study_import(self):
        (self.root/'replay_checkpoint.py').write_bytes(SNAPSHOT.read_bytes()+b'\n# changed\n')
        with mock.patch.object(self.m,'study_module',side_effect=forbidden):
            with self.assertRaisesRegex(self.m.Stop,'HELPER_EXTERNAL_PIN_DIFFERS'):
                self.authenticate()

    def test_independent_changed_control_after_auth_rejects_before_import(self):
        v=self.authenticate()
        (self.root/'REPLAY_PINS.json').write_bytes((self.root/'REPLAY_PINS.json').read_bytes()+b'\n')
        with mock.patch.object(self.m.sys,'flags',types.SimpleNamespace(isolated=True)),mock.patch.object(self.m,'study_module',side_effect=forbidden):
            with self.assertRaisesRegex(self.m.Stop,'REPLAY_PINS_EXTERNAL_PIN_DIFFERS'):
                self.m.replay(v,True)

    def test_independent_manifested_native_code_rejected(self):
        (self.root/'extra.so').write_bytes(b'MANUFACTURED_NATIVE_NAME_NOT_NATIVE_CODE')
        self.manifest['files']['extra.so']=F.pin((self.root/'extra.so').read_bytes());self.seal_fixture()
        with self.assertRaisesRegex(self.m.Stop,'CACHED_OR_NATIVE_PACKAGE_CODE_FORBIDDEN'):
            self.authenticate()

    def test_independent_isolated_mode_required_before_guard_import(self):
        v=self.authenticate()
        with mock.patch.object(self.m.sys,'flags',types.SimpleNamespace(isolated=False)),mock.patch.object(self.m,'study_module',side_effect=forbidden):
            with self.assertRaisesRegex(self.m.Stop,'USE_PYTHON_ISOLATED_MODE'):
                self.m.replay(v,True)

    def test_independent_actual_guard_precedes_reader_and_verified_context(self):
        v=self.authenticate()
        with self.fake_flow():result=self.m.replay(v,True)
        kinds=[item[0] for item in self.calls]
        self.assertEqual(kinds[:3],['import','authenticate','import'])
        self.assertEqual(self.calls[0][1],'registration_guard')
        self.assertEqual(self.calls[2][1],'readback')
        self.assertEqual(len([x for x in self.calls if x[0]=='validate']),2)
        self.assertFalse(any(x[0]=='launch' for x in self.calls))
        self.assertEqual(result['new_array_decodes'],0)
        self.assertEqual(result['new_source_callbacks'],0)

    def test_independent_package_change_during_saved_readback_blocks_fresh(self):
        v=self.authenticate()
        def change(path):
            target=self.root/'PUBLIC_GO.json'
            if target.read_bytes()==self.payload['PUBLIC_GO.json']:target.write_bytes(b'changed-invalid-GO')
        with self.fake_flow(validate_hook=change):
            with self.assertRaisesRegex(self.m.Stop,'PACKAGE_MEMBER_PIN_DIFFERS'):
                self.m.replay(v,False,Path(self.tmp.name)/'fresh')
        self.assertFalse(any(x[0]=='launch' for x in self.calls))

    def test_independent_saved_failed_custodian_must_block_saved_pass(self):
        self.change_saved_receipt('normal',lambda r:r.update({'status':'EXECUTION_FAILED','wrapper_exit_code':1}))
        v=self.authenticate()
        with self.fake_flow():
            with self.assertRaises(self.m.Stop,msg='SAVED_FAILED_CUSTODIAN_ACCEPTED'):
                self.m.replay(v,True)
        self.assertFalse(any(x[0]=='launch' for x in self.calls))

    def test_independent_saved_wrong_optimized_mode_must_block_saved_pass(self):
        self.change_saved_receipt('optimized',lambda r:r.update({'optimized':False}))
        v=self.authenticate()
        with self.fake_flow():
            with self.assertRaises(self.m.Stop,msg='SAVED_OPTIMIZED_MODE_NOT_ENFORCED'):
                self.m.replay(v,True)

    def test_independent_fresh_modes_compare_to_saved_and_each_other(self):
        v=self.authenticate();out=Path(self.tmp.name)/'fresh'
        with self.fake_flow():result=self.m.replay(v,False,out)
        self.assertEqual([x[1] for x in self.calls if x[0]=='launch'],['normal','optimized'])
        pairs=[(x[1],x[2]) for x in self.calls if x[0]=='compare']
        saved=self.root/'evidence/normal'
        self.assertIn((saved,out/'normal'),pairs)
        self.assertIn((saved,out/'optimized'),pairs)
        self.assertIn((out/'normal',out/'optimized'),pairs)
        self.assertEqual(result['exact_scientific_files'],list(self.m.SCIENCE))
        self.assertEqual(result['other_exact_file'],'OUTPUT_READBACK.json')

    def test_independent_failed_normal_never_launches_optimized(self):
        v=self.authenticate()
        def fail(run,mode):
            run.mkdir()
            raise self.m.Stop('MANUFACTURED_NORMAL_FAILURE')
        with self.fake_flow(launch_hook=fail):
            with self.assertRaisesRegex(self.m.Stop,'MANUFACTURED_NORMAL_FAILURE'):
                self.m.replay(v,False,Path(self.tmp.name)/'fresh')
        self.assertEqual([x[1] for x in self.calls if x[0]=='launch'],['normal'])

    def test_independent_first_fresh_science_mismatch_never_launches_optimized(self):
        v=self.authenticate();out=Path(self.tmp.name)/'fresh'
        def mismatch(first,second):
            if second==out/'normal':raise self.m.Stop('MANUFACTURED_SCIENCE_MISMATCH')
        with self.fake_flow(compare_hook=mismatch):
            with self.assertRaisesRegex(self.m.Stop,'MANUFACTURED_SCIENCE_MISMATCH'):
                self.m.replay(v,False,out)
        self.assertEqual([x[1] for x in self.calls if x[0]=='launch'],['normal'])

    def test_independent_fresh_output_inside_package_rejected(self):
        v=self.authenticate()
        with self.fake_flow():
            with self.assertRaisesRegex(self.m.Stop,'REPLAY_OUTPUT_MUST_BE_OUTSIDE_PACKAGE'):
                self.m.replay(v,False,self.root/'fresh')
        self.assertFalse(any(x[0]=='launch' for x in self.calls))

    def test_independent_custodian_subprocess_has_only_fixed_environment(self):
        v=self.authenticate();out=Path(self.tmp.name)/'fake-process';out.mkdir()
        run=out/'optimized';observed={};outer=self
        class FakeProcess:
            def __init__(self,command,**kwargs):
                observed.update({'command':command,'kwargs':kwargs})
                run.mkdir();(run/'EXECUTION.json').write_bytes(F.encode(outer.execution_receipt('optimized',run)))
            def wait(self,timeout=None):return 0
        with mock.patch.object(self.m.subprocess,'Popen',side_effect=FakeProcess):
            self.m.launch_registered(v,run,'optimized',out)
        self.assertEqual(observed['kwargs']['env'],{'PATH':'/usr/bin:/bin','LANG':'C.UTF-8','LC_ALL':'C.UTF-8'})
        self.assertEqual(observed['command'][1:3],['-I','-B'])
        self.assertIn('--optimized',observed['command'])
        self.assertNotIn('--fabricated-only',observed['command'])
        self.assertEqual(observed['kwargs']['stdin'],self.m.subprocess.DEVNULL)

    def test_independent_interruption_signals_custodian_cleanup(self):
        v=self.authenticate();out=Path(self.tmp.name)/'fake-process';out.mkdir();signals=[]
        class FakeProcess:
            def __init__(self,*args,**kwargs):self.calls=0
            def wait(self,timeout=None):
                self.calls+=1
                if self.calls==1:raise KeyboardInterrupt()
                return 130
            def send_signal(self,value):signals.append(value)
            def kill(self):raise AssertionError('UNNECESSARY_FAKE_KILL')
        with mock.patch.object(self.m.subprocess,'Popen',side_effect=FakeProcess):
            with self.assertRaises(KeyboardInterrupt):self.m.launch_registered(v,out/'normal','normal',out)
        self.assertEqual(signals,[signal.SIGINT])

    def test_independent_scientific_identity_names_match_frozen_readback(self):
        tree=ast.parse((W2/'execution/readback.py').read_text())
        comparison=next(node for node in tree.body if isinstance(node,ast.FunctionDef) and node.name=='compare_replays')
        stable=next(ast.literal_eval(node.value) for node in comparison.body if isinstance(node,ast.Assign) and any(isinstance(target,ast.Name) and target.id=='stable' for target in node.targets))
        self.assertEqual(self.m.STABLE,stable)

if __name__=='__main__':
    names=sorted(name for name in Independent.__dict__ if name.startswith('test_independent_'))
    result=unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(Independent(name) for name in names))
    report={'status':'PASS' if result.wasSuccessful() else 'FAIL','tests_run':result.testsRun,
            'failures':[{'test':str(test),'traceback':trace} for test,trace in result.failures],
            'errors':[{'test':str(test),'traceback':trace} for test,trace in result.errors],
            'skipped':list(result.skipped),'helper_snapshot_sha256':EXPECTED_SHA,
            'suite_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'fixture_support_sha256':hashlib.sha256((DELIVERY/'test_replay_checkpoint.py').read_bytes()).hexdigest(),
            'real_network':'FORBIDDEN','real_subprocesses':0,'real_scientific_imports':0,
            'real_array_decodes':0,'real_source_callbacks':0,
            'fixture_public_go_status':'NOT_AN_AUTHORIZATION_RECEIPT',
            'actual_authorization_and_readback':'NOT_EXECUTED_IN_SUCCESS_STUB_FLOWS',
            'orchestration_basis':'IN_MEMORY_GUARD_READER_PROCESS_STUBS_ONLY_NO_SCIENTIFIC_SUCCESS_CLAIM'}
    mode='optimized' if sys.flags.optimize else 'normal'
    (HERE/('REPLAY_INDEPENDENT_RESULTS_'+mode+'.json')).write_text(json.dumps(report,indent=2)+'\n')
    raise SystemExit(0 if result.wasSuccessful() else 1)
