"""Synthetic negative controls; no physical replay, network or publishing."""
import copy
from datetime import datetime,timezone
import hashlib,json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import monitor
import scientific_projection

def fabricated_diagnostic(route):
    """Projection-only fixture, never a claim of a valid physical experiment."""
    fields=scientific_projection.PRIMARY_TOP if route=='primary' else scientific_projection.INDEPENDENT_TOP
    value={name:None for name in fields}
    value.update(schema_version=1,fixed_cases=12,old_metric_status='FAIL',
                 producer_sha256='e'*64,input_manifest_sha256='0'*64)
    elapsed='seconds' if route=='primary' else 'elapsed_seconds'
    value['resources']={elapsed:101.0,'peak_rss_kib':1000,'scope':'fabricated entire route',
                        'limit_seconds':900,'limit_peak_rss_kib':262144}
    if route=='primary':
        value.update(status='COMPLETED_SAVED_DATA_DIAGNOSTIC',source_unchanged_during_execution=True,
                     registration_sha256='c'*64,freeze_commit='a'*40,
                     records=[{'controls':[{'profile':['0.1','0.2','0.3'],'resources':'retained scientific key'}]}])
    else:
        value.update(status='COMPLETED_DIAGNOSTIC_WITHOUT_RECLASSIFYING_OLD_FAIL',
                     resource_limits={'wall_seconds':900,'peak_rss_kib':262144},
                     rows=[{'profiles':{'R':['0.1','0.2','0.3']}}])
        value['provenance']={'freeze_commit':'a'*40,'registration_sha256':'c'*64,
            'registered_files_checked':1,'input_manifest_sha256':'0'*64,
            'post_run_frozen_verification':'PASS','full_frozen_verification':{
                'public_freeze_commit':'a'*40,'registration_sha256':'c'*64,
                'freeze_receipt_sha256':'1'*64,'frozen_files':{'inputs/INPUT_MANIFEST.json':'0'*64},
                'package_manifest_sha256':None,'input_capsule_verification':{
                    'capsules_verified':4,'members_verified':76,'array_payload_values_interpreted':False,
                    'input_manifest_sha256':'0'*64}}}
    return value

class MonitorGuards(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.target=json.loads((Path(__file__).parent/'TARGET.json').read_text())
        self.target.update(state='COMPLETED_REVIEWED_CHECKPOINT',science_commit='a'*40,freeze_commit='b'*40,
            registration_sha256='c'*64,requirements_sha256='d'*64,replay_driver_sha256='e'*64,
            expected_classification='CONSISTENCY_FAILURE',expected_scientific_sha256={
                'primary/diagnostic.json':'f'*64,'independent/diagnostic.json':'f'*64})
        self.now=datetime(2026,10,5,3,23,tzinfo=timezone.utc)
        self.context=dict(repository=monitor.REPOSITORY,branch='main',default_branch='main',attempt='1',
            event='schedule',enabled='true',start='2026-10-01T00:00:00Z',end='2026-10-31T00:00:00Z',mode='watch')
    def tearDown(self):self.tmp.cleanup()
    def pending_target(self):
        target=copy.deepcopy(self.target);target['state']='PENDING_COMPLETED_ACTIVE_CHECKPOINT'
        for name in ('science_commit','freeze_commit','registration_sha256','requirements_sha256',
                     'replay_driver_sha256','expected_classification'):target[name]=None
        target['expected_scientific_sha256']={}
        return target
    def test_pending_target_never_runs_science(self):
        pending=self.pending_target()
        result=monitor.gate(pending,self.context,self.now)
        self.assertFalse(result['ready']);self.assertFalse(result['do_replay']);self.assertFalse(result['do_watch'])
    def test_guessed_pending_pins_rejected(self):
        pending=self.pending_target();pending['science_commit']='a'*40
        with self.assertRaises(ValueError):monitor.validate_target(pending)
    def test_wrong_branch_repo_retry_enable_rejected(self):
        for field,value in [('branch','research/proposal'),('repository','other/HDblast'),('attempt','2'),('enabled','false')]:
            context=dict(self.context,**{field:value})
            with self.subTest(field=field),self.assertRaises(ValueError):monitor.gate(self.target,context,self.now)
    def test_manual_validation_does_not_require_schedule_enable(self):
        result=monitor.gate(self.target,dict(self.context,event='workflow_dispatch',enabled='false',mode='all'),self.now)
        self.assertTrue(result['do_replay']);self.assertTrue(result['do_watch'])
    def test_weekly_replay_daily_watch(self):
        monday=monitor.gate(self.target,self.context,self.now)
        tuesday=monitor.gate(self.target,self.context,self.now.replace(day=6))
        self.assertTrue(monday['do_replay']);self.assertFalse(tuesday['do_replay']);self.assertTrue(tuesday['do_watch'])
    def test_campaign_expired_future_oversize_offset_rejected(self):
        for start,end in [('2026-09-01T00:00:00Z','2026-10-01T00:00:00Z'),
                          ('2026-10-10T00:00:00Z','2026-10-20T00:00:00Z'),
                          ('2026-10-01T00:00:00Z','2026-11-01T00:00:00Z'),
                          ('2026-10-01T00:00:00+00:00','2026-10-31T00:00:00Z')]:
            with self.subTest(start=start,end=end),self.assertRaises(ValueError):monitor.campaign_ok(start,end,self.now)
    def test_unknown_schema_nonimmutable_and_historical_relabel_rejected(self):
        mutations=[{'extra':True},{'science_commit':'main'},{'registration_sha256':'0'*63},
                   {'expected_old_metric_status':'PASS'},{'expected_classification':'HYPOTHESIS_VALIDATED'}]
        for mutation in mutations:
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):monitor.validate_target(dict(self.target,**mutation))
    def test_scientific_paths_cannot_escape_or_target_runtime_receipts(self):
        for name in ('../secret.json','primary/../../secret.json','primary\\secret.json','/primary/a.json',
                     'logs/resource.json','primary/../a.json','primary/a.sh'):
            target=copy.deepcopy(self.target);target['expected_scientific_sha256'][name]='f'*64
            with self.subTest(name=name),self.assertRaises(ValueError):monitor.validate_target(target)
    def test_duplicate_json_and_nonfinite_rejected(self):
        for text in ('{"a":1,"a":2}','{"a":NaN}'):
            p=self.root/'bad.json';p.write_text(text)
            with self.assertRaises(ValueError):monitor.load_json(p)
    def test_symlink_sources_and_existing_outputs_rejected(self):
        p=self.root/'payload';p.write_text('x');(self.root/'link').symlink_to(p)
        with self.assertRaises(ValueError):monitor.safe_file(self.root,'link')
        with self.assertRaises(ValueError):monitor.fresh_output(self.root)
    def metadata(self,identifier='1',title='Fabricated title'):
        return {'hits':{'total':1,'hits':[{'id':identifier,'metadata':{'titles':[{'title':title}],
            'authors':[{'full_name':'Synthetic Author'}],'arxiv_eprints':[],'dois':[],'preprint_date':'2026-10-01'}}]}}
    def test_oversized_hits_and_unknown_schema_rejected(self):
        for data in ({'hits':{'total':0,'hits':[]}},self.metadata()):monitor.normalize_records(data)
        for data in ({},{'hits':{'total':False,'hits':[]}},
                     {'hits':{'total':21,'hits':self.metadata()['hits']['hits']*21}},
                     {'hits':{'total':2,'hits':self.metadata()['hits']['hits']*2}}):
            with self.assertRaises(ValueError):monitor.normalize_records(data)
    def test_metadata_diff_keeps_absence_distinct_from_retraction(self):
        previous=monitor.normalize_records(self.metadata('1'));current=monitor.normalize_records(self.metadata('2'))
        result=monitor.bibliographic_diff(previous,current)
        self.assertEqual(result['added'][0]['id'],'2');self.assertEqual(result['absent_from_current_first_page'][0]['id'],'1')
        self.assertNotIn('retracted',result)
    def baseline(self):
        records=monitor.normalize_records(self.metadata())
        return {'schema_version':1,'queries':{k:{'query':q,'records':records} for k,q in monitor.QUERIES.items()}}
    def test_exact_two_requests_and_no_retry_after_transport_failure(self):
        calls=[]
        def getter(url):calls.append(url);raise TimeoutError('Fabricated timeout')
        output=monitor.fresh_output(self.root/'watch')
        state=monitor.watch(self.baseline(),output,getter)
        self.assertEqual(len(calls),2);self.assertEqual(state['status'],'FAILED_METADATA_READ')
        self.assertTrue(all(r['attempts']==1 for r in state['queries'].values()))
    def test_bounded_response_and_error_details_not_saved(self):
        calls=[]
        def getter(url):calls.append(url);return b'x'*(monitor.MAX_RESPONSE_BYTES+1),{}
        state=monitor.watch(self.baseline(),monitor.fresh_output(self.root/'large'),getter)
        self.assertEqual(state['status'],'FAILED_METADATA_READ');self.assertEqual(len(calls),2)
        self.assertNotIn('error_message',next(iter(state['queries'].values())))
    def test_title_cannot_inject_job_summary(self):
        data=self.metadata(title='[Click](https://example.invalid)\n<script>bad()</script>')
        raw=json.dumps(data).encode()
        state=monitor.watch(self.baseline(),monitor.fresh_output(self.root/'summary'),lambda url:(raw,{'status':200}))
        path=self.root/'summary.md';monitor.write_summary(path,state)
        self.assertNotIn('example.invalid',path.read_text());self.assertNotIn('<script>',path.read_text())
    def test_redirects_and_unapproved_endpoints_never_followed(self):
        handler=monitor.SameOriginRedirect()
        with self.assertRaises(ValueError):handler.redirect_request(None,None,302,'redirect',{},'https://inspirehep.net/api/literature')
        for url in ('http://inspirehep.net/api/literature','https://evil.invalid/api/literature',
                    'https://inspirehep.net/account','https://user:pass@inspirehep.net/api/literature'):
            with self.assertRaises(ValueError):monitor.public_get(url)
    def fixture_replay(self):
        output=self.root/'replay';(output/'fresh'/'primary').mkdir(parents=True);(output/'fresh'/'independent').mkdir()
        target=copy.deepcopy(self.target)
        for name in ('primary/diagnostic.json','independent/diagnostic.json'):
            p=output/'fresh'/name
            value=fabricated_diagnostic(name.split('/')[0]);monitor.write_json(p,value)
            target['expected_scientific_sha256'][name]=scientific_projection.project(value,name.split('/')[0])[0]
        receipt=dict(status='COMPLETED_REPLAY',plan_only=False,classification='CONSISTENCY_FAILURE',old_metric_status='FAIL',
            physical_routes_completed=2,successful_commands=4,source_verification_before='PASS',source_verification_after='PASS')
        monitor.write_json(output/'REPLAY.json',receipt);return target,output,receipt
    def test_reproducing_negative_science_is_not_relabeling_failure(self):
        target,output,_=self.fixture_replay();result=monitor.verify_completed_replay(target,output)
        self.assertEqual(result['status'],'VERIFIED_IDENTICAL_COMPLETED_REPLAY');self.assertEqual(result['classification'],'CONSISTENCY_FAILURE')
    def test_completed_comparison_cannot_accept_empty_pending_target(self):
        target,output,_=self.fixture_replay()
        with self.assertRaises(ValueError):monitor.verify_completed_replay(self.pending_target(),output)
    def test_plan_only_incomplete_source_mutation_or_relabel_rejected(self):
        target,output,receipt=self.fixture_replay()
        for changes in ({'plan_only':True},{'physical_routes_completed':1},{'successful_commands':3},
                        {'source_verification_after':'FAIL'},{'old_metric_status':'PASS'},
                        {'classification':'LEDGER_ERROR_DEMONSTRATED'}):
            monitor.write_json(output/'REPLAY.json',dict(receipt,**changes))
            with self.subTest(changes=changes),self.assertRaises(ValueError):monitor.verify_completed_replay(target,output)
    def test_scientific_output_mutation_detected(self):
        target,output,_=self.fixture_replay();(output/'fresh'/'primary'/'diagnostic.json').write_text('changed')
        with self.assertRaises(ValueError):monitor.verify_completed_replay(target,output)
    def test_resource_only_variations_accepted_and_raw_hashes_retained(self):
        target,output,_=self.fixture_replay()
        original={name:monitor.sha(output/'fresh'/name) for name in target['expected_scientific_sha256']}
        for route in ('primary','independent'):
            path=output/'fresh'/route/'diagnostic.json';value=monitor.load_json(path)
            value['resources']['seconds' if route=='primary' else 'elapsed_seconds']=102.0
            value['resources']['peak_rss_kib']=1001;monitor.write_json(path,value)
        result=monitor.verify_completed_replay(target,output)
        for name,check in result['scientific_checks'].items():
            self.assertNotEqual(original[name],check['raw_sha256']);self.assertTrue(check['equal'])
            self.assertEqual(check['retained_variable_observations']['excluded_resource_values']['peak_rss_kib'],1001)
    def test_package_receipt_variation_retained_separately(self):
        target,output,_=self.fixture_replay();path=output/'fresh'/'independent'/'diagnostic.json'
        value=monitor.load_json(path);value['provenance']['full_frozen_verification']['package_manifest_sha256']='9'*64
        monitor.write_json(path,value);result=monitor.verify_completed_replay(target,output)
        self.assertEqual(result['scientific_checks']['independent/diagnostic.json']['retained_variable_observations']['excluded_package_manifest_sha256'],'9'*64)
    def test_resource_limits_and_scope_cannot_be_excluded(self):
        for route in ('primary','independent'):
            original=fabricated_diagnostic(route);reference=scientific_projection.project(original,route)[0]
            value=copy.deepcopy(original);value['resources']['scope']='changed scope'
            self.assertNotEqual(reference,scientific_projection.project(value,route)[0])
            value=copy.deepcopy(original);value['resources']['limit_seconds']=901
            with self.assertRaises(ValueError):scientific_projection.project(value,route)
    def test_resource_nonfinite_boolean_over_budget_or_zero_rejected(self):
        for route in ('primary','independent'):
            elapsed='seconds' if route=='primary' else 'elapsed_seconds'
            for key,replacement in [(elapsed,float('nan')),(elapsed,float('inf')),(elapsed,True),
                    (elapsed,901.0),(elapsed,0),('peak_rss_kib',True),('peak_rss_kib',262145),('peak_rss_kib',0)]:
                value=fabricated_diagnostic(route);value['resources'][key]=replacement
                with self.subTest(route=route,key=key,replacement=replacement),self.assertRaises(ValueError):
                    scientific_projection.project(value,route)
    def test_unknown_or_missing_exclusion_ancestor_fields_rejected(self):
        for route in ('primary','independent'):
            for where in ('top','resources'):
                value=fabricated_diagnostic(route);obj=value if where=='top' else value['resources'];obj['unknown']=1
                with self.assertRaises(ValueError):scientific_projection.project(value,route)
                value=fabricated_diagnostic(route);obj=value if where=='top' else value['resources'];del obj[next(iter(obj))]
                with self.assertRaises(ValueError):scientific_projection.project(value,route)
        for where in ('provenance','full_frozen_verification','input_capsule_verification'):
            value=fabricated_diagnostic('independent');obj=value['provenance']
            if where!='provenance':obj=obj['full_frozen_verification']
            if where=='input_capsule_verification':obj=obj['input_capsule_verification']
            obj['unknown']=1
            with self.assertRaises(ValueError):scientific_projection.project(value,'independent')
    def test_late_scientific_profile_and_same_named_resource_key_retained(self):
        for route in ('primary','independent'):
            value=fabricated_diagnostic(route);reference=scientific_projection.project(value,route)[0]
            if route=='primary':value['records'][0]['controls'][0]['profile'][-1]='0.3001'
            else:value['rows'][0]['profiles']['R'][-1]='0.3001'
            self.assertNotEqual(reference,scientific_projection.project(value,route)[0])
        value=fabricated_diagnostic('primary');reference=scientific_projection.project(value,'primary')[0]
        value['records'][0]['controls'][0]['resources']='changed scientific key'
        self.assertNotEqual(reference,scientific_projection.project(value,'primary')[0])
    def test_retained_freeze_producer_and_receipt_pins_change_digest(self):
        value=fabricated_diagnostic('independent');reference=scientific_projection.project(value,'independent')[0]
        for key in ('producer_sha256',):
            changed=copy.deepcopy(value);changed[key]='8'*64
            self.assertNotEqual(reference,scientific_projection.project(changed,'independent')[0])
        changed=copy.deepcopy(value);changed['provenance']['full_frozen_verification']['freeze_receipt_sha256']='8'*64
        self.assertNotEqual(reference,scientific_projection.project(changed,'independent')[0])
        changed=copy.deepcopy(value);changed['provenance']['freeze_commit']='b'*40
        with self.assertRaises(ValueError):scientific_projection.project(changed,'independent')
    def test_unknown_projection_policy_rejected(self):
        with self.assertRaises(ValueError):monitor.validate_target(dict(self.target,scientific_hash_policy='RAW_OR_GUESS'))
    def test_json_formatting_variation_preserves_exact_scientific_values(self):
        target,output,_=self.fixture_replay()
        for route in ('primary','independent'):
            path=output/'fresh'/route/'diagnostic.json';value=monitor.load_json(path)
            path.write_text(json.dumps(value,sort_keys=False,separators=(',',':')))
        result=monitor.verify_completed_replay(target,output)
        self.assertTrue(all(item['equal'] for item in result['scientific_checks'].values()))
    def test_invalid_package_receipt_or_capsule_membership_rejected(self):
        value=fabricated_diagnostic('independent')
        value['provenance']['full_frozen_verification']['package_manifest_sha256']='main'
        with self.assertRaises(ValueError):scientific_projection.project(value,'independent')
        value=fabricated_diagnostic('independent')
        value['provenance']['full_frozen_verification']['input_capsule_verification']['members_verified']=True
        with self.assertRaises(ValueError):scientific_projection.project(value,'independent')
    def source_fixture(self):
        source=self.root/'source';root=source/monitor.CHECKPOINT
        (root/'code').mkdir(parents=True)
        files={}
        for name in ('EXPERIMENT.json','requirements-replay.txt','code/replay_active_source.py','code/verify_public_freeze.py'):
            path=root/name;path.write_text('fabricated source, never executed\n');files[name]=monitor.sha(path)
        monitor.write_json(root/'FULL_REGISTRATION.json',{'schema_version':1,'files':files})
        target=copy.deepcopy(self.target);target.update(registration_sha256=monitor.sha(root/'FULL_REGISTRATION.json'),
            requirements_sha256=files['requirements-replay.txt'],replay_driver_sha256=files['code/replay_active_source.py'])
        return target,source,root
    def test_all_frozen_bytes_checked_before_any_science_code_execution(self):
        target,source,root=self.source_fixture()
        with patch('monitor.subprocess.run') as run:
            run.return_value.stdout=target['science_commit']+'\n'
            monitor.verify_checkout(target,source)
            (root/'code'/'verify_public_freeze.py').write_text('fabricated source mutation')
            with self.assertRaises(ValueError):monitor.verify_checkout(target,source)
            self.assertTrue(all(call.args[0][0]=='git' for call in run.call_args_list))
    def test_unregistered_module_shadowing_rejected_before_science_launch(self):
        target,source,root=self.source_fixture();(root/'code'/'json.py').write_text('fabricated shadow, never executed')
        with patch('monitor.subprocess.run') as run:
            run.return_value.stdout=target['science_commit']+'\n'
            with self.assertRaises(ValueError):monitor.verify_checkout(target,source)
            self.assertTrue(all(call.args[0][0]=='git' for call in run.call_args_list))

if __name__=='__main__':unittest.main(verbosity=2)
