#!/usr/bin/env python3
"""Additional manufactured continuation guards; forbidden real network."""
from pathlib import Path
import copy,importlib.util,json,sys,unittest
HERE=Path(__file__).resolve().parent
s=importlib.util.spec_from_file_location('continuation_flow_fixtures',HERE/'offline_tests.py');f=importlib.util.module_from_spec(s);sys.modules[s.name]=f;s.loader.exec_module(f)
M=f.M
class ContinuationControls(f.ControllerFixtures):
    # Run these controls separately; inherited suite runs independently.
    def test_pending_inventory_never_instantiates_write_controller(self):
        for field,value,reason in [('sealed_upload_manifest',False,'UNSEALED_OR_WRONG_FAMILY_INVENTORY'),('scientific_replay_status','PENDING','COMPLETED_SCIENTIFIC_REPLAY_REQUIRED'),('expected_total_files_final',29.0,'FINAL_INVENTORY_TOTAL_DIFFERS')]:
            original=self.inventory[field];self.inventory[field]=value
            self.stop(reason,self.make_controller)
            self.assertFalse(self.server.calls)
            self.inventory[field]=original
    def test_old_completed_state_rejected_without_remote_access(self):
        p=self.root/'evidence/STATE.json'
        p.write_text(json.dumps({'pending':None,'draft_id':'23225288','published':True,'publish_post_attempted':True}))
        self.stop('JOURNAL_BASELINE_OR_CONTROLLER_BINDING_DIFFERS',self.make_controller)
        self.assertFalse(self.server.calls)
    def test_state_binding_rejects_a_changed_baseline(self):
        self.baseline['files'][0]['sha256']='0'*64
        self.prior['rows'][0]['observed']['sha256']='0'*64
        self.inventory['inherited_files'][0]['sha256']='0'*64
        self.stop('JOURNAL_BASELINE_OR_CONTROLLER_BINDING_DIFFERS',self.make_controller)
        self.assertFalse(self.server.calls)
    def test_each_inherited_preservation_guard_rejects(self):
        candidates=[('creators',[{'person_or_org':{'type':'personal','name':'Changed author'}}],'INHERITED_METADATA_CREATORS_DIFFERS'),('rights',[{'id':'cc0-1.0'}],'INHERITED_METADATA_RIGHTS_DIFFERS'),('copyright','Changed','INHERITED_METADATA_COPYRIGHT_DIFFERS')]
        for field,value,reason in candidates:
            original=self.metadata['metadata'][field];self.metadata['metadata'][field]=value
            self.stop(reason,self.make_controller);self.metadata['metadata'][field]=original
        self.assertFalse(self.server.calls)
    def test_inherited_notes_and_references_cannot_be_replaced(self):
        notes=self.metadata['metadata']['additional_descriptions'];notes[1]['description']='<p>Changed old notes.</p>'
        self.stop('HISTORICAL_NOTES_DIFFERS',self.make_controller);notes[1]['description']='<p>Preserved old notes.</p>'
        refs=self.metadata['metadata']['related_identifiers'];saved=copy.deepcopy(refs);refs.clear()
        self.stop('INHERITED_REFERENCES_REMOVED',self.make_controller);refs.extend(saved)
        refs[0]['identifier']='10.5281/zenodo.999'
        self.stop('INHERITED_REFERENCES_DIFFERS',self.make_controller)
        self.assertFalse(self.server.calls)
    def test_unreviewed_append_rejected_before_transport(self):
        self.metadata['metadata']['related_identifiers'].append({'identifier':'10.5281/zenodo.999','scheme':'doi','relation_type':{'id':'references'},'resource_type':{'id':'software'}})
        self.stop('UNREVIEWED_REFERENCE_ADDITION',self.make_controller)
        self.assertFalse(self.server.calls)
    def test_original_prior_reference_append_is_allowed(self):
        self.metadata['metadata']['related_identifiers'].append({'identifier':'10.5281/zenodo.'+M.PRIOR_ID,'scheme':'doi','relation_type':{'id':'references'},'resource_type':{'id':'software'}})
        (self.root/'evidence').rename(self.root/'evidence-preserved-original-metadata')
        self.controller=self.make_controller();self.assertEqual(len(self.controller.metadata['metadata']['related_identifiers']),2)
        self.assertFalse(self.server.calls)
    def test_all29_streamed_before_and_after_even_when_identity_matches(self):
        result=self.prepared();self.assertEqual(len(result['content']),29)
        self.assertEqual(len(self.server.streams),29)
        result=self.controller.publish()
        self.assertEqual(len(result['content']),29)
        # publish repeats full draft verification, then all29 public streams.
        self.assertEqual(len(self.server.streams),87)
        self.assertEqual(sum(h.get('_public') is True for u,h in self.server.streams),29)
    def test_all27_imported_ids_can_change_without_stream_shortcut(self):
        self.controller.create()
        for name,row in self.server.files.items():row['file_id']='new-id-'+name;row['version_id']='new-version-'+name
        result=self.controller.prepare();self.assertEqual(len(self.server.streams),29)
        self.assertTrue(all(x['basis']=='FRESH_COMPLETE_CONTENT_STREAM' for x in result['content']))
    def test_public_prior_wrong_owner_blocks_every_write(self):
        original=self.server.document
        def altered(prior=False,public=False):
            obj=original(prior=prior,public=public)
            if prior:obj['parent']['access']['owned_by']['user']='999'
            return obj
        self.server.document=altered
        self.stop('CURRENT_PRIOR_OWNER_DIFFERS',self.controller.create)
        self.assertFalse(any(x[0]!='GET' for x in self.server.calls))
    def test_draft_wrong_owner_blocks_adoption(self):
        original=self.server.document
        def altered(prior=False,public=False):
            obj=original(prior=prior,public=public)
            if not prior:obj['parent']['access']['owned_by']['user']='999'
            return obj
        self.server.document=altered
        self.stop('DRAFT_OWNER_DIFFERS',self.controller.create)
        self.assertEqual(len(self.posts('/actions/newversion')),1)
        self.assertFalse(self.posts('/actions/publish'))
    def test_identity_only_edit_after_verify_blocks_publish_latch(self):
        for field in ('file_id','version_id'):
            self.prepared();name=self.baseline['files'][0]['filename']
            def change():self.server.files[name][field]='identity-only-change'
            self.edit_after_publish_verification(change)
            self.stop('DRAFT_CONTENT_IDENTITIES_CHANGED_SINCE_COMPLETE_VERIFICATION',self.controller.publish)
            self.assertFalse(self.posts('/actions/publish'))
            self.assertFalse(self.controller.state.get('publish_post_attempted',False))
            self.reset_synthetic_case('identity-only-'+field)
    def test_complete_stream_receipt_requires_both_identity_fields(self):
        self.controller.create();name=self.baseline['files'][0]['filename'];self.server.files[name]['version_id']=''
        self.stop('COMPLETE_REMOTE_FILE_AND_VERSION_IDS_REQUIRED',self.controller.prepare)
        self.assertFalse(self.posts('/actions/publish'))
        self.assertEqual(self.server.streams,[])
    def test_complete_identity_helper_rejects_absent_or_nonstring_ids(self):
        for field in ('file_id','version_id'):
            for value in ('',None,False,7):
                entry={'file_id':'valid-file','version_id':'valid-version'};entry[field]=value
                self.stop('COMPLETE_REMOTE_FILE_AND_VERSION_IDS_REQUIRED',lambda:M.immutable_file_identity(entry))

    def test_production_baseline_budget_pin(self):
        base=M.frozen_json(M.HERE/'read-only-observation/CURRENT_PUBLISHED_INVENTORY.json',M.BASELINE_SHA)
        self.assertEqual(base['expected_bytes'],448387915);self.assertEqual(len(base['files']),27)
        self.assertEqual(base['current_record'],23225288)
    def test_metadata_template_remains_nonsealed(self):
        template=json.loads((M.HERE/'metadata-candidates/PROSPECTIVE_METADATA_MODERN_NONSEALED.json').read_text())
        self.assertIn('PENDING',template['metadata']['version'])
        self.assertIn('UNCOMPUTED',template['metadata']['description'])
        self.assertEqual(template['metadata']['publication_date'],'2026-10-07')

if __name__=='__main__':
    names=[n for n in ContinuationControls.__dict__ if n.startswith('test_')]
    suite=unittest.TestSuite(ContinuationControls(n) for n in names)
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    report={'status':'PASS' if result.wasSuccessful() else 'FAIL','tests_run':result.testsRun,'skipped':len(result.skipped),'failures':[str(x[0]) for x in result.failures],'errors':[str(x[0]) for x in result.errors],'controller_sha256':M.sha((M.HERE/'new_edition_controller.py').read_bytes()),'suite_sha256':M.sha(Path(__file__).read_bytes()),'real_network':'FORBIDDEN_BY_SOCKET_SENTINEL','remote_mutations':0,'fixture_basis':'27smallmanufacturedfiles+2new; no retained scientific payloads'}
    (HERE/'CONTINUATION_CONTROLS_RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
    raise SystemExit(0 if result.wasSuccessful() else 1)
