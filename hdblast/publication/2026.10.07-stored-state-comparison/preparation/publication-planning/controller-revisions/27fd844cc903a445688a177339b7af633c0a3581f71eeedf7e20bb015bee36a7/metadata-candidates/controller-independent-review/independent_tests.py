#!/usr/bin/env python3
"""Independent offline controls on a pinned controller snapshot; no real sockets."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import socket
import ssl
import sys
import tempfile
import unittest
from urllib import request

HERE = Path(__file__).resolve().parent
PLANNING = HERE.parents[1]
SNAPSHOT = HERE / 'candidate_f1a381d8_snapshot.py'
EXPECTED_SHA = 'f1a381d8d7c1eba81d1f6079647f6b63a8a86d9bccc9fb75cbadee490085dc98'
if hashlib.sha256(SNAPSHOT.read_bytes()).hexdigest() != EXPECTED_SHA:
    raise RuntimeError('OFFLINE_SOURCE_SNAPSHOT_PIN_DIFFERS')

def forbidden_network(*args, **kwargs):
    raise AssertionError('REAL_NETWORK_FORBIDDEN_IN_INDEPENDENT_REVIEW')

# ssl and urllib have been imported before substituting the socket constructor.
socket.socket = forbidden_network
socket.create_connection = forbidden_network
request.urlopen = forbidden_network
request.OpenerDirector.open = forbidden_network

def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    sys.modules[name] = result
    spec.loader.exec_module(result)
    return result

M = module(SNAPSHOT, 'independent_pinned_controller')
M.Transport = forbidden_network
F = module(PLANNING / 'controller-review/offline_tests.py', 'independent_fake_fixture_support')
F.M = M
F.SOURCE = SNAPSHOT
work = HERE / 'temporary-manufactured-evidence'
work.mkdir(exist_ok=True)
tempfile.tempdir = str(work)

class IndependentControls(F.ControllerFixtures):
    def validate_changed(self, metadata=None, inventory=None, reason=None):
        if metadata is not None:
            self.controller.metadata = metadata
        if inventory is not None:
            self.controller.inventory = inventory
        self.stop(reason, self.controller.validate_inputs)
        self.assertFalse(any(method != 'GET' for method, _, _ in self.server.calls))

    def test_independent_missing_replay_marker_rejects_before_any_action(self):
        inventory = copy.deepcopy(self.inventory)
        inventory.pop('scientific_replay_status')
        self.validate_changed(inventory=inventory, reason='COMPLETED_SCIENTIFIC_REPLAY_REQUIRED')

    def test_independent_wrong_replay_marker_rejects(self):
        inventory = copy.deepcopy(self.inventory)
        inventory['scientific_replay_status'] = 'PASS_FABRICATED_ONLY'
        self.validate_changed(inventory=inventory, reason='COMPLETED_SCIENTIFIC_REPLAY_REQUIRED')

    def test_independent_creator_change_rejects(self):
        metadata = copy.deepcopy(self.metadata)
        metadata['metadata']['creators'][0]['person_or_org']['name'] = 'Changed Author'
        self.validate_changed(metadata=metadata, reason='INHERITED_METADATA_CREATORS_DIFFERS')

    def test_independent_rights_change_rejects(self):
        metadata = copy.deepcopy(self.metadata)
        metadata['metadata']['rights'][0]['id'] = 'cc0-1.0'
        self.validate_changed(metadata=metadata, reason='INHERITED_METADATA_RIGHTS_DIFFERS')

    def test_independent_old_related_identifier_change_rejects(self):
        metadata = copy.deepcopy(self.metadata)
        metadata['metadata']['related_identifiers'][0]['identifier'] = '10.5281/zenodo.999'
        self.validate_changed(metadata=metadata, reason='INHERITED_REFERENCES_DIFFERS')

    def test_independent_historical_note_change_rejects(self):
        metadata = copy.deepcopy(self.metadata)
        metadata['metadata']['additional_descriptions'][1]['description'] = '<p>Altered historical result.</p>'
        self.validate_changed(metadata=metadata, reason='HISTORICAL_NOTES_DIFFERS')

    def test_independent_first_note_language_change_rejects(self):
        metadata = copy.deepcopy(self.metadata)
        metadata['metadata']['additional_descriptions'][0]['lang']['id'] = 'spa'
        self.validate_changed(metadata=metadata, reason='CURRENT_NOTES_TYPE_OR_LANGUAGE_DIFFERS')

    def test_independent_unreviewed_related_identifier_append_rejects(self):
        metadata = copy.deepcopy(self.metadata)
        metadata['metadata']['related_identifiers'].append({'identifier':'10.5281/zenodo.999',
            'scheme':'doi','relation_type':{'id':'references'},'resource_type':{'id':'software'}})
        self.validate_changed(metadata=metadata, reason='UNREVIEWED_REFERENCE_ADDITION')

    def test_independent_exact_prior_reference_append_accepts(self):
        self.controller.metadata['metadata']['related_identifiers'].append({
            'identifier':'10.5281/zenodo.' + M.PRIOR_ID, 'scheme':'doi',
            'relation_type':{'id':'references'}, 'resource_type':{'id':'software'}})
        self.controller.validate_inputs()
        self.assertFalse(self.server.calls)

    def test_independent_access_change_rejects(self):
        metadata = copy.deepcopy(self.metadata)
        metadata['access']['files'] = 'restricted'
        self.validate_changed(metadata=metadata, reason='INHERITED_ACCESS_DIFFERS')

    def test_independent_custom_field_change_rejects(self):
        metadata = copy.deepcopy(self.metadata)
        metadata['custom_fields']['unreviewed'] = True
        self.validate_changed(metadata=metadata, reason='INHERITED_CUSTOM_FIELDS_DIFFERS')

    def test_independent_stale_state_binding_rejects_on_restart(self):
        path = self.root / 'evidence/STATE.json'
        value = json.loads(path.read_text())
        value['continuation_binding']['prior_record_id'] = '23114217'
        path.write_text(json.dumps(value))
        self.stop('JOURNAL_BASELINE_OR_CONTROLLER_BINDING_DIFFERS', self.make_controller)
        self.assertFalse(self.server.calls)

    def test_independent_orphan_journal_rejects_on_restart(self):
        (self.root / 'evidence/STATE.json').unlink()
        (self.root / 'evidence/JOURNAL.jsonl').write_text('{}\n')
        self.stop('ORPHAN_JOURNAL_REQUIRES_MANUAL_REVIEW', self.make_controller)
        self.assertFalse(self.server.calls)

    def test_independent_current_prior_owner_change_rejects_pre_action(self):
        document = self.server.document
        def wrong_owner(prior=False, public=False):
            result = document(prior=prior, public=public)
            if prior:
                result['parent']['access']['owned_by']['user'] = str(M.EXPECTED_OWNER + 1)
            return result
        self.server.document = wrong_owner
        self.stop('CURRENT_PRIOR_OWNER_DIFFERS', self.controller.create)
        self.assertFalse(any(method != 'GET' for method, _, _ in self.server.calls))

    def test_independent_wrong_draft_owner_cannot_be_adopted(self):
        document = self.server.document
        def wrong_owner(prior=False, public=False):
            result = document(prior=prior, public=public)
            if not prior:
                result['parent']['access']['owned_by']['user'] = str(M.EXPECTED_OWNER + 1)
            return result
        self.server.document = wrong_owner
        self.stop('DRAFT_OWNER_DIFFERS', self.controller.create)
        self.assertIsNone(self.controller.state['draft_id'])
        self.assertFalse(self.posts('/actions/publish'))

    def test_independent_same_identity_inherited_bad_bytes_block_publish(self):
        self.prepared()
        name = self.baseline['files'][0]['filename']
        self.server.payloads[name] = b'X' * len(self.server.payloads[name])
        self.stop('CONTENT_COMPLETE_MD5_SHA256_DIFFERS', self.controller.publish)
        self.assertFalse(self.posts('/actions/publish'))
        self.assertFalse(self.controller.state.get('publish_post_attempted', False))

    def test_independent_all29_fresh_streams_both_states(self):
        receipt = self.prepared()
        self.assertEqual(len(receipt['content']), 29)
        self.assertTrue(all(x['basis'] == 'FRESH_COMPLETE_CONTENT_STREAM' for x in receipt['content']))
        before = len(self.server.streams)
        public = self.controller.publish()
        self.assertEqual(len(public['content']), 29)
        self.assertEqual(len(self.server.streams) - before, 58)
        self.assertTrue(all(headers['_public'] is True for _, headers in self.server.streams[-29:]))

    def test_independent_raw_etag_change_stops_before_decimal_conversion(self):
        self.prepared()
        self.edit_after_publish_verification(lambda: setattr(self.server, 'draft_etag', '"07"'))
        self.stop('DRAFT_CHANGED_SINCE_COMPLETE_VERIFICATION', self.controller.publish)
        self.assertFalse(self.controller.state.get('publish_post_attempted', False))
        self.assertFalse(self.posts('/actions/publish'))

    def test_independent_final_file_identity_change_must_block_publish(self):
        """Critical guard reproduction: unchanged bytes, MD5, metadata and raw ETag."""
        self.prepared()
        name = self.baseline['files'][0]['filename']
        def change():
            self.server.files[name]['file_id'] = 'replacement-file-identity'
            self.server.files[name]['version_id'] = 'replacement-content-version'
        self.edit_after_publish_verification(change)
        with self.assertRaises(M.Stop, msg='FINAL_IDENTITY_CHANGE_WAS_NOT_BLOCKED'):
            self.controller.publish()
        self.assertFalse(self.controller.state.get('publish_post_attempted', False))
        self.assertFalse(self.posts('/actions/publish'))

if __name__ == '__main__':
    names = sorted(name for name in IndependentControls.__dict__ if name.startswith('test_independent_'))
    suite = unittest.TestSuite(IndependentControls(name) for name in names)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    report = {
        'status':'PASS' if result.wasSuccessful() else 'FAIL',
        'tests_run':result.testsRun,
        'failures':[{'test':str(test),'traceback':trace} for test,trace in result.failures],
        'errors':[{'test':str(test),'traceback':trace} for test,trace in result.errors],
        'skipped':list(result.skipped),
        'controller_snapshot_sha256':EXPECTED_SHA,
        'suite_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'fixture_support_sha256':hashlib.sha256((PLANNING/'controller-review/offline_tests.py').read_bytes()).hexdigest(),
        'real_network':'FORBIDDEN_BY_SOCKET_AND_URLLIB_SENTINELS',
        'remote_mutations':0,
        'credential_access':'NONE_FAKE_CONSTANT_ONLY',
        'scientific_callbacks':'NOT_IMPORTED_OR_CALLED',
        'real_study_status':'UNCOMPUTED',
        'release_assets_and_seal':'PENDING',
    }
    mode = 'optimized' if sys.flags.optimize else 'normal'
    (HERE / ('INDEPENDENT_RESULTS_' + mode + '.json')).write_text(json.dumps(report,indent=2)+'\n')
    raise SystemExit(0 if result.wasSuccessful() else 1)
