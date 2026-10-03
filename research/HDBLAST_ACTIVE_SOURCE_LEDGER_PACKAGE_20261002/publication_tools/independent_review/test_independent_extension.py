"""Independent fabricated extension guards; no real candidate or science execution."""
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch
import warnings
import zipfile

EXTENSION = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(EXTENSION))
sys.path.insert(0, str(EXTENSION / 'helpers'))
import build_extension as builder
import prepare_candidate as prepare
import verify_saved_record as saved
import extension_inventory as inventory
import test_extension as existing_tests


class IndependentExtensionGuards(unittest.TestCase):
    def setUp(self):
        self.fixture = existing_tests.ExtensionGuards(methodName='test_wrong_original_pin_rejected')
        self.fixture.setUp()

    def tearDown(self):
        self.fixture.tearDown()

    def test_full_build_preserves_parent_bytes_and_ordered_prefix(self):
        manifest = self.fixture.build()
        parent = json.loads((self.fixture.old / 'FILE_MANIFEST.json').read_text())
        self.assertEqual(manifest['inherited_files'], parent['inherited_files'])
        self.assertEqual(manifest['new_files'][:9], parent['new_files'])
        self.assertEqual(len(manifest['inherited_files'] + manifest['new_files']), 22)
        self.assertEqual(len({r['filename'] for r in manifest['inherited_files'] + manifest['new_files']}), 22)
        self.assertEqual(manifest['active_main_draft'], inventory.CURRENT_MAIN_DRAFT)
        self.assertEqual(manifest['latest_published_record'], inventory.CURRENT_PUBLISHED_RECORD)
        self.assertEqual(manifest['concept_doi'], '10.5281/zenodo.17088132')
        self.assertEqual(manifest['automatic_github_archiving'], 'ALL_EIGHT_OFF')
        for name in ('FILE_MANIFEST.json', 'METADATA.json'):
            self.assertEqual((self.fixture.out / 'source_free_snapshot' / name).read_bytes(),
                             (self.fixture.old / name).read_bytes())
        for row in manifest['new_files'][9:]:
            payload = (self.fixture.out / row['source']['path']).read_bytes()
            self.assertEqual(len(payload), row['bytes'])
            self.assertEqual(hashlib.sha256(payload).hexdigest(), row['sha256'])
            self.assertEqual(hashlib.md5(payload).hexdigest(), row['md5'])

    def test_corrupted_prior_candidate_source_cannot_emit_success_receipt(self):
        parent = json.loads((self.fixture.old / 'FILE_MANIFEST.json').read_text())
        row = next(row for row in parent['new_files'] if row['source']['kind'] == 'candidate')
        source = self.fixture.old / row['source']['path']
        original = source.read_bytes()
        try:
            source.write_bytes(b'X' * len(original))
            with self.assertRaises(RuntimeError):
                self.fixture.build()
            self.assertFalse((self.fixture.out / 'BUILD_RECEIPT.json').exists())
        finally:
            source.write_bytes(original)

    def test_corrupted_prior_candidate_copy_cannot_emit_success_receipt(self):
        parent = json.loads((self.fixture.old / 'FILE_MANIFEST.json').read_text())
        row = next(row for row in parent['new_files'] if row['source']['kind'] == 'candidate')
        original_copy = builder.shutil.copyfile
        def corrupted_copy(source, destination, *args, **kwargs):
            result = original_copy(source, destination, *args, **kwargs)
            if Path(source) == self.fixture.old / row['source']['path']:
                Path(destination).write_bytes(b'X' * row['bytes'])
            return result
        with patch('build_extension.shutil.copyfile', side_effect=corrupted_copy):
            with self.assertRaises(RuntimeError):
                self.fixture.build()
        self.assertFalse((self.fixture.out / 'BUILD_RECEIPT.json').exists())

    def test_every_required_fresh_receipt_execution_field_is_checked(self):
        fields = [
            ((), 'status', 'FAIL_FRESH_STANDALONE_ZIP_REPLAY'),
            ((), 'classification', 'CONSISTENCY_FAILURE'),
            ((), 'old_metric_status', 'PASS'),
            ((), 'zip_sha256', '0' * 64),
            ((), 'zip_bytes', True),
            (('replay',), 'status', 'PLANNED_REPLAY'),
            (('replay',), 'plan_only', True),
            (('replay',), 'physical_routes_completed', 1),
            (('replay',), 'successful_commands', 4.0),
            (('replay',), 'source_verification_before', 'FAIL'),
            (('replay',), 'source_verification_after', 'FAIL'),
            (('replay',), 'classification', 'CONSISTENCY_FAILURE'),
            (('replay',), 'old_metric_status', 'PASS'),
        ]
        zip_record = {'sha256': self.fixture.proof['zip_sha256'], 'bytes': self.fixture.zip.stat().st_size}
        for path, key, wrong in fields:
            for remove in (False, True):
                receipt = copy.deepcopy(self.fixture.fresh_data)
                target = receipt
                for part in path:
                    target = target[part]
                if remove:
                    del target[key]
                else:
                    target[key] = wrong
                with self.subTest(path=path, key=key, removed=remove), self.assertRaises(RuntimeError):
                    builder.completion_receipt(receipt, self.fixture.proof, zip_record)
        for checks in ({}, {'primary_all_scientific_records_equal': True},
                       {'primary_all_scientific_records_equal': 1, 'independent_all_scientific_records_equal': True},
                       {'primary_all_scientific_records_equal': True, 'independent_all_scientific_records_equal': True, 'another_check': False}):
            with self.subTest(checks=checks), self.assertRaises(RuntimeError):
                builder.completion_receipt(dict(self.fixture.fresh_data, checks=checks), self.fixture.proof, zip_record)

    def test_all_completion_evidence_pins_and_prefix_rows_are_guarded(self):
        manifest = self.fixture.build()
        for key in ('freeze_commit', 'science_commit', 'registration_sha256', 'zip_sha256',
                    'package_manifest_sha256', 'fresh_receipt_sha256', 'result_report_sha256'):
            changed = copy.deepcopy(manifest)
            changed['active_source_evidence'][key] = 'main'
            with self.subTest(key=key), self.assertRaises(RuntimeError):
                prepare.validate_extension_inventory(changed, self.fixture.out)
        for index in range(9):
            changed = copy.deepcopy(manifest)
            changed['new_files'][index]['role'] += '_changed'
            with self.subTest(prefix_row=index), self.assertRaises(RuntimeError):
                prepare.validate_extension_inventory(changed, self.fixture.out)

    def test_saved_draft_binding_and_published_selected_record(self):
        manifest, request, response = self.fixture.saved_fixture()
        clone = copy.deepcopy(response)
        clone['id'] = 99999999
        self.assertEqual(saved.verify(clone, manifest, request, 99999999)['status'], 'FAIL_SAVED_RECORD_VALIDATION')
        published = copy.deepcopy(response)
        published.update(id=99999999, submitted=True, state='done', doi='10.5281/zenodo.99999999')
        self.assertEqual(saved.verify(published, manifest, request, 99999999, published=True)['status'],
                         'PASS_VERIFIED_PUBLISHED_RECORD')
        for key, value in (('submitted', False), ('state', 'unsubmitted'), ('doi', '10.5281/zenodo.23112891')):
            changed = dict(published, **{key: value})
            with self.subTest(key=key):
                self.assertEqual(saved.verify(changed, manifest, request, 99999999, published=True)['status'],
                                 'FAIL_SAVED_RECORD_VALIDATION')

    def test_ambiguous_saved_ids_and_conflicting_size_aliases_rejected(self):
        manifest, request, response = self.fixture.saved_fixture()
        changed = copy.deepcopy(response)
        changed['id'] = float(inventory.CURRENT_MAIN_DRAFT)
        self.assertEqual(saved.verify(changed, manifest, request, inventory.CURRENT_MAIN_DRAFT)['status'], 'FAIL_SAVED_RECORD_VALIDATION')
        for index in range(22):
            changed = copy.deepcopy(response)
            changed['files'][index]['size'] = changed['files'][index]['filesize'] + 1
            with self.subTest(index=index):
                self.assertEqual(saved.verify(changed, manifest, request, inventory.CURRENT_MAIN_DRAFT)['status'], 'FAIL_SAVED_RECORD_VALIDATION')

    def test_saved_cli_rejects_duplicate_json_keys(self):
        _, _, response = self.fixture.saved_fixture()
        selected = inventory.CURRENT_MAIN_DRAFT
        serialized = json.dumps(response).replace('"id": ' + str(selected), '"id": 23111008, "id": ' + str(selected), 1)
        path = self.fixture.root / 'duplicate-saved-response.json'
        path.write_text(serialized)
        result = subprocess.run([sys.executable, '-B', str(self.fixture.out / 'verify_saved_record.py'),
            '--saved-response', str(path), '--expected-record-id', str(selected),
            '--expected-manifest-sha256', prepare.sha256(self.fixture.out / 'FILE_MANIFEST.json'),
            '--output', str(self.fixture.root / 'saved-result.json')], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn('PASS_COMPLETE_SAVED_DRAFT', result.stdout)

    def test_reused_previous_files_each_require_actual_pin(self):
        manifest = self.fixture.build()
        previous = self.fixture.root / 'previous'
        previous.mkdir()
        for name, raw in self.fixture.expected.items():
            (previous / name).write_bytes(raw)
        for source in self.fixture.repo.iterdir():
            source.unlink()
        for index, row in enumerate(manifest['new_files'][:9]):
            path = previous / row['filename']
            original = path.read_bytes()
            try:
                path.write_bytes(b'X' * len(original))
                with self.subTest(index=index), self.assertRaises(RuntimeError):
                    prepare.prepare(self.fixture.repo, self.fixture.out, self.fixture.out / 'FILE_MANIFEST.json',
                        prepare.sha256(self.fixture.out / 'FILE_MANIFEST.json'), self.fixture.root / ('bad-reuse-' + str(index)), previous)
            finally:
                path.write_bytes(original)
        receipt = prepare.prepare(self.fixture.repo, self.fixture.out, self.fixture.out / 'FILE_MANIFEST.json',
            prepare.sha256(self.fixture.out / 'FILE_MANIFEST.json'), self.fixture.root / 'good-reuse', previous)
        self.assertEqual(len(receipt['outputs']), 12)
        self.assertEqual([item['filename'] for item in receipt['outputs']], [item['filename'] for item in manifest['new_files']])

    def test_output_inside_old_candidate_rejected(self):
        self.fixture.out = self.fixture.old / 'nested-derived-output'
        with self.assertRaises(RuntimeError):
            self.fixture.build()
        self.assertFalse(self.fixture.out.exists())

    def test_each_authenticated_get_rejects_missing_or_wrong_required_fields(self):
        inherited = json.loads((self.fixture.old / 'FILE_MANIFEST.json').read_text())['inherited_files']
        for path, published in zip(self.fixture.target_paths, (True, False)):
            original = inventory.read_json(path)
            fields = [
                ((), 'method', 'POST'), ((), 'authenticated_request', False),
                ((), 'http_status', True), ((), 'requested_url', 'https://zenodo.org/api/deposit/depositions/1'),
                ((), 'final_url', 'https://zenodo.org/api/deposit/depositions/1'),
                ((), 'response_sha256', 'invalid'), ((), 'started_utc', ''),
                (('data',), 'id', 23111008), (('data',), 'owner', True),
                (('data',), 'conceptdoi', '10.5281/zenodo.22922927'),
                (('data',), 'conceptrecid', '22922927'),
                (('data',), 'submitted', not published),
                (('data',), 'state', 'unsubmitted' if published else 'done'),
                (('data', 'metadata'), 'creators', []), (('data', 'metadata'), 'license', ''),
            ]
            for parent_path, key, wrong in fields:
                for removed in (False, True):
                    wrapper = copy.deepcopy(original)
                    target = wrapper
                    for part in parent_path:
                        target = target[part]
                    if removed:
                        del target[key]
                    else:
                        target[key] = wrong
                    with self.subTest(published=published, key=key, removed=removed), self.assertRaises(RuntimeError):
                        inventory.authenticated_target(wrapper, '8' * 64, inherited, published)

    def test_authenticated_inventory_rejects_all_contradictory_aliases(self):
        inherited = json.loads((self.fixture.old / 'FILE_MANIFEST.json').read_text())['inherited_files']
        for path, published in zip(self.fixture.target_paths, (True, False)):
            original = inventory.read_json(path)
            for index in range(10):
                for alias in ('bytes', 'md5'):
                    wrapper = copy.deepcopy(original)
                    row = wrapper['data']['files'][index]
                    row[alias] = row['filesize'] + 1 if alias == 'bytes' else '0' * 32
                    with self.subTest(published=published, index=index, alias=alias), self.assertRaises(RuntimeError):
                        inventory.authenticated_target(wrapper, '8' * 64, inherited, published)

    def test_projected_retarget_receipt_rechecks_endpoint_and_timestamp(self):
        manifest = self.fixture.build()
        metadata = inventory.read_json(self.fixture.out / 'METADATA.json')['metadata']
        for route in ('published', 'draft'):
            for field, wrong in (('requested_url', 'https://zenodo.org/api/deposit/depositions/1'),
                                 ('recorded_utc', '')):
                for removed in (False, True):
                    receipt = copy.deepcopy(manifest['target_reconciliation'])
                    if removed:
                        del receipt[route][field]
                    else:
                        receipt[route][field] = wrong
                    with self.subTest(route=route, field=field, removed=removed), self.assertRaises(RuntimeError):
                        inventory.validate_target_receipt(receipt, manifest['inherited_files'], metadata)

    def test_current_published_doi_is_checked_in_wrapper_and_receipt(self):
        inherited = json.loads((self.fixture.old / 'FILE_MANIFEST.json').read_text())['inherited_files']
        wrapper = inventory.read_json(self.fixture.target_paths[0])
        for wrong in (None, '', '10.5281/zenodo.22347452', '10.5281/zenodo.23114217'):
            changed = copy.deepcopy(wrapper)
            changed['data']['doi'] = wrong
            with self.subTest(wrapper_doi=wrong), self.assertRaises(RuntimeError):
                inventory.authenticated_target(changed, '8' * 64, inherited, True)
        manifest = self.fixture.build()
        metadata = inventory.read_json(self.fixture.out / 'METADATA.json')['metadata']
        for wrong in (None, '', '10.5281/zenodo.22347452', '10.5281/zenodo.23114217'):
            receipt = copy.deepcopy(manifest['target_reconciliation'])
            receipt['published']['doi'] = wrong
            with self.subTest(projected_doi=wrong), self.assertRaises(RuntimeError):
                inventory.validate_target_receipt(receipt, inherited, metadata)

    def test_saved_twenty_two_reject_all_size_and_digest_alias_conflicts(self):
        manifest, request, response = self.fixture.saved_fixture()
        for index in range(22):
            for alias in ('bytes', 'md5'):
                changed = copy.deepcopy(response)
                row = changed['files'][index]
                row[alias] = row['filesize'] + 1 if alias == 'bytes' else '0' * 32
                with self.subTest(index=index, alias=alias):
                    self.assertEqual(saved.verify(changed, manifest, request, inventory.CURRENT_MAIN_DRAFT)['status'],
                                     'FAIL_SAVED_RECORD_VALIDATION')
        matching = copy.deepcopy(response)
        for row in matching['files']:
            row['bytes'] = row['size'] = row['filesize']
            row['md5'] = row['checksum'].removeprefix('md5:')
        self.assertEqual(saved.verify(matching, manifest, request, inventory.CURRENT_MAIN_DRAFT)['status'],
                         'PASS_COMPLETE_SAVED_DRAFT')

    def test_zip_root_full_membership_payloads_and_modes_are_authenticated(self):
        prefix = builder.CHECKPOINT_NAME + '/'
        payload = b'fabricated payload only\n'
        digest = hashlib.sha256(payload).hexdigest()
        manifest = json.dumps({'files': {'sample.txt': digest}}).encode()
        cases = [
            ('valid', [(prefix + 'MANIFEST.json', manifest, 0o100644),
                       (prefix + 'sample.txt', payload, 0o100644)], None),
            ('flat', [('MANIFEST.json', manifest, 0o100644), ('sample.txt', payload, 0o100644)], 'manifest'),
            ('extra', [(prefix + 'MANIFEST.json', manifest, 0o100644),
                       (prefix + 'sample.txt', payload, 0o100644),
                       (prefix + 'extra.txt', b'extra', 0o100644)], 'membership'),
            ('missing', [(prefix + 'MANIFEST.json', manifest, 0o100644)], 'membership'),
            ('corrupted', [(prefix + 'MANIFEST.json', manifest, 0o100644),
                           (prefix + 'sample.txt', b'changed payload', 0o100644)], 'payload bytes'),
            ('symlink', [(prefix + 'MANIFEST.json', manifest, 0o100644),
                         (prefix + 'sample.txt', payload, 0o120777)], 'nonregular'),
            ('directory', [(prefix + 'MANIFEST.json', manifest, 0o100644),
                           (prefix + 'sample.txt', payload, 0o040755)], 'nonregular'),
            ('outside', [(prefix + 'MANIFEST.json', manifest, 0o100644),
                         ('OTHER/sample.txt', payload, 0o100644)], 'Unsafe'),
            ('escape', [(prefix + 'MANIFEST.json', manifest, 0o100644),
                        (prefix + '../sample.txt', payload, 0o100644)], 'Unsafe'),
            ('duplicate', [(prefix + 'MANIFEST.json', manifest, 0o100644),
                           (prefix + 'sample.txt', payload, 0o100644),
                           (prefix + 'sample.txt', payload, 0o100644)], 'Duplicate'),
            ('duplicate_json', [(prefix + 'MANIFEST.json', b'{"files":{},"files":{}}', 0o100644)], 'Duplicate package JSON'),
        ]
        for label, entries, error in cases:
            path = self.fixture.root / ('package-' + label + '.zip')
            manifest_payload = next(raw for name, raw, mode in entries if name.endswith('MANIFEST.json'))
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', UserWarning)
                with zipfile.ZipFile(path, 'w') as archive:
                    for name, raw, mode in entries:
                        info = zipfile.ZipInfo(name)
                        info.external_attr = mode << 16
                        archive.writestr(info, raw)
            with self.subTest(label=label):
                if error is None:
                    self.assertEqual(builder.package_metadata(path, hashlib.sha256(manifest_payload).hexdigest()), 1)
                else:
                    with self.assertRaisesRegex(RuntimeError, error):
                        builder.package_metadata(path, hashlib.sha256(manifest_payload).hexdigest())
        valid_bytes = (self.fixture.root / 'package-valid.zip').read_bytes()
        nul_path = self.fixture.root / 'package-nul.zip'
        nul_path.write_bytes(valid_bytes.replace((prefix + 'sample.txt').encode(),
                                                (prefix + 'sample.\x00xt').encode()))
        with self.assertRaisesRegex(RuntimeError, 'Unsafe'):
            builder.package_metadata(nul_path, hashlib.sha256(manifest).hexdigest())
        encrypted_bytes = bytearray(valid_bytes)
        for signature, flag_offset in ((b'PK\x03\x04', 6), (b'PK\x01\x02', 8)):
            start = 0
            while True:
                offset = valid_bytes.find(signature, start)
                if offset < 0:
                    break
                encrypted_bytes[offset + flag_offset] |= 1
                start = offset + 4
        encrypted_path = self.fixture.root / 'package-encrypted.zip'
        encrypted_path.write_bytes(encrypted_bytes)
        with self.assertRaisesRegex(RuntimeError, 'Unsafe'):
            builder.package_metadata(encrypted_path, hashlib.sha256(manifest).hexdigest())


if __name__ == '__main__':
    unittest.main(verbosity=2)
