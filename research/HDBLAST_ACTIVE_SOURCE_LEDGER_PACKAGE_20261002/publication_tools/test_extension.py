"""Fabricated publication guards: no network, scientific data or model execution."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
import zipfile

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).parent / 'helpers'))
import build_extension as builder
import prepare_candidate as prepare
import verify_saved_record as saved

TEMPLATE = Path(__file__).parent / 'TEMPLATE_SOURCE_FREE_MANIFEST.json'


def dig(raw, kind='sha256'):
    return hashlib.new(kind, raw).hexdigest()


def js(path, value):
    builder.json_write(path, value)
    return dig(path.read_bytes())


class ExtensionGuards(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='hdblast-derived-candidate-')
        self.root = Path(self.temp.name); self.repo = self.root / 'repo'; self.repo.mkdir()
        self.old = self.root / 'old'; self.old.mkdir(); self.out = self.root / 'derived'
        old = json.loads(TEMPLATE.read_text()); self.expected = {}
        for i, row in enumerate(old['new_files']):
            raw = f'fabricated prior addition {i}\n'.encode(); self.expected[row['filename']] = raw
            row.update(bytes=len(raw), sha256=dig(raw), md5=dig(raw, 'md5'))
            if row['source']['kind'] == 'candidate':
                (self.old / row['filename']).write_bytes(raw)
                row['source'] = {'kind': 'candidate', 'path': row['filename'], 'bytes': len(raw), 'sha256': dig(raw)}
            else:
                path = f'prior-{i}.bin'; (self.repo / path).write_bytes(raw)
                row['source'] = {'kind': 'repository', 'path': path, 'bytes': len(raw), 'sha256': dig(raw)}
        metadata = {'metadata': {'version': '2026.10.02-v25', 'upload_type': 'software',
            'license': 'cc-by-4.0', 'title': 'Fabricated historical candidate',
            'description': 'Fabricated guard only.', 'notes': 'Historical snapshot.',
            'creators': [{'name': 'Synthetic fixture'}], 'related_identifiers': []}}
        self.old_meta_pin = js(self.old / 'METADATA.json', metadata)
        old['metadata_request'] = {'path': 'METADATA.json', 'sha256': self.old_meta_pin}
        self.old_pin = js(self.old / 'FILE_MANIFEST.json', old)
        self.old_bytes = {p.name: p.read_bytes() for p in self.old.iterdir()}
        self.zip = self.root / 'active.zip'
        manifest = (json.dumps({'files': {'synthetic.txt': dig(b'No science\n')}}) + '\n').encode()
        with zipfile.ZipFile(self.zip, 'w') as archive:
            for name, raw in [('MANIFEST.json', manifest), ('synthetic.txt', b'No science\n')]:
                info = zipfile.ZipInfo(builder.CHECKPOINT_NAME + '/' + name)
                info.external_attr = 0o100644 << 16; archive.writestr(info, raw)
        self.report = self.root / 'report.md'; self.report.write_text('Fabricated LEDGER_ERROR_DEMONSTRATED; old FAIL.\n')
        self.fresh = self.root / 'fresh.json'
        self.fresh_data = {'status': 'PASS_FRESH_STANDALONE_ZIP_REPLAY',
            'classification': 'LEDGER_ERROR_DEMONSTRATED', 'old_metric_status': 'FAIL',
            'zip_sha256': prepare.sha256(self.zip), 'zip_bytes': self.zip.stat().st_size,
            'checks': {'primary_all_scientific_records_equal': True, 'independent_all_scientific_records_equal': True},
            'replay': {'status': 'COMPLETED_REPLAY', 'plan_only': False,
                'classification': 'LEDGER_ERROR_DEMONSTRATED', 'old_metric_status': 'FAIL',
                'physical_routes_completed': 2, 'successful_commands': 4,
                'source_verification_before': 'PASS', 'source_verification_after': 'PASS'}}
        js(self.fresh, self.fresh_data)
        self.proof = {'status': 'COMPLETE_ORIGINAL_AND_FRESH_REPLAY',
            'classification': 'LEDGER_ERROR_DEMONSTRATED', 'old_metric_status': 'FAIL',
            'scientific_comparison': 'EXACT_REGISTERED_FULL_FRAMES',
            'freeze_commit': 'a' * 40, 'science_commit': 'b' * 40,
            'registration_sha256': 'c' * 64, 'zip_sha256': prepare.sha256(self.zip),
            'package_manifest_sha256': dig(manifest), 'fresh_receipt_sha256': prepare.sha256(self.fresh),
            'result_report_sha256': prepare.sha256(self.report)}
        self.proof_path = self.root / 'completion.json'; js(self.proof_path, self.proof)
        self.target_paths = []
        for identifier, published in [(23112891, True), (23114217, False)]:
            path = self.root / ('published-get.json' if published else 'draft-get.json')
            url = 'https://zenodo.org/api/deposit/depositions/' + str(identifier)
            wrapper = {'method': 'GET', 'authenticated_request': True, 'http_status': 200,
                'requested_url': url, 'final_url': url, 'started_utc': '2026-10-03T03:55:00Z',
                'response_sha256': '7' * 64,
                'data': {'id': identifier, 'owner': 1386319, 'conceptdoi': '10.5281/zenodo.17088132',
                    'conceptrecid': '17088132', 'submitted': published,
                    'doi': '10.5281/zenodo.' + str(identifier) if published else None,
                    'state': 'done' if published else 'unsubmitted',
                    'metadata': {'title': 'Nonempty inherited title', 'description': 'Nonempty inherited description',
                        'creators': metadata['metadata']['creators'], 'license': 'cc-by-4.0'},
                    'files': [{'filename': f['filename'], 'filesize': f['bytes'], 'checksum': 'md5:' + f['md5']}
                              for f in old['inherited_files']]}}
            js(path, wrapper); self.target_paths.append(path)

    def tearDown(self):
        self.assertEqual(self.old_bytes, {p.name: p.read_bytes() for p in self.old.iterdir()})
        self.temp.cleanup()

    def build(self):
        return builder.build(self.old, self.old_pin, self.old_meta_pin,
            [self.zip, self.report, self.fresh], self.proof_path, '2026-10-03', self.out, *self.target_paths)

    def test_complete_build_snapshot_and_twelve_addition_preparation(self):
        manifest = self.build(); self.assertEqual(manifest['expected_final_files_count'], 22)
        self.assertEqual(manifest['new_files'][:9], json.loads((self.old / 'FILE_MANIFEST.json').read_text())['new_files'])
        receipt = prepare.prepare(self.repo, self.out, self.out / 'FILE_MANIFEST.json',
            prepare.sha256(self.out / 'FILE_MANIFEST.json'), self.root / 'upload')
        self.assertEqual(len(receipt['outputs']), 12)
        for row in manifest['new_files']:
            self.assertEqual(prepare.sha256(self.root / 'upload' / row['filename']), row['sha256'])
        self.assertEqual(manifest['inherited_files'], json.loads((self.old / 'FILE_MANIFEST.json').read_text())['inherited_files'])

    def test_wrong_original_pin_rejected(self):
        self.old_pin = '0' * 64
        with self.assertRaises(Exception): self.build()

    def test_verified_prepared_previous_files_reused_without_reassembly(self):
        manifest = self.build(); previous = self.root / 'previous'; previous.mkdir()
        for name, raw in self.expected.items(): (previous / name).write_bytes(raw)
        # Removing repository sources proves reuse does not rebuild old ZIPs.
        for path in self.repo.iterdir(): path.unlink()
        receipt = prepare.prepare(self.repo, self.out, self.out / 'FILE_MANIFEST.json',
            prepare.sha256(self.out / 'FILE_MANIFEST.json'), self.root / 'reused-upload', previous)
        self.assertEqual(len(receipt['outputs']), 12)
        (previous / manifest['new_files'][0]['filename']).write_bytes(b'mutated')
        with self.assertRaises(Exception): prepare.prepare(self.repo, self.out, self.out / 'FILE_MANIFEST.json',
            prepare.sha256(self.out / 'FILE_MANIFEST.json'), self.root / 'wrong-upload', previous)

    def test_missing_fresh_file_rejected(self):
        self.fresh.unlink()
        with self.assertRaises(Exception): self.build()

    def test_empty_fresh_receipt_rejected(self):
        js(self.fresh, {}); self.proof['fresh_receipt_sha256'] = prepare.sha256(self.fresh); js(self.proof_path, self.proof)
        with self.assertRaises(Exception): self.build()

    def test_incomplete_replay_or_false_science_comparison_rejected(self):
        for key, value in [('successful_commands', 3), ('plan_only', True),
                           ('physical_routes_completed', True), ('source_verification_after', 'FAIL')]:
            receipt = copy.deepcopy(self.fresh_data); receipt['replay'][key] = value
            with self.subTest(key=key), self.assertRaises(Exception): builder.completion_receipt(receipt, self.proof,
                {'sha256': self.proof['zip_sha256'], 'bytes': self.zip.stat().st_size})
        for checks in [{}, {'primary_all_scientific_records_equal': True},
                       {'primary_all_scientific_records_equal': True, 'independent_all_scientific_records_equal': False}]:
            receipt = dict(self.fresh_data, checks=checks)
            with self.assertRaises(Exception): builder.completion_receipt(receipt, self.proof,
                {'sha256': self.proof['zip_sha256'], 'bytes': self.zip.stat().st_size})

    def test_wrong_package_manifest_rejected(self):
        self.proof['package_manifest_sha256'] = '0' * 64; js(self.proof_path, self.proof)
        with self.assertRaises(Exception): self.build()

    def test_attachment_mutation_rejected(self):
        self.report.write_text('changed FAIL LEDGER_ERROR_DEMONSTRATED')
        with self.assertRaises(Exception): self.build()

    def test_previous_candidate_owned_attachment_mutation_rejected(self):
        manifest = json.loads((self.old / 'FILE_MANIFEST.json').read_text())
        for row in manifest['new_files']:
            if row['source']['kind'] == 'candidate':
                path = self.old / row['source']['path']; original = path.read_bytes()
                try:
                    path.write_bytes(b'changed historical overview')
                    with self.assertRaises(Exception): self.build()
                finally:
                    path.write_bytes(original)

    def test_output_collision_or_symlink_rejected(self):
        self.out.mkdir()
        with self.assertRaises(Exception): self.build()
        self.out.rmdir(); self.out.symlink_to(self.old, target_is_directory=True)
        with self.assertRaises(Exception): self.build()

    def test_duplicate_json_or_nonfinite_completion_rejected(self):
        for raw in ['{"status":"a","status":"b"}', '{"value":NaN}']:
            self.proof_path.write_text(raw)
            with self.assertRaises(Exception): self.build()

    def test_exact_inventory_mutation_guards(self):
        manifest = self.build()
        mutations = [('new_count', lambda m: m['new_files'].pop()),
            ('old_addition', lambda m: m['new_files'][0].update(md5='0' * 32)),
            ('inherited', lambda m: m['inherited_files'][0].update(bytes=1)),
            ('draft', lambda m: m.update(active_main_draft=23111008)),
            ('family', lambda m: m.update(concept_doi='10.5281/zenodo.22922927')),
            ('archiving', lambda m: m.update(automatic_github_archiving='ON')),
            ('duplicate', lambda m: m['new_files'][9].update(filename=m['new_files'][0]['filename'])),
            ('count_bool', lambda m: m.update(expected_final_files_count=True)),
            ('metric', lambda m: m.update(old_metric_scientific_status='PASS')),
            ('role', lambda m: m['new_files'][10].update(role='other')),
            ('uncomputed', lambda m: m['active_source_evidence'].update(status='UNCOMPUTED'))]
        for label, mutate in mutations:
            changed = copy.deepcopy(manifest); mutate(changed)
            with self.subTest(label=label), self.assertRaises(Exception): prepare.validate_extension_inventory(changed, self.out)

    def test_authenticated_retarget_owner_identity_inventory_and_lifecycle_guards(self):
        original = json.loads(self.target_paths[1].read_text())
        mutations = [lambda w: w.update(authenticated_request=False),
            lambda w: w.update(http_status=403), lambda w: w.update(final_url='https://example.invalid'),
            lambda w: w['data'].update(id=23112891), lambda w: w['data'].update(owner=1),
            lambda w: w['data'].update(conceptrecid='22922927'),
            lambda w: w['data'].update(submitted=True, state='done'),
            lambda w: w['data']['metadata'].update(description=''),
            lambda w: w['data']['metadata'].update(license='mit'),
            lambda w: w['data']['metadata'].update(creators=[{'name':'Other author'}]),
            lambda w: w['data']['files'][0].update(filesize=1),
            lambda w: w['data']['files'].pop()]
        for mutate in mutations:
            changed = copy.deepcopy(original); mutate(changed); js(self.target_paths[1], changed)
            with self.assertRaises(Exception): self.build()
            # Some failures occur after a fresh folder is created; remove only
            # this synthetic candidate before the next independently bad case.
            if self.out.exists():
                import shutil; shutil.rmtree(self.out)
        js(self.target_paths[1], original)

    def saved_fixture(self):
        manifest = self.build(); request = json.loads((self.out / 'METADATA.json').read_text())
        record = {'id': 23114217, 'conceptrecid': '17088132', 'conceptdoi': '10.5281/zenodo.17088132',
            'submitted': False, 'state': 'unsubmitted', 'metadata': copy.deepcopy(request['metadata']),
            'files': [{'filename': f['filename'], 'filesize': f['bytes'], 'checksum': 'md5:' + f['md5']}
                      for f in manifest['inherited_files'] + manifest['new_files']]}
        return manifest, request, record

    def test_saved_draft_and_published_record_require_full_twenty_two(self):
        manifest, request, record = self.saved_fixture()
        self.assertEqual(saved.verify(record, manifest, request, 23114217)['status'], 'PASS_COMPLETE_SAVED_DRAFT')
        clone = dict(record, id=99999999)
        self.assertEqual(saved.verify(clone, manifest, request, 99999999)['status'], 'FAIL_SAVED_RECORD_VALIDATION')
        self.assertEqual(saved.verify(dict(record, id=23114217.0), manifest, request, 23114217)['status'], 'FAIL_SAVED_RECORD_VALIDATION')
        aliases = copy.deepcopy(record); aliases['files'][0]['size'] = aliases['files'][0]['filesize'] + 1
        self.assertEqual(saved.verify(aliases, manifest, request, 23114217)['status'], 'FAIL_SAVED_RECORD_VALIDATION')
        published = dict(record, submitted=True, state='done', doi='10.5281/zenodo.23114217')
        self.assertEqual(saved.verify(published, manifest, request, 23114217, True)['status'], 'PASS_VERIFIED_PUBLISHED_RECORD')
        for i in range(22):
            for key, replacement in [('filesize', record['files'][i]['filesize'] + 1), ('checksum', 'md5:' + '0' * 32), ('filename', 'extra.bin')]:
                changed = copy.deepcopy(record); changed['files'][i][key] = replacement
                with self.subTest(index=i,key=key): self.assertEqual(saved.verify(changed, manifest, request, 23114217)['status'], 'FAIL_SAVED_RECORD_VALIDATION')
        for key, replacement in [('conceptrecid','22922927'), ('conceptdoi','10.5281/zenodo.22922927'),
                                 ('submitted',True), ('state','done'), ('id',23111008)]:
            self.assertEqual(saved.verify(dict(record, **{key:replacement}), manifest, request, 23114217)['status'], 'FAIL_SAVED_RECORD_VALIDATION')
        for key in request['metadata']:
            changed = copy.deepcopy(record); del changed['metadata'][key]
            self.assertEqual(saved.verify(changed, manifest, request, 23114217)['status'], 'FAIL_SAVED_RECORD_VALIDATION')


if __name__ == '__main__': unittest.main(verbosity=2)
