#!/usr/bin/env python3
"""Manufactured new-record controls with authenticated dated publication seeds.

No fixture is a publication assertion. No network, numerical library, original
array, source or target is used. The expected verifier/seed hashes are external.
"""
from pathlib import Path
import argparse
import copy
import hashlib
import importlib.util
import json
import re
import socket
import sys
import tempfile

sys.dont_write_bytecode = True
HERE = Path(__file__).absolute().parent
ID = '90000001'
PREFIX = 'https://zenodo.org/api/records/' + ID


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False,
                      allow_nan=False).encode()


def no_network(*args, **kwargs):
    raise RuntimeError('Network forbidden in manufactured witness controls')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verifier-sha256', required=True)
    parser.add_argument('--seed-manifest-sha256', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    verifier = HERE / 'verify_publication.py'
    require(re.fullmatch('[0-9a-f]{64}', args.verifier_sha256 or '') and
            sha(verifier.read_bytes()) == args.verifier_sha256, 'External verifier source pin differs')
    seed_root = HERE / 'manufactured-seeds'
    seed_manifest_raw = (seed_root / 'SEED_MANIFEST.json').read_bytes()
    require(re.fullmatch('[0-9a-f]{64}', args.seed_manifest_sha256 or '') and
            sha(seed_manifest_raw) == args.seed_manifest_sha256, 'External seed manifest pin differs')
    seeds = {}
    for name, pin in json.loads(seed_manifest_raw)['files'].items():
        require(Path(name).name == name and '/' not in name and '\\' not in name, 'Unsafe seed filename')
        path = seed_root / name
        require(path.is_file() and all(not p.is_symlink() for p in (path, *path.parents)), 'Real seed file required')
        raw = path.read_bytes()
        require(len(raw) == pin['bytes'] and sha(raw) == pin['sha256'], 'Seed content pin differs')
        seeds[name] = raw
    spec = importlib.util.spec_from_file_location('manufactured_29_witness_verifier', verifier)
    v = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(v)
    socket.socket = no_network
    socket.create_connection = no_network
    prior_raw = seeds['PRIOR27_CONTENT_RECEIPT.json']
    require(sha(prior_raw) == v.PRIOR_STREAM_PIN, 'Dated prior receipt pin differs')
    prior = v.load(prior_raw)
    old = [dict(filename=row['filename'], **row['observed']) for row in prior['rows']]
    new = []
    for name, raw in [('FABRICATED_COMPARISON_CONTROL.zip', b'fabricated payload'),
                      ('FABRICATED_COMPARISON_CONTROL.md', b'fabricated overview')]:
        new.append({'filename': name, 'bytes': len(raw), 'md5': hashlib.md5(raw).hexdigest(), 'sha256': sha(raw)})
    inventory = {'sealed_upload_manifest': True, 'concept_doi': '10.5281/zenodo.' + v.PARENT,
                 'scientific_replay_status': v.SCIENTIFIC_REPLAY_STATUS,
                 'inherited_files': old, 'new_files': new, 'expected_total_files_final': 29,
                 'total_bytes_final': sum(row['bytes'] for row in old + new)}
    seed_prior = v.load(seeds['23225288-record.json'])
    metadata = {key: copy.deepcopy(seed_prior[key]) for key in ('metadata', 'custom_fields', 'access')}
    metadata['access'].pop('status', None)
    metadata['metadata']['title'] = 'Explicitly fabricated 29-file verifier fixture'
    metadata['metadata']['version'] = 'FABRICATED_ONLY_29_FILE_CONTROL'
    inventory_raw = canonical(inventory); metadata_raw = canonical(metadata)
    inventory_sha = sha(inventory_raw); metadata_sha = sha(metadata_raw)

    def fixture():
        entries = []
        for index, pin in enumerate(old + new):
            name = pin['filename']
            entries.append({'key': name, 'size': pin['bytes'], 'checksum': 'md5:' + pin['md5'],
                            'file_id': 'fabricated-new-file-' + str(index),
                            'version_id': 'fabricated-new-version-' + str(index),
                            'status': 'completed',
                            'links': {'content': PREFIX + '/files/' + v.quote(name, safe='') + '/content'}})
        # Genuine public record embeddings expose id and may omit version_id.
        embedded = {item['key']: {key: copy.deepcopy(value) for key, value in item.items()
                                 if key not in ('file_id', 'version_id')}
                    for item in entries}
        for item in entries:
            embedded[item['key']]['id'] = item['file_id']
        record = dict(copy.deepcopy(metadata), id=ID, parent={'id': v.PARENT},
                      is_published=True, is_draft=False,
                      pids={'doi': {'identifier': '10.5281/zenodo.' + ID}},
                      files={'entries': embedded}, links={'self': PREFIX, 'files': PREFIX + '/files'})
        files = {'entries': copy.deepcopy(entries)}
        verified = {'status': 'PASS_COMPLETE_29_FILE_EDITION_READBACK', 'published': True,
                    'draft_id': ID, 'parent_id': v.PARENT, 'inventory_sha256': inventory_sha,
                    'metadata_sha256': metadata_sha, 'file_count': 29,
                    'total_bytes': inventory['total_bytes_final'], 'content': []}
        journal = [{'kind': 'GET', 'url': PREFIX + '/files', 'status': 200,
                    'authenticated_request': False}]
        for pin, item in zip(old + new, entries):
            verified['content'].append({'filename': pin['filename'], 'sha256': pin['sha256'],
                'basis': 'FRESH_COMPLETE_CONTENT_STREAM', 'remote_file_id': item['file_id'],
                'remote_version_id': item['version_id']})
            journal.append({'kind': 'FULL_CONTENT_STREAM_VERIFIED', 'filename': pin['filename'],
                            'observed': {key: pin[key] for key in ('bytes', 'md5', 'sha256')}})
        preserved = {}
        for rid in v.PRESERVED:
            source_record = v.load(seeds[rid + '-record.json'])
            source_files = v.load(seeds[rid + '-files.json'])
            preserved[rid] = {'record': copy.deepcopy(source_record), 'baseline': copy.deepcopy(source_record),
                             'files': copy.deepcopy(source_files), 'baseline_files': copy.deepcopy(source_files)}
        return [copy.deepcopy(inventory), copy.deepcopy(metadata), record, files, copy.deepcopy(record),
                verified, copy.deepcopy(prior), journal, preserved, ID, inventory_sha, metadata_sha]

    positives = []
    rejected = []
    def checked(values):
        result = v.verify_values(*values)
        require(result['status'] == 'PASS_OFFLINE_29_FILE_PUBLICATION_WITNESSES'
                and result['fresh_complete_public_streams'] == 29
                and result['prior_same_immutable_stream_reuses'] == 0
                and result['historical_original_literal_Notes_guard'] == 'FAIL_PRESERVED',
                'Manufactured positive failed')
        return result
    checked(fixture()); positives.append('all29 manufactured fresh streams with distinct inherited new-edition IDs')
    decorated = fixture()
    for index in (2, 4):
        decorated[index]['metadata']['resource_type']['title'] = {'en': 'Software'}
        decorated[index]['access']['status'] = 'open'
    checked(decorated); positives.append('allowed vocabulary decoration and access status')
    explicit_versions = fixture()
    for index in (2, 4):
        for item in explicit_versions[index]['files']['entries'].values():
            item['version_id'] = next(row['version_id'] for row in explicit_versions[3]['entries']
                                      if row['key'] == item['key'])
    checked(explicit_versions); positives.append('matching optional embedded content-version IDs')

    def reject(name, mutate):
        values = fixture(); mutate(values)
        try:
            checked(values)
        except (ValueError, TypeError, KeyError, AttributeError):
            rejected.append(name)
        else:
            raise ValueError('Accepted corrupted manufactured witness: ' + name)
    def embedded(values, record_index=2):
        return next(iter(values[record_index]['files']['entries'].values()))
    def listed(values, record_id):
        return next(iter(v.file_map(values[8][record_id]['files']).values()))
    def inline_space_mutant(values):
        expected = '<b>Hello</b> <b>world</b>'
        values[1]['metadata']['description'] = expected
        values[11] = sha(canonical(values[1]))
        values[5]['metadata_sha256'] = values[11]
        values[2]['metadata']['description'] = '<b>Hello</b><b>world</b>'
        values[4]['metadata']['description'] = expected
    def mismatched_dictionary_key(values):
        entries = values[2]['files']['entries']
        name = next(iter(entries)); entries['WRONG_DICTIONARY_FILENAME'] = entries.pop(name)
    reject('wrong expected new ID', lambda x: x.__setitem__(9, '90000002'))
    reject('prior ID reused', lambda x: x.__setitem__(9, v.PRIOR))
    reject('wrong concept family', lambda x: x[2]['parent'].__setitem__('id', '22922927'))
    reject('new version unpublished', lambda x: x[2].__setitem__('is_published', False))
    reject('new version marked draft', lambda x: x[2].__setitem__('is_draft', True))
    reject('new version missing draft state', lambda x: x[2].pop('is_draft'))
    reject('new version DOI wrong', lambda x: x[2]['pids']['doi'].__setitem__('identifier', '10.5281/zenodo.1'))
    reject('latest still prior', lambda x: x[4].__setitem__('id', v.PRIOR))
    reject('new public files URL wrong', lambda x: x[2]['links'].__setitem__('files', PREFIX + '/draft/files'))
    reject('new listing omitted file', lambda x: x[3]['entries'].pop())
    reject('new listing duplicate file', lambda x: x[3]['entries'].append(copy.deepcopy(x[3]['entries'][0])))
    reject('new listing changed size', lambda x: x[3]['entries'][0].__setitem__('size', 1))
    reject('new listing size boolean', lambda x: x[3]['entries'][0].__setitem__('size', True))
    reject('new listing changed MD5', lambda x: x[3]['entries'][0].__setitem__('checksum', 'md5:' + '0' * 32))
    reject('new listing incomplete', lambda x: x[3]['entries'][0].__setitem__('status', 'pending'))
    reject('new listing missing file ID', lambda x: x[3]['entries'][0].pop('file_id'))
    reject('new listing missing version ID', lambda x: x[3]['entries'][0].pop('version_id'))
    reject('new listing file ID boolean', lambda x: x[3]['entries'][0].__setitem__('file_id', True))
    reject('new embedded file ID differs', lambda x: embedded(x).__setitem__('id', 'wrong'))
    reject('latest embedded file ID differs', lambda x: embedded(x, 4).__setitem__('id', 'wrong'))
    reject('new embedded file ID missing', lambda x: embedded(x).pop('id'))
    reject('new embedded optional version differs', lambda x: embedded(x).__setitem__('version_id', 'wrong'))
    reject('dictionary filename contradicts nested filename', mismatched_dictionary_key)
    reject('metadata title changed', lambda x: x[2]['metadata'].__setitem__('title', 'Altered'))
    reject('metadata references omitted', lambda x: x[2]['metadata'].pop('references'))
    reject('metadata description alters meaning', lambda x: x[2]['metadata'].__setitem__('description', 'Higher-dimensional origin established'))
    reject('metadata description masquerades as normalized HTML event list',
           lambda x: x[2]['metadata'].__setitem__('description', v.normalized(x[1]['metadata']['description'],
                                             x[1]['metadata']['description'], ('metadata', 'description'))))
    reject('metadata visible inline space deleted', inline_space_mutant)
    reject('metadata processing instruction inserted',
           lambda x: x[2]['metadata'].__setitem__('description', '<?fabricated instruction?>' + x[1]['metadata']['description']))
    reject('metadata declaration inserted',
           lambda x: x[2]['metadata'].__setitem__('description', '<!DOCTYPE html>' + x[1]['metadata']['description']))
    reject('metadata extra field', lambda x: x[2]['metadata'].__setitem__('extra', True))
    reject('custom fields deleted', lambda x: x[2].__setitem__('custom_fields', {}))
    reject('access restricted', lambda x: x[2]['access'].__setitem__('files', 'restricted'))
    reject('verification marked draft', lambda x: x[5].__setitem__('published', False))
    reject('verification 27-file status', lambda x: x[5].__setitem__('status', 'PASS_COMPLETE_27_FILE_EDITION_READBACK'))
    reject('verification wrong count', lambda x: x[5].__setitem__('file_count', 28))
    reject('verification wrong inventory pin', lambda x: x[5].__setitem__('inventory_sha256', '0' * 64))
    reject('verification wrong metadata pin', lambda x: x[5].__setitem__('metadata_sha256', '0' * 64))
    reject('verification wrong bytes', lambda x: x[5].__setitem__('total_bytes', 1))
    reject('verification missing basis row', lambda x: x[5]['content'].pop())
    reject('verification duplicate basis row', lambda x: x[5]['content'].append(copy.deepcopy(x[5]['content'][0])))
    reject('verification changed SHA256', lambda x: x[5]['content'][0].__setitem__('sha256', '0' * 64))
    reject('inherited asset attempts old stream reuse', lambda x: x[5]['content'][0].__setitem__('basis', 'SAME_IMMUTABLE_FILE_AND_CONTENT_VERSION_AS_PRIOR_FULL_STREAM'))
    reject('addition attempts old stream reuse', lambda x: x[5]['content'][-1].__setitem__('basis', 'SAME_IMMUTABLE_FILE_AND_CONTENT_VERSION_AS_PRIOR_FULL_STREAM'))
    reject('receipt file ID differs', lambda x: x[5]['content'][0].__setitem__('remote_file_id', 'wrong'))
    reject('receipt version ID differs', lambda x: x[5]['content'][0].__setitem__('remote_version_id', 'wrong'))
    reject('receipt file ID omitted', lambda x: x[5]['content'][0].pop('remote_file_id'))
    reject('receipt version ID omitted', lambda x: x[5]['content'][0].pop('remote_version_id'))
    reject('terminal public GET absent', lambda x: x[7].pop(0))
    reject('terminal GET authenticated', lambda x: x[7][0].__setitem__('authenticated_request', True))
    reject('stream appears before latest listing', lambda x: x[7].append(copy.deepcopy(x[7][0])))
    reject('one inherited stream missing', lambda x: x[7].pop(1))
    reject('one addition stream missing', lambda x: x[7].pop())
    reject('stream duplicate', lambda x: x[7].append(copy.deepcopy(x[7][-1])))
    reject('stream wrong byte count', lambda x: x[7][1]['observed'].__setitem__('bytes', 1))
    reject('stream wrong MD5', lambda x: x[7][1]['observed'].__setitem__('md5', '0' * 32))
    reject('stream wrong SHA256', lambda x: x[7][1]['observed'].__setitem__('sha256', '0' * 64))
    reject('wrong streamed URL binding', lambda x: x[3]['entries'][0]['links'].__setitem__('content', PREFIX + '/draft/files/wrong/content'))
    reject('prior27 attestation incomplete', lambda x: x[6]['rows'].pop())
    reject('prior27 attestation status wrong', lambda x: x[6]['rows'][0].__setitem__('status', 'PASS_COMPLETE_REMOTE_BYTES_MD5_SHA256'))
    reject('prior27 attestation checksum changed', lambda x: x[6]['rows'][0]['observed'].__setitem__('sha256', '0' * 64))
    reject('prior27 attestation dated file ID changed', lambda x: x[6]['rows'][0].__setitem__('remote_file_id', 'wrong'))
    for rid in v.PRESERVED:
        reject('preserved ' + rid + ' omitted', lambda x, rid=rid: x[8].pop(rid))
        reject('preserved ' + rid + ' metadata changed', lambda x, rid=rid: x[8][rid]['record']['metadata'].__setitem__('title', 'Altered'))
        reject('preserved ' + rid + ' embedded files omitted', lambda x, rid=rid: x[8][rid]['record'].__setitem__('files', []))
        reject('preserved ' + rid + ' file ID changed', lambda x, rid=rid: listed(x, rid).__setitem__('file_id', 'wrong'))
        reject('preserved ' + rid + ' version ID changed', lambda x, rid=rid: listed(x, rid).__setitem__('version_id', 'wrong'))
        reject('preserved ' + rid + ' dated file list omitted', lambda x, rid=rid: x[8][rid].pop('baseline_files'))
    reject('preserved empty companion custom_fields omitted', lambda x: x[8]['23111008']['record'].pop('custom_fields'))
    reject('inventory unsealed', lambda x: x[0].__setitem__('sealed_upload_manifest', False))
    reject('inventory replay incomplete', lambda x: x[0].__setitem__('scientific_replay_status', 'PENDING'))
    reject('inventory count changes type', lambda x: x[0].__setitem__('expected_total_files_final', 29.0))
    reject('inventory byte total wrong', lambda x: x[0].__setitem__('total_bytes_final', 1))
    reject('inventory inherited pin omitted', lambda x: x[0]['inherited_files'].pop())
    reject('inventory addition duplicated', lambda x: x[0]['new_files'].append(copy.deepcopy(x[0]['new_files'][0])))
    reject('inventory boolean bytes', lambda x: x[0]['new_files'][0].__setitem__('bytes', True))

    with tempfile.TemporaryDirectory(prefix='manufactured-29-publication-', dir=HERE) as temporary:
        root = Path(temporary); values = fixture(); files = {}; roles = {}; groups = {}
        names = ['inventory', 'metadata', 'new_public_record', 'new_public_files', 'latest_public_record',
                 'verified_public', 'prior_full_stream']
        for role, value in zip(names, values[:7]):
            name = role + '.json'
            raw = prior_raw if role == 'prior_full_stream' else canonical(value)
            (root / name).write_bytes(raw); roles[role] = name
            files[name] = {'bytes': len(raw), 'sha256': sha(raw)}
        name = 'journal.jsonl'; raw = b''.join(canonical(event) + b'\n' for event in values[7])
        (root / name).write_bytes(raw); roles['controller_journal'] = name
        files[name] = {'bytes': len(raw), 'sha256': sha(raw)}
        for rid, group in values[8].items():
            groups[rid] = {}
            for role, value in group.items():
                name = rid + '-' + role + '.json'
                raw = (seeds[rid + '-record.json'] if role == 'baseline' else
                       seeds[rid + '-files.json'] if role == 'baseline_files' else canonical(value))
                (root / name).write_bytes(raw); groups[rid][role] = name
                files[name] = {'bytes': len(raw), 'sha256': sha(raw)}
        manifest = {'schema_version': 1, 'files': files, 'witnesses': roles, 'preserved_records': groups,
                    'historical_original_literal_Notes_guard': 'FAIL_PRESERVED'}
        manifest_raw = canonical(manifest); manifest_path = root / 'WITNESS_MANIFEST.json'
        manifest_path.write_bytes(manifest_raw); manifest_sha = sha(manifest_raw)
        def directory(manifest_pin=manifest_sha, inventory_pin=inventory_sha, metadata_pin=metadata_sha):
            return v.verify_directory(root, manifest_path, manifest_pin, ID, inventory_pin, metadata_pin)
        require(directory()['status'] == 'PASS_OFFLINE_29_FILE_PUBLICATION_WITNESSES', 'Directory positive failed')
        positives.append('complete externally pinned manufactured witness directory')
        def directory_reject(name, operation):
            try:
                operation()
            except (ValueError, TypeError, KeyError, AttributeError):
                rejected.append(name)
            else:
                raise ValueError('Accepted directory mutant: ' + name)
        directory_reject('wrong external witness manifest pin', lambda: directory(manifest_pin='0' * 64))
        directory_reject('wrong external inventory pin', lambda: directory(inventory_pin='0' * 64))
        directory_reject('wrong external metadata pin', lambda: directory(metadata_pin='0' * 64))
        target = root / roles['verified_public']; original = target.read_bytes()
        target.write_bytes(original + b' ')
        directory_reject('saved witness bytes changed', directory); target.write_bytes(original)
        real = root / 'real.json'; real.write_bytes(original); target.unlink(); target.symlink_to(real)
        directory_reject('symlink witness', directory); target.unlink(); target.write_bytes(original)
        altered = copy.deepcopy(manifest); altered['preserved_records'][v.PRIOR].pop('baseline_files')
        raw = canonical(altered); manifest_path.write_bytes(raw)
        directory_reject('dated baseline file role omitted', lambda: directory(manifest_pin=sha(raw)))
        altered = copy.deepcopy(manifest); altered['historical_original_literal_Notes_guard'] = 'PASS'
        raw = canonical(altered); manifest_path.write_bytes(raw)
        directory_reject('historic Notes failure promoted', lambda: directory(manifest_pin=sha(raw)))
        manifest_path.write_bytes(manifest_raw)
        for rid in v.PRESERVED:
            for role in ('baseline', 'baseline_files'):
                changed = root / groups[rid][role]; original = changed.read_bytes()
                changed.write_bytes(original + b' ')
                altered = copy.deepcopy(manifest)
                altered['files'][groups[rid][role]] = {'bytes': changed.stat().st_size, 'sha256': sha(changed.read_bytes())}
                raw = canonical(altered); manifest_path.write_bytes(raw)
                directory_reject('re-pinned dated ' + rid + '/' + role + ' replacement',
                                 lambda raw=raw: directory(manifest_pin=sha(raw)))
                changed.write_bytes(original); manifest_path.write_bytes(manifest_raw)
        directory_reject('duplicate JSON keys', lambda: v.load('{"a":1,"a":2}'))
        directory_reject('nonfinite JSON value', lambda: v.load('{"a":NaN}'))
        directory_reject('overflowing finite JSON syntax', lambda: v.load('{"a":1e9999}'))
        directory_reject('unsafe relative witness path', lambda: v.safe_file(root, '../escape.json'))
    checked(fixture()); positives.append('controls preserve the original manufactured fixture')
    receipt = {'schema_version': 1, 'status': 'PASS_MANUFACTURED_OFFLINE_29_FILE_PUBLICATION_CONTROLS',
               'positive_checks': len(positives), 'rejected_mutants': len(rejected),
               'positives': positives, 'rejections': rejected,
               'fixture_scope': 'FABRICATED_NEW_RECORD_WITNESSES_WITH_PINNED_DATED_PRIOR_SEEDS_ONLY',
               'fabricated_record_id': ID, 'verifier_sha256': args.verifier_sha256,
               'seed_manifest_sha256': args.seed_manifest_sha256,
               'test_sha256': sha(Path(__file__).read_bytes()),
               'optimized_python': bool(sys.flags.optimize),
               'network_calls': 0, 'numeric_imports': 0, 'source_callbacks': 0,
               'array_decodes': 0, 'target_evaluations': 0, 'remote_writes': 0,
               'actual_new_publication_readback_claim': False}
    output = Path(args.output).absolute()
    require(not output.exists(), 'Fresh controls output required')
    with output.open('x') as stream:
        stream.write(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps({key: receipt[key] for key in ('status', 'positive_checks', 'rejected_mutants', 'network_calls')}))


if __name__ == '__main__':
    main()
