#!/usr/bin/env python3
"""Nonauthor, offline manufactured controls; this never reads publication payloads."""
import argparse
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import socket
import sys
import tempfile

sys.dont_write_bytecode = True
HERE = Path(__file__).absolute().parent
ROOT = HERE.parent
SOURCE_PIN = '4a030499aab3771819a7f52431ae7415d9e82edde6ca996d7f63ae62e089364d'
SEED_PIN = 'b778c6c6583365db2892a722eadb4719f8246572dc20e688e6e1b8417aee2db7'
ID = '90000007'
PREFIX = 'https://zenodo.org/api/records/' + ID


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'),
                      ensure_ascii=False, allow_nan=False).encode()


def forbidden_network(*args, **kwargs):
    raise RuntimeError('Independent manufactured controls prohibit network access')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True)
    parser.add_argument('--source', default=str(ROOT / 'verify_publication.py'))
    parser.add_argument('--verifier-sha256', default=SOURCE_PIN)
    args = parser.parse_args()
    source = Path(args.source).absolute()
    require(sha(source.read_bytes()) == args.verifier_sha256, 'Externally provided verifier pin differs')
    seed_root = ROOT / 'manufactured-seeds'
    seed_raw = (seed_root / 'SEED_MANIFEST.json').read_bytes()
    require(sha(seed_raw) == SEED_PIN, 'Externally provided seed manifest pin differs')
    seeds = {}
    for name, pin in json.loads(seed_raw)['files'].items():
        require(Path(name).name == name, 'Leaf seed name required')
        raw = (seed_root / name).read_bytes()
        require(len(raw) == pin['bytes'] and sha(raw) == pin['sha256'], 'Seed bytes differ')
        seeds[name] = raw
    socket.socket = forbidden_network
    socket.create_connection = forbidden_network
    spec = importlib.util.spec_from_file_location('independent_publication_verifier', source)
    verifier = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(verifier)
    prior = verifier.load(seeds['PRIOR27_CONTENT_RECEIPT.json'])
    require(sha(seeds['PRIOR27_CONTENT_RECEIPT.json']) == verifier.PRIOR_STREAM_PIN,
            'Historical prior composite receipt pin differs')
    inherited = [dict(filename=row['filename'], **row['observed']) for row in prior['rows']]
    additions = []
    for name, payload in [('INDEPENDENT_MANUFACTURED_ONLY.md', b'Independent review fixture text'),
                          ('INDEPENDENT_MANUFACTURED_ONLY.zip', b'Independent review fixture bytes')]:
        additions.append({'filename': name, 'bytes': len(payload),
                          'md5': hashlib.md5(payload).hexdigest(), 'sha256': sha(payload)})
    inventory = {'sealed_upload_manifest': True, 'concept_doi': '10.5281/zenodo.17088132',
                 'scientific_replay_status': 'PASS_COMPLETE_STORED_STATE_COMPARISON_REPLAY',
                 'inherited_files': inherited, 'new_files': additions,
                 'expected_total_files_final': 29,
                 'total_bytes_final': sum(row['bytes'] for row in inherited + additions)}
    metadata = {key: copy.deepcopy(verifier.load(seeds['23225288-record.json'])[key])
                for key in ('metadata', 'custom_fields', 'access')}
    metadata['access'].pop('status', None)
    metadata['metadata']['title'] = 'Independent manufactured fixture only'
    metadata['metadata']['description'] = '<b>Independent</b> <b>manufactured metadata.</b>'
    metadata['metadata']['version'] = 'INDEPENDENT_MANUFACTURED_ONLY'
    inventory_hash = sha(canonical(inventory))
    metadata_hash = sha(canonical(metadata))

    def fixture():
        entries = []
        for number, pin in enumerate(inherited + additions):
            entries.append({'key': pin['filename'], 'size': pin['bytes'], 'checksum': 'md5:' + pin['md5'],
                            'file_id': 'independent-file-' + str(number),
                            'version_id': 'independent-version-' + str(number), 'status': 'completed',
                            'links': {'content': PREFIX + '/files/' + verifier.quote(pin['filename'], safe='') + '/content'}})
        embedded = {row['key']: {'key': row['key'], 'size': row['size'], 'checksum': row['checksum'],
                                 'id': row['file_id'], 'status': row['status']}
                    for row in entries}
        record = dict(copy.deepcopy(metadata), id=ID, parent={'id': '17088132'}, is_published=True,
                      is_draft=False, pids={'doi': {'identifier': '10.5281/zenodo.' + ID}},
                      files={'entries': embedded}, links={'self': PREFIX, 'files': PREFIX + '/files'})
        receipt = {'status': 'PASS_COMPLETE_29_FILE_EDITION_READBACK', 'published': True,
                   'draft_id': ID, 'parent_id': '17088132', 'inventory_sha256': inventory_hash,
                   'metadata_sha256': metadata_hash, 'file_count': 29,
                   'total_bytes': inventory['total_bytes_final'], 'content': []}
        journal = [{'kind': 'GET', 'url': PREFIX + '/files', 'status': 200, 'authenticated_request': False}]
        for pin, item in zip(inherited + additions, entries):
            receipt['content'].append({'filename': pin['filename'], 'sha256': pin['sha256'],
                                      'basis': 'FRESH_COMPLETE_CONTENT_STREAM',
                                      'remote_file_id': item['file_id'], 'remote_version_id': item['version_id']})
            journal.append({'kind': 'FULL_CONTENT_STREAM_VERIFIED', 'filename': pin['filename'],
                            'observed': {key: pin[key] for key in ('bytes', 'md5', 'sha256')}})
        preserved = {}
        for rid in ('23225288', '23114217', '22347452', '23111008'):
            record_seed = verifier.load(seeds[rid + '-record.json'])
            files_seed = verifier.load(seeds[rid + '-files.json'])
            preserved[rid] = {'record': copy.deepcopy(record_seed), 'baseline': copy.deepcopy(record_seed),
                              'files': copy.deepcopy(files_seed), 'baseline_files': copy.deepcopy(files_seed)}
        return [copy.deepcopy(inventory), copy.deepcopy(metadata), record, {'entries': entries},
                copy.deepcopy(record), receipt, copy.deepcopy(prior), journal, preserved,
                ID, inventory_hash, metadata_hash]

    positives, rejected, accepted_corruptions = [], [], []

    def check(values):
        result = verifier.verify_values(*values)
        require(result['status'] == 'PASS_OFFLINE_29_FILE_PUBLICATION_WITNESSES'
                and result['fresh_complete_public_streams'] == 29
                and result['prior_same_immutable_stream_reuses'] == 0
                and result['original_binary80_state_error'] == 'ENCLOSED_RETAINED_NODES'
                and result['full_twelve_case_certificate'] == 'UNRESOLVED'
                and result['higher_dimensional_Big_Bang_origin'] == 'NOT_ESTABLISHED'
                and result['metric_calibration'] == 'FAIL'
                and result['external_novelty'] == 'NOT_ASSESSED', 'Result contract differs')
        return result

    def positive(name, mutation=lambda values: None):
        values = fixture(); mutation(values); check(values); positives.append(name)

    def negative(name, mutation):
        values = fixture(); mutation(values)
        try:
            check(values)
        except (ValueError, TypeError, KeyError, AttributeError):
            rejected.append(name)
        else:
            accepted_corruptions.append(name)

    def embedded(values, index=2):
        return next(iter(values[index]['files']['entries'].values()))

    def current_preserved(values, rid, role='files'):
        return next(iter(verifier.file_map(values[8][rid][role]).values()))

    positive('Complete independent 29-file fixture with record id and no embedded version IDs')
    def use_file_id(values):
        for index in (2, 4):
            for item in values[index]['files']['entries'].values():
                item['file_id'] = item.pop('id')
    positive('Embedded immutable identity accepted through file_id alias', use_file_id)
    def use_both_and_versions(values):
        by_name = {item['key']: item for item in values[3]['entries']}
        for index in (2, 4):
            for item in values[index]['files']['entries'].values():
                item['file_id'] = item['id']
                item['version_id'] = by_name[item['key']]['version_id']
    positive('Matching aliases and optional embedded versions accepted', use_both_and_versions)
    def historical_ids_with_fresh_streams(values):
        by_name = {row['filename']: row for row in values[6]['rows']}
        for item in values[3]['entries']:
            if item['key'] in by_name:
                item['file_id'] = by_name[item['key']]['remote_file_id']
                item['version_id'] = by_name[item['key']]['remote_version_id']
        listed = {item['key']: item for item in values[3]['entries']}
        for index in (2, 4):
            for item in values[index]['files']['entries'].values():
                item['id'] = listed[item['key']]['file_id']
        for row in values[5]['content']:
            row['remote_file_id'] = listed[row['filename']]['file_id']
            row['remote_version_id'] = listed[row['filename']]['version_id']
    positive('All29 fresh streams accepted even when inherited immutable IDs equal historical IDs', historical_ids_with_fresh_streams)
    def semantically_same_html(values):
        values[2]['metadata']['description'] = '<b>Independent</b> <b>manufactured metadata&#46;</b>'
    positive('Semantic HTML character reference accepted', semantically_same_html)

    for index in (2, 4):
        negative('Conflicting record file_id and id aliases index ' + str(index),
                 lambda x, index=index: embedded(x, index).__setitem__('file_id', 'different'))
        negative('Optional embedded content version boolean index ' + str(index),
                 lambda x, index=index: embedded(x, index).__setitem__('version_id', True))
        negative('Empty embedded file identity index ' + str(index),
                 lambda x, index=index: embedded(x, index).__setitem__('id', ''))
        negative('Numeric embedded immutable file identity index ' + str(index),
                 lambda x, index=index: embedded(x, index).__setitem__('id', 123))
    negative('Public list conflicting identity aliases', lambda x: x[3]['entries'][0].__setitem__('id', 'other'))
    negative('Public list content-version identity is list', lambda x: x[3]['entries'][0].__setitem__('version_id', ['bad']))
    negative('Receipt fresh basis wrong case', lambda x: x[5]['content'][0].__setitem__('basis', 'fresh_complete_content_stream'))
    negative('Receipt fresh basis trailing whitespace', lambda x: x[5]['content'][0].__setitem__('basis', 'FRESH_COMPLETE_CONTENT_STREAM '))
    negative('Receipt historical basis masked by terminal fresh event', lambda x: x[5]['content'][0].__setitem__('basis', 'PASS_PRIOR_COMPLETE_BYTES_WITH_FRESH_IMMUTABLE_IDENTITY'))
    negative('Receipt fresh file identity boolean', lambda x: x[5]['content'][0].__setitem__('remote_file_id', True))
    negative('Receipt fresh version identity omitted', lambda x: x[5]['content'][0].pop('remote_version_id'))
    negative('Terminal stream observed bytes float', lambda x: x[7][1]['observed'].__setitem__('bytes', float(x[7][1]['observed']['bytes'])))
    negative('Terminal stream extra checksum claim', lambda x: x[7][1]['observed'].__setitem__('extra', 'claim'))
    negative('Terminal stream changed to historical kind', lambda x: x[7][1].__setitem__('kind', 'PRIOR_FULL_CONTENT_STREAM_VERIFIED'))
    negative('Second terminal listing invalidates preceding streams', lambda x: x[7].append(copy.deepcopy(x[7][0])))
    negative('Metadata creator added', lambda x: x[2]['metadata']['creators'].append(copy.deepcopy(x[2]['metadata']['creators'][0])))
    negative('Metadata extra custom field', lambda x: x[2]['custom_fields'].__setitem__('independent:extra', 'claim'))
    negative('Metadata vocab identity changed', lambda x: x[2]['metadata']['resource_type'].__setitem__('id', 'wrong'))
    negative('Metadata nonvocabulary title object', lambda x: x[2]['metadata'].__setitem__('title', {'id': 'wrong', 'title': 'decoration'}))
    negative('Metadata description text changed', lambda x: x[2]['metadata'].__setitem__('description', '<p>Established origin.</p>'))
    negative('Metadata description HTML event-list type collision',
             lambda x: x[2]['metadata'].__setitem__('description', [['start', 'b', []], ['text', 'Independent'], ['end', 'b'], ['start', 'b', []], ['text', 'manufactured metadata.'], ['end', 'b']]))
    negative('Metadata description exact preserved-space event-list type collision',
             lambda x: x[2]['metadata'].__setitem__('description', [['start', 'b', []], ['text', 'Independent'], ['end', 'b'], ['text', ' '], ['start', 'b', []], ['text', 'manufactured metadata.'], ['end', 'b']]))
    negative('Metadata description visible interword HTML whitespace lost',
             lambda x: x[2]['metadata'].__setitem__('description', '<b>Independent</b><b>manufactured metadata.</b>'))
    negative('Metadata description processing instruction added',
             lambda x: x[2]['metadata'].__setitem__('description', '<?independent mutated?><b>Independent</b> <b>manufactured metadata.</b>'))
    negative('Metadata description declaration added',
             lambda x: x[2]['metadata'].__setitem__('description', '<!DOCTYPE invented><b>Independent</b> <b>manufactured metadata.</b>'))
    negative('Metadata malformed empty closing tag silently deleted',
             lambda x: x[2]['metadata'].__setitem__('description', '<b>Independent</b></> <b>manufactured metadata.</b>'))
    saved_description = metadata['metadata']['description']
    saved_metadata_hash = metadata_hash
    metadata['metadata']['description'] = '<p style="display:inline">Independent</p> <p style="display:inline">manufactured metadata.</p>'
    metadata_hash = sha(canonical(metadata))
    negative('Metadata inline-styled block interword whitespace deleted',
             lambda x: x[2]['metadata'].__setitem__('description', '<p style="display:inline">Independent</p><p style="display:inline">manufactured metadata.</p>'))
    metadata['metadata']['description'] = saved_description
    metadata_hash = saved_metadata_hash
    negative('Preserved companion empty custom_fields omitted',
             lambda x: x[8]['23111008']['record'].pop('custom_fields'))
    def rename_embedded_map_key(values):
        mapping = values[2]['files']['entries']; first = next(iter(mapping))
        mapping['CONTRADICTORY_MAP_FILENAME'] = mapping.pop(first)
    negative('Record entries map key contradicts inner filename', rename_embedded_map_key)
    for rid in ('23225288', '23114217', '22347452', '23111008'):
        negative('Preserved ' + rid + ' current file conflicting aliases',
                 lambda x, rid=rid: current_preserved(x, rid).__setitem__('id', 'different'))
        negative('Preserved ' + rid + ' dated content version absent',
                 lambda x, rid=rid: current_preserved(x, rid, 'baseline_files').pop('version_id'))
        negative('Preserved ' + rid + ' embedded id changed',
                 lambda x, rid=rid: next(iter(verifier.file_map(x[8][rid]['record']).values())).__setitem__('id', 'different'))

    directory_rejected = []
    with tempfile.TemporaryDirectory(prefix='independent-witness-', dir=HERE) as temp:
        path = Path(temp)
        values = fixture()
        pinned_files, roles, groups = {}, {}, {}
        def save(name, raw):
            (path / name).write_bytes(raw)
            pinned_files[name] = {'bytes': len(raw), 'sha256': sha(raw)}
        names = ('inventory', 'metadata', 'new_public_record', 'new_public_files',
                 'latest_public_record', 'verified_public', 'prior_full_stream')
        for role, value in zip(names, values[:7]):
            name = role + '.json'
            save(name, seeds['PRIOR27_CONTENT_RECEIPT.json'] if role == 'prior_full_stream' else canonical(value))
            roles[role] = name
        save('journal.jsonl', b''.join(canonical(row) + b'\n' for row in values[7]))
        roles['controller_journal'] = 'journal.jsonl'
        for rid, group in values[8].items():
            groups[rid] = {}
            for role, value in group.items():
                name = rid + '-' + role + '.json'
                raw = seeds[rid + ('-record.json' if role == 'baseline' else '-files.json')] if role in ('baseline', 'baseline_files') else canonical(value)
                save(name, raw); groups[rid][role] = name
        manifest = {'schema_version': 1, 'files': pinned_files, 'witnesses': roles,
                    'preserved_records': groups, 'historical_original_literal_Notes_guard': 'FAIL_PRESERVED'}
        manifest_path = path / 'manifest.json'
        def directory(value):
            raw = canonical(value); manifest_path.write_bytes(raw)
            return verifier.verify_directory(path, manifest_path, sha(raw), ID, inventory_hash, metadata_hash)
        check(values)
        require(directory(manifest)['status'] == 'PASS_OFFLINE_29_FILE_PUBLICATION_WITNESSES', 'Directory fixture fails')
        positives.append('Independent externally pinned manufactured directory')
        def reject_schema(name, mutation):
            altered = copy.deepcopy(manifest); mutation(altered)
            try:
                directory(altered)
            except (ValueError, TypeError, KeyError, AttributeError):
                directory_rejected.append(name)
            else:
                accepted_corruptions.append(name)
        reject_schema('Manifest schema version boolean', lambda x: x.__setitem__('schema_version', True))
        reject_schema('Manifest unknown top-level field', lambda x: x.__setitem__('unknown', True))
        reject_schema('Manifest role omitted', lambda x: x['witnesses'].pop('latest_public_record'))
        reject_schema('Manifest role extra', lambda x: x['witnesses'].__setitem__('unknown', roles['inventory']))
        reject_schema('Manifest preserved group extra key', lambda x: x['preserved_records']['23225288'].__setitem__('extra', roles['inventory']))
        reject_schema('Manifest preserved dated file role omitted', lambda x: x['preserved_records']['23111008'].pop('baseline_files'))
        reject_schema('Manifest pin bytes float', lambda x: x['files'][roles['inventory']].__setitem__('bytes', float(x['files'][roles['inventory']]['bytes'])))
        reject_schema('Manifest historical Notes guard improved', lambda x: x.__setitem__('historical_original_literal_Notes_guard', 'PASS'))
        for rid in ('23114217', '22347452', '23111008'):
            altered = copy.deepcopy(manifest)
            for role in ('record', 'baseline'):
                name = groups[rid][role]
                body = verifier.load((path / name).read_bytes())
                body['metadata']['title'] = 'Substituted dated history despite a new witness pin'
                raw = canonical(body)
                (path / name).write_bytes(raw)
                altered['files'][name] = {'bytes': len(raw), 'sha256': sha(raw)}
            try:
                directory(altered)
            except (ValueError, TypeError, KeyError, AttributeError):
                directory_rejected.append('Re-pinned dated/current baseline metadata substitution ' + rid)
            else:
                accepted_corruptions.append('Re-pinned dated/current baseline metadata substitution ' + rid)
            for role in ('record', 'baseline'):
                name = groups[rid][role]
                (path / name).write_bytes(seeds[rid + '-record.json'] if role == 'baseline' else canonical(values[8][rid]['record']))
        for rid in ('23114217', '22347452', '23111008'):
            altered = copy.deepcopy(manifest)
            filename = next(iter(verifier.file_map(values[8][rid]['baseline_files'])))
            for role in ('record', 'baseline', 'files', 'baseline_files'):
                name = groups[rid][role]
                body = verifier.load((path / name).read_bytes())
                item = verifier.file_map(body)[filename]
                for key in ('file_id', 'id'):
                    if key in item:
                        item[key] = 'substituted-dated-file-identity'
                if 'version_id' in item:
                    item['version_id'] = 'substituted-dated-content-version'
                raw = canonical(body)
                (path / name).write_bytes(raw)
                altered['files'][name] = {'bytes': len(raw), 'sha256': sha(raw)}
            try:
                directory(altered)
            except (ValueError, TypeError, KeyError, AttributeError):
                directory_rejected.append('Re-pinned dated/current listing and embedded immutable-identity substitution ' + rid)
            else:
                accepted_corruptions.append('Re-pinned dated/current listing and embedded immutable-identity substitution ' + rid)
            for role in ('record', 'baseline', 'files', 'baseline_files'):
                name = groups[rid][role]
                raw = seeds[rid + ('-record.json' if role == 'baseline' else '-files.json')] if role in ('baseline', 'baseline_files') else canonical(values[8][rid][role])
                (path / name).write_bytes(raw)

    receipt = {'schema_version': 1,
               'status': 'BLOCKED_INDEPENDENT_MANUFACTURED_PUBLICATION_CORRUPTIONS' if accepted_corruptions else 'PASS_INDEPENDENT_MANUFACTURED_OFFLINE_PUBLICATION_CONTROLS',
               'verifier_sha256': args.verifier_sha256, 'seed_manifest_sha256': SEED_PIN,
               'independent_controls_sha256': sha(Path(__file__).read_bytes()),
               'positive_checks': len(positives), 'rejected_mutants': len(rejected) + len(directory_rejected),
               'accepted_corruptions': accepted_corruptions, 'positives': positives,
               'rejections': rejected + directory_rejected, 'optimized_python': bool(sys.flags.optimize),
               'fixture_scope': 'INDEPENDENT_MANUFACTURED_NEW_RECORD_WITH_EXTERNALLY_PINNED_DATED_SEEDS_ONLY',
               'fabricated_record_id': ID, 'network_calls': 0, 'remote_writes': 0, 'numeric_imports': 0,
               'array_decodes': 0, 'source_callbacks': 0, 'target_evaluations': 0,
               'actual_new_publication_readback_claim': False}
    output = Path(args.output).absolute()
    require(not output.exists(), 'Fresh independent receipt required')
    with output.open('x') as stream:
        stream.write(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps({key: receipt[key] for key in ('status', 'positive_checks', 'rejected_mutants', 'accepted_corruptions')}))


if __name__ == '__main__':
    main()
