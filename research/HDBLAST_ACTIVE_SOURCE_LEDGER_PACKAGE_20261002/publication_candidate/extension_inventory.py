"""Validate the exact same-family 10+9+3 publication inventory, locally only."""
from datetime import date
import hashlib
import json
from pathlib import Path
import re

ROLES = ('active_source_complete_package', 'active_source_final_result',
         'active_source_external_fresh_replay')
CURRENT_PUBLISHED_RECORD = 23112891
CURRENT_MAIN_DRAFT = 23114217
AUTHORIZED_OWNER = 1386319


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    def pairs(items):
        result = {}
        for key, value in items:
            require(key not in result, 'Duplicate JSON key')
            result[key] = value
        return result
    def constant(_value):
        raise RuntimeError('Nonfinite JSON literal')
    return json.loads(Path(path).read_text(), object_pairs_hook=pairs, parse_constant=constant)


def checked_file(root, name):
    root = Path(root).absolute()
    require(root.is_dir() and not root.is_symlink()
            and all(not p.is_symlink() for p in root.parents), 'Invalid source directory')
    require(type(name) is str and name and '\\' not in name and ':' not in name,
            'Unsafe candidate path')
    parts = name.split('/')
    require(not Path(name).is_absolute() and all(p not in ('', '.', '..') for p in parts),
            'Unsafe candidate path')
    path = root
    for part in parts:
        path = path / part
        require(not path.is_symlink(), 'Symlink candidate path')
    require(path.is_file(), 'Candidate input missing')
    return path


def valid_hash(value, length=64):
    return type(value) is str and re.fullmatch('[0-9a-f]{%d}' % length, value) is not None


def file_inventory(files):
    require(type(files) is list and len(files) == 10, 'Expected exact ten inherited files')
    result = {}
    for row in files:
        require(type(row) is dict, 'Malformed remote file')
        name = row.get('filename')
        require(type(name) is str and name not in result, 'Malformed or duplicate remote filename')
        sizes = [row[key] for key in ('filesize', 'size', 'bytes') if key in row]
        digests = [row[key] for key in ('checksum', 'md5') if key in row]
        require(sizes and all(type(value) is int and value > 0 for value in sizes)
                and len(set(sizes)) == 1, 'Conflicting or invalid remote byte counts')
        require(digests and all(type(value) is str for value in digests),
                'Missing remote size or checksum')
        digests = [value.lower().removeprefix('md5:') for value in digests]
        require(all(valid_hash(value, 32) for value in digests)
                and len(set(digests)) == 1, 'Conflicting or invalid remote MD5')
        result[name] = {'filename': name, 'bytes': sizes[0], 'md5': digests[0]}
    return [result[name] for name in sorted(result)]


def authenticated_target(wrapper, wrapper_sha256, inherited, published):
    require(type(wrapper) is dict and wrapper.get('method') == 'GET'
            and wrapper.get('authenticated_request') is True
            and type(wrapper.get('http_status')) is int and wrapper['http_status'] == 200,
            'Target proof is not an authenticated successful GET')
    data = wrapper.get('data')
    expected_id = CURRENT_PUBLISHED_RECORD if published else CURRENT_MAIN_DRAFT
    require(type(data) is dict and type(data.get('id')) is int and data['id'] == expected_id
            and data.get('conceptdoi') == '10.5281/zenodo.17088132'
            and str(data.get('conceptrecid')) == '17088132'
            and type(data.get('owner')) is int and data['owner'] == AUTHORIZED_OWNER,
            'Wrong authenticated target id, family or owner')
    url = 'https://zenodo.org/api/deposit/depositions/' + str(expected_id)
    require(wrapper.get('requested_url') == url and wrapper.get('final_url') == url,
            'Authenticated target GET used the wrong endpoint')
    require(data.get('submitted') is published and data.get('state') == ('done' if published else 'unsubmitted'),
            'Authenticated target publication state differs')
    if published:
        require(data.get('doi') == '10.5281/zenodo.' + str(expected_id),
                'Authenticated published record DOI differs')
    metadata = data.get('metadata')
    require(type(metadata) is dict and metadata.get('title') and metadata.get('description')
            and type(metadata.get('creators')) is list and metadata['creators']
            and metadata.get('license'), 'Authenticated target metadata is incomplete')
    remote_files = file_inventory(data.get('files'))
    require(remote_files == file_inventory(inherited), 'Authenticated target inherited inventory differs')
    require(valid_hash(wrapper_sha256) and valid_hash(wrapper.get('response_sha256'))
            and type(wrapper.get('started_utc')) is str and wrapper['started_utc'],
            'Missing authenticated target response provenance')
    return {'id': expected_id, 'owner': AUTHORIZED_OWNER, 'concept_doi': data['conceptdoi'],
        'submitted': published, 'state': data['state'], 'creators': metadata['creators'],
        'doi': data.get('doi'),
        'license': metadata['license'], 'files': remote_files, 'method': 'GET',
        'http_status': 200, 'authenticated_request': True, 'requested_url': url,
        'recorded_utc': wrapper['started_utc'], 'source_wrapper_sha256': wrapper_sha256,
        'source_response_sha256': wrapper['response_sha256']}


def validate_target_receipt(receipt, inherited, metadata):
    require(type(receipt) is dict and receipt.get('status') == 'VERIFIED_RECONCILED_MAIN_TARGET',
            'Target is a historical snapshot, not reconciled current readiness')
    published, draft = receipt.get('published'), receipt.get('draft')
    require(type(published) is dict and type(draft) is dict, 'Missing reconciled target evidence')
    for row, expected_id, expected_state in [(published, CURRENT_PUBLISHED_RECORD, 'done'),
                                             (draft, CURRENT_MAIN_DRAFT, 'unsubmitted')]:
        require(type(row.get('id')) is int and row['id'] == expected_id
                and type(row.get('owner')) is int and row['owner'] == AUTHORIZED_OWNER
                and row.get('state') == expected_state and row.get('submitted') is (expected_state == 'done')
                and row.get('concept_doi') == '10.5281/zenodo.17088132'
                and row.get('authenticated_request') is True and row.get('method') == 'GET'
                and type(row.get('http_status')) is int and row['http_status'] == 200,
                'Reconciled target owner/family/lifecycle differs')
        require(file_inventory(row.get('files')) == file_inventory(inherited),
                'Reconciled target inherited bytes differ')
        require(valid_hash(row.get('source_wrapper_sha256')) and valid_hash(row.get('source_response_sha256')),
                'Reconciled target readback provenance missing')
        require(row.get('requested_url') == 'https://zenodo.org/api/deposit/depositions/' + str(expected_id)
                and type(row.get('recorded_utc')) is str and row['recorded_utc'],
                'Reconciled target endpoint or timestamp missing')
        if expected_state == 'done':
            require(row.get('doi') == '10.5281/zenodo.' + str(expected_id),
                    'Reconciled published record DOI differs')
    require(published.get('creators') == draft.get('creators') == metadata.get('creators')
            and published.get('license') == draft.get('license') == metadata.get('license'),
            'Owner identity, creators or license changed during retarget')
    return receipt


def validate_extension_inventory(manifest, candidate):
    parent = manifest.get('parent_candidate')
    require(type(parent) is dict and set(parent) == {'manifest_path', 'manifest_sha256',
            'metadata_path', 'metadata_sha256'}, 'Missing exact parent candidate provenance')
    original_path = checked_file(candidate, parent['manifest_path'])
    require(valid_hash(parent['manifest_sha256']) and sha256(original_path) == parent['manifest_sha256'],
            'Original candidate snapshot changed')
    original = read_json(original_path)
    require(original.get('concept_doi') == '10.5281/zenodo.17088132'
            and original.get('latest_published_record') == 22347452
            and original.get('candidate_version') == '2026.10.02-v25'
            and len(original.get('inherited_files', [])) == 10
            and len(original.get('new_files', [])) == 9
            and original.get('expected_final_files_count') == 19,
            'Wrong original candidate family, version or inventory')
    metadata_path = checked_file(candidate, parent['metadata_path'])
    require(valid_hash(parent['metadata_sha256']) and sha256(metadata_path) == parent['metadata_sha256']
            == original['metadata_request']['sha256'], 'Original metadata snapshot changed')
    require(manifest.get('concept_doi') == original['concept_doi']
            and manifest.get('latest_published_record') == CURRENT_PUBLISHED_RECORD
            and manifest.get('active_main_draft') == CURRENT_MAIN_DRAFT
            and type(manifest.get('active_main_draft')) is int,
            'Wrong existing main-family target')
    current_metadata = read_json(checked_file(candidate, manifest['metadata_request']['path']))['metadata']
    validate_target_receipt(manifest.get('target_reconciliation'), original['inherited_files'], current_metadata)
    require(manifest.get('automatic_github_archiving') == 'ALL_EIGHT_OFF',
            'Automatic GitHub archiving must remain OFF')
    require(manifest.get('inherited_files') == original['inherited_files'],
            'An inherited historical file was changed')
    files = manifest.get('new_files')
    require(type(files) is list and len(files) == 12
            and files[:9] == original['new_files'],
            'Nine reviewed prior additions must remain byte-pinned and unchanged')
    require(tuple(f.get('role') for f in files[9:]) == ROLES,
            'Missing or reordered active-source extension roles')
    require(type(manifest.get('expected_final_files_count')) is int
            and manifest['expected_final_files_count'] == 22,
            'Expected exactly ten inherited and twelve new files')
    names = [f['filename'] for f in manifest['inherited_files'] + files]
    require(len(set(names)) == 22, 'Duplicate publication filename')
    require(manifest.get('old_metric_scientific_status') == 'FAIL', 'Historical FAIL was changed')
    version = manifest.get('candidate_version')
    require(type(version) is str and re.fullmatch(r'[0-9]{4}\.[0-9]{2}\.[0-9]{2}-v25', version),
            'Expected pending v25 date/version')
    date.fromisoformat(version[:-4].replace('.', '-'))
    require(manifest.get('active_source_classification') in
            {'LEDGER_ERROR_DEMONSTRATED', 'NO_GATE_SCALE_ATTRIBUTION', 'CONSISTENCY_FAILURE'},
            'Missing accepted active-source classification')
    proof = manifest.get('active_source_evidence')
    require(type(proof) is dict and proof.get('status') == 'COMPLETE_ORIGINAL_AND_FRESH_REPLAY',
            'Completed original/fresh evidence absent')
    require(proof.get('classification') == manifest['active_source_classification']
            and proof.get('old_metric_status') == 'FAIL'
            and proof.get('scientific_comparison') == 'EXACT_REGISTERED_FULL_FRAMES',
            'Invalid active-source completion summary')
    for key in ('freeze_commit', 'science_commit'):
        require(valid_hash(proof.get(key), 40), 'Missing full scientific Git pin')
    for key in ('registration_sha256', 'zip_sha256', 'package_manifest_sha256',
                'fresh_receipt_sha256', 'result_report_sha256'):
        require(valid_hash(proof.get(key)), 'Missing completed evidence hash')
    require(proof['zip_sha256'] == files[9]['sha256']
            and proof['result_report_sha256'] == files[10]['sha256']
            and proof['fresh_receipt_sha256'] == files[11]['sha256'],
            'Completion summary disagrees with attachment pins')
    require(manifest.get('status') == 'PREPARED_DERIVED_CANDIDATE_NOT_PROOF_OF_PUBLICATION',
            'Unrecognized candidate preparation state')
    return original
