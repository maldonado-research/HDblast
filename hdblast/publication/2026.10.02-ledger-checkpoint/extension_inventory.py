"""Validate the exact same-family 10+9+3 publication inventory, locally only."""
from datetime import date
import hashlib
import json
from pathlib import Path
import re

ROLES = ('active_source_complete_package', 'active_source_final_result',
         'active_source_external_fresh_replay')
CURRENT_PUBLISHED_RECORD = 22347452
CURRENT_MAIN_DRAFT = 23114217
AUTHORIZED_OWNER = 1386319
CURRENT_CANDIDATE_VERSION = '2026.10.02-ledger-checkpoint'
COMPANION_FILES = [
    ('HDBLAST_SCALAR_BRANCH_SUPPLEMENT_20260923.zip', 964953, 'efe197fdcee11e93a4678efa6636703e'),
    ('HDBLAST_SCALAR_JUNCTION_STATIC_BRANCH_20260923.pdf', 1426918, '907cb21c11f817a466e06e3b11a0ff8e')]


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


def normalize_public_files(files):
    require(type(files) is list, 'Public file inventory missing')
    return [{'filename': row.get('key'), 'bytes': row.get('size'),
             'md5': row.get('checksum')} for row in files if type(row) is dict]


def proof_wrapper(candidate, proof, key, url, status, auth):
    relative = proof.get(key)
    require(type(relative) is str, 'Missing post-removal proof path')
    path = checked_file(candidate, relative)
    require(valid_hash(proof.get('proof_pins', {}).get(relative))
            and sha256(path) == proof['proof_pins'][relative], 'Post-removal proof bytes changed')
    wrapper = read_json(path)
    require(wrapper.get('method') == 'GET' and wrapper.get('http_status') == status
            and wrapper.get('authenticated_request') is auth
            and wrapper.get('requested_url') == url and 'error' not in wrapper
            and valid_hash(wrapper.get('response_sha256'))
            and type(wrapper.get('started_utc')) is str, 'Invalid post-removal GET proof')
    return wrapper


def validate_post_removal_proof(manifest, candidate, metadata):
    proof = manifest.get('post_removal_proof')
    require(type(proof) is dict and proof.get('removed_record') == 23112891
            and proof.get('preserved_latest_published_record') == CURRENT_PUBLISHED_RECORD
            and proof.get('existing_main_draft') == CURRENT_MAIN_DRAFT,
            'Wrong post-removal scope')
    for key, role, record_id, published in [
            ('authenticated_published', 'published', CURRENT_PUBLISHED_RECORD, True),
            ('authenticated_draft', 'draft', CURRENT_MAIN_DRAFT, False)]:
        relative = proof.get(key)
        wrapper = proof_wrapper(candidate, proof, key,
            'https://zenodo.org/api/deposit/depositions/' + str(record_id), 200, True)
        recomputed = authenticated_target(wrapper, proof['proof_pins'][relative],
                                          manifest['inherited_files'], published)
        require(recomputed == manifest['target_reconciliation'][role],
                'Authenticated target summary differs from pinned GET bytes')
        if published:
            require(wrapper['data']['metadata'].get('version') == '2026.09.05-v24',
                    'Authenticated latest semantic version differs')
    removed = proof_wrapper(candidate, proof, 'public_removed',
        'https://zenodo.org/api/records/23112891', 410, False)
    require(removed.get('final_url') == removed['requested_url']
            and removed.get('redirects') == [], 'Removed record proof redirected')
    tombstone = removed.get('data', {}).get('tombstone', {})
    require(tombstone.get('is_visible') is True
            and tombstone.get('removed_by', {}).get('user') == str(AUTHORIZED_OWNER)
            and tombstone.get('removal_reason', {}).get('id') == 'test-record'
            and type(tombstone.get('removal_date')) is str
            and tombstone['removal_date'].startswith('2026-10-03T')
            and '10.5281/zenodo.23112891' in tombstone.get('citation_text', ''),
            'Missing visible owner removal tombstone')
    for key, url in [('public_latest', 'https://zenodo.org/api/records/23112891/versions/latest'),
                     ('public_concept', 'https://zenodo.org/api/records/17088132')]:
        wrapper = proof_wrapper(candidate, proof, key, url, 200, False)
        data = wrapper.get('data', {})
        require(wrapper.get('final_url') == 'https://zenodo.org/api/records/22347452'
                and type(data.get('id')) is int and data['id'] == CURRENT_PUBLISHED_RECORD
                and data.get('conceptdoi') == '10.5281/zenodo.17088132'
                and str(data.get('conceptrecid')) == '17088132'
                and data.get('doi') == '10.5281/zenodo.22347452'
                and data.get('submitted') is True and data.get('state') == 'done'
                and not data.get('tombstone')
                and data.get('metadata', {}).get('version') == '2026.09.05-v24',
                'Preserved latest public record identity/lifecycle differs')
        links = data.get('links', {})
        require(links.get('self') == 'https://zenodo.org/api/records/22347452'
                and links.get('parent') == 'https://zenodo.org/api/records/17088132'
                and links.get('parent_doi') == 'https://doi.org/10.5281/zenodo.17088132'
                and links.get('latest') == 'https://zenodo.org/api/records/22347452/versions/latest',
                'Latest record family links differ')
        require(file_inventory(normalize_public_files(data.get('files')))
                == file_inventory(manifest['inherited_files']), 'Public inherited files differ')
        require(data['metadata'].get('creators') == metadata['creators']
                and data['metadata'].get('license', {}).get('id') == metadata['license'],
                'Latest public creators/license differ')
    modern = proof_wrapper(candidate, proof, 'modern_draft',
        'https://zenodo.org/api/records/23114217/draft', 200, True)
    data = modern.get('data', {})
    require(modern.get('final_url') == modern['requested_url'] and modern.get('redirects') == []
            and type(data.get('id')) is int and data['id'] == CURRENT_MAIN_DRAFT
            and data.get('recid') == str(CURRENT_MAIN_DRAFT)
            and data.get('conceptdoi') == '10.5281/zenodo.17088132'
            and str(data.get('conceptrecid')) == '17088132'
            and data.get('owners') == [{'id': str(AUTHORIZED_OWNER)}]
            and data.get('submitted') is False and data.get('state') == 'unsubmitted'
            and data.get('status') == 'draft', 'Modern draft owner/family/lifecycle differs')
    require(data.get('metadata', {}).get('relations', {}).get('version')
            == [proof.get('modern_draft_observed_relation')]
            and proof.get('modern_draft_observed_relation')
            == {'index': 25, 'is_last': False, 'parent': {'pid_type': 'recid', 'pid_value': '17088132'}}
            and proof.get('semantic_version_is_distinct_from_ui_ordinal') is True,
            'Observed draft relation changed; no UI ordinal may be inferred')
    require(data.get('metadata', {}).get('creators') == metadata['creators']
            and data.get('metadata', {}).get('license', {}).get('id') == metadata['license']
            and data.get('links', {}).get('parent_doi') == 'https://doi.org/10.5281/zenodo.17088132'
            and file_inventory(normalize_public_files(data.get('files')))
                == file_inventory(manifest['inherited_files']), 'Modern draft metadata/files differ')
    companion = proof_wrapper(candidate, proof, 'protected_companion',
        'https://zenodo.org/api/deposit/depositions/23111008', 200, True)
    data = companion.get('data', {})
    require(data.get('id') == 23111008 and data.get('owner') == AUTHORIZED_OWNER
            and data.get('submitted') is True and data.get('state') == 'done'
            and data.get('doi') == '10.5281/zenodo.23111008'
            and data.get('conceptdoi') == '10.5281/zenodo.22922927'
            and str(data.get('conceptrecid')) == '22922927'
            and len(data.get('files', [])) == 2, 'Protected companion scope differs')
    require(sorted((row.get('filename'), row.get('filesize'), row.get('checksum'))
                   for row in data['files']) == COMPANION_FILES,
            'Protected companion original filenames/bytes/MD5 differ')
    companion_public = proof_wrapper(candidate, proof, 'protected_companion_concept',
        'https://zenodo.org/api/records/22922927', 200, False)
    data = companion_public.get('data', {})
    require(companion_public.get('final_url') == 'https://zenodo.org/api/records/23111008'
            and data.get('id') == 23111008 and data.get('submitted') is True
            and data.get('state') == 'done' and not data.get('tombstone')
            and len(data.get('files', [])) == 2, 'Companion public witness differs')
    require(sorted((row.get('key'), row.get('size'), row.get('checksum', '').removeprefix('md5:'))
                   for row in data['files']) == COMPANION_FILES,
            'Public companion original filenames/bytes/MD5 differ')
    return proof


def validate_extension_inventory(manifest, candidate):
    prior = manifest.get('prior_reviewed_candidate')
    require(type(prior) is dict and set(prior) == {'manifest_path', 'manifest_sha256',
            'metadata_path', 'metadata_sha256'}, 'Missing reviewed complete candidate provenance')
    prior_manifest_path = checked_file(candidate, prior['manifest_path'])
    prior_metadata_path = checked_file(candidate, prior['metadata_path'])
    require(valid_hash(prior['manifest_sha256']) and sha256(prior_manifest_path) == prior['manifest_sha256']
            and valid_hash(prior['metadata_sha256']) and sha256(prior_metadata_path) == prior['metadata_sha256'],
            'Reviewed complete candidate snapshot changed')
    prior_manifest = read_json(prior_manifest_path)
    prior_metadata = read_json(prior_metadata_path)['metadata']
    require(prior_manifest.get('new_files') == manifest.get('new_files')
            and prior_manifest.get('inherited_files') == manifest.get('inherited_files')
            and prior_manifest.get('active_source_evidence') == manifest.get('active_source_evidence'),
            'Approved files or scientific completion evidence changed')
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
    changed_metadata = {'title', 'version', 'publication_date', 'notes'}
    require({k: v for k, v in current_metadata.items() if k not in changed_metadata}
            == {k: v for k, v in prior_metadata.items() if k not in changed_metadata},
            'Scientific scope, creators, license or related evidence changed')
    validate_target_receipt(manifest.get('target_reconciliation'), original['inherited_files'], current_metadata)
    validate_post_removal_proof(manifest, candidate, current_metadata)
    require(manifest.get('latest_published_version') == '2026.09.05-v24',
            'Preserved latest semantic version differs')
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
    require(version == CURRENT_CANDIDATE_VERSION
            and current_metadata.get('version') == CURRENT_CANDIDATE_VERSION
            and current_metadata.get('publication_date') == '2026-10-02',
            'Expected distinct Pacific-dated semantic checkpoint version')
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
