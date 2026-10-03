#!/usr/bin/env python3
"""Check a saved legacy Zenodo deposition response; never contacts or publishes."""
import argparse
import json
from pathlib import Path
import hashlib
import re
AUTHORIZED_OWNER = 1386319

class CandidateError(RuntimeError):
    pass

def require(condition, message):
    if not condition:
        raise CandidateError(message)

def read_json(path):
    def unique(items):
        result = {}
        for key, value in items:
            require(key not in result, 'Duplicate JSON key')
            result[key] = value
        return result
    return json.loads(path.read_text(), object_pairs_hook=unique)

def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def safe_file(root, relative):
    require(isinstance(relative, str) and relative != '', 'Invalid relative file path')
    require(not Path(relative).is_absolute() and all(x not in ('', '.', '..') for x in relative.split('/')) and '\\' not in relative, 'Unsafe relative file path')
    require(root.is_dir() and not root.is_symlink(), 'Invalid candidate root')
    target = root
    for part in Path(relative).parts:
        target = target / part
        require(not target.is_symlink(), 'Symlink candidate input')
    require(target.is_file() and target.resolve().is_relative_to(root.resolve()), 'Missing or escaped candidate input')
    return target

def validate_candidate(manifest, candidate):
    require(isinstance(manifest, dict) and manifest.get('schema_version') == 2, 'Wrong candidate schema')
    require(manifest.get('concept_doi') == '10.5281/zenodo.17088132', 'Wrong main family')
    require(type(manifest.get('active_main_draft')) is int and manifest['active_main_draft'] == 23114217, 'Wrong bound existing draft')
    require(type(manifest.get('latest_published_record')) is int and manifest['latest_published_record'] == 22347452, 'Wrong preserved published record')
    require(manifest.get('candidate_version') == '2026.10.03-verified-source-operator', 'Wrong candidate version')
    require(manifest.get('automatic_github_archiving') == 'ALL_EIGHT_OFF', 'Automatic integration state changed')
    require(manifest.get('old_metric_scientific_status') == 'FAIL', 'Historical metric FAIL was changed')
    require(manifest.get('ledger_classification') == 'LEDGER_ERROR_DEMONSTRATED', 'Historical ledger diagnosis changed')
    require(type(manifest.get('expected_final_files_count')) is int and manifest['expected_final_files_count'] == 23, 'Expected23completefiles')
    inherited, additions = manifest.get('inherited_files'), manifest.get('new_files')
    require(isinstance(inherited, list) and len(inherited) == 10 and isinstance(additions, list) and len(additions) == 13, 'Wrong inherited/addition cardinalities')
    names = set()
    for item in inherited + additions:
        require(isinstance(item, dict), 'Invalid file entry')
        name = item.get('filename')
        require(isinstance(name, str) and re.fullmatch('[A-Za-z0-9_.-]+', name) and name not in ('.', '..') and name not in names, 'Invalid/duplicate file name')
        names.add(name)
        require(type(item.get('bytes')) is int and item['bytes'] > 0, 'Invalid file size')
        require(isinstance(item.get('md5'), str) and re.fullmatch('[0-9a-f]{32}', item['md5']), 'Invalid file MD5')
        require(isinstance(item.get('sha256'), str) and re.fullmatch('[0-9a-f]{64}', item['sha256']), 'Invalid local SHA256')
    require(all(x.get('policy') == 'RETAIN_UNCHANGED_HISTORICAL_V24' for x in inherited), 'Inherited-file preservation policy differs')
    require(additions[-1]['filename'] == 'HDBLAST_CHECKPOINT_20261003_VERIFIED_SOURCE_OPERATOR.zip' and additions[-1]['bytes'] == 5331624 and additions[-1]['sha256'] == 'dfc8c6c6144fa1d9039390bcbaa7a9bd062355863c0fe1494d63824721193cd7', 'New registered package differs')
    require(manifest.get('expected_total_bytes') == sum(x['bytes'] for x in inherited + additions), 'Complete total byte count differs')
    require(manifest.get('expected_addition_bytes') == sum(x['bytes'] for x in additions), 'Addition byte count differs')
    proof = manifest.get('source_operator_evidence', {})
    require(proof.get('scientific_status') == 'PASS_UNIFORM_MODEL_AND_REGISTERED_PROBE_CERTIFICATE' and proof.get('fresh_status') == 'PASS_EXACT_FRESH_REPLAY', 'Scientific or replay status differs')
    require(proof.get('full_twelve_case_pressure_contact_certificate') == 'UNRESOLVED' and proof.get('external_mathematical_novelty') == 'NOT_ASSESSED', 'Scientific scope differs')
    require(manifest.get('separate_companion', {}).get('concept_doi') == '10.5281/zenodo.22922927' and manifest['separate_companion'].get('published_record') == 23111008, 'Separate companion family differs')
    request = manifest.get('metadata_request', {})
    require(request.get('path') == 'METADATA.json' and isinstance(request.get('sha256'), str) and re.fullmatch('[0-9a-f]{64}', request['sha256']), 'Metadata pin is absent')
    metadata_path = safe_file(candidate, request['path'])
    require(sha256(metadata_path) == request['sha256'], 'Metadata bytes changed')
    wrapper = read_json(metadata_path)
    require(isinstance(wrapper, dict) and set(wrapper) == {'metadata'} and isinstance(wrapper['metadata'], dict), 'Expected documented metadata wrapper')
    metadata = wrapper['metadata']
    require(metadata.get('version') == manifest['candidate_version'] and metadata.get('publication_date') == '2026-10-03' and metadata.get('upload_type') == 'software' and metadata.get('license') == 'cc-by-4.0', 'Metadata identity differs')
    require(all(metadata.get(x) for x in ('title', 'description', 'creators', 'notes')), 'Required metadata is empty')


def comparable(key, value):
    if key in ('keywords', 'references') and isinstance(value, list):
        return sorted(value)
    if key == 'related_identifiers' and isinstance(value, list):
        return sorted(json.dumps(item, sort_keys=True) for item in value)
    return value

def check_wrapper(saved, expected_record_id, public=False):
    failures = []
    endpoint = ('https://zenodo.org/api/records/' if public else
                'https://zenodo.org/api/deposit/depositions/') + str(expected_record_id)
    if not isinstance(saved, dict):
        return ['Response is not a JSON object'], {}
    if (saved.get('method') != 'GET' or saved.get('http_status') != 200
            or saved.get('requested_url') != endpoint or saved.get('final_url') != endpoint
            or saved.get('redirects') != [] or 'error' in saved):
        failures.append('Expected direct successful saved GET at the selected endpoint')
    if not public and saved.get('authenticated_request') is not True:
        failures.append('Deposition proof is not an authenticated GET')
    data = saved.get('data')
    if not isinstance(data, dict):
        failures.append('Response has no deposition/record object')
        data = {}
    if data.get('tombstone') or data.get('status') == 410:
        failures.append('Record has a removal tombstone')
    return failures, data


def check_files(actual, expected_files, public=False):
    failures = []
    if not isinstance(actual, list):
        return ['Remote file list is absent']
    by_name = {}
    for item in actual:
        if not isinstance(item, dict):
            failures.append('Malformed remote file record')
            continue
        name = item.get('key') if public else item.get('filename')
        if not isinstance(name, str) or name in by_name:
            failures.append('Malformed or duplicate remote filename')
            continue
        by_name[name] = item
    wanted = {item['filename'] for item in expected_files}
    if len(expected_files) != len(wanted):
        failures.append('Duplicate expected publication filename')
    if set(by_name) != wanted:
        failures.append('Remote inventory differs: missing=' + ','.join(sorted(wanted - set(by_name)))
                        + '; extra=' + ','.join(sorted(set(by_name) - wanted)))
    for item in expected_files:
        remote = by_name.get(item['filename'])
        if remote is None:
            continue
        sizes = [remote[key] for key in ('filesize', 'size', 'bytes') if key in remote]
        if not sizes or any(type(value) is not int or value != item['bytes'] for value in sizes):
            failures.append('Remote byte count differs: ' + item['filename'])
        checksums = [remote[key] for key in ('checksum', 'md5') if key in remote]
        if not checksums or any(type(value) is not str or value.lower().removeprefix('md5:') != item['md5'] for value in checksums):
            failures.append('Remote MD5 differs: ' + item['filename'])
    return failures


def public_metadata(metadata):
    result = dict(metadata) if isinstance(metadata, dict) else {}
    if isinstance(result.get('license'), dict):
        result['license'] = result['license'].get('id')
    if 'upload_type' not in result and isinstance(result.get('resource_type'), dict):
        result['upload_type'] = result['resource_type'].get('type')
    return result


def verify(saved, manifest, request, expected_record_id, published=False,
           public_saved=None, metadata_only=False):
    if not isinstance(saved, dict):
        return {'status': 'FAIL_SAVED_RECORD_VALIDATION', 'failures': ['Response is not a JSON object']}
    failures, data = check_wrapper(saved, expected_record_id)
    if type(data.get('id')) is not int or type(expected_record_id) is not int or data.get('id') != expected_record_id:
        failures.append('Record ID differs from the explicitly selected existing draft/record')
    if type(expected_record_id) is not int or expected_record_id != manifest.get('active_main_draft'):
        failures.append('Selected target differs from the bound existing main-family draft')
    if not published and (type(expected_record_id) is not int or
            expected_record_id != manifest.get('active_main_draft')):
        failures.append('Unpublished target differs from the bound existing main-family draft')
    if str(data.get('conceptrecid')) != '17088132' or data.get('conceptdoi') != '10.5281/zenodo.17088132':
        failures.append('Record belongs to the wrong concept DOI family')
    if type(data.get('owner')) is not int or data.get('owner') != AUTHORIZED_OWNER:
        failures.append('Record owner differs from the authorized owner')
    if published and metadata_only:
        failures.append('Metadata-only checking cannot verify publication')
    if published:
        if data.get('submitted') is not True or data.get('state') != 'done':
            failures.append('Record is not submitted and done; a reserved DOI is insufficient')
        if data.get('doi') != '10.5281/zenodo.' + str(expected_record_id):
            failures.append('Published record DOI is absent or differs from its selected record ID')
    elif data.get('submitted') is not False or data.get('state') != 'unsubmitted':
        failures.append('Expected the selected unsubmitted draft')
    metadata = data.get('metadata', {})
    if not isinstance(metadata, dict):
        metadata = {}
        failures.append('Remote metadata is not an object')
    mismatches = []
    for key, expected in request['metadata'].items():
        if comparable(key, metadata.get(key)) != comparable(key, expected):
            mismatches.append(key)
    if mismatches:
        failures.append('Saved metadata differs in: ' + ', '.join(sorted(mismatches)))
    complete_files = manifest['inherited_files'] + manifest['new_files']
    if len(complete_files) != 23 or len({f['filename'] for f in complete_files}) != 23:
        failures.append('Expected inventory is not twenty-three unique files')
    expected_files = manifest['inherited_files'] if metadata_only else complete_files
    failures.extend(check_files(data.get('files'), expected_files))
    public_count = None
    if published:
        public_failures, public_data = check_wrapper(public_saved, expected_record_id, public=True)
        failures.extend(['Public witness: ' + value for value in public_failures])
        if (type(public_data.get('id')) is not int or public_data.get('id') != expected_record_id
                or public_data.get('submitted') is not True or public_data.get('state') != 'done'
                or public_data.get('doi') != '10.5281/zenodo.' + str(expected_record_id)
                or public_data.get('conceptdoi') != '10.5281/zenodo.17088132'
                or str(public_data.get('conceptrecid')) != '17088132'):
            failures.append('Public witness identity/family/publication state differs')
        links = public_data.get('links', {})
        if (not isinstance(links, dict)
                or links.get('self') != 'https://zenodo.org/api/records/' + str(expected_record_id)
                or links.get('parent') != 'https://zenodo.org/api/records/17088132'
                or links.get('parent_doi') != 'https://doi.org/10.5281/zenodo.17088132'
                or links.get('latest') != 'https://zenodo.org/api/records/' + str(expected_record_id) + '/versions/latest'):
            failures.append('Public witness parent/latest record links differ')
        normalized = public_metadata(public_data.get('metadata'))
        for key, expected in request['metadata'].items():
            if comparable(key, normalized.get(key)) != comparable(key, expected):
                failures.append('Public witness metadata differs: ' + key)
        failures.extend(['Public witness: ' + value for value in check_files(public_data.get('files'), complete_files, public=True)])
        if isinstance(public_data.get('files'), list):
            public_count = len(public_data['files'])
    actual = data.get('files') if isinstance(data.get('files'), list) else []
    if not failures:
        status = ('PASS_VERIFIED_PUBLISHED_RECORD' if published else
                  'PASS_EXACT_SAVED_DRAFT_METADATA' if metadata_only else 'PASS_COMPLETE_SAVED_DRAFT')
    else:
        status = 'FAIL_SAVED_RECORD_VALIDATION'
    return {'schema_version': 2, 'status': status, 'record_id': expected_record_id, 'concept_doi': '10.5281/zenodo.17088132', 'expected_files': len(expected_files), 'complete_expected_files': 23, 'observed_files': len(actual), 'public_observed_files': public_count, 'metadata_only': metadata_only, 'publication_ready': status == 'PASS_COMPLETE_SAVED_DRAFT', 'metadata_mismatch_fields': sorted(mismatches), 'failures': failures, 'scope': 'Saved response comparison only; no network, upload or publish action.', 'publication_performed_by_this_tool': False}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--saved-response', type=Path, required=True)
    parser.add_argument('--expected-record-id', type=int, required=True)
    parser.add_argument('--expected-manifest-sha256', required=True)
    parser.add_argument('--published', action='store_true')
    parser.add_argument('--public-response', type=Path,
                        help='Direct public GET200/no tombstone witness required with --published')
    parser.add_argument('--metadata-only', action='store_true',
                        help='Require exact saved metadata and the unchanged ten inherited files before uploads')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    candidate = Path(__file__).resolve().parent
    manifest_path = safe_file(candidate, 'FILE_MANIFEST.json')
    require(re.fullmatch('[0-9a-f]{64}', args.expected_manifest_sha256) is not None, 'External manifest SHA256 pin is required')
    require(sha256(manifest_path) == args.expected_manifest_sha256, 'Candidate manifest hash changed')
    manifest = read_json(manifest_path)
    require(manifest['concept_doi'] == '10.5281/zenodo.17088132', 'Wrong expected DOI family')
    validate_candidate(manifest, candidate)
    metadata_path = safe_file(candidate, manifest['metadata_request']['path'])
    require(sha256(metadata_path) == manifest['metadata_request']['sha256'], 'Pinned metadata request changed')
    report = verify(read_json(args.saved_response), manifest, read_json(metadata_path), args.expected_record_id, args.published,
                    read_json(args.public_response) if args.public_response else None, args.metadata_only)
    args.output.open('x').write(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': report['status'], 'observed_files': report['observed_files'], 'metadata_mismatch_fields': report['metadata_mismatch_fields']}))
    raise SystemExit(0 if report['status'].startswith('PASS_') else 1)

if __name__ == '__main__':
    main()
