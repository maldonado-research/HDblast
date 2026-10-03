#!/usr/bin/env python3
"""Check a saved legacy Zenodo deposition response; never contacts or publishes."""
import argparse
import json
from pathlib import Path
from prepare_candidate import CandidateError, require, safe_file, sha256
from extension_inventory import validate_extension_inventory, read_json

def comparable(key, value):
    if key in ('keywords', 'references') and isinstance(value, list):
        return sorted(value)
    if key == 'related_identifiers' and isinstance(value, list):
        return sorted(json.dumps(item, sort_keys=True) for item in value)
    return value

def verify(saved, manifest, request, expected_record_id, published=False):
    if not isinstance(saved, dict):
        return {'status': 'FAIL_SAVED_RECORD_VALIDATION', 'failures': ['Response is not a JSON object']}
    data = saved.get('data', saved)
    failures = []
    if not isinstance(data, dict):
        return {'status': 'FAIL_SAVED_RECORD_VALIDATION', 'failures': ['Response has no deposition object']}
    if type(data.get('id')) is not int or type(expected_record_id) is not int or data.get('id') != expected_record_id:
        failures.append('Record ID differs from the explicitly selected existing draft/record')
    if not published and (type(expected_record_id) is not int or
            expected_record_id != manifest.get('active_main_draft')):
        failures.append('Unpublished target differs from the bound existing main-family draft')
    if str(data.get('conceptrecid')) != '17088132' or data.get('conceptdoi') != '10.5281/zenodo.17088132':
        failures.append('Record belongs to the wrong concept DOI family')
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
    expected_files = manifest['inherited_files'] + manifest['new_files']
    if len(expected_files) != 22 or len({f['filename'] for f in expected_files}) != 22:
        failures.append('Expected inventory is not twenty-two unique files')
    actual = data.get('files', [])
    if not isinstance(actual, list):
        actual = []
        failures.append('Remote file list is absent')
    by_name = {}
    for item in actual:
        if not isinstance(item, dict) or not isinstance(item.get('filename'), str):
            failures.append('Malformed remote file record')
            continue
        if item['filename'] in by_name:
            failures.append('Duplicate remote filename: ' + item['filename'])
        by_name[item['filename']] = item
    wanted = {item['filename'] for item in expected_files}
    if set(by_name) != wanted:
        failures.append('Remote inventory differs: missing=' + ','.join(sorted(wanted - set(by_name))) + '; extra=' + ','.join(sorted(set(by_name) - wanted)))
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
    if not failures:
        status = 'PASS_VERIFIED_PUBLISHED_RECORD' if published else 'PASS_COMPLETE_SAVED_DRAFT'
    else:
        status = 'FAIL_SAVED_RECORD_VALIDATION'
    return {'schema_version': 1, 'status': status, 'record_id': expected_record_id, 'concept_doi': '10.5281/zenodo.17088132', 'expected_files': 22, 'observed_files': len(actual), 'metadata_mismatch_fields': sorted(mismatches), 'failures': failures, 'scope': 'Saved response comparison only; no network, upload or publish action.', 'publication_performed_by_this_tool': False}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--saved-response', type=Path, required=True)
    parser.add_argument('--expected-record-id', type=int, required=True)
    parser.add_argument('--expected-manifest-sha256', required=True)
    parser.add_argument('--published', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    candidate = Path(__file__).absolute().parent
    manifest_path = safe_file(candidate, 'FILE_MANIFEST.json')
    require(sha256(manifest_path) == args.expected_manifest_sha256, 'Candidate manifest hash changed')
    manifest = read_json(manifest_path)
    require(manifest['concept_doi'] == '10.5281/zenodo.17088132', 'Wrong expected DOI family')
    validate_extension_inventory(manifest, candidate)
    metadata_path = safe_file(candidate, manifest['metadata_request']['path'])
    require(sha256(metadata_path) == manifest['metadata_request']['sha256'], 'Pinned metadata request changed')
    report = verify(read_json(args.saved_response), manifest, read_json(metadata_path), args.expected_record_id, args.published)
    args.output.open('x').write(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': report['status'], 'observed_files': report['observed_files'], 'metadata_mismatch_fields': report['metadata_mismatch_fields']}))
    raise SystemExit(0 if report['status'].startswith('PASS_') else 1)

if __name__ == '__main__':
    main()
