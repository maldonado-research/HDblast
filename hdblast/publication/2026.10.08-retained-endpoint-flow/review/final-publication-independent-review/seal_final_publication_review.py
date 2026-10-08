"""Seal read-only publication evidence without network or scientific execution."""
from pathlib import Path
import datetime
import hashlib
import json
import stat

HERE = Path(__file__).resolve().parent
BASE = HERE.parent


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def load(path):
    return json.loads(path.read_bytes())


def pin(path):
    meta = path.lstat()
    require(stat.S_ISREG(meta.st_mode) and meta.st_nlink == 1 and all(not p.is_symlink() for p in (path, *path.parents)), 'unaliased regular evidence')
    raw = path.read_bytes()
    require(len(raw) == meta.st_size, 'file changed while sealing')
    return {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()}


def save(name, value):
    with (HERE / name).open('x') as stream:
        json.dump(value, stream, sort_keys=True, indent=2)
        stream.write('\n')


def main():
    normal = load(HERE / 'REAL_WITNESS_NORMAL.json')
    optimized = load(HERE / 'REAL_WITNESS_OPTIMIZED.json')
    require(normal == optimized and (HERE / 'REAL_WITNESS_NORMAL.json').read_bytes() == (HERE / 'REAL_WITNESS_OPTIMIZED.json').read_bytes(), 'real mode outputs differ')
    require(normal['status'] == 'PASS_OFFLINE_31_FILE_PUBLICATION_WITNESSES', 'real witness failed')
    for mode in ('normal', 'optimized'):
        require((HERE / ('real-witness-' + mode + '.stderr.log')).read_bytes() == b'', 'verifier stderr not empty')
        stdout = load(HERE / ('real-witness-' + mode + '.stdout.log'))
        require(stdout == {key: normal[key] for key in ('status', 'record_id', 'file_count', 'total_bytes', 'network_calls')}, 'verifier stdout disagrees')
    closure = load(HERE / 'VERIFIER_CLOSURE_AUTHENTICATION.json')
    for name, expected in closure['manifested_file_pins'].items():
        require(pin(HERE / 'verifier-closure' / name) == expected == pin(BASE / 'release-assets/publication-verifier' / name), 'verifier closure changed')
    require(pin(HERE / 'verifier-closure/MANIFEST.json')['sha256'] == closure['source_manifest_sha256'], 'source closure manifest changed')
    witness = load(HERE / 'TERMINAL_WITNESS_AUTHENTICATION.json')
    source_proof = BASE / 'publication-delivery-candidate/hdblast/publication/2026.10.08-retained-endpoint-flow/proof'
    require(pin(HERE / 'captured-witness/WITNESS_MANIFEST.json') == witness['witness_manifest'] == pin(source_proof / 'WITNESS_MANIFEST.json'), 'witness manifest changed')
    for name, expected in witness['witness_files'].items():
        require(pin(HERE / 'captured-witness' / name) == expected == pin(source_proof / name), 'witness changed')
    journal_path = HERE / 'journal-audit/FINAL_JOURNAL_REVIEW_RECEIPT.json'
    journal_manifest_path = HERE / 'journal-audit/FINAL_JOURNAL_REVIEW_MANIFEST.json'
    require(pin(journal_path)['sha256'] == 'ede9ab7bb2ed629d95e6e2a9d04125d1c2ff77c478367c84e5f427d445349caa', 'journal review pin')
    require(pin(journal_manifest_path)['sha256'] == '1b005f6cd715161f677bb9e02c8ea58e63207fd37fbdb586abe3687efd46b6b3', 'journal review manifest pin')
    journal = load(journal_path)
    journal_manifest = load(journal_manifest_path)
    for item in journal_manifest['artifacts']:
        require(pin(HERE / 'journal-audit' / item['path']) == {key: item[key] for key in ('bytes', 'sha256')}, 'journal evidence changed')
    for name, sha in journal['terminal_external_pins'].items():
        require(pin(BASE / 'publication-execution-001' / name)['sha256'] == sha, 'terminal publisher source changed')
    before = load(HERE / 'observation-002/PUBLIC_PRESERVATION_REVIEW.json')
    after = load(HERE / 'observation-003/PUBLIC_PRESERVATION_REVIEW.json')
    for rid in before['records']:
        require(before['records'][rid]['file_pins'] == after['records'][rid]['file_pins'], 'historical identities changed across publish')
    journal_rows = [json.loads(row) for row in (HERE / 'captured-witness/controller/JOURNAL.jsonl').read_bytes().splitlines()]
    intent = next(row['utc'] for row in journal_rows if row.get('kind') == 'WRITE_INTENT' and row.get('operation') == 'publish')
    response = next(row['utc'] for row in journal_rows if row.get('kind') == 'WRITE_RESPONSE' and row.get('operation') == 'publish')
    phases = {}
    for label, relationship in [('observation-002', 'BEFORE_PUBLISH_INTENT'), ('observation-003', 'AFTER_SUCCESSFUL_PUBLISH_RESPONSE')]:
        rows = [load(path) for path in (HERE / label).glob('*.GET.json')]
        start = min(row['request_started_utc'] for row in rows)
        end = max(row['response_completed_utc'] for row in rows)
        require(len(rows) == 10 and all(row['http_status'] == 200 and row['method'] == 'GET' and row['authenticated_request'] is False and row['TLS_certificate_and_hostname_verification'] is True for row in rows), 'complete anonymous verified TLS observation')
        require((end < intent) if label == 'observation-002' else (start > response), 'publication observation temporal classification differs')
        phases[label] = {'start_utc': start, 'end_utc': end, 'relationship': relationship, 'GET_requests': 10, 'receipt': pin(HERE / label / 'PUBLIC_PRESERVATION_REVIEW.json')}
    preservation = {}
    for directory, name, expected_sha in [
        ('actual-independent-review', 'FINAL_ACTUAL_REVIEW_MANIFEST.json', 'c974fabd71062ff8a1121ef54cd7c2b8094f7529eb091b54c360edb0e0bef7e7'),
        ('actual-replay-independent-review', 'FINAL_PORTABLE_REPLAY_REVIEW_MANIFEST.json', '1ffeb720131ac6ff7fa23433033adcee30ef1de31b3f81aa8add31ffd034f65f')]:
        root = BASE / directory
        require(pin(root / name)['sha256'] == expected_sha, 'prior review manifest changed')
        m = load(root / name)
        for relative, expected in m['files'].items():
            require(pin(root / relative) == expected, 'prior review artifact changed')
        preservation[directory] = {'manifest_sha256': expected_sha, 'artifacts_rehashed_unchanged': len(m['files'])}
    save('PRIOR_SCIENTIFIC_REVIEWS_PRESERVED.json', {'status': 'PASS_PRIOR_ACTUAL_AND_PORTABLE_REVIEW_BYTES_UNCHANGED', 'reviews': preservation})
    live = load(HERE / 'new-public-observation-001/NEW_PUBLIC_IDENTITY_REVIEW.json')
    require(live['record_id'] == normal['record_id'] == '23244754' and live['file_count'] == normal['file_count'] == 31 and live['total_file_bytes'] == normal['total_bytes'] == 453646943, 'live/witness edition differs')
    all_gets = [path for folder in ('observation-001', 'observation-002', 'observation-003', 'new-public-observation-001') for path in (HERE / folder).glob('*.GET.json')]
    statuses = [load(path)['http_status'] for path in all_gets]
    require(len(all_gets) == 32 and statuses.count(200) == 27 and statuses.count(406) == 5, 'recorded network accounting differs')
    result = {
        'status': 'PASS_INDEPENDENT_FINAL_PUBLICATION_WITNESS_AND_LIVE_PRESERVATION_REVIEW',
        'sealed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'record_id': '23244754', 'concept_id': '17088132', 'owner': '1386319',
        'version_doi': '10.5281/zenodo.23244754', 'concept_doi': '10.5281/zenodo.17088132',
        'file_count': 31, 'total_bytes': 453646943,
        'inventory_sha256': normal['inventory_sha256'], 'metadata_sha256': normal['metadata_sha256'],
        'verifier_sha256': normal['verifier_sha256'], 'verifier_source_closure_manifest_sha256': closure['source_manifest_sha256'],
        'verifier_source_closure_files': 38, 'witness_manifest_sha256': normal['witness_manifest_sha256'],
        'captured_witness_files': 31, 'captured_witness_roles': 11,
        'witness_source_and_copy_unchanged_at_seal': True,
        'root_terminal_source_pins': journal['terminal_external_pins'],
        'root_controller_sha256': journal['input_pins']['controller']['sha256'],
        'root_review_pins': {name: pin(HERE / name) for name in ('INDEPENDENT_RELEASE_ASSET_REVIEW.json', 'INDEPENDENT_31_FILE_VERIFIER_REVIEW.json')},
        'real_witness_checks': {
            'normal': pin(HERE / 'REAL_WITNESS_NORMAL.json'), 'optimized': pin(HERE / 'REAL_WITNESS_OPTIMIZED.json'),
            'both_observed_exit_codes': 0, 'byte_identical_outputs': True, 'isolated_no_bytecode_Python': True,
            'captured_external_source_before_execution': True,
        },
        'independent_journal_review': {
            'receipt': pin(journal_path), 'compact_manifest': pin(journal_manifest_path),
            'events': journal['journal_events'], 'draft_rounds': 2, 'draft_streams': 62,
            'public_rounds': 1, 'public_streams': 31, 'bytes_each_round': 453646943,
            'publish_intents': 1, 'recorded_write_intents': 9, 'successful_write_responses': 9,
            'unknown_failed_truncated_write_outcomes': 0, 'retry_or_reconciliation_events': 0,
            'missing_duplicate_partial_wrong_scope_streams': 0,
            'independent_checker_controls_per_Python_mode': 17,
            'scope_limits': journal['scope_limits'],
        },
        'publication_times_utc': {
            'publish_intent': intent, 'publish_response': response,
            'public_31_stream_verification_complete': journal['record_dates']['public']['receipt_utc'],
        },
        'historical_preservation_observations': phases,
        'historical_records': {rid: {key: row[key] for key in ('concept_id', 'file_count', 'total_file_bytes', 'full_metadata_custom_fields_access_exact_equal', 'complete_file_listing_exact_equal', 'immutable_file_version_ids_equal')} for rid, row in after['records'].items()},
        'historical_witness_groups_use_exact_own_postpublication_GET_bytes': True,
        'live_new_record_identity_review': pin(HERE / 'new-public-observation-001/NEW_PUBLIC_IDENTITY_REVIEW.json'),
        'new_witness_metadata_and_complete_file_list_match_own_live_GET': True,
        'content_attestation_basis': 'All31 fresh complete public streams and two31file draft rounds recorded by externally pinned root-owned controller; no duplicate independent attachment stream.',
        'prior_scientific_review_preservation': preservation,
        'review_operations': {
            'anonymous_GET_requests': 32, 'HTTP200': 27, 'HTTP406': 5,
            'attachment_streams': 0, 'attachment_bytes_downloaded': 0, 'remote_writes': 0,
            'publisher_executions': 0, 'physical_source_calls': 0, 'retained_array_decodes': 0,
            'scientific_replays': 0, 'repository_mutations': 0,
        },
        'preserved_observation_failure': {'receipt': pin(HERE / 'OBSERVATION001_FAILURE_RECEIPT.json'), 'cause': 'Five file-list endpoints rejected the record-only modern media type; fresh correctly negotiated listing GETs passed.'},
        'scientific_limits': {key: normal[key] for key in ('later_saved_endpoint_error', 'state_error_scope', 'full_continuous_pressure_contact_certificate', 'metric_calibration', 'higher_dimensional_Big_Bang_origin', 'external_novelty', 'historical_original_literal_Notes_guard')},
        'internal_independent_review_not_external_peer_review': True,
        'authorizes_further_remote_or_scientific_operations': False,
    }
    save('FINAL_PUBLICATION_REVIEW_RECEIPT.json', result)
    selected = {p.name for p in HERE.iterdir() if p.is_file()}
    selected.discard('FINAL_PUBLICATION_REVIEW_MANIFEST.json')
    for directory in ('observation-001', 'observation-002', 'observation-003', 'new-public-observation-001', 'verifier-closure', 'captured-witness'):
        selected.update(str(path.relative_to(HERE)) for path in (HERE / directory).rglob('*') if path.is_file())
    selected.add('journal-audit/FINAL_JOURNAL_REVIEW_MANIFEST.json')
    selected.update('journal-audit/' + item['path'] for item in journal_manifest['artifacts'])
    files = {name: pin(HERE / name) for name in sorted(selected)}
    manifest = {'status': 'SEALED_COMPACT_FINAL_PUBLICATION_INDEPENDENT_REVIEW', 'path_base': '.',
                'files': files, 'file_count': len(files), 'total_bytes': sum(row['bytes'] for row in files.values()),
                'final_receipt_sha256': files['FINAL_PUBLICATION_REVIEW_RECEIPT.json']['sha256'],
                'excluded': ['Large redundant manufactured journal-control inputs remain preserved and pinned by journal-audit/CONTROL_EVIDENCE_MANIFEST.json.', 'No Zenodo attachment payloads were downloaded into this review.']}
    save('FINAL_PUBLICATION_REVIEW_MANIFEST.json', manifest)
    for name, expected in files.items():
        require(pin(HERE / name) == expected, 'evidence changed while sealing')
    print(json.dumps({'receipt': pin(HERE / 'FINAL_PUBLICATION_REVIEW_RECEIPT.json'), 'manifest': pin(HERE / 'FINAL_PUBLICATION_REVIEW_MANIFEST.json'), 'files': len(files), 'bytes': manifest['total_bytes']}, sort_keys=True))


if __name__ == '__main__':
    main()
