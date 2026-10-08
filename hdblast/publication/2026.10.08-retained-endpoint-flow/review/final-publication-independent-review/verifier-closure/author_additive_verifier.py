#!/usr/bin/env python3
"""Create a separate verifier; do not edit the previous publication packet."""
from pathlib import Path
import difflib
import hashlib
import json

HERE = Path(__file__).absolute().parent
W3 = HERE.parents[1]
OLD = Path('/workspace/HDblast/hdblast/publication/2026.10.07-stored-state-comparison')
OBS = W3 / 'publication-planning/read-only-observation'
source = (OLD / 'verify_publication.py').read_text()
old_sha = hashlib.sha256(source.encode()).hexdigest()
require_old = 'e2f8be2'
if not old_sha.startswith(require_old):
    raise ValueError('Previous verifier external source prefix differs')

replacements = [
    ("PRIOR = '23225288'", "PRIOR = '23228395'"),
    ("PRIOR_STREAM_PIN = '2664ec6a867d806e6146847fcf677d2e5c20061a400262acb8563d3ca552529c'", "PRIOR_STREAM_PIN = '8577ff0f006a2ef5ad75d4132066c2b0fa784a569d221c47a723fd1e21081e99'"),
    ("SCIENTIFIC_REPLAY_STATUS = 'PASS_COMPLETE_STORED_STATE_COMPARISON_REPLAY'", "SCIENTIFIC_REPLAY_STATUS = 'PASS_COMPLETE_LATER_STATE_AND_BULK_RESPONSE_REPLAY'\nEXPECTED_VERSION = '2026.10.08-later-state-and-bulk-response'\nEXPECTED_OWNER = '1386319'\nNOTES_HISTORY_MARKER = '\\n<h3>Preserved Notes from published record 23228395 (2026.10.07-stored-state-comparison)</h3>\\n'"),
    ("'latest_public_record', 'verified_public', 'controller_journal', 'prior_full_stream'}", "'latest_public_record', 'verified_public', 'controller_journal', 'prior_full_stream',\n         'new_draft_record', 'new_draft_files', 'verified_draft'}"),
    ("PRESERVED = {PRIOR: PARENT, '23114217': PARENT", "PRESERVED = {PRIOR: PARENT, '23225288': PARENT, '23114217': PARENT"),
    ("BASELINE_PINS = {", "BASELINE_PINS = {\n    '23228395': {\n        'baseline': '1f446fd3be07ea67ae2847dd0b09de6d844f62ca23ccf7523881cb65f2375c1c',\n        'baseline_files': '" + hashlib.sha256((OBS/'CURRENT_FILES.raw.json').read_bytes()).hexdigest() + "'},"),
    ("prior_stream, journal, preserved, expected_id, inventory_sha, metadata_sha):", "prior_stream, journal, preserved, expected_id, inventory_sha, metadata_sha,\n                  draft_record, draft_files, verified_draft):"),
    ("Completed stored-state scientific replay contract required", "Completed later-state and bulk-response scientific replay contract required"),
    ("pin_rows(inventory.get('inherited_files'), 27)", "pin_rows(inventory.get('inherited_files'), 29)"),
    ("448387915", "451598148"),
    ("inventory['expected_total_files_final'] == 29", "inventory['expected_total_files_final'] == 31"),
    ("'PASS_COMPLETE_29_FILE_EDITION_READBACK'", "'PASS_COMPLETE_31_FILE_EDITION_READBACK'"),
    ("verified['file_count'] == 29", "verified['file_count'] == 31"),
    ("len(rows) == 27", "len(rows) == 29"),
    ("prior27", "prior29"),
    ("PRIOR_FULL_STREAM_ATTESTATION_WITH_FRESH_GET_IDENTITY", "DATED_PRIOR29_COMPLETE_CONTENT_PLUS_FRESH_IMMUTABLE_IDENTITY"),
    ("len(content) == 29", "len(content) == 31"),
    ("Complete29 SHA256", "Complete31 SHA256"),
    ("all29", "all31"),
    ("All four preserved", "All five preserved"),
    ("'PASS_OFFLINE_29_FILE_PUBLICATION_WITNESSES'", "'PASS_OFFLINE_31_FILE_PUBLICATION_WITNESSES'"),
    ("'file_count': 29", "'file_count': 31"),
    ("'fresh_complete_public_streams': 29", "'fresh_complete_public_streams': 31"),
    ("'original_binary80_state_error': 'ENCLOSED_RETAINED_NODES'", "'later_saved_endpoint_error': 'ENCLOSED_RETAINED_NODES'"),
    ("'49,152 saved capsule-node occurrences; conditional fixed-weight transport under identical subsequent forcing/contacts'", "'49,152 saved capsule-node occurrences,98,304 endpoint comparisons and24 finite weighted cases; exact flow begins at the same represented incoming state; historical continuous path and full continuum pressure/contact/time/UV certificate remain unresolved'"),
    ("'full_twelve_case_certificate'", "'full_continuous_pressure_contact_certificate'"),
    ("expected_id, inventory_sha, metadata_sha)\n    result.update", "expected_id, inventory_sha, metadata_sha, data['new_draft_record'], data['new_draft_files'],\n        data['verified_draft'])\n    result.update"),
]
for before, after in replacements:
    if before not in source:
        raise ValueError('Missing exact adaptation needle: ' + before)
    source = source.replace(before, after)

extra = '''\ndef new_identity(record, expected_id, *, draft=False):
    require(type(record) is dict and record.get('id') == expected_id
            and record.get('parent', {}).get('id') == PARENT, 'Strict modern new-record identity differs')
    require(record.get('is_draft') is draft and record.get('is_published') is (not draft),
            'Strict modern new-record state differs')
    parent = record.get('parent', {})
    owner = parent.get('access', {}).get('owned_by', {}).get('user')
    require(type(owner) in (str, int) and str(owner) == EXPECTED_OWNER, 'Strict new-record owner differs')
    require(parent.get('pids', {}).get('doi', {}).get('identifier') == '10.5281/zenodo.' + PARENT,
            'Strict new-record concept DOI differs')
    pids = record.get('pids')
    require(type(pids) is dict, 'Strict new-record PIDs object required')
    if not draft or 'doi' in pids:
        require(type(pids.get('doi')) is dict and
                pids['doi'].get('identifier') == '10.5281/zenodo.' + expected_id,
                'Strict new-record version DOI differs')
    require(record.get('versions', {}).get('is_latest_draft' if draft else 'is_latest') is True,
            'Strict new-record latest state differs')


def preserved_edition_contract(metadata, prior):
    require(type(metadata) is dict and set(metadata) == {'metadata', 'custom_fields', 'access'},
            'Full frozen edition contract required')
    old = prior['metadata']; new = metadata['metadata']
    require(type(new) is dict and set(new) == set(old) and new.get('version') == EXPECTED_VERSION,
            'Exact new semantic version or inherited metadata fields differ')
    for key in set(old) - {'title', 'description', 'additional_descriptions', 'publication_date', 'version', 'related_identifiers'}:
        require(same(new[key], old[key], ('metadata', key)), 'Protected inherited metadata differs: ' + key)
    for key in ('custom_fields', 'access'):
        expected = prior[key]
        if key == 'access': expected = {k: v for k, v in expected.items() if k != 'status'}
        require(same(metadata[key], expected, (key,)), 'Protected inherited ' + key + ' differs')
    before = old['additional_descriptions']; after = new['additional_descriptions']
    require(type(after) is list and len(after) == len(before) > 0,
            'Preserved Notes list shape differs')
    require(same(after[1:], before[1:], ('metadata', 'additional_descriptions')),
            'Older Notes descriptions differ')
    require(type(after[0].get('description')) is str and
            after[0]['description'].endswith(NOTES_HISTORY_MARKER + before[0]['description']),
            'Complete prior Notes exact history suffix differs')
    require(set(after[0]) == set(before[0]), 'Current Notes field set differs')
    for key in set(before[0]) - {'description'}:
        require(same(after[0][key], before[0][key], ('metadata', 'additional_descriptions', 0, key)),
                'Current Notes type/language differs')
    previous = old.get('related_identifiers', []); related = new.get('related_identifiers', [])
    require(type(related) is list and len(related) in (len(previous), len(previous) + 1) and
            same(related[:len(previous)], previous, ('metadata', 'related_identifiers')),
            'Preserved DOI relations differ')
    if len(related) > len(previous):
        allowed = {'identifier': '10.5281/zenodo.' + PRIOR, 'relation_type': {'id': 'references'},
                   'resource_type': {'id': 'software'}, 'scheme': 'doi'}
        require(same(related[-1], allowed, ('metadata', 'related_identifiers', len(previous))),
                'Unreviewed additional DOI relation')


def draft_evidence(record, listing, receipt, pins, journal, expected_id, inventory_sha, metadata_sha, metadata):
    new_identity(record, expected_id, draft=True)
    metadata_equal(record, metadata)
    prefix = 'https://zenodo.org/api/records/' + expected_id + '/draft'
    require(record.get('links', {}).get('self') == prefix and
            record.get('links', {}).get('files') == prefix + '/files', 'Draft URL binding differs')
    files = file_map(listing); assert_pins(files, pins); assert_pins(file_map(record), pins)
    identities = {name: (immutable_file_id(row), immutable_version_id(row)) for name, row in files.items()}
    for name, row in file_map(record).items():
        require(immutable_file_id(row) == identities[name][0], 'Draft embedded file identity differs')
        if 'version_id' in row:
            require(immutable_version_id(row) == identities[name][1], 'Draft embedded version identity differs')
    require(receipt.get('status') == 'PASS_COMPLETE_31_FILE_EDITION_READBACK' and
            receipt.get('published') is False and receipt.get('draft_id') == expected_id and
            receipt.get('parent_id') == PARENT and receipt.get('inventory_sha256') == inventory_sha and
            receipt.get('metadata_sha256') == metadata_sha and type(receipt.get('file_count')) is int and
            receipt['file_count'] == 31 and type(receipt.get('total_bytes')) is int and
            receipt['total_bytes'] == sum(row['bytes'] for row in pins.values()), 'Draft complete verification receipt differs')
    etag = receipt.get('etag'); revision = record.get('revision_id')
    require(type(etag) is str and re.fullmatch(r'\"(0|[1-9][0-9]*)\"', etag) is not None and
            type(revision) is int and revision >= 0 and etag == '\"' + str(revision) + '\"',
            'Canonical draft ETag/integer revision differs')
    bases = {}
    require(type(receipt.get('content')) is list and len(receipt['content']) == 31,
            'Complete31 fresh draft basis rows required')
    for row in receipt['content']:
        name = row.get('filename')
        require(name in pins and name not in bases and row.get('sha256') == pins[name]['sha256'] and
                row.get('basis') == 'FRESH_COMPLETE_CONTENT_STREAM' and
                (row.get('remote_file_id'), row.get('remote_version_id')) == identities[name],
                'Draft fresh stream basis or identities differ')
        require(files[name].get('links', {}).get('content') == prefix + '/files/' + quote(name, safe='') + '/content',
                'Draft content URL differs')
        bases[name] = row
    publish = [i for i, event in enumerate(journal) if event.get('kind') == 'WRITE_INTENT' and
               event.get('operation') == 'publish']
    require(len(publish) == 1 and journal[publish[0]].get('url') == prefix + '/actions/publish',
            'Exactly one bound publish intent required')
    completed = []; current = None; last_stream_index = None; last_partial_index = -1
    for index, event in enumerate(journal[:publish[0]]):
        if event.get('kind') == 'GET' and event.get('url') == prefix + '/files':
            if current:
                if set(current) == set(pins): completed.append((current, last_stream_index))
                else: last_partial_index = last_stream_index
            require(event.get('status') == 200 and event.get('authenticated_request') is True,
                    'Draft file-list GET authentication/status differs')
            current = {}
        elif event.get('kind') == 'FULL_CONTENT_STREAM_VERIFIED' and current is not None:
            name = event.get('filename')
            require(name in pins and name not in current, 'Draft stream membership duplicate/unknown')
            require(canonical(event.get('observed')) == canonical({k: pins[name][k] for k in ('bytes', 'md5', 'sha256')}),
                    'Fresh draft content stream differs')
            current[name] = event
            last_stream_index = index
    if current:
        if set(current) == set(pins): completed.append((current, last_stream_index))
        else: last_partial_index = last_stream_index
    require(len(completed) >= 2, 'Preparation and immediately-before-publish complete31 draft stream rounds required')
    last_complete_index = completed[-1][1]
    require(last_partial_index < last_complete_index, 'Partial fresh draft round followed the last complete round')
    require(not any(event.get('kind') == 'WRITE_INTENT' for event in
                    journal[last_complete_index + 1:publish[0]]),
            'Mutation followed the final complete draft stream round')
    return len(completed)

'''
source = source.replace('\ndef verify_values(', extra + '\ndef verify_values(', 1)
source = source.replace("    identity(public_record, expected_id, PARENT)\n    identity(latest, expected_id, PARENT)",
                        "    new_identity(public_record, expected_id)\n    new_identity(latest, expected_id)\n    require(set(preserved) == set(PRESERVED), 'All five preserved record groups required')\n    preserved_edition_contract(metadata, preserved[PRIOR]['baseline'])\n    complete_draft_rounds = draft_evidence(draft_record, draft_files, verified_draft, pins, journal,\n        expected_id, inventory_sha, metadata_sha, metadata)")
source = source.replace("'fresh_complete_public_streams': 31,", "'fresh_complete_public_streams': 31, 'fresh_complete_draft_stream_rounds': complete_draft_rounds,\n            'fresh_complete_draft_streams_per_round': 31, 'publish_intents': 1,")
compile(source, str(HERE/'verify_publication.py'), 'exec')
(HERE/'verify_publication.py').write_text(source)
(HERE/'VERIFIER_ADAPTATION.diff').write_text(''.join(difflib.unified_diff(
    (OLD/'verify_publication.py').read_text().splitlines(True), source.splitlines(True),
    fromfile='previous/verify_publication.py', tofile='additive31/verify_publication.py')))
seeds = HERE/'manufactured-seeds'; seeds.mkdir(exist_ok=True)
seed_files = {name.name: name.read_bytes() for name in (OLD/'manufactured-seeds').iterdir() if name.is_file() and name.name != 'SEED_MANIFEST.json'}
seed_files.pop('PRIOR27_CONTENT_RECEIPT.json', None)
seed_files.update({'PRIOR29_CONTENT_RECEIPT.json': (OBS/'PRIOR29_CONTENT_RECEIPT.json').read_bytes(),
                   '23228395-record.json': (OBS/'CURRENT_PRIOR.raw.json').read_bytes(),
                   '23228395-files.json': (OBS/'CURRENT_FILES.raw.json').read_bytes()})
for name, raw in seed_files.items(): (seeds/name).write_bytes(raw)
(seeds/'SEED_MANIFEST.json').write_text(json.dumps({'schema_version': 1, 'scope': 'DATED_PUBLIC_SEEDS_ONLY_NO_NEW_PUBLICATION_CLAIM',
    'files': {name: {'bytes': len(raw), 'sha256': hashlib.sha256(raw).hexdigest()} for name, raw in sorted(seed_files.items())}}, indent=2, sort_keys=True) + '\n')
(HERE/'SOURCE_PREPARATION_RECEIPT.json').write_text(json.dumps({'schema_version': 1,
    'scope': 'ADDITIVE_OFFLINE_VERIFIER_PREPARATION_UNTESTED', 'old_source_sha256': old_sha,
    'source_sha256': hashlib.sha256(source.encode()).hexdigest(),
    'diff_sha256': hashlib.sha256((HERE/'VERIFIER_ADAPTATION.diff').read_bytes()).hexdigest(),
    'new_published_readback_claim': False, 'network_calls': 0, 'remote_writes': 0}, indent=2, sort_keys=True) + '\n')
