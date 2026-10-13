#!/usr/bin/env python3
"""Authenticate saved publication witnesses offline, not present server state.

Adapted from the 11 October public checker (SHA256
635adb7beff6047ed5d4e32cdd5337effb0adc8c162cc003d7a2d9c9b74cbfcf).
No networking, third-party imports, assertions, or source execution.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import re

CONCEPT = '17088132'
PRIOR = '23294979'
COUNT, TOTAL = 35, 454964066
INHERITED_COUNT, INHERITED_BYTES = 33, 454129865
INVENTORY_SHA = '98f7db2f79c8fa29e0552a6b97cd923d2f25a2cd907d85d5893c370aabcc0454'
METADATA_SHA = 'debbaf5a0fbe50d886d303940c348cfe4aa52d23041d519b68804c03d6bb7e80'
PRIOR_INVENTORY_SHA = '2395097106b2d22a091ff9fedc9342b86eea3cfd559943f5fac0edd133f41446'
SCIENCE_SHA = 'd1d47c66cef1ed79979c22ae7145dc66d35fa35be69dba900f9912ac01f31e2f'
REPLAY_SHA = '228a305cc2d84321cf861ff5d685e11ccd16107cae6f0e77a342db7c7ab28953'
COMMIT = '401efdfaaef657bcbdaec588b21b82789daddd0e'
CAPSULE = 'research/HDBLAST_REGULAR_SHELL_SPECTRA_20261013'
REPOSITORY = 'maldonado-research/HDblast'
SOURCE_URL = f'https://github.com/{REPOSITORY}/tree/{COMMIT}/{CAPSULE}'
VERSION = '2026.10.13-regular-shell-spectra'
TITLE = 'HDBLAST: Conditional Scalar Spectral Positivity and Tensor Mass Thresholds on Regular de Sitter Shells'
GROUPS = {rid: CONCEPT for rid in ('23294979','23244754','23228395','23225288','23114217','22347452')}
GROUPS['23111008'] = '22922927'
NEW_ASSETS = [
    {'filename':'HDBLAST_REGULAR_SHELL_SPECTRA_20261013.zip','bytes':830422,
     'md5':'8c0e82a903a220cca6e3921228f51abb','sha256':'e15bf5bf6ba4432e230cf411350314f3122bbec61778d2d9c34b81ed27dc5e44'},
    {'filename':'HDBLAST_REGULAR_SHELL_SPECTRA_OVERVIEW_20261013.md','bytes':3779,
     'md5':'7299777e115acd77c268bdd596bb1031','sha256':'4b73fcc429fa254dc2ce29beb6009c1abac111dc427ba21f2d79b08d4946eab3'},
]
VOCABULARY_PATHS = {
    ('metadata','resource_type'),('metadata','rights'),('metadata','languages'),
    ('metadata','creators','role'),('metadata','contributors','role'),
    ('metadata','related_identifiers','relation_type'),
    ('metadata','related_identifiers','resource_type'),
    ('metadata','additional_descriptions','type'),('metadata','additional_descriptions','lang'),
    ('custom_fields','code:programmingLanguage'),('custom_fields','code:developmentStatus'),
}

def require(test, message):
    if not test:
        raise ValueError(message)

def sha(data):
    return hashlib.sha256(data).hexdigest()

def no_duplicates(pairs):
    obj = {}
    for key, value in pairs:
        require(key not in obj, 'DUPLICATE_JSON_KEY_'+key)
        obj[key] = value
    return obj

def parse(raw):
    return json.loads(raw, object_pairs_hook=no_duplicates,
                      parse_constant=lambda token: (_ for _ in ()).throw(ValueError('NONFINITE_JSON_'+token)))

def read(path):
    return parse(path.read_bytes())

def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False)

def exact(a, b):
    return canonical(a) == canonical(b)

def hexadecimal(value, length=64):
    return type(value) is str and re.fullmatch('[a-f0-9]{'+str(length)+'}', value) is not None

def safe_relative(name):
    return (type(name) is str and bool(name) and '\\' not in name
            and not Path(name).is_absolute() and '..' not in Path(name).parts
            and str(Path(name)) == name)

def timestamp(value):
    require(type(value) is str, 'EVENT_TIMESTAMP_TYPE')
    value = datetime.datetime.fromisoformat(value)
    require(value.tzinfo is not None, 'EVENT_TIMEZONE_REQUIRED')
    return value

def filemap(obj):
    rows = obj['entries'] if 'entries' in obj else obj['files']['entries']
    if type(rows) is dict:
        require(all(type(key) is str and type(value) is dict and value.get('key') == key
                    for key, value in rows.items()), 'FILE_MAPPING_KEY_DIFFERS')
        rows = list(rows.values())
    require(type(rows) is list, 'INVALID_FILE_LIST')
    names = [r['key'] for r in rows]
    require(all(type(n) is str for n in names) and len(names) == len(set(names)), 'DUPLICATE_FILE_NAMES')
    return dict(zip(names, rows))

def normalize(value, path=()):
    # Only display decorations of explicitly enumerated vocabulary objects vary.
    # Scientific rich text is compared byte-for-byte as JSON strings.
    if type(value) is dict:
        drop = set()
        words = tuple(x for x in path if type(x) is not int)
        if words in VOCABULARY_PATHS:
            require(type(value.get('id')) is str and set(value) <= {'id','title','description','icon','props'},
                    'UNRECOGNIZED_VOCABULARY_SHAPE')
            drop = {'title','description','icon','props'}
        return {k: normalize(v, path+(k,)) for k, v in value.items() if k not in drop}
    if type(value) is list:
        return [normalize(v, path+(i,)) for i, v in enumerate(value)]
    return value

def compare_metadata(record, expected):
    require(type(expected) is dict and set(expected) == {'metadata','custom_fields','access'}, 'METADATA_BODY_FIELDS_DIFFER')
    for key in ('metadata','custom_fields'):
        require(exact(normalize(record[key], (key,)), normalize(expected[key], (key,))), 'COMPLETE_'+key.upper()+'_DIFFERS')
    access = dict(record['access'])
    status = access.pop('status', None)
    require(status in (None, 'open') and exact(access, expected['access']), 'ACCESS_DIFFERS')

def verify_source(source, manifest_raw, replay_raw):
    require(sha(manifest_raw) == SCIENCE_SHA, 'SCIENTIFIC_MANIFEST_PIN_DIFFERS')
    require(sha(replay_raw) == REPLAY_SHA, 'SCIENTIFIC_REPLAY_PIN_DIFFERS')
    manifest, replay = parse(manifest_raw), parse(replay_raw)
    require(source['status'] == 'PASS_IMMUTABLE_PUBLIC_GIT_FULL_BYTE_READBACK', 'SOURCE_SCHEMA_STATUS_DIFFERS')
    require(source['repository'] == REPOSITORY and source['commit'] == COMMIT
            and source['capsule_path'] == CAPSULE and source['public_source_url'] == SOURCE_URL, 'IMMUTABLE_SOURCE_IDENTITY_DIFFERS')
    require(source['complete_file_count'] == 130 and type(source['complete_file_count']) is int
            and source['complete_bytes'] == 1670878 and type(source['complete_bytes']) is int, 'SOURCE_COUNT_BYTES_DIFFERS')
    for key in ('all_files_full_stream_sha256_and_direct_byte_matches','manifest_byte_match',
                'public_git_tree_inventory_exact','workflow_only_update_preserved_scientific_subtree',
                'ssl_certificate_verification_enabled','inherited_proxy_and_CA_preserved'):
        require(source[key] is True, 'SOURCE_REQUIRED_TRUE_'+key)
    require(source['explicit_authorization_headers_or_credentials_used'] is False, 'SOURCE_NOT_ANONYMOUS')
    require(source['manifest_sha256'] == SCIENCE_SHA and source['frozen_science_replay_receipt_sha256'] == REPLAY_SHA, 'SOURCE_SCIENCE_BINDING_DIFFERS')
    require(source['prior_scientific_commit'] == 'ab391f8f018c5aa3d2d6a7a6201fe93dbf75a9af'
            and source['prior_capsule_tree_sha'] == source['public_capsule_tree_sha'] == '62c232212e4cb6681c241930efa494da2b2d6f8f', 'SOURCE_TREE_PRESERVATION_DIFFERS')
    timestamp(source['verified_utc'])
    require(exact(source['local_zip'], {'bytes':NEW_ASSETS[0]['bytes'],'file':NEW_ASSETS[0]['filename'],
                                      'sha256':NEW_ASSETS[0]['sha256']}), 'SOURCE_ZIP_BINDING_DIFFERS')
    leaves = manifest['files']
    require(type(leaves) is list and len(leaves) == manifest['file_count_excluding_this_manifest'] == 129, 'SCIENCE_MANIFEST_COUNT_DIFFERS')
    pins = {row['path']: row for row in leaves}
    require(len(pins) == len(leaves) and 'MANIFEST.json' not in pins, 'SCIENCE_MANIFEST_DUPLICATE')
    require(sum(row['bytes'] for row in leaves) == manifest['bytes_excluding_this_manifest'] == 1647459, 'SCIENCE_MANIFEST_TOTAL_DIFFERS')
    pins['MANIFEST.json'] = {'path':'MANIFEST.json','bytes':len(manifest_raw),'sha256':SCIENCE_SHA}
    rows = source['files']
    require(type(rows) is list and len(rows) == 130, 'SOURCE_FILE_LIST_COUNT_DIFFERS')
    observed = {row['path']:row for row in rows}
    require(len(observed) == len(rows) and set(observed) == set(pins), 'SOURCE_FILE_MEMBERSHIP_DIFFERS')
    for name, pin in pins.items():
        require(safe_relative(name), 'UNSAFE_SCIENCE_PATH')
        row = observed[name]
        url = f'https://raw.githubusercontent.com/{REPOSITORY}/{COMMIT}/{CAPSULE}/{name}'
        require(row['url'] == row['final_url'] == url and type(row['http_status']) is int and row['http_status'] == 200, 'SOURCE_URL_OR_HTTP_DIFFERS_'+name)
        require(type(row['bytes']) is int and type(row['expected_bytes']) is int
                and row['bytes'] == row['expected_bytes'] == pin['bytes']
                and row['sha256'] == row['expected_sha256'] == pin['sha256']
                and row['exact_bytes_match'] is True, 'SOURCE_FILE_BYTES_DIFFERS_'+name)
    require(sum(row['bytes'] for row in rows) == source['complete_bytes'], 'SOURCE_TOTAL_DIFFERS')
    require(replay['status'] == 'PASS_REVIEWED_REGULAR_SHELL_SPECTRA_20261013'
            and replay['source_manifest_sha256'] == SCIENCE_SHA
            and replay['new_positive_checks_per_mode'] == 103 and replay['new_wrong_formula_controls_per_mode'] == 30,
            'SCIENCE_REPLAY_SUMMARY_DIFFERS')
    for mode in ('normal','optimized'):
        selected = [row for row in replay['checks'] if row['mode'] == mode]
        require(len(selected) == 5 and sum(row['positive_checks'] for row in selected) == 103
                and sum(row['wrong_formula_controls'] for row in selected) == 30, 'SCIENCE_REPLAY_MODE_COUNTS_DIFFER')
    return len(rows)

def verify_inventory(inventory):
    require(set(inventory) == {'schema_version','concept_id','record_id','derived_from_frozen_private_inventory_sha256',
            'derivation','file_count','total_bytes','inherited_file_count','inherited_bytes','files'}, 'PUBLIC_INVENTORY_FIELDS_DIFFER')
    require(inventory['schema_version'] == 1 and inventory['concept_id'] == CONCEPT
            and inventory['derived_from_frozen_private_inventory_sha256'] == INVENTORY_SHA, 'PUBLIC_INVENTORY_PROVENANCE_DIFFERS')
    rid = inventory['record_id']
    require(type(rid) is str and rid.isdigit() and rid not in GROUPS and rid != CONCEPT, 'NEW_RECORD_ID_INVALID')
    rows = inventory['files']
    require(type(rows) is list, 'INVENTORY_FILES_TYPE')
    pins = {r['filename']:r for r in rows}
    require(len(rows) == len(pins) == COUNT and type(inventory['file_count']) is int and inventory['file_count'] == COUNT, 'FILE_COUNT_OR_DUPLICATE_DIFFERS')
    for name, pin in pins.items():
        require(set(pin) == {'filename','bytes','md5','sha256'} and safe_relative(name) and '/' not in name,
                'PUBLIC_FILE_PIN_FIELDS_OR_NAME_DIFFERS')
        require(type(pin['bytes']) is int and pin['bytes'] > 0 and hexadecimal(pin['md5'],32)
                and hexadecimal(pin['sha256']), 'FILE_PIN_FORMAT_DIFFERS_'+name)
    require(type(inventory['total_bytes']) is int and inventory['total_bytes'] == sum(p['bytes'] for p in rows) == TOTAL,
            'TOTAL_BYTES_DIFFERS')
    require(type(inventory['inherited_file_count']) is int and inventory['inherited_file_count'] == INHERITED_COUNT
            and type(inventory['inherited_bytes']) is int and inventory['inherited_bytes'] == INHERITED_BYTES, 'INHERITED_COUNT_BYTES_DIFFERS')
    for pin in NEW_ASSETS:
        require(exact(pins.get(pin['filename']),pin), 'NEW_ASSET_PIN_DIFFERS')
    return rid, pins

def verify_writes(writes, rid, expected):
    require(writes['draft_id'] == rid and hexadecimal(writes['source_journal_sha256']), 'WRITE_JOURNAL_IDENTITY_DIFFERS')
    intents, responses, unknown = writes['write_intents'], writes['write_responses'], writes['unknown_outcomes']
    events = intents+responses+unknown
    indexes = [row['source_journal_index'] for row in events]
    require(all(type(i) is int and 0 <= i < writes['source_journal_event_count'] for i in indexes)
            and len(indexes) == len(set(indexes)), 'WRITE_EVENT_INDEX_DIFFERS')
    for rows, kind in ((intents,'WRITE_INTENT'),(responses,'WRITE_RESPONSE'),(unknown,'WRITE_OUTCOME_UNKNOWN')):
        require(all(row['kind'] == kind for row in rows), 'WRITE_EVENT_KIND_DIFFERS')
        for row in rows:
            timestamp(row['utc'])
    publish = [row for row in intents if row['operation'] == 'publish']
    publish_responses = [row for row in responses if row['operation'] == 'publish']
    require(len(publish) == writes['publish_intent_count'] == 1 and len(publish_responses) == writes['publish_response_count'] <= 1,
            'PUBLICATION_INTENT_COUNTS_DIFFER')
    require(writes['write_intent_count'] == len(intents) and writes['unknown_write_outcome_count'] == len(unknown)
            and writes['successful_write_response_count'] == sum(type(r['status']) is int and 200 <= r['status'] < 300 for r in responses),
            'WRITE_SUMMARY_COUNTS_DIFFER')
    p = publish[0]
    require(p['url'] == f'https://zenodo.org/api/records/{rid}/draft/actions/publish'
            and p['body_bytes'] == 0 and p['body_sha256'] == sha(b''), 'PUBLISH_INTENT_DIFFERS')
    metadata = [row for row in intents if row['operation'] == 'save_metadata']
    require(len(metadata) == 1, 'EXACT_METADATA_WRITE_REQUIRED')
    serialized = json.dumps(expected, ensure_ascii=False).encode()
    m = metadata[0]
    require(m['url'] == f'https://zenodo.org/api/records/{rid}/draft' and m['body_bytes'] == len(serialized)
            and m['body_sha256'] == sha(serialized) and m['source_journal_index'] < p['source_journal_index'], 'METADATA_WRITE_INTENT_DIFFERS')
    for row in publish_responses:
        require(row['source_journal_index'] > p['source_journal_index'], 'PUBLISH_RESPONSE_PRECEDES_INTENT')
    return p

def verify(root, manifest_sha):
    root = Path(root).resolve()
    raw = (root/'WITNESS_MANIFEST.json').read_bytes()
    require(hexadecimal(manifest_sha) and sha(raw) == manifest_sha, 'WITNESS_MANIFEST_PIN_DIFFERS')
    manifest = parse(raw)
    require(manifest['schema_version'] == 1 and type(manifest['files']) is dict, 'WITNESS_MANIFEST_SCHEMA_DIFFERS')
    expectedpaths = set(manifest['files'])
    require(all(safe_relative(name) for name in expectedpaths), 'UNSAFE_WITNESS_PATH')
    present = {str(x.relative_to(root)) for x in root.rglob('*') if x.is_file() or x.is_symlink()}
    require(present == expectedpaths | {'WITNESS_MANIFEST.json'}, 'WITNESS_FILE_MEMBERSHIP_DIFFERS')
    for name, pin in manifest['files'].items():
        path = root/name
        require(not path.is_symlink() and not any(p.is_symlink() for p in path.parents if p != root.parent)
                and path.resolve().is_relative_to(root), 'UNSAFE_WITNESS_PATH')
        data = path.read_bytes()
        require(type(pin['bytes']) is int and len(data) == pin['bytes'] and sha(data) == pin['sha256'], 'WITNESS_BYTES_DIFFER_'+name)
    expected_raw = (root/'FINAL_METADATA_MODERN.json').read_bytes()
    require(sha(expected_raw) == METADATA_SHA, 'FROZEN_METADATA_PIN_DIFFERS')
    expected = parse(expected_raw)
    require(expected['metadata']['title'] == TITLE and expected['metadata']['version'] == VERSION
            and expected['metadata']['publication_date'] == '2026-10-13', 'SCIENTIFIC_INTENT_DIFFERS')
    require(SOURCE_URL in [r['identifier'] for r in expected['metadata']['related_identifiers']], 'SCIENCE_SOURCE_REFERENCE_MISSING')
    record, listing = read(root/'NEW_PUBLIC_RECORD.json'), read(root/'NEW_PUBLIC_FILES.json')
    inventory, verified = read(root/'PUBLIC_INVENTORY.json'), read(root/'VERIFIED_PUBLIC.json')
    streams, writes = read(root/'PUBLIC_STREAM_OBSERVATIONS.json'), read(root/'WRITE_SUMMARY.json')
    rid, pins = verify_inventory(inventory)
    require(record['id'] == rid and record['parent']['id'] == CONCEPT and record['is_published'] is True
            and record['is_draft'] is False, 'RECORD_IDENTITY_STATE_DIFFERS')
    require(record['pids']['doi']['identifier'] == '10.5281/zenodo.'+rid
            and record['parent']['pids']['doi']['identifier'] == '10.5281/zenodo.'+CONCEPT, 'DOI_FAMILY_DIFFERS')
    require(record['versions']['is_latest'] is True and str(record['parent']['access']['owned_by']['user']) == '1386319', 'LATEST_OR_OWNER_DIFFERS')
    compare_metadata(record, expected)
    publication_intent = verify_writes(writes, rid, expected)
    preceding = streams['preceding_public_file_list_event']
    require(streams['record_id'] == rid and preceding['kind'] == 'GET'
            and preceding['url'] == f'https://zenodo.org/api/records/{rid}/files'
            and type(preceding['status']) is int and preceding['status'] == 200
            and preceding['authenticated_request'] is False, 'PUBLIC_STREAM_ROUND_SCOPE_DIFFERS')
    require(streams['source_journal_sha256'] == writes['source_journal_sha256'], 'PUBLIC_JOURNAL_BINDING_DIFFERS')
    index = preceding['source_journal_index']
    require(type(index) is int and index >= 0 and index+COUNT < writes['source_journal_event_count'], 'PUBLIC_LIST_EVENT_INDEX_INVALID')
    previous_time = timestamp(preceding['utc'])
    require(index > publication_intent['source_journal_index'] and previous_time >= timestamp(publication_intent['utc']),
            'PUBLIC_STREAM_ROUND_PRECEDES_PUBLISH_INTENT')
    write_indexes = {r['source_journal_index'] for group in ('write_intents','write_responses','unknown_outcomes') for r in writes[group]}
    require(write_indexes.isdisjoint(range(index,index+COUNT+1)), 'WRITE_AND_PUBLIC_READBACK_EVENT_INDEX_COLLISION')
    for offset, event in enumerate(streams['files'], 1):
        require(event['kind'] == 'FULL_CONTENT_STREAM_VERIFIED' and type(event['source_journal_index']) is int
                and event['source_journal_index'] == index+offset, 'PUBLIC_STREAM_EVENT_ORDER_DIFFERS')
        observed_time = timestamp(event['utc'])
        require(observed_time >= previous_time, 'PUBLIC_STREAM_TIME_ORDER_DIFFERS')
        previous_time = observed_time
    require(timestamp(verified['utc']) >= previous_time, 'PUBLIC_RECEIPT_PRECEDES_STREAMS')
    files = filemap(listing)
    observed = {r['filename']:r for r in streams['files']}
    vrows = {r['filename']:r for r in verified['content']}
    require(len(streams['files']) == len(observed) == len(verified['content']) == len(vrows) == COUNT, 'PUBLIC_STREAM_COUNT_DIFFERS')
    require(set(files) == set(pins) == set(observed) == set(vrows), 'FILE_MEMBERSHIP_DIFFERS')
    for name, pin in pins.items():
        e, obs, v = files[name], observed[name], vrows[name]
        require(type(e['size']) is int and e['size'] == pin['bytes'] and e['checksum'] == 'md5:'+pin['md5']
                and e['status'] == 'completed', 'FILE_PIN_DIFFERS_'+name)
        require(exact(obs['observed'], {k:pin[k] for k in ('bytes','md5','sha256')}), 'COMPLETE_STREAM_OBSERVATION_DIFFERS_'+name)
        require(v['sha256'] == pin['sha256'] and v['basis'] == 'FRESH_COMPLETE_CONTENT_STREAM'
                and e['file_id'] == v['remote_file_id'] and e['version_id'] == v['remote_version_id'], 'VERIFIED_IDENTITY_DIFFERS_'+name)
    require(verified['status'] == 'PASS_COMPLETE_35_FILE_EDITION_READBACK' and verified['published'] is True
            and verified['draft_id'] == rid and verified['parent_id'] == CONCEPT
            and type(verified['file_count']) is int and verified['file_count'] == COUNT
            and type(verified['total_bytes']) is int and verified['total_bytes'] == TOTAL, 'VERIFICATION_RECEIPT_IDENTITY_DIFFERS')
    require(verified['metadata_sha256'] == METADATA_SHA and verified['inventory_sha256'] == INVENTORY_SHA, 'APPROVED_INPUT_PIN_DIFFERS')
    preserved = read(root/'PRESERVED_RECORDS.json')
    require(preserved['status'] == 'PASS_ALL_SEVEN_PRIOR_PUBLIC_RECORDS_PRESERVED' and len(preserved['records']) == 7
            and {g['record_id']:g['concept_id'] for g in preserved['records']} == GROUPS, 'EXACT_SEVEN_PRESERVED_RECORDS_REQUIRED')
    timestamp(preserved['utc'])
    prior_record, prior_files = None, None
    for group in preserved['records']:
        for key in ('before_record','after_record','before_files','after_files'):
            relative = group[key]
            require(safe_relative(relative) and relative.startswith('preserved/') and relative in expectedpaths, 'UNSAFE_PRESERVATION_WITNESS_PATH')
        before, after = read(root/group['before_record']), read(root/group['after_record'])
        bf, af = read(root/group['before_files']), read(root/group['after_files'])
        require(all(exact(before[k], after[k]) for k in ('metadata','custom_fields','access')) and exact(filemap(bf),filemap(af)), 'PRESERVED_RECORD_CHANGED_'+group['record_id'])
        require(all(str(r['id']) == group['record_id'] and str(r['parent']['id']) == group['concept_id']
                    and r['is_published'] is True and r['is_draft'] is False for r in (before,after)), 'PRESERVED_IDENTITY_CHANGED')
        if group['record_id'] == PRIOR:
            prior_record, prior_files = before, filemap(bf)
    inherited = set(pins)-{r['filename'] for r in NEW_ASSETS}
    require(set(prior_files) == inherited and len(inherited) == INHERITED_COUNT
            and sum(pins[name]['bytes'] for name in inherited) == INHERITED_BYTES, 'INHERITED_FILE_MEMBERSHIP_DIFFERS')
    for name in inherited:
        require(type(prior_files[name]['size']) is int and prior_files[name]['size'] == pins[name]['bytes']
                and prior_files[name]['checksum'] == 'md5:'+pins[name]['md5'], 'INHERITED_FILE_PIN_DIFFERS_'+name)
    prior_inventory_raw = (root/'PRIOR_PUBLIC_INVENTORY.json').read_bytes()
    require(sha(prior_inventory_raw) == PRIOR_INVENTORY_SHA, 'PRIOR_PUBLIC_INVENTORY_PIN_DIFFERS')
    prior_inventory = parse(prior_inventory_raw)
    require(prior_inventory['record_id'] == PRIOR and prior_inventory['concept_id'] == CONCEPT
            and prior_inventory['file_count'] == INHERITED_COUNT and prior_inventory['total_bytes'] == INHERITED_BYTES,
            'PRIOR_PUBLIC_INVENTORY_IDENTITY_DIFFERS')
    prior_pins = {r['filename']:r for r in prior_inventory['files']}
    require(len(prior_inventory['files']) == len(prior_pins) == INHERITED_COUNT and set(prior_pins) == inherited
            and all(exact(prior_pins[name],pins[name]) for name in inherited), 'INHERITED_PUBLIC_SHA256_PINS_DIFFER')
    em, pm = expected['metadata'], prior_record['metadata']
    require(em['description'].endswith('\n<h3>Preserved Description from published record 23294979</h3>\n'+pm['description']), 'PRIOR_DESCRIPTION_HISTORY_SUFFIX_DIFFERS')
    marker = '\n<h3>Preserved Notes from published record 23294979 (2026.10.11-scalar-response-monotonicity)</h3>\n'
    require(em['additional_descriptions'][0]['description'].endswith(marker+pm['additional_descriptions'][0]['description'])
            and exact(em['additional_descriptions'][1:],pm['additional_descriptions'][1:]), 'PRIOR_NOTES_HISTORY_SUFFIX_DIFFERS')
    source_count = verify_source(read(root/'IMMUTABLE_PUBLIC_SOURCE_READBACK.json'),
                                 (root/'SCIENTIFIC_MANIFEST.json').read_bytes(), (root/'SCIENTIFIC_REPLAY_RECEIPT.json').read_bytes())
    return {'status':'PASS_OFFLINE_SAVED_35_FILE_PUBLICATION_WITNESSES','record_id':rid,'concept_id':CONCEPT,
            'file_count':COUNT,'total_bytes':TOTAL,'inherited_files':INHERITED_COUNT,'inherited_bytes':INHERITED_BYTES,
            'witness_files':len(expectedpaths),'preserved_records':7,'immutable_source_files':source_count,
            'immutable_source_bytes':1670878,'scientific_source_commit':COMMIT,'manifest_sha256':manifest_sha,'network_calls':0,
            'scope':'Authenticates pinned saved client observations and scientific intent; does not re-download attachments, rerun the mathematics, establish current server state, or prove server-wide atomicity.'}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',required=True)
    p.add_argument('--manifest-sha256',required=True)
    p.add_argument('--output')
    args = p.parse_args()
    result = verify(args.root,args.manifest_sha256)
    encoded = json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output:
        output = Path(args.output).resolve()
        require(not output.is_relative_to(Path(args.root).resolve()), 'OUTPUT_MUST_BE_OUTSIDE_WITNESS_ROOT')
        output.write_text(encoded)
    print(encoded,end='')

if __name__ == '__main__':
    main()
