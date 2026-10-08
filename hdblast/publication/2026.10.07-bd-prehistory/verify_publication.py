#!/usr/bin/env python3
"""Verify externally pinned saved publication witnesses; no network or source calls.

This verifies captured evidence. It does not establish a fresh remote state or
turn self-produced hash pins into authenticated public readback.
"""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import html.parser
import json
import re
import sys
from urllib.parse import quote

sys.dont_write_bytecode = True
MAX_INPUT = 20 * 1024 * 1024
PARENT = '17088132'
PRIOR = '23114217'
PRIOR_STREAM_PIN = 'f3cb8d3c5280fe90307c98988bb2bf97de252bc2bebfbc20aef296612deb96ff'
ROLES = {'inventory', 'metadata', 'new_public_record', 'new_public_files',
         'latest_public_record', 'verified_public', 'controller_journal', 'prior_full_stream'}
PRESERVED = {PRIOR: PARENT, '22347452': PARENT, '23111008': '22922927'}


def require(value, message):
    if not value:
        raise ValueError(message)


def load(raw):
    def unique(items):
        result = {}
        for key, value in items:
            require(key not in result, 'Duplicate JSON key')
            result[key] = value
        return result
    def invalid(value):
        raise ValueError('Nonfinite JSON value: ' + value)
    return json.loads(raw, object_pairs_hook=unique, parse_constant=invalid)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'),
                      ensure_ascii=False, allow_nan=False).encode()


def safe_file(root, name):
    require(type(name) is str and name and '\\' not in name, 'Relative witness filename required')
    parts = PurePosixPath(name)
    require(not parts.is_absolute() and parts.as_posix() == name
            and all(p not in ('', '.', '..') for p in name.split('/')), 'Unsafe witness path')
    path = root / name
    require(path.is_file() and all(not p.is_symlink() for p in (path, *path.parents)),
            'Real witness file required: ' + name)
    require(path.stat().st_size <= MAX_INPUT, 'Witness too large')
    return path


class HtmlEvents(html.parser.HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.events = []
        self.depth = 0
    def handle_starttag(self, tag, attrs):
        self.events.append(('start', tag, tuple(sorted(attrs))))
        self.depth += 1
    def handle_endtag(self, tag):
        self.events.append(('end', tag))
        self.depth = max(0, self.depth - 1)
    def handle_startendtag(self, tag, attrs):
        self.events.append(('empty', tag, tuple(sorted(attrs))))
    def handle_data(self, data):
        if self.depth == 0 and not data.strip():
            return
        if self.events and self.events[-1][0] == 'text':
            self.events[-1] = ('text', self.events[-1][1] + data)
        else:
            self.events.append(('text', data))
    def handle_comment(self, data):
        self.events.append(('comment', data))


def vocabulary_path(path):
    words = tuple(p for p in path if not isinstance(p, int))
    return words in {
        ('metadata', 'resource_type'), ('metadata', 'rights'), ('metadata', 'languages'),
        ('metadata', 'creators', 'role'), ('metadata', 'contributors', 'role'),
        ('metadata', 'related_identifiers', 'relation_type'),
        ('metadata', 'related_identifiers', 'resource_type'),
        ('metadata', 'additional_descriptions', 'type'),
        ('metadata', 'additional_descriptions', 'lang'),
        ('custom_fields', 'code:programmingLanguage'), ('custom_fields', 'code:developmentStatus'),
    }


def normalized(actual, expected, path=()):
    if type(actual) is dict:
        drop = ({'title', 'description', 'icon', 'props'}
                if type(expected) is dict and 'id' in expected and 'id' in actual
                and vocabulary_path(path) else set())
        return {k: normalized(v, expected.get(k) if type(expected) is dict else None, path + (k,))
                for k, v in actual.items() if k not in drop}
    if type(actual) is list:
        return [normalized(v, expected[i] if type(expected) is list and i < len(expected) else None,
                           path + (i,)) for i, v in enumerate(actual)]
    if type(actual) is str and type(expected) is str and path and path[-1] == 'description':
        parser = HtmlEvents()
        parser.feed(actual)
        parser.close()
        return parser.events
    return actual


def same(actual, expected, path):
    return canonical(normalized(actual, expected, path)) == canonical(normalized(expected, expected, path))


def metadata_equal(actual, expected):
    require(type(expected) is dict and set(expected) == {'metadata', 'custom_fields', 'access'},
            'Complete frozen editable metadata body required')
    for key in ('metadata', 'custom_fields'):
        require(same(actual.get(key, {}), expected[key], (key,)), 'Complete ' + key + ' differs')
    access = actual.get('access', {})
    require(type(access) is dict and set(access) <= set(expected['access']) | {'status'},
            'Unexpected access fields')
    require(same({k: access.get(k) for k in expected['access']}, expected['access'], ('access',)),
            'Complete access differs')


def identity(record, expected_id, parent):
    require(type(record) is dict, 'Saved record object required')
    value = record.get('id')
    require(type(value) in (str, int) and str(value) == expected_id, 'Record identity differs')
    if 'is_published' in record:
        require(record['is_published'] is True and record.get('is_draft') is False
                and str(record.get('parent', {}).get('id')) == parent, 'Published family/state differs')
    else:
        require(record.get('submitted') is True and record.get('state') == 'done'
                and str(record.get('conceptrecid')) == parent, 'Legacy published family/state differs')
    doi = record.get('pids', {}).get('doi', {}).get('identifier', record.get('metadata', {}).get('doi'))
    require(doi == '10.5281/zenodo.' + expected_id, 'Published version DOI differs')


def file_map(record):
    files = record.get('files')
    entries = files.get('entries', {}) if type(files) is dict else record.get('entries', files)
    if type(entries) is dict:
        entries = list(entries.values())
    require(type(entries) is list, 'Complete saved file list required')
    result = {}
    for item in entries:
        require(type(item) is dict, 'File entry object required')
        name = item.get('key', item.get('filename'))
        require(type(name) is str and name and name not in result, 'Duplicate/missing filename')
        result[name] = item
    return result


def pin_rows(rows, count):
    require(type(rows) is list and len(rows) == count, 'Exact inventory count required')
    result = {}
    for row in rows:
        require(type(row) is dict, 'Inventory object required')
        name = row.get('filename')
        require(type(name) is str and name and name not in result and '/' not in name and '\\' not in name,
                'Unique inventory leaf filename required')
        require(type(row.get('bytes')) is int and row['bytes'] > 0
                and type(row.get('md5')) is str and re.fullmatch('[0-9a-f]{32}', row['md5'])
                and type(row.get('sha256')) is str and re.fullmatch('[0-9a-f]{64}', row['sha256']),
                'Canonical complete checksum pin required')
        result[name] = row
    return result


def assert_pins(files, pins):
    require(set(files) == set(pins), 'Exact file membership differs')
    for name, row in pins.items():
        remote = files[name]
        size = remote.get('size', remote.get('filesize'))
        require(type(size) is int and size == row['bytes'], 'File size differs: ' + name)
        require(remote.get('checksum') in (row['md5'], 'md5:' + row['md5']), 'File MD5 differs: ' + name)
        if 'status' in remote:
            require(remote['status'] == 'completed', 'File not completed: ' + name)


def preserve(actual, baseline, current_files):
    require(type(baseline.get('metadata')) is dict and baseline['metadata'],
            'Full preserved baseline metadata required')
    if all(key in baseline for key in ('metadata', 'custom_fields', 'access')):
        metadata_equal(actual, {k: baseline[k] for k in ('metadata', 'custom_fields', 'access')})
    else:
        require(same(actual.get('metadata'), baseline.get('metadata'), ('metadata',)),
                'Preserved legacy metadata differs')
    old = file_map(baseline)
    embedded = file_map(actual)
    require(set(current_files) == set(old), 'Preserved record file membership differs')
    require(set(embedded) == set(old), 'Preserved embedded file membership differs')
    for name, f in old.items():
        pin = {'bytes': f.get('size', f.get('filesize')), 'md5': f.get('checksum', '').removeprefix('md5:')}
        assert_pins({name: current_files[name]}, {name: pin})
        assert_pins({name: embedded[name]}, {name: pin})
        # Where both captures expose immutable identity, preserve it too.
        for key in ('file_id', 'version_id'):
            if key in f:
                require(current_files[name].get(key) == f[key], 'Preserved immutable identity differs')
                require(embedded[name].get(key) == f[key], 'Preserved embedded immutable identity differs')
            if key in embedded[name]:
                require(embedded[name][key] == current_files[name].get(key),
                        'Preserved current file captures disagree')


def verify_values(inventory, metadata, public_record, public_files, latest, verified,
                  prior_stream, journal, preserved, expected_id, inventory_sha, metadata_sha):
    require(re.fullmatch('[1-9][0-9]*', expected_id or '') and expected_id not in PRESERVED,
            'Explicit distinct expected new record id required')
    require(inventory.get('sealed_upload_manifest') is True
            and inventory.get('concept_doi') == '10.5281/zenodo.' + PARENT, 'Frozen edition scope differs')
    old = pin_rows(inventory.get('inherited_files'), 23)
    new = pin_rows(inventory.get('new_files'), 4)
    require(not (set(old) & set(new)), 'Inherited/new filename overlap')
    pins = {**old, **new}
    require(sum(p['bytes'] for p in old.values()) == 440222994
            and inventory.get('expected_total_files_if_four_additions_finalized') == 27
            and type(inventory.get('total_bytes_final')) is int
            and inventory['total_bytes_final'] == sum(p['bytes'] for p in pins.values()),
            'Frozen count/byte totals differ')
    identity(public_record, expected_id, PARENT)
    identity(latest, expected_id, PARENT)
    for record in (public_record, latest):
        metadata_equal(record, metadata)
        assert_pins(file_map(record), pins)
    files = file_map(public_files)
    assert_pins(files, pins)
    for captured in (file_map(public_record), file_map(latest)):
        for name, item in captured.items():
            for key in ('file_id', 'version_id'):
                if key in item:
                    require(item[key] == files[name].get(key), 'New public file identities disagree')
    prefix = 'https://zenodo.org/api/records/' + expected_id
    require(public_record.get('links', {}).get('self') == prefix
            and public_record.get('links', {}).get('files') == prefix + '/files', 'Direct public URL binding differs')
    require(verified.get('status') == 'PASS_COMPLETE_27_FILE_EDITION_READBACK'
            and verified.get('published') is True and str(verified.get('draft_id')) == expected_id
            and str(verified.get('parent_id')) == PARENT
            and verified.get('inventory_sha256') == inventory_sha
            and verified.get('metadata_sha256') == metadata_sha
            and type(verified.get('file_count')) is int and verified['file_count'] == 27
            and type(verified.get('total_bytes')) is int and verified['total_bytes'] == inventory['total_bytes_final'],
            'Published verification receipt differs')
    rows = prior_stream.get('rows')
    require(type(rows) is list and len(rows) == 23, 'Complete prior23 content receipt required')
    prior = {}
    for row in rows:
        name = row.get('filename')
        require(name in old and name not in prior, 'Prior full-stream membership differs')
        expected = {k: old[name][k] for k in ('bytes', 'md5', 'sha256')}
        require(row.get('status') == 'PASS_COMPLETE_REMOTE_BYTES_MD5_SHA256'
                and canonical(row.get('observed')) == canonical(expected)
                and row.get('remote_file_id') and row.get('remote_version_id'), 'Prior full stream differs')
        prior[name] = row
    content = verified.get('content')
    require(type(content) is list and len(content) == 27, 'Complete27 SHA256 basis rows required')
    bases = {}
    markers = [i for i, item in enumerate(journal)
               if item.get('kind') == 'GET' and item.get('url') == prefix + '/files'
               and item.get('status') == 200 and item.get('authenticated_request') is False]
    require(markers, 'Terminal direct public file-list GET absent from journal')
    terminal = journal[markers[-1] + 1:]
    streams = {}
    for event in terminal:
        if event.get('kind') == 'FULL_CONTENT_STREAM_VERIFIED':
            name = event.get('filename')
            require(name in pins and name not in streams, 'Terminal content stream membership differs')
            observed = {k: pins[name][k] for k in ('bytes', 'md5', 'sha256')}
            require(canonical(event.get('observed')) == canonical(observed), 'Terminal full stream differs')
            streams[name] = event
    for row in content:
        name = row.get('filename')
        require(name in pins and name not in bases and row.get('sha256') == pins[name]['sha256'],
                'Published SHA256 basis identity differs')
        basis = row.get('basis')
        if basis == 'FRESH_COMPLETE_CONTENT_STREAM':
            require(name in streams, 'Fresh public stream evidence absent: ' + name)
            require(files[name].get('links', {}).get('content') == prefix + '/files/' + quote(name, safe='') + '/content',
                    'Stream content URL binding differs')
        elif basis == 'SAME_IMMUTABLE_FILE_AND_CONTENT_VERSION_AS_PRIOR_FULL_STREAM':
            require(name in old and files[name].get('file_id', files[name].get('id')) == prior[name]['remote_file_id']
                    and files[name].get('version_id') == prior[name]['remote_version_id'],
                    'Inherited immutable SHA256 reuse differs')
        else:
            raise ValueError('Unrecognized SHA256 evidence basis')
        require(name not in new or basis == 'FRESH_COMPLETE_CONTENT_STREAM', 'Every addition needs fresh stream')
        bases[name] = basis
    require(set(streams) == {name for name, basis in bases.items() if basis == 'FRESH_COMPLETE_CONTENT_STREAM'},
            'Terminal streams and published basis disagree')
    require(set(preserved) == set(PRESERVED), 'All three preserved records required')
    for record_id, group in preserved.items():
        identity(group['record'], record_id, PRESERVED[record_id])
        identity(group['baseline'], record_id, PRESERVED[record_id])
        preserve(group['record'], group['baseline'], file_map(group['files']))
    assert_pins(file_map(preserved[PRIOR]['files']), old)
    return {'status': 'PASS_OFFLINE_27_FILE_PUBLICATION_WITNESSES', 'record_id': expected_id,
            'version_doi': '10.5281/zenodo.' + expected_id, 'concept_id': PARENT,
            'file_count': 27, 'total_bytes': inventory['total_bytes_final'],
            'fresh_complete_public_streams': sum(x == 'FRESH_COMPLETE_CONTENT_STREAM' for x in bases.values()),
            'prior_same_immutable_stream_reuses': sum(x != 'FRESH_COMPLETE_CONTENT_STREAM' for x in bases.values()),
            'preserved_records': sorted(PRESERVED),
            'full_metadata_custom_fields_access': 'MATCH_SEPARATELY_PINNED_EDITION_CONTRACT',
            'historical_original_literal_Notes_guard': 'FAIL_PRESERVED',
            'metadata_normalization': 'Only whitelisted RDM vocabulary decorations and semantic HTML character references; complete editable values and types retained',
            'trust_boundary': 'Externally supplied witness pins authenticate saved evidence; no live remote authentication occurs in this verifier',
            'network_calls': 0, 'numeric_imports': 0, 'source_callbacks': 0, 'array_decodes': 0,
            'original_binary80_state_error': 'NOT_ENCLOSED', 'full_twelve_case_certificate': 'UNRESOLVED',
            'metric_calibration': 'FAIL', 'higher_dimensional_Big_Bang_origin': 'NOT_ESTABLISHED',
            'external_novelty': 'NOT_ASSESSED'}


def verify_directory(root, manifest_path, manifest_sha, expected_id, inventory_sha, metadata_sha):
    root = Path(root).absolute()
    require(root.is_dir() and root.resolve() == root and all(not p.is_symlink() for p in (root, *root.parents)),
            'Real saved-witness root required')
    manifest_path = Path(manifest_path).absolute()
    require(manifest_path.is_file() and manifest_path.stat().st_size <= MAX_INPUT
            and all(not p.is_symlink() for p in (manifest_path, *manifest_path.parents)), 'Real witness manifest required')
    raw = manifest_path.read_bytes()
    require(re.fullmatch('[0-9a-f]{64}', manifest_sha or '') and sha(raw) == manifest_sha,
            'External witness manifest hash differs')
    manifest = load(raw)
    require(type(manifest) is dict and set(manifest) == {'schema_version', 'files', 'witnesses', 'preserved_records',
            'historical_original_literal_Notes_guard'} and type(manifest['schema_version']) is int
            and manifest['schema_version'] == 1 and manifest['historical_original_literal_Notes_guard'] == 'FAIL_PRESERVED',
            'Exact witness manifest schema/history required')
    values = {}
    require(type(manifest['files']) is dict and manifest['files'], 'Witness file pins required')
    for name, pin in manifest['files'].items():
        require(type(pin) is dict and set(pin) == {'bytes', 'sha256'} and type(pin['bytes']) is int
                and 0 <= pin['bytes'] <= MAX_INPUT and re.fullmatch('[0-9a-f]{64}', pin['sha256']), 'Witness pin differs')
        path = safe_file(root, name)
        content = path.read_bytes()
        require(len(content) == pin['bytes'] and sha(content) == pin['sha256'], 'Saved witness bytes differ: ' + name)
        values[name] = content
    roles = manifest['witnesses']
    require(type(roles) is dict and set(roles) == ROLES and all(type(v) is str and v in values for v in roles.values()),
            'Complete witness roles required')
    require(re.fullmatch('[0-9a-f]{64}', inventory_sha or '') and sha(values[roles['inventory']]) == inventory_sha
            and re.fullmatch('[0-9a-f]{64}', metadata_sha or '') and sha(values[roles['metadata']]) == metadata_sha,
            'Independent edition input hashes differ')
    require(sha(values[roles['prior_full_stream']]) == PRIOR_STREAM_PIN, 'Preserved prior full-stream receipt pin differs')
    groups = manifest['preserved_records']
    require(type(groups) is dict and set(groups) == set(PRESERVED), 'Preserved record roles required')
    preserved = {}
    for record_id, group in groups.items():
        require(type(group) is dict and set(group) == {'record', 'files', 'baseline'}
                and all(type(v) is str and v in values for v in group.values()), 'Preserved witnesses incomplete')
        preserved[record_id] = {k: load(values[v]) for k, v in group.items()}
    data = {k: load(values[v]) for k, v in roles.items() if k != 'controller_journal'}
    journal = [load(line) for line in values[roles['controller_journal']].splitlines()]
    require(all(type(item) is dict for item in journal), 'Controller journal object rows required')
    result = verify_values(data['inventory'], data['metadata'], data['new_public_record'], data['new_public_files'],
        data['latest_public_record'], data['verified_public'], data['prior_full_stream'], journal, preserved,
        expected_id, inventory_sha, metadata_sha)
    result.update(witness_manifest_sha256=manifest_sha, inventory_sha256=inventory_sha,
                  metadata_sha256=metadata_sha, verifier_sha256=sha(Path(__file__).read_bytes()))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True)
    parser.add_argument('--witness-manifest', required=True)
    parser.add_argument('--witness-manifest-sha256', required=True)
    parser.add_argument('--expected-record-id', required=True)
    parser.add_argument('--inventory-sha256', required=True)
    parser.add_argument('--metadata-sha256', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    result = verify_directory(args.root, args.witness_manifest, args.witness_manifest_sha256,
                              args.expected_record_id, args.inventory_sha256, args.metadata_sha256)
    output = Path(args.output).absolute()
    require(not output.exists() and all(not p.is_symlink() for p in (output, *output.parents)), 'Fresh output file required')
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('x') as stream:
        stream.write(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: result[k] for k in ('status', 'record_id', 'file_count', 'total_bytes', 'network_calls')}))


if __name__ == '__main__':
    main()
