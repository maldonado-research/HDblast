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
PRIOR = '23225288'
PRIOR_STREAM_PIN = '2664ec6a867d806e6146847fcf677d2e5c20061a400262acb8563d3ca552529c'
SCIENTIFIC_REPLAY_STATUS = 'PASS_COMPLETE_STORED_STATE_COMPARISON_REPLAY'
ROLES = {'inventory', 'metadata', 'new_public_record', 'new_public_files',
         'latest_public_record', 'verified_public', 'controller_journal', 'prior_full_stream'}
PRESERVED = {PRIOR: PARENT, '23114217': PARENT, '22347452': PARENT, '23111008': '22922927'}
BASELINE_PINS = {
    '23225288': {
        'baseline': '261d81302396c1683f5c0e1185b071189cf83bbb5fde5fb4a098ca64e8fe47ca',
        'baseline_files': 'b9355b33cb3ba6e88dd15db52c53b85f51d68a122faccbacc5153880fc8b8276'},
    '23114217': {
        'baseline': '2e090396eb95ecd2f3a5bd5c1fb8305923827d4206600458f72f282b2b873623',
        'baseline_files': 'd7d90ef99d524647875ad3403b5eb608cd09e4088d819741f6063df637678dbb'},
    '22347452': {
        'baseline': 'd3b310f3c5138c40720fdcc0ab6ce97da99c5e6395fb33fe98f09910053c8289',
        'baseline_files': 'ded3b02e0c3523d01ade4e266eac7405bdaa02a24f94e666280b6b9460f95d8d'},
    '23111008': {
        'baseline': '95b6ecf9b3f57af3857365712167ae2562aa9e5dd56c893b2e5cdacd2ecc5005',
        'baseline_files': '23482f413f563a4aca79f2df3f17ac062a6d9dc02223d290ba1caedde11fd3a7'},
}


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
    result = json.loads(raw, object_pairs_hook=unique, parse_constant=invalid)
    # The JSON decoder also maps overflowing exponents to infinity.
    json.dumps(result, allow_nan=False)
    return result


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'),
                      ensure_ascii=False, allow_nan=False).encode()


def safe_file(root, name):
    require(type(name) is str and name and '\\' not in name and '\0' not in name,
            'Relative witness filename required')
    parts = PurePosixPath(name)
    require(not parts.is_absolute() and parts.as_posix() == name
            and all(p not in ('', '.', '..') for p in name.split('/')), 'Unsafe witness path')
    path = root / name
    require(path.is_file() and all(not p.is_symlink() for p in (path, *path.parents)),
            'Real witness file required: ' + name)
    require(path.stat().st_size <= MAX_INPUT, 'Witness too large')
    return path


class HtmlEvents(html.parser.HTMLParser):
    BLOCK_TAGS = {'address', 'article', 'aside', 'blockquote', 'br', 'div', 'dl', 'dt', 'dd',
                  'fieldset', 'figcaption', 'figure', 'footer', 'form', 'h1', 'h2', 'h3',
                  'h4', 'h5', 'h6', 'header', 'hr', 'li', 'main', 'nav', 'ol', 'p',
                  'pre', 'section', 'table', 'tbody', 'thead', 'tfoot', 'tr', 'td', 'th', 'ul'}
    VOID_TAGS = {'area', 'base', 'br', 'col', 'embed', 'hr', 'img', 'input', 'link',
                 'meta', 'param', 'source', 'track', 'wbr'}
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.events = []
        self.depth = 0
    def handle_starttag(self, tag, attrs):
        self.events.append(('start', tag, tuple(sorted(attrs))))
        if tag not in self.VOID_TAGS:
            self.depth += 1
    def handle_endtag(self, tag):
        self.events.append(('end', tag))
        self.depth = max(0, self.depth - 1)
    def handle_startendtag(self, tag, attrs):
        self.events.append(('empty', tag, tuple(sorted(attrs))))
    def handle_data(self, data):
        if self.depth == 0 and data and all(character in '\t\n\f\r ' for character in data):
            self.events.append(('top_whitespace', data))
            return
        if self.events and self.events[-1][0] == 'text':
            self.events[-1] = ('text', self.events[-1][1] + data)
        else:
            self.events.append(('text', data))
    def handle_comment(self, data):
        self.events.append(('comment', data))
    def handle_decl(self, decl):
        self.events.append(('declaration', decl))
    def handle_pi(self, data):
        self.events.append(('processing_instruction', data))
    def unknown_decl(self, data):
        self.events.append(('unknown_declaration', data))
    def semantic_events(self):
        result = []
        for index, event in enumerate(self.events):
            if event[0] == 'top_whitespace':
                previous = self.events[index - 1] if index else None
                following = self.events[index + 1] if index + 1 < len(self.events) else None
                def block(value):
                    return value is not None and value[0] in ('start', 'end', 'empty') and value[1] in self.BLOCK_TAGS
                if previous is None or following is None or block(previous) or block(following):
                    continue
                event = ('text', event[1])
            if result and result[-1][0] == event[0] == 'text':
                result[-1] = ('text', result[-1][1] + event[1])
            else:
                result.append(event)
        return result


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
    require(type(actual) is type(expected), 'Complete editable value type differs: ' + repr(path))
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
        return parser.semantic_events()
    return actual


def same(actual, expected, path):
    return canonical(normalized(actual, expected, path)) == canonical(normalized(expected, expected, path))


def metadata_equal(actual, expected):
    require(type(expected) is dict and set(expected) == {'metadata', 'custom_fields', 'access'},
            'Complete frozen editable metadata body required')
    require(type(actual) is dict and all(key in actual for key in ('metadata', 'custom_fields', 'access')),
            'Complete editable metadata/custom_fields/access presence required')
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
        for name, item in entries.items():
            require(type(item) is dict and item.get('key', item.get('filename')) == name,
                    'Dictionary file key and nested filename disagree')
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


def immutable_file_id(item):
    value = item.get('file_id', item.get('id'))
    require(type(value) is str and value and len(value) <= 256,
            'Immutable file identity required')
    if 'file_id' in item and 'id' in item:
        require(item['file_id'] == item['id'], 'Conflicting immutable file identities')
    return value


def immutable_version_id(item):
    value = item.get('version_id')
    require(type(value) is str and value and len(value) <= 256,
            'Immutable content-version identity required')
    return value


def preserved_embedded_identity(item, baseline_item):
    for key in ('file_id', 'id'):
        if key in baseline_item:
            require(immutable_file_id(item) == immutable_file_id(baseline_item),
                    'Preserved embedded immutable file identity differs')
    if 'version_id' in baseline_item:
        require(immutable_version_id(item) == immutable_version_id(baseline_item),
                'Preserved embedded immutable version identity differs')


def preserve(actual, baseline, current_files, baseline_files):
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
    require(set(baseline_files) == set(old), 'Preserved dated file-list membership differs')
    require(set(embedded) == set(old), 'Preserved embedded file membership differs')
    for name, f in old.items():
        pin = {'bytes': f.get('size', f.get('filesize')), 'md5': f.get('checksum', '').removeprefix('md5:')}
        assert_pins({name: current_files[name]}, {name: pin})
        assert_pins({name: baseline_files[name]}, {name: pin})
        assert_pins({name: embedded[name]}, {name: pin})
        # The dated /files endpoint provides both identities even when public
        # record embeddings expose only id. Require the complete dated pair.
        require(immutable_file_id(current_files[name]) == immutable_file_id(baseline_files[name])
                and immutable_version_id(current_files[name]) == immutable_version_id(baseline_files[name]),
                'Preserved dated file/version identity differs')
        preserved_embedded_identity(embedded[name], f)
        require(immutable_file_id(embedded[name]) == immutable_file_id(current_files[name]),
                'Preserved current embedded file identity differs')
        if 'version_id' in embedded[name]:
            require(immutable_version_id(embedded[name]) == immutable_version_id(current_files[name]),
                    'Preserved current embedded version identity differs')


def verify_values(inventory, metadata, public_record, public_files, latest, verified,
                  prior_stream, journal, preserved, expected_id, inventory_sha, metadata_sha):
    require(re.fullmatch('[1-9][0-9]*', expected_id or '') and expected_id not in PRESERVED,
            'Explicit distinct expected new record id required')
    require(inventory.get('sealed_upload_manifest') is True
            and inventory.get('concept_doi') == '10.5281/zenodo.' + PARENT, 'Frozen edition scope differs')
    require(inventory.get('scientific_replay_status') == SCIENTIFIC_REPLAY_STATUS,
            'Completed stored-state scientific replay contract required')
    old = pin_rows(inventory.get('inherited_files'), 27)
    new = pin_rows(inventory.get('new_files'), 2)
    require(not (set(old) & set(new)), 'Inherited/new filename overlap')
    pins = {**old, **new}
    require(sum(p['bytes'] for p in old.values()) == 448387915
            and type(inventory.get('expected_total_files_final')) is int
            and inventory['expected_total_files_final'] == 29
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
    public_identities = {}
    for name, item in files.items():
        public_identities[name] = (immutable_file_id(item), immutable_version_id(item))
    for captured in (file_map(public_record), file_map(latest)):
        for name, item in captured.items():
            require(immutable_file_id(item) == public_identities[name][0],
                    'New public file identities disagree')
            if 'version_id' in item:
                require(immutable_version_id(item) == public_identities[name][1],
                        'New public version identities disagree')
    prefix = 'https://zenodo.org/api/records/' + expected_id
    require(public_record.get('links', {}).get('self') == prefix
            and public_record.get('links', {}).get('files') == prefix + '/files', 'Direct public URL binding differs')
    require(verified.get('status') == 'PASS_COMPLETE_29_FILE_EDITION_READBACK'
            and verified.get('published') is True and str(verified.get('draft_id')) == expected_id
            and str(verified.get('parent_id')) == PARENT
            and verified.get('inventory_sha256') == inventory_sha
            and verified.get('metadata_sha256') == metadata_sha
            and type(verified.get('file_count')) is int and verified['file_count'] == 29
            and type(verified.get('total_bytes')) is int and verified['total_bytes'] == inventory['total_bytes_final'],
            'Published verification receipt differs')
    rows = prior_stream.get('rows')
    require(prior_stream.get('status') == 'PRIOR_FULL_STREAM_ATTESTATION_WITH_FRESH_GET_IDENTITY',
            'Pinned dated prior27 attestation required')
    require(type(rows) is list and len(rows) == 27, 'Complete prior27 content receipt required')
    prior = {}
    for row in rows:
        name = row.get('filename')
        require(name in old and name not in prior, 'Prior full-stream membership differs')
        expected = {k: old[name][k] for k in ('bytes', 'md5', 'sha256')}
        require(row.get('status') == 'PASS_PRIOR_COMPLETE_BYTES_WITH_FRESH_IMMUTABLE_IDENTITY'
                and canonical(row.get('observed')) == canonical(expected)
                and type(row.get('remote_file_id')) is str and row['remote_file_id']
                and type(row.get('remote_version_id')) is str and row['remote_version_id'],
                'Prior dated content attestation differs')
        prior[name] = row
    content = verified.get('content')
    require(type(content) is list and len(content) == 29, 'Complete29 SHA256 basis rows required')
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
        require(basis == 'FRESH_COMPLETE_CONTENT_STREAM', 'Every one of all29 assets needs a fresh complete stream')
        require(name in streams, 'Fresh public stream evidence absent: ' + name)
        require(files[name].get('links', {}).get('content') == prefix + '/files/' + quote(name, safe='') + '/content',
                'Stream content URL binding differs')
        require((row.get('remote_file_id'), row.get('remote_version_id')) == public_identities[name],
                'Verified fresh stream immutable identities differ: ' + name)
        bases[name] = basis
    require(set(streams) == set(pins) == set(bases),
            'Terminal streams and published basis disagree')
    require(set(preserved) == set(PRESERVED), 'All four preserved records required')
    for record_id, group in preserved.items():
        identity(group['record'], record_id, PRESERVED[record_id])
        identity(group['baseline'], record_id, PRESERVED[record_id])
        preserve(group['record'], group['baseline'], file_map(group['files']), file_map(group['baseline_files']))
    assert_pins(file_map(preserved[PRIOR]['files']), old)
    dated_prior_files = file_map(preserved[PRIOR]['baseline_files'])
    for name, row in prior.items():
        require(immutable_file_id(dated_prior_files[name]) == row['remote_file_id']
                and immutable_version_id(dated_prior_files[name]) == row['remote_version_id'],
                'Dated prior27 receipt and preserved dated identities disagree')
    return {'status': 'PASS_OFFLINE_29_FILE_PUBLICATION_WITNESSES', 'record_id': expected_id,
            'version_doi': '10.5281/zenodo.' + expected_id, 'concept_id': PARENT,
            'file_count': 29, 'total_bytes': inventory['total_bytes_final'],
            'fresh_complete_public_streams': 29, 'prior_same_immutable_stream_reuses': 0,
            'all29_content_policy': 'MUST_FRESH_COMPLETE_CONTENT_STREAM',
            'preserved_records': sorted(PRESERVED),
            'full_metadata_custom_fields_access': 'MATCH_SEPARATELY_PINNED_EDITION_CONTRACT',
            'historical_original_literal_Notes_guard': 'FAIL_PRESERVED',
            'metadata_normalization': 'Whitelisted RDM vocabulary decorations, HTML character references, and harmless ASCII whitespace at document/block boundaries; exact JSON value types, inline spaces, declarations and processing instructions retained',
            'trust_boundary': 'Externally supplied witness pins authenticate saved evidence; no live remote authentication occurs in this verifier',
            'network_calls': 0, 'numeric_imports': 0, 'source_callbacks': 0, 'array_decodes': 0,
            'original_binary80_state_error': 'ENCLOSED_RETAINED_NODES',
            'state_error_scope': '49,152 saved capsule-node occurrences; conditional fixed-weight transport under identical subsequent forcing/contacts',
            'scientific_status_basis': 'Externally pinned edition contract; this publication verifier does not reprove the numerical enclosure',
            'full_twelve_case_certificate': 'UNRESOLVED',
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
        require(type(group) is dict and set(group) == {'record', 'files', 'baseline', 'baseline_files'}
                and all(type(v) is str and v in values for v in group.values()), 'Preserved witnesses incomplete')
        for role, pin in BASELINE_PINS[record_id].items():
            require(sha(values[group[role]]) == pin, 'Dated preserved baseline source pin differs: ' + record_id + '/' + role)
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
