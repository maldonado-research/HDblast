#!/usr/bin/env python3
"""Independent local-data journal audit; never import or invoke publisher code."""
import argparse
import collections
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import sys
from urllib.parse import quote

BASE = Path('/workspace/hdblast-research-work/continuation-trajectory-20261008')
OWN = BASE / 'final-publication-independent-review/journal-audit'
EVIDENCE = BASE / 'publication-execution-001'
INVENTORY = BASE / 'release-assets/sealed-candidate-001/SEALED_UPLOAD_INVENTORY.json'
METADATA = BASE / 'release-assets/sealed-candidate-001/SEALED_METADATA_MODERN.json'
CONTROLLER = BASE / 'publication-planning/new_edition_controller_v2.py'
INVENTORY_SHA = '2c6d0a3d4c4e227c3a8aa5f6e48dd3ff8494ac4e2f1b9045c603cf35617ced3c'
METADATA_SHA = '43feec5ca6bd9ab8c5c6a2cd925f6f42ff13418bd90b334de7ab604367dc385c'
CONTROLLER_SHA = '4d6dfa5e9ca2add8654e5c1803f5b65fa25ecce87d71c3a7590fbb4d29f9f117'
RID, PRIOR, FAMILY, OWNER = '23244754', '23228395', '17088132', '1386319'
ORIGIN = 'https://zenodo.org/api/records/'
PUBLIC, DRAFT = ORIGIN + RID, ORIGIN + RID + '/draft'
VERSION = '2026.10.08-later-state-and-bulk-response'
SHA = re.compile('[0-9a-f]{64}')

def require(condition, message):
    if not condition:
        raise ValueError(message)

def digest(raw):
    return hashlib.sha256(raw).hexdigest()

def capture(path):
    require(all(not p.is_symlink() for p in (path, *path.parents)), 'Symlink: ' + str(path))
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    try:
        first = os.fstat(fd)
        require(stat.S_ISREG(first.st_mode) and first.st_nlink == 1 and first.st_size <= 20*1024*1024,
                'Unsafe/noncompact evidence: ' + str(path))
        raw = b''
        while chunk := os.read(fd, 65536):
            raw += chunk
            require(len(raw) <= 20*1024*1024, 'Evidence grew too large')
        last = os.fstat(fd)
        named = path.stat(follow_symlinks=False)
        ident = lambda s: (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns, s.st_nlink)
        require(ident(first) == ident(last) == ident(named), 'Evidence changed while captured: ' + str(path))
        return raw
    finally:
        os.close(fd)

def load(raw):
    def unique(items):
        answer = {}
        for key, value in items:
            require(key not in answer, 'Duplicate JSON key')
            answer[key] = value
        return answer
    def reject(value):
        raise ValueError('Nonfinite JSON constant: ' + value)
    obj = json.loads(raw, object_pairs_hook=unique, parse_constant=reject)
    json.dumps(obj, allow_nan=False)
    return obj

def pin(path, raw):
    return {'path': str(path), 'bytes': len(raw), 'sha256': digest(raw)}

def file_map(obj):
    entries = obj.get('entries', obj.get('files', {}).get('entries'))
    require(type(entries) in (dict, list), 'File collection missing')
    rows = list(entries.values()) if type(entries) is dict else entries
    answer = {}
    for row in rows:
        name = row.get('key')
        require(type(name) is str and name not in answer, 'Duplicate file entry')
        answer[name] = row
    if type(entries) is dict:
        require(set(answer) == set(entries), 'Embedded key differs')
    return answer

def check_files(obj, pins, prefix, embedded=False):
    rows = file_map(obj)
    require(set(rows) == set(pins), 'File membership differs')
    for name, row in rows.items():
        p = pins[name]
        require(type(row.get('size')) is int and row['size'] == p['bytes'] and
                row.get('checksum') == 'md5:' + p['md5'], 'Remote size/MD5 differs: ' + name)
        require(type(row.get('file_id', row.get('id'))) is str, 'File identity missing')
        if not embedded:
            require(row.get('status') == 'completed', 'File incomplete: ' + name)
            require(row.get('links', {}).get('content') == prefix + '/files/' + quote(name, safe='') + '/content',
                    'Content endpoint differs: ' + name)
            require(type(row.get('version_id')) is str and row['version_id'], 'Version identity missing')
    return rows

def identity(record, draft):
    parent = record.get('parent', {})
    require(record.get('id') == RID and parent.get('id') == FAMILY, 'Record/family identity differs')
    require(str(parent.get('access', {}).get('owned_by', {}).get('user')) == OWNER, 'Owner differs')
    require(parent.get('pids', {}).get('doi', {}).get('identifier') == '10.5281/zenodo.' + FAMILY,
            'Concept DOI differs')
    require(record.get('is_draft') is draft and record.get('is_published') is (not draft), 'Record state differs')
    if not draft or 'doi' in record.get('pids', {}):
        require(record.get('pids', {}).get('doi', {}).get('identifier') == '10.5281/zenodo.' + RID,
                'Version DOI differs')
    require(record.get('versions', {}).get('is_latest_draft' if draft else 'is_latest') is True, 'Latest flag differs')
    require(record.get('metadata', {}).get('version') == VERSION, 'Semantic version differs')
    return {'id': RID, 'parent_id': FAMILY, 'owner': OWNER,
            'created': record.get('created'), 'updated': record.get('updated'),
            'publication_date': record.get('metadata', {}).get('publication_date'),
            'revision_id': record.get('revision_id')}

def main():
    require(sys.flags.isolated and sys.flags.dont_write_bytecode, 'Use Python -I -B')
    parser = argparse.ArgumentParser()
    parser.add_argument('--snapshot', required=True)
    parser.add_argument('--partial', action='store_true')
    parser.add_argument('--journal-sha256')
    args = parser.parse_args()
    require(re.fullmatch('[a-z0-9-]+', args.snapshot), 'Unsafe snapshot name')
    out = OWN / args.snapshot
    out.mkdir()
    inputs = {}
    for label, path, expected in [('inventory', INVENTORY, INVENTORY_SHA), ('metadata', METADATA, METADATA_SHA),
                                   ('controller', CONTROLLER, CONTROLLER_SHA)]:
        raw = capture(path)
        require(digest(raw) == expected, 'External source/input pin differs: ' + label)
        inputs[label] = pin(path, raw)
        (out / path.name).write_bytes(raw)
    inventory = load((out / INVENTORY.name).read_bytes())
    metadata = load((out / METADATA.name).read_bytes())
    rows = inventory['inherited_files'] + inventory['new_files']
    pins = {x['filename']: x for x in rows}
    require(len(inventory['inherited_files']) == 29 and len(inventory['new_files']) == 2 and len(pins) == 31,
            'Expected 29 plus 2 unique files')
    require(sum(p['bytes'] for p in pins.values()) == inventory['total_bytes_final'] == 453646943,
            'Total bytes differ')
    raw = capture(EVIDENCE / 'JOURNAL.jsonl')
    if not args.partial:
        require(type(args.journal_sha256) is str and SHA.fullmatch(args.journal_sha256) and
                digest(raw) == args.journal_sha256, 'Root terminal journal external pin required/differs')
    require(raw.endswith(b'\n'), 'Truncated final JSONL line')
    partial_raw = capture(OWN / 'JOURNAL_PARTIAL_OBSERVATION_001.jsonl')
    require(raw.startswith(partial_raw), 'Previously observed journal prefix changed')
    (out / 'JOURNAL.jsonl').write_bytes(raw)
    evidence_pins = {'JOURNAL.jsonl': pin(EVIDENCE / 'JOURNAL.jsonl', raw)}
    journal = [load(line) for line in raw.splitlines()]
    stamps = [datetime.datetime.fromisoformat(e['utc']) for e in journal]
    require(all(x.utcoffset() == datetime.timedelta(0) and x.date().isoformat() == '2026-10-08' for x in stamps),
            'Unexpected journal date/timezone')
    require(stamps == sorted(stamps), 'Journal UTC timestamps regressed')
    allowed = {'GET', 'LATEST_LINK_GET', 'EXPLICIT_PUBLIC_LATEST_RESOLUTION', 'OWNER_PAGE_GET',
               'WRITE_INTENT', 'WRITE_RESPONSE', 'OWNED_LATEST_ONLY_JOURNALED_RECORD_MATCH',
               'OWNED_DRAFT_BOUND_FROM_GETS', 'FULL_CONTENT_STREAM_VERIFIED', 'PUBLICATION_GET_POLL'}
    require(all(e.get('kind') in allowed for e in journal), 'Unknown/failure/reconciliation event preserved; requires review')
    names = [p['filename'] for p in rows]
    empty_sha = digest(b'')
    expected_writes = [('create_version', 'https://zenodo.org/api/deposit/depositions/' + PRIOR + '/actions/newversion',
                        empty_sha, 0, {})]
    for p in inventory['new_files']:
        filename = p['filename']
        initializer = json.dumps([{'key': filename}], ensure_ascii=False).encode()
        endpoint = DRAFT + '/files/' + quote(filename, safe='')
        expected_writes.extend([
            ('initialize_file', DRAFT + '/files', digest(initializer), len(initializer), {'filename': filename}),
            ('put_content', endpoint + '/content', p['sha256'], p['bytes'], {'filename': filename}),
            ('commit_file', endpoint + '/commit', empty_sha, 0, {'filename': filename})])
    metaraw = json.dumps(metadata, ensure_ascii=False).encode()
    expected_writes += [('save_metadata', DRAFT, digest(metaraw), len(metaraw), {}),
                        ('publish', DRAFT + '/actions/publish', empty_sha, 0, {})]
    intents, responses, rounds = [], [], []
    pending = None
    group = None
    def finish_group():
        nonlocal group
        if group and group['filenames']:
            group['complete'] = group['filenames'] == names
            require(group['complete'] or args.partial, 'Partial/wrong-order content round')
            rounds.append(group)
        group = None
    for index, e in enumerate(journal):
        kind = e['kind']
        if kind != 'FULL_CONTENT_STREAM_VERIFIED' and group is not None:
            finish_group()
        if kind == 'GET':
            require(type(e.get('status')) is int and e['status'] == 200 and SHA.fullmatch(e.get('response_sha256', '')),
                    'GET failed or hash malformed')
            allowed_urls = {ORIGIN + PRIOR, ORIGIN + PRIOR + '/files', PUBLIC, PUBLIC + '/files', DRAFT,
                            DRAFT + '/files', 'https://zenodo.org/api/deposit/depositions/' + PRIOR}
            require(e.get('url') in allowed_urls, 'Unexpected GET endpoint')
            authenticated = e['url'].startswith(DRAFT) or '/api/deposit/' in e['url']
            require(e.get('authenticated_request') is authenticated, 'GET auth mode differs')
            if e['url'] in (DRAFT + '/files', PUBLIC + '/files'):
                group = {'scope': 'draft' if e['url'].startswith(DRAFT) else 'public',
                         'list_event_index': index, 'list_response_sha256': e['response_sha256'],
                         'filenames': [], 'first_stream_index': None, 'last_stream_index': None}
        elif kind == 'FULL_CONTENT_STREAM_VERIFIED':
            require(group is not None, 'Unscoped content stream')
            name = e.get('filename')
            require(name in pins and name not in group['filenames'], 'Unknown/duplicate stream filename')
            observed = {k: pins[name][k] for k in ('bytes', 'md5', 'sha256')}
            require(type(e.get('observed', {}).get('bytes')) is int and e['observed'] == observed,
                    'Incomplete/wrong bytes or checksum: ' + name)
            if group['first_stream_index'] is None:
                group['first_stream_index'] = index
            group['last_stream_index'] = index
            group['filenames'].append(name)
        elif kind == 'WRITE_INTENT':
            require(pending is None and len(intents) < len(expected_writes), 'Concurrent/extra remote write')
            actual = tuple(e.get(k) for k in ('operation', 'url', 'body_sha256', 'body_bytes', 'context'))
            require(type(e.get('body_bytes')) is int and actual == expected_writes[len(intents)],
                    'Write operation, target, body or order differs')
            pending = e['operation']
            intents.append({'index': index, **e})
        elif kind == 'WRITE_RESPONSE':
            require(pending is not None and e.get('operation') == pending, 'Unmatched/duplicate write response')
            require(type(e.get('status')) is int and e['status'] in (200, 201, 202, 204) and
                    e.get('body_truncated') is False and e.get('body_read_failed') is False and
                    SHA.fullmatch(e.get('response_sha256', '')), 'Unsuccessful/unknown/truncated write outcome')
            responses.append({'index': index, **e})
            pending = None
        elif kind == 'LATEST_LINK_GET':
            require(e.get('url') == ORIGIN + PRIOR + '/versions/latest' and e.get('authenticated_request') is False
                    and e.get('status') in (200, 301, 302, 303, 307, 308), 'Latest lookup differs')
        elif kind == 'EXPLICIT_PUBLIC_LATEST_RESOLUTION':
            require(e.get('source') == ORIGIN + PRIOR + '/versions/latest' and
                    e.get('target') in (ORIGIN + PRIOR, PUBLIC), 'Wrong latest resolution')
        elif kind == 'OWNER_PAGE_GET':
            require(e.get('status') == 200 and type(e.get('page')) is int and e['page'] >= 1, 'Owner page failed')
        elif kind in ('OWNED_LATEST_ONLY_JOURNALED_RECORD_MATCH', 'OWNED_DRAFT_BOUND_FROM_GETS'):
            require(e.get('record_id') == RID, 'Owner/draft binding differs')
        elif kind == 'PUBLICATION_GET_POLL':
            require(e.get('status') in (200, 404, 503) and any(x['operation'] == 'publish' for x in responses),
                    'Unscoped public poll')
    finish_group()
    draft_rounds = [x for x in rounds if x['scope'] == 'draft']
    public_rounds = [x for x in rounds if x['scope'] == 'public']
    pub = [x for x in intents if x['operation'] == 'publish']
    record_checks = {}
    if not args.partial:
        require(pending is None and len(intents) == len(responses) == 9 and len(pub) == 1, 'Publication operation count differs')
        require(len(draft_rounds) >= 2 and all(x['complete'] and x['last_stream_index'] < pub[0]['index'] for x in draft_rounds),
                'Two complete draft rounds before publication required')
        require(len(public_rounds) == 1 and public_rounds[0]['complete'] and public_rounds[0]['list_event_index'] > pub[0]['index'],
                'Exactly one complete public round after publication required')
        require(not any(x['operation'] != 'publish' for x in intents if x['index'] > draft_rounds[-1]['last_stream_index']),
                'Mutation follows last complete draft round')
        files_to_capture = sorted(p for p in EVIDENCE.iterdir() if p.suffix == '.json')
        docs = {}
        for path in files_to_capture:
            data = capture(path)
            docs[path.name] = load(data)
            (out / path.name).write_bytes(data)
            evidence_pins[path.name] = pin(path, data)
        state = docs['STATE.json']
        require(state.get('pending') is None and state.get('draft_id') == RID and state.get('published') is True and
                state.get('create_post_attempted') is True and state.get('publish_post_attempted') is True,
                'Final persisted publication state differs')
        require(state.get('inventory_sha256') == INVENTORY_SHA and state.get('metadata_sha256') == METADATA_SHA and
                state.get('continuation_binding', {}).get('controller_sha256') == CONTROLLER_SHA, 'State source/input pins differ')
        require(state.get('verified_receipt_sha256') == evidence_pins['VERIFIED_DRAFT.json']['sha256'],
                'State final draft receipt pin differs')
        for scope, recordname, listname, receiptname, prefix in [
            ('draft', 'DRAFT_' + RID + '.json', 'DRAFT_FILES_' + RID + '.json', 'VERIFIED_DRAFT.json', DRAFT),
            ('public', 'NEW_PUBLIC_RECORD.json', 'NEW_PUBLIC_FILES.json', 'VERIFIED_PUBLIC.json', PUBLIC)]:
            record, listing, receipt = docs[recordname], docs[listname], docs[receiptname]
            draft = scope == 'draft'
            record_checks[scope] = identity(record, draft)
            remote = check_files(listing, pins, prefix)
            embedded = check_files(record, pins, prefix, embedded=True)
            for name in pins:
                require(embedded[name]['id'] == remote[name]['file_id'], 'Embedded file identity differs')
            require(receipt.get('status') == 'PASS_COMPLETE_31_FILE_EDITION_READBACK' and
                    receipt.get('draft_id') == RID and receipt.get('parent_id') == FAMILY and
                    receipt.get('published') is (not draft) and receipt.get('file_count') == 31 and
                    receipt.get('total_bytes') == 453646943 and receipt.get('inventory_sha256') == INVENTORY_SHA and
                    receipt.get('metadata_sha256') == METADATA_SHA, 'Verification receipt scope differs')
            content = receipt.get('content')
            require(type(content) is list and [x.get('filename') for x in content] == names, 'Receipt content roster differs')
            for row in content:
                name = row['filename']
                require(row.get('sha256') == pins[name]['sha256'] and row.get('basis') == 'FRESH_COMPLETE_CONTENT_STREAM' and
                        row.get('remote_file_id') == remote[name]['file_id'] and
                        row.get('remote_version_id') == remote[name]['version_id'], 'Receipt identity/content differs')
            record_checks[scope]['receipt_utc'] = receipt['utc']
            matching = draft_rounds[-1] if draft else public_rounds[-1]
            require(stamps[matching['last_stream_index']] <= datetime.datetime.fromisoformat(receipt['utc']),
                    'Receipt predates final stream')
            if draft:
                require(datetime.datetime.fromisoformat(receipt['utc']) < stamps[pub[0]['index']], 'Draft receipt after publish')
                require(type(record.get('revision_id')) is int and receipt.get('etag') == '"' + str(record['revision_id']) + '"',
                        'Draft ETag/revision differs')
        record_checks['latest'] = identity(docs['CURRENT_LATEST.json'], False)
        require([(x['remote_file_id'], x['remote_version_id']) for x in docs['VERIFIED_DRAFT.json']['content']] ==
                [(x['remote_file_id'], x['remote_version_id']) for x in docs['VERIFIED_PUBLIC.json']['content']],
                'Draft/public content identities changed')
        require(capture(EVIDENCE / 'JOURNAL.jsonl') == raw, 'Final journal changed during evidence capture')
        require(capture(EVIDENCE / 'STATE.json') == (out / 'STATE.json').read_bytes(), 'State changed during capture')
    result = {
        'status': 'PARTIAL_OBSERVATION_ONLY' if args.partial else 'PASS_INDEPENDENT_SAVED_PUBLICATION_JOURNAL_AUDIT',
        'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'input_pins': inputs, 'evidence_pins': evidence_pins,
        'record_id': RID, 'parent_id': FAMILY, 'owner': OWNER,
        'journal_event_count': len(journal), 'event_counts': dict(collections.Counter(x['kind'] for x in journal)),
        'first_event_utc': journal[0]['utc'], 'last_event_utc': journal[-1]['utc'],
        'write_intents': intents, 'write_responses': responses, 'stream_rounds': rounds,
        'complete_draft_rounds': sum(x['complete'] for x in draft_rounds),
        'complete_public_rounds': sum(x['complete'] for x in public_rounds),
        'streams_verified': sum(len(x['filenames']) for x in rounds),
        'bytes_per_complete_round': 453646943, 'record_checks': record_checks,
        'prior_partial_prefix_preserved': pin(OWN / 'JOURNAL_PARTIAL_OBSERVATION_001.jsonl', partial_raw),
        'auditor_activity': {'network_calls': 0, 'remote_mutations': 0, 'attachment_streams': 0,
                             'publisher_imports_or_execution': 0, 'scientific_source_calls': 0},
        'source_bound_stream_semantics': 'Pinned controller requires full HTTP200 at exact content endpoint; no Content-Range; identity content encoding; complete expected byte count plus MD5 and SHA256 before logging each stream. Record and content identities checked by saved listings/receipts.',
        'limits': [
            'This is local saved-evidence verification, not an independent fresh attachment fetch or remote attestation.',
            'Stream events omit URL and HTTP headers; their endpoint/status semantics rely on the externally pinned controller and scoped listing sequence.',
            'Journal is fsynced append-only by accepted source, but has no cryptographic hash chain or general forged/missing-history detection.',
            'JSON captures are normalized/redacted reserializations, so their byte hashes do not generally equal raw HTTP response_sha256 journal fields.',
            'Saved response filenames preserve per-write objects; repeated GET labels and verification receipts retain only their latest values, with earlier raw response hashes/times in journal.',
            'No account-wide single-publisher exclusion is established; exact recorded writes are checked.',
            'Complete editable metadata normalization and historical-record preservation are separate parent review obligations.'
        ]
    }
    (out / 'AUDIT_RECEIPT.json').write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': result['status'], 'receipt': str(out / 'AUDIT_RECEIPT.json'),
                      'sha256': digest((out / 'AUDIT_RECEIPT.json').read_bytes()),
                      'events': len(journal), 'streams': result['streams_verified']}))

if __name__ == '__main__':
    main()
