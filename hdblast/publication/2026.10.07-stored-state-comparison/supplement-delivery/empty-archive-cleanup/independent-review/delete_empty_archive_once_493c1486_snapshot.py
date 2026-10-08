#!/usr/bin/env python3
"""One reviewed DELETE of our exact unavailable initialized draft archive only.

Default is GET-only. Original controller state, unknown PUT and both failed
recovery attempts remain unchanged. Root owns the explicit mutation invocation.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path('/workspace/hdblast-research-work/continuation-state-comparison-20261007')
HERE = ROOT / 'publication-final/empty-archive-cleanup'
SIDE = HERE / 'live-delete-001'
REVIEWED = ROOT / 'publication-final/unknown-put-recovery/recover_put_chunked_once.py'
REVIEWED_SHA = 'f76aeb162e6d1cad6f9d53afbf20bf68f7f4ee0f4b1504789f2880f97a4315eb'
SECOND_FAILURE = REVIEWED.parent / 'live-chunked-attempt-001'
SECOND_FAILURE_PINS = {
    'ONE_PUT_ATTEMPT_LATCH.json': 'e194ea7422f0f1e36b6aa8a4ce5d86ec741a46e4d966272a0bc49e95dca46081',
    'RECOVERY_RESULT.json': '93b432ec3ad3ef04bb9e04a0c52245af171d0647bab881aabf9bfd24d25109e5',
    'curl_body.raw': 'f8ac594b632063665f0991ca2c37e9b79e3bda1218c2404ad1d36a23134472c7',
    'curl_headers.raw': '0124837ddaee9cb65b10864e1c30e14aa6c39a6f3c440a38c4e1f2111d3c2bec',
    'curl_stderr.raw': '3eb77d2fec058bf3791641be32d61a5157f3f017cc7d0e78492618f802e7ee98'}
EVIDENCE = ROOT / 'publication-final/new-edition-execution'
ORIGINAL_STATE_SHA = '9927d046277e4c808e43fff49d3c6add2495df8795cff7b166941694830132ad'
ORIGINAL_JOURNAL_SHA = 'a413a83242ab9bfa68557b0a1f8230f3e9507f0525f48b349ae92cba1ddfb162'
DRAFT_BASELINE = ROOT / 'publication-final/unknown-put-read-only/draft_vendor.body'
DRAFT_BASELINE_SHA = 'fc5fe357c2171567de78d5baa874a420e3c1a5782508097aebafe1997cb86035'
DELETE_URL = 'https://zenodo.org/api/records/23228395/draft/files/HDBLAST_CHECKPOINT_20261007_STORED_STATE_COMPARISON.zip'
NATIVE_URL = 'https://zenodo.org/api/files/b183fa23-92b9-4f4c-8b18-815ebdc4da12/HDBLAST_CHECKPOINT_20261007_STORED_STATE_COMPARISON.zip?version_id=21a587cc-fcd9-4255-9423-f1acfba25013'

class Stop(RuntimeError):
    pass

def require(ok, code):
    if not ok:
        raise Stop(code)

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def load_reviewed():
    raw = REVIEWED.read_bytes()
    require(sha(raw) == REVIEWED_SHA, 'REVIEWED_CHUNKED_SOURCE_PIN_DIFFERS')
    name = '_reviewed_empty_archive_cleanup_primitives'
    module = importlib.util.module_from_spec(importlib.util.spec_from_loader(name, loader=None, origin=str(REVIEWED)))
    module.__file__ = str(REVIEWED)
    sys.modules[name] = module
    exec(compile(raw, str(REVIEWED), 'exec'), module.__dict__)
    return module

def verify_two_failures(h):
    first = h.verify_prior_failure()
    raw = {name: h.read_pinned(SECOND_FAILURE / name, pin) for name, pin in SECOND_FAILURE_PINS.items()}
    intent = json.loads(raw['ONE_PUT_ATTEMPT_LATCH.json'])
    require(intent.get('url') == h.URL and intent.get('asset') == h.ASSET_PIN
        and intent.get('controller_sha256') == h.CONTROLLER_SHA
        and intent.get('whole_body_transfer_encoding') == 'chunked', 'SECOND_FAILURE_INTENT_DIFFERS')
    result = json.loads(raw['RECOVERY_RESULT.json'])
    metrics = result.get('result', {}).get('metrics', {})
    require(result.get('original_state_journal_unchanged') is True and result.get('pending_cleared') is False
        and type(result.get('commits')) is int and result['commits'] == 0
        and type(result.get('publishes')) is int and result['publishes'] == 0, 'SECOND_FAILURE_SCOPE_DIFFERS')
    require(type(metrics.get('http_status')) is int and metrics['http_status'] == 401
        and type(metrics.get('uploaded_bytes')) is int and metrics['uploaded_bytes'] == 30210945,
        'SECOND_FAILURE_NOT_KNOWN_401')
    require(raw['curl_body.raw'] == b'Unauthorized: authentication failed for integration codex-secret-zenodo.org\n'
        and b'x-at-upstream-error: false\r\n' in raw['curl_headers.raw'].lower()
        and b'HTTP/1.1 401 Unauthorized\r\n' in raw['curl_headers.raw'], 'SECOND_FAILURE_BRIDGE_INDICATORS_DIFFERS')
    return {'first_failure': first, 'second_source_sha256': REVIEWED_SHA, 'second_failure_pins': SECOND_FAILURE_PINS,
        'both_http_statuses': [401, 401], 'both_upstream_error': False, 'histories_preserved': True,
        'upload_cap': 'NOT_ESTABLISHED'}

def make_view(h, module, token, side):
    class View(module.Controller):
        capture_number = 0
        def save(self):
            pass
        def write(self, *args, **kwargs):
            raise Stop('ORIGINAL_CONTROLLER_WRITES_FORBIDDEN')
        def event(self, kind, **values):
            row = h.redact_value({'utc': h.utc(), 'kind': kind, **values}, token)
            with (side / 'GUARD_JOURNAL.jsonl').open('a') as out:
                out.write(json.dumps(row) + '\n')
                out.flush()
                os.fsync(out.fileno())
        def capture(self, label, obj):
            self.capture_number += 1
            h.write_json(side / (str(self.capture_number) + '_' + label + '.json'), obj, token)
    c = View(EVIDENCE, module.Transport(token), token,
        json.loads(h.read_pinned(h.INVENTORY, h.INVENTORY_SHA)),
        json.loads(h.read_pinned(h.METADATA, h.METADATA_SHA)), h.INVENTORY_SHA, h.METADATA_SHA)
    c.module = module
    return c

def guard(c, h, initialized, expected_editable, target_present):
    h.validate_state(c.state)
    c.current()
    require(c.state.get('published') is False, 'NEW_EDITION_ALREADY_PUBLISHED')
    owner = c.owner_list()
    require(owner and str(owner.get('id')) == h.DRAFT, 'FRESH_OWNED_DRAFT_DIFFERS')
    draft, headers, files = c.draft(h.DRAFT)
    c.module.metadata_equal(draft, expected_editable)
    expected = {p['filename'] for p in c.old} | ({h.NAME} if target_present else set())
    require(set(files) == expected, 'EXACT_INHERITED_AND_TARGET_MEMBERSHIP_DIFFERS')
    c.module.assert_pins(files, c.old, exact=False)
    for pin in c.old:
        f = files[pin['filename']]
        require(c.module.immutable_file_identity(f) == c.module.immutable_file_identity(initialized[pin['filename']])
            and f.get('bucket_id') == h.BUCKET_ID, 'INHERITED_DRAFT_IDENTITY_OR_BUCKET_DIFFERS')
    if target_present:
        h.validate_remote(c, files, initialized)
        target = files[h.NAME]
        require(target.get('links', {}).get('self') == DELETE_URL and target.get('key') == h.NAME,
            'TARGET_SELF_LINK_OR_KEY_DIFFERS')
        # Actual initialized pending metadata omits size/checksum; explicit null
        # is also unavailable. Any non-null value, including zero, forbids DELETE.
        require(target.get('size') is None and target.get('checksum') is None, 'TARGET_HAS_STORAGE_METADATA_NO_DELETE')
        identity = c.module.immutable_file_identity(target)
        require(all(identity[0] != c.module.immutable_file_identity(files[p['filename']])[0]
            and identity[1] != c.module.immutable_file_identity(files[p['filename']])[1] for p in c.old),
            'TARGET_IDENTITY_COLLIDES_WITH_INHERITED')
    return draft, files

def unavailable(c, h, side, token, label, url):
    require(url in (h.URL, NATIVE_URL), 'UNREVIEWED_UNAVAILABLE_PROBE_URL')
    response = c.transport.request('GET', url, headers={'Accept': 'application/json', 'Accept-Encoding': 'identity'})
    h.write_json(side / (label + '_UNAVAILABLE.json'), {'utc': h.utc(), 'method': 'GET', 'url': url,
        'http_status': response.status, 'body_bytes': len(response.body), 'body_sha256': sha(response.body),
        'body': response.body.decode(errors='replace')}, token)
    h.validate_missing(response)
    return {'http_status': 400, 'body_sha256': sha(response.body)}

def original_unchanged(before):
    require(before == {n: sha((EVIDENCE / n).read_bytes()) for n in before}, 'ORIGINAL_STATE_OR_JOURNAL_CHANGED')

def delete_once(c, h, side, token, initialized, expected_editable, failure_proof, before):
    guard(c, h, initialized, expected_editable, True)
    modern = unavailable(c, h, side, token, 'MODERN', h.URL)
    native = unavailable(c, h, side, token, 'NATIVE_EXACT_VERSION', NATIVE_URL)
    # Re-read all owned draft/inherited/empty target identities after both probes.
    guard(c, h, initialized, expected_editable, True)
    original_unchanged(before)
    require(sha(REVIEWED.read_bytes()) == REVIEWED_SHA, 'REVIEWED_SOURCE_CHANGED_BEFORE_DELETE')
    intent = {'utc': h.utc(), 'operation': 'ONE_DELETE_EXACT_OWN_EMPTY_INITIALIZED_ARCHIVE', 'method': 'DELETE',
        'url': DELETE_URL, 'draft_id': h.DRAFT, 'filename': h.NAME, 'file_id': h.FILE_ID,
        'version_id': h.VERSION_ID, 'bucket_id': h.BUCKET_ID, 'modern_unavailable': modern,
        'native_unavailable': native, 'failure_proof': failure_proof, 'original_state_journal': before,
        'old_files_deleted': 0, 'other_user_files_deleted': 0, 'original_pending_cleared': False}
    h.write_json(side / 'ONE_DELETE_ATTEMPT_LATCH.json', intent, token)
    try:
        response = c.transport.request('DELETE', DELETE_URL, headers={'Accept': 'application/json'})
    except Exception as exception:
        h.write_json(side / 'DELETE_UNKNOWN.json', {'status': 'OUTCOME_UNKNOWN_GET_ONLY_NO_REPEAT',
            'exception_type': type(exception).__name__, 'reason': str(exception)}, token)
        raise Stop('DELETE_OUTCOME_UNKNOWN_GET_ONLY_NO_REPEAT') from None
    h.write_json(side / 'DELETE_RESPONSE.json', {'http_status': response.status, 'body_bytes': len(response.body),
        'body_sha256': sha(response.body), 'body': response.body.decode(errors='replace'),
        'body_truncated': response.body_truncated, 'body_read_failed': response.body_read_failed}, token)
    require(type(response.status) is int and response.status in (200, 202, 204)
        and not response.body_truncated and not response.body_read_failed, 'DELETE_NOT_CONFIRMED_GET_ONLY_NO_REPEAT')
    guard(c, h, initialized, expected_editable, False)
    original_unchanged(before)
    return {'status': 'PASS_EXACT_EMPTY_OWN_ARCHIVE_ABSENT_27_INHERITED_PRESERVED', 'deleted_file': h.NAME,
        'original_pending_cleared': False, 'original_state_journal_unchanged': True,
        'inherited_files': len(c.old), 'metadata_writes': 0, 'publishes': 0, 'new_drafts': 0}

def reconcile(c, h, initialized, expected_editable):
    # Absence proof is the complete fresh owned draft/list, not a content error.
    guard(c, h, initialized, expected_editable, False)
    return {'status': 'PASS_GET_ONLY_EXACT_TARGET_ABSENT_27_INHERITED_PRESERVED', 'writes': 0,
        'original_pending_cleared': False, 'repeat_delete_authorized': False}

def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--helper-sha256', required=True)
    parser.add_argument('--delete-empty-owned-archive', action='store_true')
    args = parser.parse_args(argv)
    require(sha(Path(__file__).read_bytes()) == args.helper_sha256, 'EXTERNAL_CLEANUP_HELPER_PIN_DIFFERS')
    token = os.environ.get('ZENODO_ACCESS_TOKEN', '')
    h = load_reviewed()
    h.token_config(token)
    proof = verify_two_failures(h)
    module = h.load_controller()
    if args.delete_empty_owned_archive:
        require(not SIDE.exists(), 'FIXED_CLEANUP_SIDE_ALREADY_EXISTS_GET_ONLY_NO_REPEAT')
        SIDE.mkdir(mode=0o700)
        parent_fd = os.open(SIDE.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        try:
            os.fsync(parent_fd)
        finally:
            os.close(parent_fd)
        side = SIDE
    else:
        side = HERE / ('get-only-' + str(time.time_ns()))
        side.mkdir(mode=0o700)
    before = {n: sha((EVIDENCE / n).read_bytes()) for n in ('STATE.json', 'JOURNAL.jsonl')}
    try:
        with module.EvidenceLock(EVIDENCE):
            require(before == {'STATE.json': ORIGINAL_STATE_SHA, 'JOURNAL.jsonl': ORIGINAL_JOURNAL_SHA}, 'ORIGINAL_EVIDENCE_PINS_DIFFERS')
            c = make_view(h, module, token, side)
            initialized = module.file_map(json.loads(h.read_pinned(h.INITIALIZATION, h.INITIALIZATION_SHA)))
            baseline = json.loads(h.read_pinned(DRAFT_BASELINE, DRAFT_BASELINE_SHA))
            editable = {k: baseline[k] for k in ('metadata', 'custom_fields', 'access')}
            result = delete_once(c, h, side, token, initialized, editable, proof, before) if args.delete_empty_owned_archive else reconcile(c, h, initialized, editable)
            original_unchanged(before)
            h.write_json(side / 'RESULT.json', result, token)
            print(json.dumps(result))
    except Exception as exception:
        original_unchanged(before)
        h.write_json(side / 'STOPPED.json', {'status': 'STOPPED_GET_ONLY_NO_REPEAT', 'exception_type': type(exception).__name__,
            'reason': str(exception), 'original_pending_cleared': False}, token)
        raise
    return 0

if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except Exception as exception:
        print(json.dumps({'status': 'STOPPED_GET_ONLY_NO_REPEAT', 'reason': str(exception) if isinstance(exception, Stop) else 'RECORDED_UNEXPECTED_FAILURE'}))
        raise SystemExit(2)
