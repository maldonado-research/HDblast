#!/usr/bin/env python3
"""One explicitly authorized identical PUT; all original controller files remain intact.

This preparation revision has no accepted live invocation yet. Its transport is
deliberately a single curl PUT with an empty Expect header: libcurl's HTTP417
fallback must not turn a subprocess attempt into two wire requests.
"""
from __future__ import annotations
import argparse
import datetime
import fcntl
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import selectors
import stat
import subprocess
import sys
import time

ROOT = Path('/workspace/hdblast-research-work/continuation-state-comparison-20261007')
HERE = ROOT / 'publication-final/unknown-put-recovery'
EVIDENCE = ROOT / 'publication-final/new-edition-execution'
CONTROLLER = ROOT / 'publication-planning/new_edition_controller.py'
CONTROLLER_SHA = 'aa7e6cfe271559731ed3b45c5339f8eb612ebe6f017ddafa5fdadad9579e084b'
INVENTORY = ROOT / 'publication-final/FINAL_UPLOAD_INVENTORY.json'
INVENTORY_SHA = 'f6b7db33e1cefc22e9868569700f90969368f058d5c4fdc59251cdc810e45b21'
METADATA = ROOT / 'publication-final/FINAL_METADATA_MODERN.json'
METADATA_SHA = '0d0896174064c4e5c796d2a14becd817467d457774a24e82f5ce61772b271ef8'
INITIALIZATION = EVIDENCE / 'response_initialize_file_1791431996909506944.json'
INITIALIZATION_SHA = '7b4e819dd237601efcbf4d6c68de6b968585f9173ad324f6dd15bf09376fcb88'
ORIGINAL_STATE_SHA = '9927d046277e4c808e43fff49d3c6add2495df8795cff7b166941694830132ad'
ORIGINAL_JOURNAL_SHA = 'ef1ba41a215b203356efd46564919dff7df30d5c9f9c8825c346fdcb7c6b5f60'
DRAFT = '23228395'
NAME = 'HDBLAST_CHECKPOINT_20261007_STORED_STATE_COMPARISON.zip'
URL = 'https://zenodo.org/api/records/' + DRAFT + '/draft/files/' + NAME + '/content'
FILE_ID = '4d09e310-3b8c-401c-bfea-585bea657e25'
VERSION_ID = '21a587cc-fcd9-4255-9423-f1acfba25013'
BUCKET_ID = 'b183fa23-92b9-4f4c-8b18-815ebdc4da12'
ASSET_PIN = {'filename': NAME, 'bytes': 240135519,
    'sha256': 'ead5ffa3987c760d99da7a79189ff5ad766fb9bef6e032c489f4efb2316cd189',
    'md5': 'e424502f12dbaacd36abaaabdd3fb74e'}
SIDE = HERE / 'live-attempt-001'
MAX_JSON = 8 * 1024 * 1024
CURL = '/usr/bin/curl'
# Linux UAPI linux/fcntl.h: F_LINUX_SPECIFIC_BASE=1024; ADD_SEALS=+9,
# GET_SEALS=+10 and the four immutable memfd bits below. Python 3.12 in
# this reviewed execution environment does not expose these constants.
F_ADD_SEALS = getattr(fcntl, 'F_ADD_SEALS', 1033)
F_GET_SEALS = getattr(fcntl, 'F_GET_SEALS', 1034)
IMMUTABLE_SEALS = 0x0001 | 0x0002 | 0x0004 | 0x0008

class Stop(RuntimeError):
    pass

def require(ok, code):
    if not ok:
        raise Stop(code)

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def utc():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def redact_text(text, token):
    if not token:
        return text
    escaped = json.dumps(token, ensure_ascii=False)[1:-1]
    return text.replace(token, '[REDACTED]').replace(escaped, '[REDACTED]')

def redact_value(value, token):
    if type(value) is str:
        return redact_text(value, token)
    if type(value) is list:
        return [redact_value(x, token) for x in value]
    if type(value) is dict:
        return {redact_text(k, token): redact_value(v, token) for k, v in value.items()}
    return value

def redact_bytes(raw, token):
    escaped = json.dumps(token, ensure_ascii=False)[1:-1].encode()
    return raw.replace(token.encode(), b'[REDACTED]').replace(escaped, b'[REDACTED]')

def write_json(path, obj, token=''):
    raw = (json.dumps(redact_value(obj, token), indent=2, sort_keys=True) + '\n').encode()
    fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, 'wb') as f:
        f.write(raw)
        f.flush()
        os.fsync(f.fileno())
    d = os.open(Path(path).parent, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(d)
    finally:
        os.close(d)

def read_pinned(path, pin):
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    with os.fdopen(fd, 'rb') as f:
        s = os.fstat(f.fileno())
        require(stat.S_ISREG(s.st_mode) and s.st_nlink == 1, 'PINNED_FILE_NOT_SINGLE_REGULAR_LEAF')
        raw = f.read(MAX_JSON + 1)
        require(len(raw) <= MAX_JSON and sha(raw) == pin, 'PINNED_FILE_HASH_DIFFERS')
        return raw

def stat_identity(s):
    return (s.st_dev, s.st_ino, s.st_mode, s.st_nlink, s.st_size, s.st_mtime_ns, s.st_ctime_ns)

def capture_asset(path, pin):
    require(sys.platform == 'linux' and hasattr(os, 'memfd_create'), 'LINUX_SEALED_MEMFD_RUNTIME_REQUIRED')
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    captured = os.memfd_create('HDBLAST-pinned-identical-upload', os.MFD_ALLOW_SEALING)
    try:
        before = os.fstat(fd)
        require(stat.S_ISREG(before.st_mode) and before.st_nlink == 1 and before.st_size == pin['bytes'], 'ASSET_STAT_DIFFERS')
        digest = hashlib.sha256()
        md5 = hashlib.md5()
        count = 0
        while chunk := os.read(fd, 1024 * 1024):
            count += len(chunk)
            require(count <= pin['bytes'], 'ASSET_SIZE_EXCEEDED')
            digest.update(chunk)
            md5.update(chunk)
            view = memoryview(chunk)
            while view:
                written = os.write(captured, view)
                require(written > 0, 'CAPTURE_WRITE_FAILED')
                view = view[written:]
        require(stat_identity(before) == stat_identity(os.fstat(fd)), 'ASSET_CHANGED_DURING_CAPTURE')
        require(stat_identity(before) == stat_identity(os.stat(path, follow_symlinks=False)), 'ASSET_PATH_CHANGED_DURING_CAPTURE')
        require(count == pin['bytes'] and digest.hexdigest() == pin['sha256'] and md5.hexdigest() == pin['md5'], 'ASSET_COMPLETE_PINS_DIFFER')
        seals = IMMUTABLE_SEALS
        fcntl.fcntl(captured, F_ADD_SEALS, seals)
        require(fcntl.fcntl(captured, F_GET_SEALS) == seals, 'CAPTURE_NOT_IMMUTABLE')
        os.lseek(captured, 0, os.SEEK_SET)
        return captured, {'bytes': count, 'sha256': digest.hexdigest(), 'md5': md5.hexdigest(), 'stat_identity': list(stat_identity(before)), 'memfd_seals': seals}
    except Exception:
        os.close(captured)
        raise
    finally:
        os.close(fd)

def validate_state(state):
    require(type(state) is dict and state.get('draft_id') == DRAFT and state.get('published') is False, 'ORIGINAL_DRAFT_STATE_DIFFERS')
    require(state.get('continuation_binding', {}).get('controller_sha256') == CONTROLLER_SHA, 'ORIGINAL_CONTROLLER_BINDING_DIFFERS')
    require(state.get('inventory_sha256') == INVENTORY_SHA and state.get('metadata_sha256') == METADATA_SHA, 'ORIGINAL_INPUT_BINDING_DIFFERS')
    p = state.get('pending')
    require(type(p) is dict and set(p) == {'kind', 'url', 'body_sha256', 'body_bytes', 'context', 'utc'}, 'ORIGINAL_PENDING_SHAPE_DIFFERS')
    require(p['kind'] == 'put_content' and p['url'] == URL and p['body_sha256'] == ASSET_PIN['sha256'] and type(p['body_bytes']) is int and p['body_bytes'] == ASSET_PIN['bytes'] and p['context'] == {'filename': NAME}, 'ORIGINAL_PENDING_EXACT_INTENT_DIFFERS')
    require(state.get('create_post_attempted') is True and ('publish_post_attempted' not in state or state['publish_post_attempted'] is False) and NAME not in state.get('content_verified_files', []), 'ORIGINAL_EXECUTION_STAGE_DIFFERS')

def validate_remote(c, files, initialized):
    require(set(files) == {p['filename'] for p in c.old} | {NAME}, 'FRESH_DRAFT_FILE_MEMBERSHIP_DIFFERS')
    c.module.assert_pins(files, c.old, exact=False)
    for pin in c.old:
        fresh = files[pin['filename']]
        original = initialized[pin['filename']]
        require(c.module.immutable_file_identity(fresh) == c.module.immutable_file_identity(original), 'INHERITED_DRAFT_FILE_IDENTITIES_DIFFERS')
        require(fresh.get('bucket_id') == BUCKET_ID, 'INHERITED_DRAFT_BUCKET_DIFFERS')
    f = files[NAME]
    require(c.module.immutable_file_identity(f) == c.module.immutable_file_identity(initialized[NAME]), 'INITIALIZED_NEW_FILE_IDENTITIES_DIFFERS')
    require(f.get('key') == NAME and f.get('status') == 'pending' and f.get('size') is None and f.get('checksum') is None, 'NEW_FILE_NOT_SAME_UNAVAILABLE_PENDING_ENTRY')
    require(f.get('file_id') == FILE_ID and f.get('version_id') == VERSION_ID and f.get('bucket_id') == BUCKET_ID, 'NEW_FILE_IMMUTABLE_IDENTITIES_DIFFERS')
    require(f.get('links', {}).get('content') == URL and f.get('links', {}).get('commit') == URL.removesuffix('/content') + '/commit', 'NEW_FILE_LINKS_DIFFERS')

def validate_missing(response):
    require(response.status == 400 and not response.body_truncated and not response.body_read_failed, 'MISSING_CONTENT_GET_NOT_DEFINITIVE_400')
    expected = {'status': 400, 'message': 'File with key "' + NAME + '" is not available.'}
    require(type(response.json()) is dict and response.json() == expected, 'MISSING_CONTENT_ERROR_DIFFERS')

def load_controller():
    raw = read_pinned(CONTROLLER, CONTROLLER_SHA)
    module_name = '_reviewed_hdblast_identical_put_controller'
    spec = importlib.util.spec_from_loader(module_name, loader=None, origin=str(CONTROLLER))
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(CONTROLLER)
    sys.modules[module_name] = module
    exec(compile(raw, str(CONTROLLER), 'exec'), module.__dict__)
    return module

def controller_view(module, token, side):
    class View(module.Controller):
        def save(self):
            pass
        def event(self, kind, **values):
            row = module.redact({'utc': utc(), 'kind': kind, **values}, token)
            with (side / 'GUARD_JOURNAL.jsonl').open('a') as out:
                out.write(json.dumps(row) + '\n')
                out.flush()
                os.fsync(out.fileno())
        def capture(self, label, obj):
            write_json(side / (label + '.json'), obj, token)
        def write(self, *args, **kwargs):
            raise Stop('ORIGINAL_CONTROLLER_WRITES_FORBIDDEN')
    inventory = json.loads(read_pinned(INVENTORY, INVENTORY_SHA))
    metadata = json.loads(read_pinned(METADATA, METADATA_SHA))
    c = View(EVIDENCE, module.Transport(token), token, inventory, metadata, INVENTORY_SHA, METADATA_SHA)
    c.module = module
    return c

def durable_latch(side, intent):
    write_json(side / 'ONE_PUT_ATTEMPT_LATCH.json', intent)

def curl_command(asset_fd, header_fd):
    # Empty Expect avoids libcurl's implicit HTTP417 resending. No -L or retry.
    return [CURL, '-q', '--config', '-', '--http1.1', '--request', 'PUT',
        '--upload-file', '/proc/self/fd/' + str(asset_fd), '--header', 'Expect:',
        '--header', 'Accept: application/json', '--header', 'Content-Type: application/octet-stream',
        '--header', 'Content-Length: ' + str(ASSET_PIN['bytes']), '--proto', '=https',
        '--proto-redir', '=https', '--max-redirs', '0', '--retry', '0',
        '--connect-timeout', '60', '--max-time', '180', '--max-filesize', str(MAX_JSON),
        '--dump-header', '/proc/self/fd/' + str(header_fd), '--output', '-',
        '--silent', '--show-error', '--write-out',
        '%{stderr}\nHDBLAST_CURL_METRICS|%{http_code}|%{size_upload}|%{size_download}|%{time_total}|%{ssl_verify_result}|%{num_redirects}\n', URL]

def token_config(token):
    require(bool(token) and all(32 <= ord(x) < 127 for x in token), 'TOKEN_NOT_SINGLE_ASCII_HEADER_VALUE')
    return ('header = "Authorization: Bearer ' + token.replace('\\', '\\\\').replace('"', '\\"') + '"\n').encode()

def parse_metrics(stderr):
    require(stderr.count(b'HDBLAST_CURL_METRICS|') == 1, 'CURL_METRICS_MISSING_OR_AMBIGUOUS')
    match = re.search(rb'\nHDBLAST_CURL_METRICS\|([0-9]{3})\|([0-9]+)\|([0-9]+)\|([0-9]+(?:\.[0-9]+)?)\|([0-9]+)\|([0-9]+)\n\Z', stderr)
    require(match is not None, 'CURL_FINAL_METRICS_UNAVAILABLE')
    status, uploaded, downloaded, duration, tls, redirects = match.groups()
    require(int(status) == 0 or 100 <= int(status) <= 599, 'CURL_HTTP_STATUS_INVALID')
    return {'http_status': int(status), 'uploaded_bytes': int(uploaded), 'downloaded_bytes': int(downloaded), 'seconds': duration.decode(), 'ssl_verify_result': int(tls), 'redirects': int(redirects)}

def run_curl(asset_fd, token):
    read_fd, write_fd = os.pipe()
    process = None
    selector = selectors.DefaultSelector()
    buffers = {'body': bytearray(), 'stderr': bytearray(), 'headers': bytearray()}
    outcome = {'response_bound_exceeded': False, 'wall_deadline_exceeded': False}
    try:
        command = curl_command(asset_fd, write_fd)
        process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, pass_fds=(asset_fd, write_fd), close_fds=True)
        os.close(write_fd)
        write_fd = None
        process.stdin.write(token_config(token))
        process.stdin.close()
        for label, fd in [('body', process.stdout.fileno()), ('stderr', process.stderr.fileno()), ('headers', read_fd)]:
            os.set_blocking(fd, False)
            selector.register(fd, selectors.EVENT_READ, label)
        deadline = time.monotonic() + 190
        killed_at = None
        while selector.get_map():
            now = time.monotonic()
            if now > deadline and killed_at is None:
                outcome['wall_deadline_exceeded'] = True
                if process.poll() is None:
                    process.kill()
                killed_at = now
            if killed_at is not None and now > killed_at + 5:
                outcome['pipe_drain_deadline_exceeded'] = True
                break
            for key, _ in selector.select(0.25):
                raw = os.read(key.fd, 65536)
                if not raw:
                    selector.unregister(key.fd)
                    continue
                buf = buffers[key.data]
                room = MAX_JSON - len(buf)
                if len(raw) > room:
                    outcome['response_bound_exceeded'] = True
                    if killed_at is None:
                        if process.poll() is None:
                            process.kill()
                        killed_at = time.monotonic()
                buf.extend(raw[:max(room, 0)])
        outcome['curl_exit'] = process.wait(timeout=5)
        try:
            outcome['metrics'] = parse_metrics(bytes(buffers['stderr']))
        except Stop as stopped:
            outcome['metrics_unavailable'] = str(stopped)
        return outcome, {k: redact_bytes(bytes(v), token)[:MAX_JSON] for k, v in buffers.items()}
    finally:
        if process is not None and process.poll() is None:
            process.kill()
            process.wait(timeout=5)
        selector.close()
        os.close(read_fd)
        if write_fd is not None:
            os.close(write_fd)

def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--helper-sha256', required=True)
    parser.add_argument('--attempt-identical-put', action='store_true')
    args = parser.parse_args(argv)
    require(sha(Path(__file__).read_bytes()) == args.helper_sha256, 'EXTERNAL_HELPER_SOURCE_PIN_DIFFERS')
    require(args.attempt_identical_put, 'EXPLICIT_IDENTICAL_PUT_FLAG_REQUIRED')
    token = os.environ.get('ZENODO_ACCESS_TOKEN', '')
    token_config(token)
    module = load_controller()
    require(not SIDE.exists(), 'FIXED_RECOVERY_DIRECTORY_ALREADY_EXISTS_ROOT_REVIEW_REQUIRED')
    SIDE.mkdir(mode=0o700)
    parent_fd = os.open(SIDE.parent, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(parent_fd)
    finally:
        os.close(parent_fd)
    before = {name: sha((EVIDENCE / name).read_bytes()) for name in ('STATE.json', 'JOURNAL.jsonl')}
    captured = None
    try:
        with module.EvidenceLock(EVIDENCE):
            require(before == {'STATE.json': ORIGINAL_STATE_SHA, 'JOURNAL.jsonl': ORIGINAL_JOURNAL_SHA}, 'ORIGINAL_EXECUTION_EVIDENCE_PINS_DIFFERS')
            validate_state(json.loads((EVIDENCE / 'STATE.json').read_bytes()))
            c = controller_view(module, token, SIDE)
            validate_state(c.state)
            asset = next(p for p in c.add if p['filename'] == NAME)
            require({k: asset[k] for k in ASSET_PIN} == ASSET_PIN, 'INVENTORY_ASSET_PINS_DIFFERS')
            captured, capture_receipt = capture_asset(asset['local_path'], ASSET_PIN)
            c.current()
            require(c.state.get('published') is False, 'PUBLICATION_ALREADY_COMPLETED')
            owner = c.owner_list()
            require(owner and str(owner['id']) == DRAFT, 'FRESH_OWNED_DRAFT_DIFFERS')
            draft, headers, files = c.draft(DRAFT)
            initialized = module.file_map(json.loads(read_pinned(INITIALIZATION, INITIALIZATION_SHA)))
            validate_remote(c, files, initialized)
            response = c.transport.request('GET', URL, headers={'Accept': '*/*', 'Accept-Encoding': 'identity'})
            write_json(SIDE / 'FRESH_MISSING_CONTENT.json', {'status': response.status, 'body_sha256': sha(response.body), 'body': response.body.decode(errors='replace')}, token)
            validate_missing(response)
            require(sha(CONTROLLER.read_bytes()) == CONTROLLER_SHA, 'CONTROLLER_CHANGED_BEFORE_UPLOAD')
            require(before == {name: sha((EVIDENCE / name).read_bytes()) for name in before}, 'ORIGINAL_STATE_OR_JOURNAL_CHANGED_BEFORE_UPLOAD')
            intent = {'utc': utc(), 'operation': 'ONE_IDENTICAL_PUT_RECOVERY', 'url': URL, 'asset': ASSET_PIN,
                'captured_asset': capture_receipt, 'controller_sha256': CONTROLLER_SHA,
                'original_state_journal': before, 'proxy_environment_preserved': True,
                'original_pending_remains': True, 'reinitializations': 0, 'commits': 0, 'publishes': 0,
                'expect_header': 'EMPTY_TO_FORBID_IMPLICIT_417_RESEND'}
            durable_latch(SIDE, intent)
            result, outputs = run_curl(captured, token)
            for label, raw in outputs.items():
                fd = os.open(SIDE / ('curl_' + label + '.raw'), os.O_CREAT | os.O_EXCL | os.O_WRONLY | os.O_NOFOLLOW, 0o600)
                with os.fdopen(fd, 'wb') as out:
                    out.write(raw)
                    out.flush()
                    os.fsync(out.fileno())
            require(before == {name: sha((EVIDENCE / name).read_bytes()) for name in before}, 'ORIGINAL_STATE_OR_JOURNAL_CHANGED_DURING_UPLOAD')
            receipt = {'utc': utc(), 'status': 'ONE_IDENTICAL_PUT_ATTEMPT_RECORDED_RECONCILE_REQUIRED', 'result': result,
                'outputs': {k: {'bytes': len(v), 'sha256': sha(v)} for k, v in outputs.items()},
                'original_state_journal_unchanged': True, 'pending_cleared': False, 'commits': 0, 'publishes': 0}
            write_json(SIDE / 'RECOVERY_RESULT.json', receipt, token)
            print(json.dumps(receipt))
    except Exception as exception:
        write_json(SIDE / 'STOPPED.json', {'utc': utc(), 'status': 'STOPPED_OR_UNKNOWN_NO_REPEAT',
            'reason': str(exception).replace(token, '[REDACTED]')[:1024], 'exception_type': type(exception).__name__,
            'put_latch_exists': (SIDE / 'ONE_PUT_ATTEMPT_LATCH.json').exists(), 'pending_cleared': False}, token)
        raise
    finally:
        if captured is not None:
            os.close(captured)
    return 0

if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except Exception as exception:
        print(json.dumps({'status': 'STOPPED_OR_UNKNOWN_NO_REPEAT', 'reason': str(exception) if isinstance(exception, Stop) else 'UNEXPECTED_FAILURE_RECORDED_IN_SIDECAR'}))
        raise SystemExit(2)
