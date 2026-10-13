#!/usr/bin/env python3
"""Reviewed GET-only recovery of one published edition; no publication methods."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import sys
import time
import types
from urllib import error

HERE = Path(__file__).resolve().parent
BASE_PATH = HERE / 'regular_shell_edition_controller.py'
BASE_SHA = '2c53e32a24fdc922d07f711cb67b8451da0be3325f740658d91d1edc0d97ca34'
INVENTORY_PATH = HERE / 'frozen-inputs-001/FINAL_UPLOAD_INVENTORY.json'
INVENTORY_SHA = '98f7db2f79c8fa29e0552a6b97cd923d2f25a2cd907d85d5893c370aabcc0454'
METADATA_PATH = HERE / 'frozen-inputs-001/FINAL_METADATA_MODERN.json'
METADATA_SHA = 'debbaf5a0fbe50d886d303940c348cfe4aa52d23041d519b68804c03d6bb7e80'
EVIDENCE = HERE / 'live-edition-001'
RECORD_ID = '23335485'
MAX_ATTEMPTS = 3
RETRYABLE = frozenset(('CONTENT_STREAM_FAILED', 'CONTENT_COMPLETE_MD5_SHA256_DIFFERS'))
SAFE_HEADERS = ('Content-Length', 'Content-Type', 'Content-Encoding', 'Content-Range', 'ETag')

def digest(data):
    return hashlib.sha256(data).hexdigest()

def load_base():
    raw = BASE_PATH.read_bytes()
    if digest(raw) != BASE_SHA:
        raise RuntimeError('BYTE_PINNED_BASE_CONTROLLER_DIFFERS')
    module = types.ModuleType('_hdblast_byte_pinned_readonly_base')
    module.__file__ = str(BASE_PATH)
    sys.modules[module.__name__] = module
    exec(compile(raw, str(BASE_PATH), 'exec'), module.__dict__)
    return module

BASE = load_base()

def public_get_guard(method, data, headers):
    BASE.require(method == 'GET' and data is None, 'READONLY_RECOVERY_NON_GET_OR_BODY_DENIED')
    BASE.require(headers.get('_public') is True, 'READONLY_RECOVERY_AUTHENTICATED_REQUEST_DENIED')
    BASE.require(not any(k.lower() == 'authorization' for k in headers), 'READONLY_RECOVERY_AUTHORIZATION_DENIED')

class ObservedOpener:
    def __init__(self, inner, owner):
        self.inner, self.owner = inner, owner
    def open(self, req, *args, **kwargs):
        BASE.require(req.get_method() == 'GET' and req.data is None, 'READONLY_RECOVERY_NON_GET_OR_BODY_DENIED')
        BASE.require(not any(k.lower() == 'authorization' for k, _ in req.header_items()), 'READONLY_RECOVERY_AUTHORIZATION_DENIED')
        BASE.safe_url(req.full_url, query=True)
        try:
            response = self.inner.open(req, *args, **kwargs)
        except error.HTTPError as exc:
            self.note(exc.code, exc.headers)
            exc.close()
            raise
        self.note(response.status, response.headers)
        return response
    def note(self, status, headers):
        if self.owner.observation is not None:
            self.owner.observation['http_status'] = status
            self.owner.observation['headers'] = {k: str(headers.get(k))[:500] if headers.get(k) is not None else None for k in SAFE_HEADERS}

class ReadOnlyTransport(BASE.Transport):
    def __init__(self, token):
        super().__init__(token)
        self.observation = None
        self.opener = ObservedOpener(self.opener, self)
    def request(self, method, url, data=None, headers=None):
        public_get_guard(method, data, headers or {})
        return super().request(method, url, data=data, headers=headers)
    def begin_attempt(self):
        self.observation = {'bytes': 0, 'ended_at_eof': False}
        self.observed_md5 = hashlib.md5()
        self.observed_sha256 = hashlib.sha256()
    def stream(self, url, headers=None):
        public_get_guard('GET', None, headers or {})
        for chunk in super().stream(url, headers=headers):
            self.observation['bytes'] += len(chunk)
            self.observed_md5.update(chunk)
            self.observed_sha256.update(chunk)
            yield chunk
        self.observation['ended_at_eof'] = True
    def observed(self):
        return {**self.observation, 'md5': self.observed_md5.hexdigest(), 'sha256': self.observed_sha256.hexdigest()}

class Sidecar:
    def __init__(self, path, runner_sha, token):
        self.path, self.runner_sha, self.token = path, runner_sha, token
        self.run_id = datetime.datetime.now(datetime.timezone.utc).isoformat()
    def event(self, kind, **values):
        row = {'utc': BASE.utc(), 'kind': kind, 'run_id': self.run_id,
               'runner_sha256': self.runner_sha, 'base_controller_sha256': BASE_SHA,
               'inventory_sha256': INVENTORY_SHA, 'metadata_sha256': METADATA_SHA,
               'record_id': RECORD_ID, **values}
        fd = os.open(self.path, os.O_WRONLY | os.O_APPEND | os.O_CREAT | os.O_NOFOLLOW, 0o600)
        with os.fdopen(fd, 'a', encoding='utf-8') as output:
            output.write(json.dumps(BASE.redact(row, self.token), ensure_ascii=False) + '\n')
            output.flush()
            os.fsync(output.fileno())

class ReadOnlyController(BASE.Controller):
    def denied(self, *args, **kwargs):
        raise BASE.Stop('READONLY_RECOVERY_MUTATION_OR_DRAFT_ENTRYPOINT_DENIED')
    write = create = prepare = publish = reconcile = owner_list = draft = require_draft_adoption = denied
    def stream_check(self, url, pin, public=False):
        BASE.require(public is True and self.state.get('published') is True
                     and self.state.get('draft_id') == RECORD_ID and self.state.get('pending') is None,
                     'READONLY_RECOVERY_PUBLISHED_STATE_REQUIRED')
        for attempt in range(1, MAX_ATTEMPTS + 1):
            self.transport.begin_attempt()
            self.recovery.event('STREAM_ATTEMPT_BEGIN', filename=pin['filename'], attempt=attempt)
            try:
                # Every call initializes the original count/MD5/SHA256 from zero.
                result = BASE.Controller.stream_check(self, url, pin, public=True)
            except BASE.Stop as exc:
                reason = str(exc)
                self.recovery.event('STREAM_ATTEMPT_FAILED', filename=pin['filename'], attempt=attempt,
                                    reason=reason, retryable=reason in RETRYABLE,
                                    observed=self.transport.observed())
                if reason not in RETRYABLE or attempt == MAX_ATTEMPTS:
                    raise
                time.sleep(1)
            else:
                self.recovery.event('STREAM_ATTEMPT_PASSED', filename=pin['filename'], attempt=attempt,
                                    observed=self.transport.observed())
                self.transport.observation = None
                return result

def published_state_guard(state):
    BASE.require(state.get('published') is True and state.get('draft_id') == RECORD_ID
                 and state.get('continuation_binding', {}).get('parent_id') == '17088132'
                 and state.get('pending') is None,
                 'READONLY_RECOVERY_EXISTING_PUBLISHED_STATE_REQUIRED')
    BASE.require(state.get('inventory_sha256') == INVENTORY_SHA and state.get('metadata_sha256') == METADATA_SHA
                 and state.get('continuation_binding', {}).get('controller_sha256') == BASE_SHA,
                 'READONLY_RECOVERY_STATE_PINS_DIFFER')

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runner-sha256', required=True)
    args = parser.parse_args(argv)
    runner_sha = digest(Path(__file__).read_bytes())
    BASE.require(args.runner_sha256 == runner_sha, 'READONLY_RECOVERY_REVIEWED_RUNNER_PIN_DIFFERS')
    BASE.require(EVIDENCE.is_dir() and not EVIDENCE.is_symlink() and (EVIDENCE / 'STATE.json').is_file(), 'READONLY_RECOVERY_EXISTING_EVIDENCE_REQUIRED')
    with BASE.EvidenceLock(EVIDENCE):
        published_state_guard(json.loads((EVIDENCE / 'STATE.json').read_text()))
        inventory = BASE.frozen_json(INVENTORY_PATH, INVENTORY_SHA)
        metadata = BASE.frozen_json(METADATA_PATH, METADATA_SHA)
        token = os.environ.get('ZENODO_ACCESS_TOKEN', '')
        controller = ReadOnlyController(EVIDENCE, ReadOnlyTransport(token), token,
                                        inventory, metadata, INVENTORY_SHA, METADATA_SHA)
        published_state_guard(controller.state)
        audit = Sidecar(EVIDENCE / 'READONLY_RECOVERY_JOURNAL.jsonl', runner_sha, token)
        controller.recovery = audit
        audit.event('READONLY_RUN_BEGIN', max_attempts_per_file=MAX_ATTEMPTS, retryable_reasons=sorted(RETRYABLE))
        try:
            # The original full 35-file verify loop, identity/metadata/file guards,
            # receipt capture and successful journal events remain unchanged.
            result = BASE.Controller.verify(controller, published=True)
            published_state_guard(controller.state)
        except BASE.Stop as exc:
            audit.event('READONLY_RUN_STOPPED', reason=str(exc), remote_mutations=0)
            raise
        except Exception as exc:
            audit.event('READONLY_RUN_STOPPED', exception_class=type(exc).__name__, remote_mutations=0)
            raise
        audit.event('READONLY_RUN_COMPLETE', verified_public_sha256=digest((EVIDENCE / 'VERIFIED_PUBLIC.json').read_bytes()),
                    file_count=result['file_count'], total_bytes=result['total_bytes'], remote_mutations=0)
        print(json.dumps({'status': 'COMPLETE_READONLY_PUBLIC_RECOVERY', 'record_id': RECORD_ID,
                          'file_count': result['file_count'], 'total_bytes': result['total_bytes'], 'remote_mutations': 0}))
    return 0

if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except BASE.Stop as exc:
        print(json.dumps({'status': 'STOPPED', 'reason': str(exc)}))
        raise SystemExit(2)
    except Exception:
        print(json.dumps({'status': 'STOPPED', 'reason': 'UNEXPECTED_FAILURE_DETAILS_NOT_OUTPUT'}))
        raise SystemExit(3)
