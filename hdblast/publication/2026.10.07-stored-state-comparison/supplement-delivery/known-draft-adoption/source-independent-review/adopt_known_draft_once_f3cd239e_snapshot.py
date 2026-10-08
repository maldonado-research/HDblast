#!/usr/bin/env python3
"""One local migration of a proved existing draft into a fresh publisher state.

No networking or controller import. The unchanged controller must perform its
fresh remote guards after this local adoption. An interrupted intent is never
automatically repeated. The original pending state and journal stay untouched.
"""
import argparse
import copy
import fcntl
import hashlib
import json
import os
from pathlib import Path
import stat
import time
from contextlib import ExitStack, contextmanager

ROOT = Path('/workspace/hdblast-research-work/continuation-state-comparison-20261007')
HERE = ROOT / 'publication-final/known-draft-adoption'
PINS_SHA256 = '4c3331184089e7eb4cc546cec3c5a0e61d24ac91670ab3cfe47e842ac4071b74'
DRAFT_ID = '23228395'
CONTROLLER_SHA256 = 'aa7e6cfe271559731ed3b45c5339f8eb612ebe6f017ddafa5fdadad9579e084b'
EXPECTED_ROLES = {
    'old_state', 'old_journal', 'new_state', 'new_journal', 'inventory',
    'metadata', 'creation', 'controller', 'cleanup_source',
    'cleanup_acceptance', 'cleanup_workflow_review',
    'cleanup_1_CURRENT_PRIOR.json', 'cleanup_2_CURRENT_LATEST.json',
    'cleanup_3_CURRENT_PRIOR_FILES.json', 'cleanup_4_owner_page_1.json',
    'cleanup_5_DRAFT_23228395.json', 'cleanup_6_DRAFT_FILES_23228395.json',
    'cleanup_GUARD_JOURNAL.jsonl', 'cleanup_RESULT.json',
}


class Refused(RuntimeError):
    pass


def require(condition, reason):
    if not condition:
        raise Refused(reason)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def regular(path):
    path = Path(path).absolute()
    for component in (path, *path.parents):
        require(not component.is_symlink(), 'SYMLINK_PATH')
    info = path.stat()
    require(stat.S_ISREG(info.st_mode), 'NOT_REGULAR_FILE')
    return path, info


def capture(spec):
    require(type(spec) is dict and set(spec) == EXPECTED_ROLES, 'ROLE_SET_DIFFERS')
    raw, seen = {}, set()
    for role, entry in spec.items():
        require(type(entry) is dict and set(entry) == {'path', 'sha256', 'bytes'},
                'PIN_ENTRY_DIFFERS')
        require(type(entry['bytes']) is int and entry['bytes'] >= 0, 'PIN_SIZE_TYPE')
        path, before = regular(entry['path'])
        identity = (before.st_dev, before.st_ino)
        require(identity not in seen, 'INPUT_ALIAS')
        seen.add(identity)
        value = path.read_bytes()
        after = path.stat()
        require((before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) ==
                (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns),
                'INPUT_CHANGED_DURING_CAPTURE')
        require(len(value) == entry['bytes'] and digest(value) == entry['sha256'],
                'PIN_MISMATCH:' + role)
        raw[role] = value
    return raw


def validate(spec):
    raw = capture(spec)
    old = json.loads(raw['old_state'])
    new = json.loads(raw['new_state'])
    require(type(old) is dict and type(new) is dict, 'STATE_TYPE')
    require(old.get('draft_id') == DRAFT_ID and old.get('create_post_attempted') is True
            and old.get('published') is False and type(old.get('pending')) is dict
            and old['pending'].get('kind') == 'put_content', 'ORIGINAL_CREATION_STATE_DIFFERS')
    require(set(new) == {'pending', 'draft_id', 'published', 'continuation_binding',
                         'inventory_sha256', 'metadata_sha256'}, 'NEW_STATE_FIELDS_DIFFERS')
    require(new['pending'] is None and new['draft_id'] is None and
            new['published'] is False, 'NEW_STATE_NOT_FRESH')
    require(type(new['continuation_binding']) is dict and
            new['continuation_binding'] == old.get('continuation_binding') and
            new['continuation_binding'].get('controller_sha256') == CONTROLLER_SHA256
            and new['continuation_binding'].get('prior_record_id') == '23225288'
            and new['continuation_binding'].get('parent_id') == '17088132',
            'CONTINUATION_BINDING_DIFFERS')
    require(new['inventory_sha256'] == digest(raw['inventory']) and
            new['metadata_sha256'] == digest(raw['metadata']) and
            digest(raw['controller']) == CONTROLLER_SHA256, 'NEW_INPUT_BINDING_DIFFERS')
    creation = json.loads(raw['creation'])
    require(type(creation.get('id')) is int and creation['id'] == int(DRAFT_ID)
            and creation.get('conceptrecid') == '17088132'
            and type(creation.get('owner')) is int and creation['owner'] == 1386319
            and creation.get('state') == 'unsubmitted'
            and creation.get('submitted') is False, 'CREATION_RESPONSE_DIFFERS')
    cleanup = json.loads(raw['cleanup_RESULT.json'])
    require(cleanup == {'original_pending_cleared': False,
                        'repeat_delete_authorized': False,
                        'status': 'PASS_GET_ONLY_EXACT_TARGET_ABSENT_27_INHERITED_PRESERVED',
                        'writes': 0} and type(cleanup['writes']) is int,
            'CLEANUP_ABSENCE_PROOF_DIFFERS')
    require(raw['new_journal'].endswith(b'\n'), 'NEW_JOURNAL_NOT_COMPLETE')
    return raw, new


def fsync_directory(path):
    fd = os.open(path, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def exclusive_json(path, value):
    data = (json.dumps(value, indent=2, sort_keys=True) + '\n').encode()
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, 0o600)
    with os.fdopen(fd, 'wb') as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())
    fsync_directory(Path(path).parent)
    return data


@contextmanager
def execution_lock(path):
    path, _ = regular(path)
    fd = os.open(path, os.O_RDWR | os.O_NOFOLLOW)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        yield
    finally:
        os.close(fd)


def append_event(path, value):
    data = (json.dumps(value, sort_keys=True) + '\n').encode()
    fd = os.open(path, os.O_WRONLY | os.O_APPEND | os.O_NOFOLLOW)
    with os.fdopen(fd, 'ab') as stream:
        stream.write(data)
        stream.flush()
        os.fsync(stream.fileno())


def adopt(spec, recovery, *, apply=False):
    """Fabricated controls may supply temporary-file specs; main fixes live pins."""
    require(type(apply) is bool, 'APPLY_TYPE')
    new_path = Path(spec['new_state']['path'])
    old_path = Path(spec['old_state']['path'])
    recovery = Path(recovery).absolute()
    require(not recovery.exists() and not recovery.is_symlink(), 'ADOPTION_ALREADY_ATTEMPTED')
    require(recovery.parent.exists() and not recovery.parent.is_symlink(), 'RECOVERY_PARENT_DIFFERS')
    require(old_path.parent != new_path.parent and
            recovery != old_path.parent and recovery != new_path.parent and
            old_path.parent not in recovery.parents and new_path.parent not in recovery.parents,
            'RECOVERY_OR_EXECUTION_ALIAS')
    with ExitStack() as stack:
        for directory in sorted({old_path.parent, new_path.parent}):
            stack.enter_context(execution_lock(directory / 'CONTROLLER.lock'))
        raw, new = validate(spec)
        if not apply:
            return {'status': 'PASS_LOCAL_ADOPTION_INPUTS_ONLY', 'writes': 0}
        recovery.mkdir(mode=0o700)
        fsync_directory(recovery.parent)
        intent = {'kind': 'LOCAL_KNOWN_DRAFT_ADOPTION_INTENT',
                  'utc_epoch_ns': time.time_ns(), 'draft_id': DRAFT_ID,
                  'input_pins': spec, 'network_writes': 0,
                  'original_pending_cleared': False,
                  'next_step': 'UNCHANGED_CONTROLLER_FRESH_REMOTE_GUARDS'}
        exclusive_json(recovery / 'ONE_LOCAL_ADOPTION_LATCH.json', intent)
        exclusive_json(recovery / 'NEW_STATE_BEFORE.json', new)
        # Recheck all inputs after the durable latch and before either local mutation.
        require(capture(spec) == raw, 'INPUT_CHANGED_AFTER_INTENT')
        append_event(spec['new_journal']['path'], intent)
        adopted = copy.deepcopy(new)
        adopted['draft_id'] = DRAFT_ID
        temp = new_path.parent / '.KNOWN_DRAFT_ADOPTION_STATE.tmp'
        encoded = exclusive_json(temp, adopted)
        os.replace(temp, new_path)
        fsync_directory(new_path.parent)
        for role in ('old_state', 'old_journal'):
            require(Path(spec[role]['path']).read_bytes() == raw[role], 'ORIGINAL_EVIDENCE_CHANGED')
        require(json.loads(new_path.read_bytes()) == adopted, 'ADOPTED_STATE_READBACK_DIFFERS')
        result = {'status': 'PASS_LOCAL_KNOWN_DRAFT_ADOPTED_FRESH_REMOTE_GUARDS_REQUIRED',
                  'draft_id': DRAFT_ID, 'new_state_sha256': digest(encoded),
                  'only_state_field_changed': 'draft_id', 'network_writes': 0,
                  'original_pending_cleared': False}
        append_event(spec['new_journal']['path'], {'kind': 'LOCAL_KNOWN_DRAFT_ADOPTED', **result})
        exclusive_json(recovery / 'RESULT.json', result)
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--helper-sha256', required=True)
    parser.add_argument('--adopt-known-draft', action='store_true')
    args = parser.parse_args()
    require(digest(Path(__file__).read_bytes()) == args.helper_sha256, 'HELPER_PIN_DIFFERS')
    pins_path, _ = regular(HERE / 'INPUT_PINS.json')
    pin_bytes = pins_path.read_bytes()
    require(digest(pin_bytes) == PINS_SHA256, 'INPUT_PINS_MANIFEST_DIFFERS')
    result = adopt(json.loads(pin_bytes), HERE / 'live-local-adoption-001',
                   apply=args.adopt_known_draft)
    print(json.dumps(result, sort_keys=True))


if __name__ == '__main__':
    main()
