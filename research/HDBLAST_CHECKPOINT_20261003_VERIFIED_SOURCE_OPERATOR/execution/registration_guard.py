"""Stdlib pre-import integrity and prospective-source chronology guard.

The public API readback is performed separately by the registering custodian.
This offline guard verifies the complete archived receipt and all source bytes;
it does not cryptographically authenticate an arbitrary supplied receipt.
"""
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import re

SOURCES = ('positive_B', 'signed_uB')
REQUIRED = {'REGISTRATION_CONTRACT.json', 'PROTOCOL.md',
            'execution/registration_guard.py', 'execution/run_registered.py',
            'execution/execute_bounded.py', 'primary/verified_moments.py',
            'execution/review_completed.py',
            'primary/route.py', 'independent/generic_operator.py',
            'independent/source_models.py', 'independent/route.py',
            'protocol/operator_probe_contract.py',
            'evidence/preparation/FINAL_PRE_FREEZE_REVIEW.json'}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def zero_counts(value):
    return (type(value) is dict and set(value) == {'primary', 'independent'} and
            all(type(v) is int and v == 0 for v in value.values()))


def load(raw):
    def unique(pairs):
        out = {}
        for k, v in pairs:
            require(k not in out, 'Duplicate JSON key')
            out[k] = v
        return out
    def invalid(value):
        raise ValueError('Nonfinite JSON constant')
    return json.loads(raw, object_pairs_hook=unique, parse_constant=invalid)


def safe_file(root, name):
    require(type(name) is str and name and '\\' not in name,
            'Canonical relative filename required')
    parts = name.split('/')
    require(not any(p in ('', '.', '..') for p in parts),
            'Relative traversal forbidden')
    path = root
    for part in parts:
        path = path / part
        require(not path.is_symlink(), 'Symlink forbidden')
    require(path.is_file(), 'Registered file absent: ' + name)
    require(path.stat().st_size <= 32 * 1024 * 1024,
            'Unexpected oversized registration member')
    return path


def authenticate(root, receipt_path, receipt_sha256, registration_sha256):
    root = Path(root).absolute()
    require(root.is_dir() and root.resolve() == root, 'Real checkpoint root required')
    require(re.fullmatch('[0-9a-f]{64}', receipt_sha256 or '') and
            re.fullmatch('[0-9a-f]{64}', registration_sha256 or ''),
            'Explicit exact receipt and registration pins required')
    raw = Path(receipt_path).read_bytes()
    require(sha(raw) == receipt_sha256, 'Readback receipt hash differs')
    receipt = load(raw)
    require(receipt['status'] == 'PASS_REMOTE_REGISTERED_SOURCE_GO' and
            receipt['public_repository'] == 'https://github.com/maldonado-research/HDblast' and
            receipt['input_scope'] == 'NO_RETAINED_ARRAYS_SOURCE_OPERATOR_ONLY' and
            receipt['all_registered_public_blobs_independently_downloaded'] is True and
            zero_counts(receipt['new_target_evaluations_before_public_verification']),
            'Unverified prospective receipt')
    require(re.fullmatch('[0-9a-f]{40}', receipt['freeze_commit']) and
            receipt['registration_sha256'] == registration_sha256,
            'Receipt freeze/registration pins differ')
    regraw = safe_file(root, 'FULL_REGISTRATION.json').read_bytes()
    require(sha(regraw) == registration_sha256, 'Registration hash differs')
    registration = load(regraw)
    require(registration['scientific_outcome_before_freeze'] == 'UNCOMPUTED_SOURCE_OPERATOR' and
            zero_counts(registration['new_target_evaluations_before_freeze']),
            'Nonprospective registration')
    files = registration['files']
    require(type(files) is dict and REQUIRED <= set(files), 'Incomplete source inventory')
    require('FULL_REGISTRATION.json' not in files, 'Circular registration inventory')
    remote = receipt['remote_registered_files']
    require(set(remote) == set(files) | {'FULL_REGISTRATION.json'} and
            receipt['remote_files_verified'] == len(remote), 'Incomplete remote readback')
    for name, item in files.items():
        require(type(item) is dict and set(item) == {'sha256', 'bytes'},
                'Canonical file inventory required')
        data = safe_file(root, name).read_bytes()
        require(sha(data) == item['sha256'] and len(data) == item['bytes'],
                'Frozen member changed: ' + name)
        verify_remote(remote[name], data)
    verify_remote(remote['FULL_REGISTRATION.json'], regraw)
    # Prevent unregistered Python modules from shadowing frozen imports.
    for path in root.rglob('*.py'):
        require(not path.is_symlink() and path.relative_to(root).as_posix() in files,
                'Unregistered Python source')
    for path in root.rglob('*'):
        require(not path.is_symlink(), 'Checkpoint symlink forbidden')
        if path.is_file():
            require(path.suffix not in ('.pyc', '.pyo', '.so', '.pyd', '.dll') and
                    '__pycache__' not in path.parts,
                    'Cached bytecode/native module forbidden in checkpoint')
    contract = load(safe_file(root, 'REGISTRATION_CONTRACT.json').read_bytes())
    require(contract == registration['frozen_configuration'], 'Frozen configuration differs')
    require(contract['freeze_authorization'] == 'NOT_AUTHORIZED_BEFORE_REMOTE_BYTE_VERIFICATION' and
            contract['include_Lg'] is True and contract['panel_count'] == 64 and
            contract['probe_momenta'] == ['0/1','1/1099511627776','1/4096','1/4',
                                         '1/1','16/1','64/1','128/1','256/1'] and
            contract['sources'] == list(SOURCES) and
            contract['archive_array_decoding_allowed'] is False,
            'Registered source universe differs')
    review = load(safe_file(root, 'evidence/preparation/FINAL_PRE_FREEZE_REVIEW.json').read_bytes())
    require(review['status'] == 'GO_FOR_PROSPECTIVE_PUBLIC_FREEZE',
            'Independent prefreeze GO missing')
    require((REQUIRED-{'evidence/preparation/FINAL_PRE_FREEZE_REVIEW.json'}) <=
            set(review['reviewed_source_sha256']), 'Incomplete reviewed source coverage')
    for name, pin in review['reviewed_source_sha256'].items():
        require(name in files and files[name]['sha256'] == pin, 'Reviewed source changed')
    return {'receipt': receipt, 'registration': registration, 'contract': contract,
            'receipt_sha256': receipt_sha256, 'registration_sha256': registration_sha256}


def verify_remote(item, data):
    require(item['sha256'] == sha(data) and item['bytes'] == len(data) and
            item['independent_public_blob_download_sha256_verified'] is True and
            item['git_blob'] == hashlib.sha1(
                b'blob ' + str(len(data)).encode() + b'\0' + data).hexdigest(),
            'Remote blob binding differs')


class SourceAuthorization:
    def __init__(self, verified, route, journal=None):
        require(route in ('primary', 'independent'), 'Known route required')
        self.verified = verified
        self.route = route
        self.events = []
        self.seen = set()
        self.journal = Path(journal) if journal is not None else None
        if self.journal is not None:
            self.journal.open('x').close()
        receipt = verified['receipt']
        self.authorization = {k: receipt[k] for k in
            ('status', 'freeze_commit', 'registration_sha256', 'input_scope')}

    def __call__(self, event, detail):
        require(event == 'real_source_taylor' and set(detail) == {'source', 'center'},
                'Unexpected source callback')
        source, center = detail['source'], Fraction(detail['center'])
        centers = {Fraction(-9, 2)+Fraction(2*j+1, 128) for j in range(64)}
        require(source in SOURCES and center in centers,
                'Source callback outside frozen universe')
        key = (source, center)
        require(key not in self.seen, 'Repeated real-source construction')
        record = {'source': source, 'center': str(center)}
        if self.journal is not None:
            # Persist the permitted attempt before the actual source callback,
            # so failure/timeout cannot erase known-outcome chronology.
            with self.journal.open('a') as handle:
                handle.write(json.dumps({'index':len(self.events)+1,
                    'route':self.route,**record},sort_keys=True)+'\n')
                handle.flush()
                os.fsync(handle.fileno())
        self.seen.add(key)
        self.events.append(record)

    def complete(self):
        require(len(self.seen) == 128, 'Incomplete real-source construction')
        return {'route': self.route, 'source_bundle_constructions': 128,
                'archive_arrays_decoded': 0, 'events': self.events}
