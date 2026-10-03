#!/usr/bin/env python3
"""Byte and path guards; this module never imports or interprets physical arrays."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
import re

CHECKPOINT_NAME = 'HDBLAST_CHECKPOINT_20261002_SOURCE_FREE_LEDGER'
ZIP_NAME = CHECKPOINT_NAME + '.zip'
ZIP_SHA_NAME = ZIP_NAME + '.sha256'
EXCLUDED = frozenset(('MANIFEST.json', ZIP_NAME, ZIP_SHA_NAME))
VERSIONS = {'numpy': '2.2.6', 'scipy': '1.15.3', 'sympy': '1.14.0',
            'mpmath': '1.3.0', 'matplotlib': '3.10.1'}
TIME_LIMIT_SECONDS = 900
MEMORY_LIMIT_KIB = 262144
POSTFREEZE_DATA_SUFFIXES = frozenset(('.json', '.jsonl', '.log', '.txt', '.md', '.csv', '.tsv',
                                     '.npy', '.npz', '.png', '.svg', '.pdf'))


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def valid_hex(value, length):
    return isinstance(value, str) and re.fullmatch('[0-9a-f]{%d}' % length, value) is not None


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def file_sha(path):
    hasher = hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            hasher.update(block)
    return hasher.hexdigest()


def safe_relative(name):
    require(isinstance(name, str) and bool(name) and '\\' not in name and ':' not in name,
            'Invalid relative path: ' + repr(name))
    parts = name.split('/')
    require(not PurePosixPath(name).is_absolute() and all(p not in ('', '.', '..') for p in parts),
            'Unsafe relative path: ' + name)
    return parts


def safe_file(root, name):
    root = Path(root).absolute()
    require(root.is_dir() and not root.is_symlink(), 'Missing or symlink checkpoint root')
    require(all(not parent.is_symlink() for parent in root.parents), 'Symlink checkpoint ancestor forbidden')
    target = root
    for component in safe_relative(name):
        target = target / component
        require(not target.is_symlink(), 'Symlink forbidden: ' + name)
    require(target.resolve().is_relative_to(root.resolve()) and target.is_file(),
            'Missing or escaped file: ' + name)
    return target


def fresh_external(root, output):
    root, output = Path(root).absolute(), Path(output).absolute()
    require(not output.exists() and not output.is_symlink()
            and all(not parent.is_symlink() for parent in output.parents)
            and not output.resolve().is_relative_to(root.resolve()),
            'Choose a fresh output directory outside checkpoint, without symlink ancestors')
    return output


def read_json(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, 'Duplicate JSON key: ' + key)
            result[key] = value
        return result
    def invalid(value):
        raise RuntimeError('Nonfinite JSON literal: ' + value)
    return json.loads(Path(path).read_bytes(), object_pairs_hook=unique, parse_constant=invalid)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')


def payload_files(root):
    root = Path(root).absolute()
    files = {}
    for path in sorted(root.rglob('*')):
        require(not path.is_symlink(), 'Symlink forbidden: ' + str(path))
        name = path.relative_to(root).as_posix()
        if '__pycache__' in path.parts or path.suffix in ('.pyc', '.pyo'):
            continue
        if path.is_file() and name not in EXCLUDED:
            files[name] = file_sha(safe_file(root, name))
    return files


def verify_package_manifest(root):
    path = Path(root) / 'MANIFEST.json'
    if not path.exists():
        return None
    manifest = read_json(safe_file(root, 'MANIFEST.json'))
    require(manifest.get('schema_version') == 1, 'Unsupported package manifest schema')
    require(manifest.get('excluded') == sorted(EXCLUDED), 'Package exclusion list differs')
    require(isinstance(manifest.get('files'), dict) and 'FREEZE_RECEIPT.json' in manifest['files'],
            'Package manifest must protect freeze receipt')
    require(payload_files(root) == manifest['files'], 'Package payload membership or bytes changed')
    return file_sha(path)


def verify_frozen(root, registration_sha256, freeze_commit):
    """Authenticate frozen bytes without opening NPY members or executing source."""
    root = Path(root).absolute()
    require(valid_hex(registration_sha256, 64), 'Expected explicit registration SHA256')
    require(valid_hex(freeze_commit, 40), 'Expected explicit full public freeze commit')
    registration_path = safe_file(root, 'FULL_REGISTRATION.json')
    require(file_sha(registration_path) == registration_sha256, 'Full registration SHA256 mismatch')
    receipt_path = safe_file(root, 'FREEZE_RECEIPT.json')
    receipt = read_json(receipt_path)
    require(receipt.get('public_freeze_commit') == freeze_commit, 'Freeze receipt commit mismatch')
    require(receipt.get('registration_sha256') == registration_sha256, 'Freeze receipt registration mismatch')
    registration = read_json(registration_path)
    require(registration.get('schema_version') == 1, 'Unsupported full registration schema')
    files = registration.get('files')
    require(isinstance(files, dict) and bool(files), 'Missing registered file table')
    require('EXPERIMENT.json' in files, 'Experiment omitted from registration')
    require('FULL_REGISTRATION.json' not in files and 'FREEZE_RECEIPT.json' not in files,
            'Self-referential registration is forbidden')
    for name, pin in files.items():
        require(name not in EXCLUDED and valid_hex(pin, 64), 'Invalid registered SHA/path: ' + str(name))
        require(file_sha(safe_file(root, name)) == pin, 'Frozen input changed: ' + name)
    experiment = read_json(safe_file(root, 'EXPERIMENT.json'))
    if 'experiment_sha256' in registration:
        require(registration['experiment_sha256'] == files['EXPERIMENT.json'], 'Experiment pin differs')
    if 'input_manifest_sha256' in registration:
        require(files.get('inputs/INPUT_MANIFEST.json') == registration['input_manifest_sha256'],
                'Input manifest pin differs')
    if 'frozen_configuration' in registration:
        require(registration['frozen_configuration'] == experiment, 'Frozen experiment differs')
    # Unregistered files must not shadow imports before a package manifest exists.
    actual = payload_files(root)
    directories = experiment.get('reproduction', {}).get('postfreeze_output_directories', [])
    require(isinstance(directories, list), 'Malformed postfreeze output namespaces')
    for directory in directories:
        safe_relative(directory)
        require(directory.split('/')[0] not in ('code', 'independent', 'inputs', 'theory', 'literature'),
                'Output namespace overlaps frozen source/input namespace')
    for name in set(actual) - set(files) - {'FULL_REGISTRATION.json', 'FREEZE_RECEIPT.json'}:
        require(any(name.startswith(directory + '/') for directory in directories)
                and Path(name).suffix in POSTFREEZE_DATA_SUFFIXES,
                'Unregistered payload file: ' + name)
    package_pin = verify_package_manifest(root)
    return {'public_freeze_commit': freeze_commit, 'registration_sha256': registration_sha256,
            'freeze_receipt_sha256': file_sha(receipt_path), 'frozen_files': files,
            'package_manifest_sha256': package_pin}
