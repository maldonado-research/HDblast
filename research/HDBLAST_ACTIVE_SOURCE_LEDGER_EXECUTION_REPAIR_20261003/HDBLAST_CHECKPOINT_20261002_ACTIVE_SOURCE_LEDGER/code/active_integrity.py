#!/usr/bin/env python3
"""Byte and path guards; this module never imports or interprets physical arrays."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import ast
import struct
import zipfile

CHECKPOINT_NAME = 'HDBLAST_CHECKPOINT_20261002_ACTIVE_SOURCE_LEDGER'
ZIP_NAME = CHECKPOINT_NAME + '.zip'
ZIP_SHA_NAME = ZIP_NAME + '.sha256'
EXCLUDED = frozenset(('MANIFEST.json', ZIP_NAME, ZIP_SHA_NAME))
VERSIONS = {'numpy': '2.2.6', 'scipy': '1.15.3', 'sympy': '1.14.0',
            'mpmath': '1.3.0', 'matplotlib': '3.10.1'}
TIME_LIMIT_SECONDS = 900
MEMORY_LIMIT_KIB = 262144
POSTFREEZE_DATA_SUFFIXES = frozenset(('.json', '.jsonl', '.log', '.txt', '.md', '.csv', '.tsv',
                                     '.npy', '.npz', '.png', '.svg', '.pdf'))
CAPSULE_NAMES = frozenset('metric_modes_' + source + '_' + setting + '.npz'
                         for source in ('positive_B', 'signed_uB')
                         for setting in ('coarse', 'fine'))
ACTIVE_MEMBERS = frozenset(name + '.npy' for name in (
    'k', 'momentum_weights', 'u_1', 'w_1', 'u_2', 'w_2', 'u_3', 'w_3',
    'observation_eta', 'history_eta', 'history_values', 'history_ledger_integrand',
    'history_baseline_contact', 'history_source_jet', 'history_forcing_jet',
    'history_geometry', 'history_cutoffs', 'history_quantity_names',
    'history_geometry_names'))


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


def json_bytes(raw):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, 'Duplicate JSON key: ' + key)
            result[key] = value
        return result
    def invalid(value):
        raise RuntimeError('Nonfinite JSON literal: ' + value)
    return json.loads(raw, object_pairs_hook=unique, parse_constant=invalid)


def read_json(path):
    return json_bytes(Path(path).read_bytes())


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


def npy_header(raw):
    """Read only NPY format metadata, never the array's numerical payload."""
    require(raw[:6] == b'\x93NUMPY' and len(raw) >= 10, 'Invalid NPY magic')
    version = tuple(raw[6:8])
    require(version in ((1, 0), (2, 0), (3, 0)), 'Unsupported NPY version')
    width = 2 if version == (1, 0) else 4
    require(len(raw) >= 8 + width, 'Truncated NPY header length')
    length = struct.unpack('<H' if width == 2 else '<I', raw[8:8+width])[0]
    end = 8 + width + length
    require(end <= len(raw), 'Truncated NPY header')
    header = ast.literal_eval(raw[8+width:end].decode('utf-8' if version == (3, 0) else 'latin1'))
    require(isinstance(header, dict) and set(header) == {'descr', 'fortran_order', 'shape'},
            'Unexpected NPY header fields')
    descriptor, shape = header['descr'], header['shape']
    require(isinstance(descriptor, str) and re.fullmatch(r'[<>=|][fciuUb][0-9]+', descriptor),
            'Nonprimitive/object NPY descriptor forbidden')
    require(isinstance(header['fortran_order'], bool) and isinstance(shape, tuple)
            and all(type(n) is int and n >= 0 for n in shape), 'Malformed NPY layout')
    elements = 1
    for dimension in shape:
        elements *= dimension
    itemsize = int(descriptor[2:]) * (4 if descriptor[1] == 'U' else 1)
    require(end + elements*itemsize == len(raw), 'NPY header/payload byte size mismatch')
    return {'npy_version': list(version), 'dtype_descriptor': descriptor,
            'shape': list(shape), 'fortran_order': header['fortran_order'],
            'header_bytes_including_magic': end, 'header_sha256': digest(raw[:end]),
            'payload_bytes': len(raw)-end, 'array_payload_interpreted': False}


def verify_input_capsules(root, registered):
    """Authenticate four complete nineteen-member capsules as opaque bytes."""
    name = 'inputs/INPUT_MANIFEST.json'
    require(name in registered, 'Input manifest omitted from registration')
    manifest = read_json(safe_file(root, name))
    require(manifest.get('schema_version') == 1, 'Unsupported input manifest schema')
    require(manifest.get('array_payload_values_loaded_interpreted_compared_or_evaluated') is False,
            'Preparation input manifest does not establish byte-only extraction')
    capsules = manifest.get('capsules')
    require(isinstance(capsules, list) and len(capsules) == 4, 'Exactly four input capsules required')
    names = [capsule.get('capsule_path') for capsule in capsules]
    require(len(names) == len(set(names)) and set(names) == CAPSULE_NAMES,
            'Input capsule case membership differs')
    for capsule in capsules:
        relative = 'inputs/' + capsule['capsule_path']
        path = safe_file(root, relative)
        require(relative in registered and registered[relative] == capsule.get('capsule_sha256')
                and file_sha(path) == capsule['capsule_sha256']
                and path.stat().st_size == capsule.get('capsule_bytes'), 'Capsule bytes/pin differ: ' + relative)
        require(capsule.get('original_members') == 631 and capsule.get('selected_members') == 19,
                'Input archive/member scope differs')
        members = capsule.get('members')
        require(isinstance(members, list) and len(members) == 19, 'Nineteen member receipts required')
        member_names = [member.get('name') for member in members]
        require(len(member_names) == len(set(member_names)) and set(member_names) == ACTIVE_MEMBERS,
                'Active input member membership differs')
        pins = {member['name']: member for member in members}
        with zipfile.ZipFile(path) as archive:
            info = archive.infolist()
            require(len(info) == 19 and len(set(item.filename for item in info)) == 19
                    and set(item.filename for item in info) == ACTIVE_MEMBERS,
                    'Capsule ZIP membership differs')
            for item in info:
                pin = pins[item.filename]
                require(item.compress_type == zipfile.ZIP_STORED and item.compress_type == pin.get('capsule_compression')
                        and item.file_size == pin.get('bytes') and format(item.CRC, '08x') == pin.get('crc32_hex'),
                        'Capsule member inventory differs: ' + item.filename)
                raw = archive.read(item)  # ZIP CRC validation and opaque byte hashing only.
                require(digest(raw) == pin.get('sha256'), 'Capsule member SHA256 differs: ' + item.filename)
                require(npy_header(raw) == pin.get('header'), 'Capsule NPY metadata differs: ' + item.filename)
    return {'capsules_verified': 4, 'members_verified': 76,
            'array_payload_values_interpreted': False, 'input_manifest_sha256': file_sha(safe_file(root, name))}


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
    capsule_receipt = (verify_input_capsules(root, files)
                       if 'inputs/INPUT_MANIFEST.json' in files else None)
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
            'package_manifest_sha256': package_pin, 'input_capsule_verification': capsule_receipt}
