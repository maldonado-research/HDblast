#!/usr/bin/env python3
"""Authenticate the whole portable package before importing any study code.

External SHA256 values for this helper, PAYLOAD_MANIFEST and REPLAY_PINS must
come from the independently sealed release, not from the extracted package.
No network, extraction, native-array loading or public-GO creation occurs here.
Original input dependencies are supplied separately at their immutable Git paths.
Manufactured controls require an explicit separate flag and never create GO.
"""
import sys
if not sys.flags.isolated or not sys.flags.dont_write_bytecode:
    raise SystemExit('Invoke the captured authenticated helper with Python -I -B')
sys.dont_write_bytecode = True
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import json
import math
import os
import re
import stat
import subprocess
import signal
import importlib.metadata

sys.dont_write_bytecode = True
PREFIX = 'research/HDBLAST_CHECKPOINT_20261008_RETAINED_ENDPOINT_FLOW'
RESERVED = {'replay_checkpoint.py', 'PAYLOAD_MANIFEST.json', 'REPLAY_PINS.json'}
SCIENCE = ('SOURCE_CERTIFICATE.json', 'TARGET_BUDGETS.json', 'NODE_CERTIFICATES.jsonl.gz',
           'DATA.json', 'SCIENCE_SUMMARY.json', 'SOURCE_ATTEMPTS.jsonl',
           'DECODE_ATTEMPTS.jsonl', 'INPUT_AUTHENTICATION.json')
STABLE = SCIENCE + ('OUTPUT_READBACK.json',)
CORE = PREFIX + '/core'
REGISTRATION = PREFIX + '/FULL_REGISTRATION.json'
DEPENDENCY_COMMIT = '13a30ef46f5c90c1b01ae83b609db65fd0f8a709'
ACTUAL_SCOPE = 'ACTUAL_RETAINED_ENDPOINT_FLOW'
CONTROL_SCOPE = 'MANUFACTURED_DELIVERY_CONTROL'
LIMITS = {'wall_seconds': 900, 'cpu_seconds': 900, 'address_space_bytes': 536870912,
          'aggregate_output_bytes': 134217728, 'file_size_bytes': 134217728,
          'nonroot': True, 'single_threads': True, 'nproc': 0, 'core': 0,
          'network_allowed': False, 'child_processes_allowed': False}
RUNTIME = {'python_version': '3.12.14', 'platform': 'linux',
           'packages': {'python-flint': '0.9.0', 'sympy': '1.14.0', 'mpmath': '1.3.0'}}
MAX_FILE = 256 * 1024 * 1024
MAX_PACKAGE = 1024 * 1024 * 1024
MAX_ENTRIES = 8192
MAX_DEPTH = 24

class Stop(ValueError):
    pass

def require(ok, reason):
    if not ok:
        raise Stop(reason)

def digest(value):
    require(type(value) is str and re.fullmatch('[0-9a-f]{64}', value) is not None,
            'EXTERNAL_SHA256_REQUIRED')
    return value

def load(raw):
    def unique(items):
        result = {}
        for key, value in items:
            require(key not in result, 'DUPLICATE_JSON_KEY')
            result[key] = value
        return result
    def bad(value):
        raise Stop('NONFINITE_JSON')
    result = json.loads(raw, object_pairs_hook=unique, parse_constant=bad)
    # json.loads also maps an overflowing exponent to infinity.
    json.dumps(result, allow_nan=False)
    return result

def relative(name):
    require(type(name) is str and name and '\\' not in name and '\0' not in name,
            'SAFE_RELATIVE_PATH_REQUIRED')
    path = PurePosixPath(name)
    require(not path.is_absolute() and path.as_posix() == name and
            all(part not in ('', '.', '..') for part in name.split('/')) and
            len(path.parts) <= MAX_DEPTH, 'UNSAFE_PACKAGE_PATH')
    return path

def directory(path):
    path = Path(os.path.abspath(os.fspath(path)))
    require(path.is_dir() and path.resolve() == path and
            all(not parent.is_symlink() for parent in (path, *path.parents)),
            'REAL_DIRECTORY_WITHOUT_SYMLINK_ANCESTORS_REQUIRED')
    return path

def identity(st):
    return (st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns, st.st_nlink)

def read_member(root_fd, name, expected=None, capture=False):
    """Read a regular single-link file through no-follow directory descriptors."""
    parts = relative(name).parts
    fd = os.dup(root_fd)
    try:
        for part in parts[:-1]:
            new = os.open(part, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC,
                          dir_fd=fd)
            os.close(fd); fd = new
        member = os.open(parts[-1], os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC,
                         dir_fd=fd)
        try:
            before = os.fstat(member)
            require(stat.S_ISREG(before.st_mode) and before.st_nlink == 1,
                    'REGULAR_UNALIASED_PACKAGE_FILE_REQUIRED:' + name)
            require(before.st_size <= MAX_FILE, 'PACKAGE_FILE_TOO_LARGE:' + name)
            if capture:
                require(before.st_size <= 20 * 1024 * 1024, 'PACKAGE_JSON_TOO_LARGE:' + name)
            total = 0; h = hashlib.sha256(); chunks = []
            while True:
                data = os.read(member, 1024 * 1024)
                if not data:
                    break
                total += len(data)
                require(total <= MAX_FILE, 'GROWING_PACKAGE_FILE:' + name)
                h.update(data)
                if capture:
                    chunks.append(data)
            after = os.fstat(member)
            named = os.stat(parts[-1], dir_fd=fd, follow_symlinks=False)
            require(identity(before) == identity(after) == identity(named) and total == before.st_size,
                    'PACKAGE_FILE_CHANGED_WHILE_HASHING:' + name)
            pin = {'bytes': total, 'sha256': h.hexdigest()}
            if expected is not None:
                require(pin == expected, 'PACKAGE_MEMBER_PIN_DIFFERS:' + name)
            return pin, b''.join(chunks) if capture else None
        finally:
            os.close(member)
    finally:
        os.close(fd)

def inventory(root_fd):
    files = set(); directories = set(); count = 0
    def walk(fd, prefix, depth):
        nonlocal count
        require(depth <= MAX_DEPTH, 'PACKAGE_DEPTH_LIMIT')
        before = os.fstat(fd)
        with os.scandir(fd) as entries:
            for item in entries:
                count += 1; require(count <= MAX_ENTRIES, 'PACKAGE_ENTRY_LIMIT')
                name = prefix + item.name; relative(name)
                st = os.stat(item.name, dir_fd=fd, follow_symlinks=False)
                if stat.S_ISDIR(st.st_mode):
                    directories.add(name)
                    child = os.open(item.name, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC,
                                    dir_fd=fd)
                    try:
                        require((os.fstat(child).st_dev, os.fstat(child).st_ino) == (st.st_dev, st.st_ino),
                                'PACKAGE_DIRECTORY_CHANGED')
                        walk(child, name + '/', depth + 1)
                    finally:
                        os.close(child)
                else:
                    require(stat.S_ISREG(st.st_mode) and st.st_nlink == 1,
                            'PACKAGE_SYMLINK_SPECIAL_OR_ALIASED_FILE:' + name)
                    require('__pycache__' not in relative(name).parts and
                            Path(name).suffix not in ('.pyc', '.pyo', '.so', '.pyd', '.dll', '.dylib'),
                            'CACHED_OR_NATIVE_PACKAGE_CODE_FORBIDDEN')
                    files.add(name)
        after = os.fstat(fd)
        require(identity(before) == identity(after), 'PACKAGE_DIRECTORY_CHANGED_DURING_SCAN')
    walk(root_fd, '', 0)
    return files, directories

def expected_directories(files):
    result = set()
    for name in files:
        path = relative(name)
        result.update(parent.as_posix() for parent in path.parents if parent.as_posix() != '.')
    return result

def authenticate_package(package, manifest_sha256, helper_sha256, replay_pins_sha256):
    """Stdlib-only complete byte closure. This does not authorize scientific work."""
    root = directory(package)
    for pin in (manifest_sha256, helper_sha256, replay_pins_sha256):
        digest(pin)
    require(Path(__file__).absolute() == root / 'replay_checkpoint.py',
            'EXECUTED_HELPER_MUST_BE_PACKAGE_HELPER')
    fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC)
    try:
        helper, unused = read_member(fd, 'replay_checkpoint.py')
        require(helper['sha256'] == helper_sha256, 'HELPER_EXTERNAL_PIN_DIFFERS')
        manifest_pin, manifest_raw = read_member(fd, 'PAYLOAD_MANIFEST.json', capture=True)
        require(manifest_pin['sha256'] == manifest_sha256, 'PAYLOAD_MANIFEST_EXTERNAL_PIN_DIFFERS')
        pins_pin, pins_raw = read_member(fd, 'REPLAY_PINS.json', capture=True)
        require(pins_pin['sha256'] == replay_pins_sha256, 'REPLAY_PINS_EXTERNAL_PIN_DIFFERS')
        manifest = load(manifest_raw); pins = load(pins_raw)
        require(type(manifest) is dict and set(manifest) == {'schema_version', 'files'} and
                type(manifest['schema_version']) is int and manifest['schema_version'] == 1 and
                type(manifest['files']) is dict and manifest['files'], 'PAYLOAD_MANIFEST_SCHEMA_DIFFERS')
        rows = manifest['files']
        require(not (set(rows) & RESERVED), 'CIRCULAR_PAYLOAD_MANIFEST')
        total = 0
        for name, item in rows.items():
            relative(name)
            require(type(item) is dict and set(item) == {'bytes', 'sha256'} and
                    type(item['bytes']) is int and 0 <= item['bytes'] <= MAX_FILE,
                    'PACKAGE_MEMBER_SCHEMA_DIFFERS')
            digest(item['sha256']); total += item['bytes']
        require(total <= MAX_PACKAGE, 'PACKAGE_BYTE_LIMIT')
        actual, dirs = inventory(fd)
        require(actual == set(rows) | RESERVED, 'COMPLETE_PACKAGE_MEMBERSHIP_DIFFERS')
        require(dirs == expected_directories(actual), 'UNREGISTERED_EMPTY_PACKAGE_DIRECTORY')
        for name, item in rows.items():
            read_member(fd, name, item)
        # Hash control files again after the whole potentially large closure.
        read_member(fd, 'replay_checkpoint.py', helper)
        read_member(fd, 'PAYLOAD_MANIFEST.json', manifest_pin)
        read_member(fd, 'REPLAY_PINS.json', pins_pin)
        require(inventory(fd) == (actual, dirs), 'PACKAGE_MEMBERSHIP_CHANGED_DURING_AUTHENTICATION')
        named = directory(root); nfd = os.open(named, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        try:
            require((os.fstat(fd).st_dev, os.fstat(fd).st_ino) == (os.fstat(nfd).st_dev, os.fstat(nfd).st_ino),
                    'PACKAGE_ROOT_REPLACED')
        finally:
            os.close(nfd)
    finally:
        os.close(fd)
    validate_pins(root, rows, pins, manifest_sha256, helper_sha256)
    return {'root': root, 'manifest': manifest, 'pins': pins,
            'external': (manifest_sha256, helper_sha256, replay_pins_sha256)}

NUMERICAL_CONTRACT = {
    'arithmetic_bits': 1024, 'momentum_degree': 2048, 'source_degree': 24,
    'source_coefficient_bits': 512, 'source_cells': 64, 'export_bits': 96,
    'export_gate': '1/1000000000000000000', 'anchor': '-9/2', 'endpoints': ['-4', '-7/2'],
    'selected_decode_members': ['k.npy', 'momentum_weights.npy', 'observation_eta.npy',
        'u_1.npy', 'w_1.npy', 'u_2.npy', 'w_2.npy', 'u_3.npy', 'w_3.npy'],
    'wall_seconds': 900, 'rss_kib': 524288, 'output_bytes': 134217728,
}
CUSTODIAN_LIMITS = {
    'wall_seconds': 900, 'cpu_seconds': 900, 'address_space_bytes': 536870912,
    'aggregate_output_bytes': 134217728, 'file_size_bytes': 134217728,
    'process_limit': 0, 'core_bytes': 0, 'custodian_receipt_reserve_bytes': 65536,
}
REQUIRED_KERNEL_DENIALS = {'socket', 'socketpair', 'connect', 'clone', 'execve',
                         'io_uring_setup', 'io_uring_enter', 'io_uring_register'}
ALLOWED_KERNEL_DENIALS = {'socket', 'socketpair', 'connect', 'bind', 'listen', 'accept',
    'accept4', 'fork', 'vfork', 'clone', 'clone3', 'execve', 'execveat',
    'io_uring_setup', 'io_uring_enter', 'io_uring_register'}

def equal_json(left, right):
    # JSON distinguishes bool from integer, unlike Python's True == 1.
    return json.dumps(left, sort_keys=True, allow_nan=False, separators=(',', ':')) == \
           json.dumps(right, sort_keys=True, allow_nan=False, separators=(',', ':'))

def package_json(root, rows, name):
    fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC)
    try:
        unused, raw = read_member(fd, name, rows[name], capture=True)
        return load(raw)
    finally:
        os.close(fd)

def dependency_roster(spec):
    require(type(spec) is dict and type(spec['capsules']) is list and len(spec['capsules']) == 4,
            'FOUR_CAPSULE_INPUT_SPEC_REQUIRED')
    required = [(spec['manifest_repository_path'], spec['manifest_sha256']),
                (spec['producer_repository_path'], spec['producer_sha256']),
                (spec['prior_static_audit']['repository_path'], spec['prior_static_audit']['sha256'])]
    sizes = {spec['prior_static_audit']['repository_path']: spec['prior_static_audit']['bytes']}
    for capsule in spec['capsules']:
        required.extend([(capsule['repository_path'], capsule['reader_inputspec']['file_sha256']),
                         (capsule['original']['repository_path'], capsule['original']['sha256'])])
        sizes[capsule['repository_path']] = capsule['reader_inputspec']['file_bytes']
        sizes[capsule['original']['repository_path']] = capsule['original']['bytes']
    require(len(required) == len({name for name, pin in required}) == 11,
            'EXACT_ELEVEN_ORIGINAL_DEPENDENCIES_REQUIRED')
    for name, pin in required:
        relative(name); digest(pin)
    return dict(required), sizes

def validate_pins(root, rows, pins, manifest_sha256, helper_sha256):
    keys = {'schema_version', 'sealed', 'checkpoint_path', 'core_path', 'registration_path',
        'public_go_path', 'registration_sha256', 'public_go_sha256', 'payload_manifest_sha256',
        'helper_sha256', 'runtime', 'resource_limits', 'saved_runs', 'stable_files',
        'retained_inputs_commit', 'retained_input_dependencies', 'execution_scope'}
    require(type(pins) is dict and set(pins) == keys and type(pins['schema_version']) is int and
            pins['schema_version'] == 1 and pins['sealed'] is True, 'SEALED_SOURCE_REPLAY_PINS_REQUIRED')
    require(pins['checkpoint_path'] == PREFIX and pins['core_path'] == CORE and
            pins['registration_path'] == REGISTRATION and
            pins['payload_manifest_sha256'] == manifest_sha256 and pins['helper_sha256'] == helper_sha256 and
            equal_json(pins['runtime'], RUNTIME) and equal_json(pins['resource_limits'], LIMITS) and
            pins['stable_files'] == list(STABLE), 'SOURCE_REPLAY_BINDING_OR_FIXED_POLICY_DIFFERS')
    digest(pins['registration_sha256'])
    require(REGISTRATION in rows and rows[REGISTRATION]['sha256'] == pins['registration_sha256'],
            'REGISTRATION_PAYLOAD_BINDING_DIFFERS')
    registration = package_json(root, rows, REGISTRATION)
    require(type(registration) is dict and type(registration.get('schema_version')) is int and
            registration['schema_version'] == 1 and
            registration.get('scope') == 'LATER_RETAINED_ENDPOINT_FLOW_CERTIFICATE' and
            equal_json(registration.get('contract'), NUMERICAL_CONTRACT), 'REGISTERED_NUMERICAL_CONTRACT_DIFFERS')
    files = registration.get('files')
    require(type(files) is dict and files and
            {name for name in rows if name.startswith(CORE + '/')} == {CORE + '/' + name for name in files},
            'EXACT_FROZEN_CORE_CLOSURE_REQUIRED')
    for name, item in files.items():
        relative(name)
        require(rows[CORE + '/' + name] == item, 'FROZEN_CORE_PIN_DIFFERS:' + name)
    require(type(pins['saved_runs']) is dict and set(pins['saved_runs']) == {'normal', 'optimized'},
            'DISTINCT_BOTH_SAVED_MODES_REQUIRED')
    paths = [relative(pins['saved_runs'][mode]).as_posix() for mode in ('normal', 'optimized')]
    require(paths[0] != paths[1] and all(path != CORE and not path.startswith(CORE + '/') for path in paths),
            'SAVED_EVIDENCE_MUST_BE_OUTSIDE_FROZEN_CORE')
    for path in paths:
        for name in STABLE + ('RESOURCE.json', 'ENTRY_RECEIPT.json', 'EXECUTION.json', 'child.log'):
            require(path + '/' + name in rows, 'COMPLETE_SAVED_EXECUTION_REQUIRED:' + path + '/' + name)
    for name in STABLE:
        require(rows[paths[0] + '/' + name] == rows[paths[1] + '/' + name],
                'SAVED_NORMAL_OPTIMIZED_STABLE_PINS_DIFFER:' + name)
    require(pins['execution_scope'] in (ACTUAL_SCOPE, CONTROL_SCOPE), 'EXPLICIT_REPLAY_SCOPE_REQUIRED')
    if pins['execution_scope'] == CONTROL_SCOPE:
        require(pins['public_go_path'] is None and pins['public_go_sha256'] is None and
                pins['retained_inputs_commit'] is None and pins['retained_input_dependencies'] == {},
                'MANUFACTURED_CONTROLS_CANNOT_CARRY_PHYSICAL_GO_OR_ORIGINAL_INPUT_ROOT')
    else:
        digest(pins['public_go_sha256'])
        go = relative(pins['public_go_path']).as_posix()
        require(go in rows and go != CORE and not go.startswith(CORE + '/') and
                rows[go]['sha256'] == pins['public_go_sha256'], 'PUBLIC_GO_OUTSIDE_CORE_EXACT_PIN_REQUIRED')
        require(pins['retained_inputs_commit'] == DEPENDENCY_COMMIT, 'IMMUTABLE_ORIGINAL_INPUT_COMMIT_REQUIRED')
        dependencies = pins['retained_input_dependencies']
        spec = package_json(root, rows, CORE + '/provenance_decoder/INPUT_SPEC.json')
        hashes, sizes = dependency_roster(spec)
        require(type(dependencies) is dict and set(dependencies) == set(hashes), 'EXACT_EXTERNAL_INPUT_ROSTER_REQUIRED')
        for name, item in dependencies.items():
            require(type(item) is dict and set(item) == {'bytes', 'sha256'} and
                    type(item['bytes']) is int and 0 <= item['bytes'] <= MAX_FILE and
                    item['sha256'] == hashes[name] and (name not in sizes or item['bytes'] == sizes[name]),
                    'EXTERNAL_INPUT_BYTE_PIN_DIFFERS:' + name)
        require(sum(item['bytes'] for item in dependencies.values()) <= MAX_PACKAGE, 'ORIGINAL_INPUT_DEPENDENCY_BOUND')
        require(not (set(dependencies) & set(rows)), 'SOURCE_BUNDLE_MUST_NOT_DUPLICATE_ORIGINAL_INPUTS')

def reauthenticate_package(v):
    return authenticate_package(v['root'], *v['external'])

def authenticate_dependencies(v, repository_root):
    require(v['pins']['execution_scope'] == ACTUAL_SCOPE and repository_root is not None,
            'ACTUAL_REPLAY_REQUIRES_EXTERNAL_ORIGINAL_INPUT_ROOT')
    root = directory(repository_root)
    require(root != v['root'] and v['root'] not in root.parents, 'ORIGINAL_INPUT_ROOT_MUST_BE_EXTERNAL_TO_BUNDLE')
    fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC)
    try:
        before = os.fstat(fd)
        rows = {name: read_member(fd, name, item)[0]
                for name, item in v['pins']['retained_input_dependencies'].items()}
        nfd = os.open(directory(root), os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC)
        try:
            require((before.st_dev, before.st_ino) == (os.fstat(nfd).st_dev, os.fstat(nfd).st_ino),
                    'ORIGINAL_INPUT_ROOT_REPLACED')
        finally:
            os.close(nfd)
    finally:
        os.close(fd)
    return {'root': root, 'opaque_byte_pins': rows, 'commit': DEPENDENCY_COMMIT,
            'new_array_decodes': 0, 'new_source_callbacks': 0}

def runtime_check():
    require(sys.flags.isolated and sys.flags.dont_write_bytecode, 'USE_PYTHON_-I_-B')
    require(sys.platform == 'linux' and sys.version_info[:3] == (3, 12, 14), 'REGISTERED_LINUX_PYTHON3_12_14_REQUIRED')
    for name, version in RUNTIME['packages'].items():
        require(importlib.metadata.version(name) == version, 'REGISTERED_PACKAGE_VERSION_DIFFERS:' + name)

def source_module(v, name, relative_path):
    require(name not in sys.modules, 'FOREIGN_OR_PREEXISTING_STUDY_MODULE:' + name)
    path = v['root'] / relative_path
    fd = os.open(v['root'], os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC)
    try:
        unused, raw = read_member(fd, relative_path, v['manifest']['files'][relative_path], capture=True)
    finally:
        os.close(fd)
    module = __import__('types').ModuleType(name)
    module.__file__ = str(path)
    sys.modules[name] = module
    try:
        exec(compile(raw, str(path), 'exec'), module.__dict__)
    except BaseException:
        del sys.modules[name]
        raise
    return module

def stable_pins(directory_path):
    root = directory(directory_path)
    fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC)
    try:
        return {name: read_member(fd, name)[0] for name in STABLE}
    finally:
        os.close(fd)

def compare_science(left, right):
    left_pins, right_pins = stable_pins(left), stable_pins(right)
    require(left_pins == right_pins, 'NINE_STABLE_SCIENTIFIC_BYTES_DIFFER')
    return {'status': 'PASS_EXACT_NINE_STABLE_SCIENTIFIC_FILES', 'files': left_pins}

def output_snapshot(root, root_fd):
    files, dirs = inventory(root_fd)
    require('EXECUTION.json' in files, 'AUTHORITATIVE_CUSTODIAN_RECEIPT_REQUIRED')
    pins = {name: read_member(root_fd, name)[0] for name in sorted(files)}
    require(sum(item['bytes'] for item in pins.values()) <= LIMITS['aggregate_output_bytes'],
            'TOTAL_ACCEPTED_OUTPUT_LIMIT_EXCEEDED')
    tree = hashlib.sha256()
    for name, item in pins.items():
        if name == 'EXECUTION.json':
            continue
        tree.update(name.encode('utf-8', 'surrogateescape') + b'\0' + str(item['bytes']).encode() + b'\0' +
                    bytes.fromhex(item['sha256']))
    summary = {'bytes': sum(item['bytes'] for name, item in pins.items() if name != 'EXECUTION.json'),
               'file_count': len(pins) - 1, 'directory_count': len(dirs), 'tree_sha256': tree.hexdigest()}
    return pins, summary

def validate_custodian(output, mode, pins, manufactured=False):
    require(mode in ('normal', 'optimized') and type(manufactured) is bool, 'EXPLICIT_MODE_AND_SCOPE_REQUIRED')
    root = directory(output)
    fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC)
    try:
        before = os.fstat(fd)
        files, summary = output_snapshot(root, fd)
        receipt_pin, raw = read_member(fd, 'EXECUTION.json', capture=True)
        require(receipt_pin['bytes'] <= 65536, 'BOUNDED_CUSTODIAN_RECEIPT_REQUIRED')
        receipt = load(raw)
        unused, entry_raw = read_member(fd, 'ENTRY_RECEIPT.json', capture=True)
        unused, readback_raw = read_member(fd, 'OUTPUT_READBACK.json', capture=True)
        unused, resource_raw = read_member(fd, 'RESOURCE.json', capture=True)
        entry, readback, resources = load(entry_raw), load(readback_raw), load(resource_raw)
        require(identity(before) == identity(os.fstat(fd)), 'OUTPUT_DIRECTORY_CHANGED_DURING_CUSTODY_CHECK')
        nfd = os.open(directory(root), os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
        try:
            require((before.st_dev, before.st_ino) == (os.fstat(nfd).st_dev, os.fstat(nfd).st_ino),
                    'OUTPUT_ROOT_REPLACED')
        finally:
            os.close(nfd)
    finally:
        os.close(fd)
    expected_status = 'PASS_BOUNDED_MANUFACTURED_EXECUTION' if manufactured else 'PASS_BOUNDED_REGISTERED_EXECUTION'
    require(type(receipt) is dict and receipt.get('status') == expected_status and
            receipt.get('manufactured_only') is manufactured and receipt.get('optimized') is (mode == 'optimized'),
            'CUSTODIAN_SCOPE_MODE_OR_SUCCESS_DIFFERS:' + mode)
    require(all(type(receipt.get(key)) is int and receipt[key] == 0
                for key in ('exit_code', 'wrapper_exit_code', 'wait_status')) and
            receipt.get('stop_reason') is None and receipt.get('monitor_error') is None and
            receipt.get('preserved_child_log') is None, 'AUTHORITATIVE_WAIT4_NOT_SUCCESS:' + mode)
    require(equal_json(receipt.get('limits'), CUSTODIAN_LIMITS) and
            receipt.get('aggregate_limit_enforcement') == 'POLL_AND_FINAL_ACCEPTANCE_CHECK_WITH_RECEIPT_RESERVE',
            'CUSTODIAN_FIXED_LIMITS_DIFFER:' + mode)
    require(receipt.get('entry_receipt') == files['ENTRY_RECEIPT.json'] and
            receipt.get('captured_child_log_sha256') == files['child.log']['sha256'] and
            receipt.get('observed_child_log_sha256') == files['child.log']['sha256'],
            'CUSTODIAN_ENTRY_OR_LOG_LINK_DIFFERS:' + mode)
    require(type(receipt.get('uid')) is int and receipt['uid'] > 0 and
            type(receipt.get('peak_rss_kib_wait4')) is int and
            0 < receipt['peak_rss_kib_wait4'] * 1024 <= LIMITS['address_space_bytes'], 'CUSTODIAN_UID_OR_RSS_DIFFERS')
    for key in ('wall_seconds', 'cpu_seconds_wait4'):
        value = receipt.get(key)
        require(type(value) in (int, float) and math.isfinite(value) and 0 <= value <= 900,
                'CUSTODIAN_MEASURED_RESOURCE_DIFFERS:' + key)
    require(equal_json(receipt.get('worker_output'), summary) and summary['bytes'] <= 134217728 - 65536,
            'CUSTODIAN_COMPLETE_OUTPUT_TREE_OR_BOUND_DIFFERS')
    expected_entry = 'PASS_MANUFACTURED_ENDPOINT_ENTRY' if manufactured else 'PASS_REGISTERED_ENDPOINT_ENTRY'
    require(type(entry) is dict and entry.get('status') == expected_entry and
            entry.get('scope') == 'LATER_RETAINED_ENDPOINT_FLOW_CERTIFICATE' and
            entry.get('manufactured') is manufactured and
            entry.get('registration_sha256') == pins['registration_sha256'] and
            entry.get('go_sha256') == pins['public_go_sha256'] and
            type(entry.get('nodes')) is int and entry['nodes'] == 49152 and
            type(entry.get('endpoint_comparisons')) is int and entry['endpoint_comparisons'] == 98304,
            'ENTRY_SCIENTIFIC_SCOPE_AUTHORIZATION_OR_COVERAGE_DIFFERS')
    require(equal_json(entry.get('outputs'), {name: files[name] for name in STABLE + ('RESOURCE.json',)}),
            'ENTRY_COMPLETE_SCIENTIFIC_RESOURCE_PIN_LINK_DIFFERS')
    denied = entry.get('kernel_denied_syscalls')
    require(type(denied) is list and all(type(x) is str for x in denied) and len(denied) == len(set(denied)) and
            REQUIRED_KERNEL_DENIALS <= set(denied) <= ALLOWED_KERNEL_DENIALS,
            'MANDATORY_KERNEL_DENIAL_EVIDENCE_MISSING')
    require(type(readback) is dict and readback.get('status') == 'PASS_STANDALONE_ENDPOINT_READBACK' and
            readback.get('manufactured') is manufactured and
            all(type(readback.get(k)) is int and readback[k] == n for k, n in
                (('nodes', 49152), ('endpoint_comparisons', 98304), ('prefix_cases', 24),
                 ('source_evaluations', 0), ('retained_arrays_opened', 0))), 'READBACK_SCOPE_OR_COVERAGE_DIFFERS')
    require(type(resources) is dict and type(resources.get('output_bytes')) is int and
            0 <= resources['output_bytes'] <= summary['bytes'] and
            type(resources.get('peak_rss_kib')) is int and
            0 < resources['peak_rss_kib'] * 1024 <= LIMITS['address_space_bytes'] and
            type(resources.get('wall_seconds')) in (int, float) and math.isfinite(resources['wall_seconds']) and
            0 <= resources['wall_seconds'] <= 900, 'ENTRY_RESOURCE_MEASUREMENTS_DIFFERS')
    # Command is host provenance; require the actual mode/registered pin tokens,
    # without demanding the saved host's absolute source/output paths still exist.
    command = receipt.get('command')
    require(type(command) is list and command and all(type(x) is str for x in command) and
            '-I' in command and '-B' in command and (('-O' in command) is (mode == 'optimized')) and
            '--registration-sha256' in command and command.count('--registration-sha256') == 1 and
            command[command.index('--registration-sha256') + 1] == pins['registration_sha256'] and
            (('--manufactured-only' in command) is manufactured), 'CUSTODIAN_COMMAND_BINDING_DIFFERS')
    if not manufactured:
        require(command.count('--go-sha256') == 1 and command[command.index('--go-sha256') + 1] == pins['public_go_sha256'] and
                '--repository-root' in command, 'CUSTODIAN_PUBLIC_GO_COMMAND_BINDING_DIFFERS')
    return {'status': 'PASS_AUTHENTICATED_CUSTODY_AND_ENTRY', 'mode': mode,
            'manufactured_control': manufactured, 'execution': receipt_pin,
            'entry': files['ENTRY_RECEIPT.json'], 'output': summary}

def registered_command(v, output, mode, repository_root=None, manufactured=False):
    require(mode in ('normal', 'optimized') and type(manufactured) is bool, 'REGISTERED_MODE_SCOPE_REQUIRED')
    core = v['root'] / CORE
    command = [sys.executable, '-I', '-B', str(core / 'execution/bounded_launcher.py'),
               '--candidate', str(core), '--output', str(output)]
    if mode == 'optimized':
        command.append('--optimized')
    if manufactured:
        command.append('--manufactured-only')
    command += ['--', '--registration', str(v['root'] / REGISTRATION),
                '--registration-sha256', v['pins']['registration_sha256']]
    if not manufactured:
        require(repository_root is not None, 'ACTUAL_REPLAY_REQUIRES_ORIGINAL_ROOT')
        command += ['--go', str(v['root'] / v['pins']['public_go_path']),
                    '--go-sha256', v['pins']['public_go_sha256'], '--repository-root', str(repository_root)]
    return command

def launch_registered(v, output, mode, logs, repository_root=None, manufactured=False):
    command = registered_command(v, output, mode, repository_root, manufactured)
    with (logs / (mode + '-custodian.stdout.log')).open('xb') as stdout, \
         (logs / (mode + '-custodian.stderr.log')).open('xb') as stderr:
        process = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=stdout, stderr=stderr,
                                   env={'PATH': '/usr/bin:/bin', 'LANG': 'C.UTF-8', 'LC_ALL': 'C.UTF-8'})
        try:
            code = process.wait()
        except BaseException:
            process.send_signal(signal.SIGINT)
            try:
                process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                process.kill(); process.wait()
            raise
    require(code == 0, 'BOUNDED_REPLAY_PROCESS_FAILED:' + mode)
    return validate_custodian(output, mode, v['pins'], manufactured)

def replay(v, verify_only, repository_root=None, output_root=None, manufactured_control=False):
    require(type(verify_only) is bool and type(manufactured_control) is bool, 'EXPLICIT_OPERATION_SCOPE_REQUIRED')
    # Consume a newly authenticated disk context, never mutable cached caller
    # maps. The three externally supplied control pins remain the trust root.
    v = reauthenticate_package(v)
    require((v['pins']['execution_scope'] == CONTROL_SCOPE) is manufactured_control,
            'MANUFACTURED_CONTROL_FLAG_MUST_MATCH_SEPARATE_PINNED_SCOPE')
    require(not manufactured_control or repository_root is None, 'CONTROL_CANNOT_RECEIVE_REAL_INPUT_ROOT')
    runtime_check()
    dependencies = None if manufactured_control else authenticate_dependencies(v, repository_root)
    pins = v['pins']; core = v['root'] / CORE
    # Capture and execute only authenticated stdlib modules. Do not add study
    # paths or read local aliases before the complete source/input authentication.
    module_names = ('registration_guard', 'endpoint_aggregate', 'validate_outputs')
    require(all(name not in sys.modules for name in module_names), 'PREEXISTING_STUDY_ALIAS_FORBIDDEN')
    original_path = list(sys.path)
    try:
        guard = source_module(v, 'registration_guard', CORE + '/execution/registration_guard.py')
        verified = guard.authenticate(core, v['root'] / REGISTRATION, pins['registration_sha256'],
            manufactured_control, None if manufactured_control else v['root'] / pins['public_go_path'],
            pins['public_go_sha256'])
        source_module(v, 'endpoint_aggregate', CORE + '/engine/endpoint_aggregate.py')
        reader = source_module(v, 'validate_outputs', CORE + '/execution/validate_outputs.py')
        saved = {mode: v['root'] / pins['saved_runs'][mode] for mode in ('normal', 'optimized')}
        custody = {mode: validate_custodian(path, mode, pins, manufactured_control) for mode, path in saved.items()}
        checks = {mode: reader.validate(path) for mode, path in saved.items()}
        for mode in checks:
            require(checks[mode]['manufactured'] is manufactured_control, 'SERIALIZED_READBACK_MODE_DIFFERS')
            require(equal_json(checks[mode], load((saved[mode] / 'OUTPUT_READBACK.json').read_bytes())),
                    'RECOMPUTED_SAVED_READBACK_DIFFERS:' + mode)
        comparison = compare_science(saved['normal'], saved['optimized'])
        guard.reauthenticate(verified); reauthenticate_package(v)
        if not manufactured_control:
            authenticate_dependencies(v, repository_root)
        if verify_only:
            require(output_root is None, 'VERIFY_ONLY_CANNOT_CREATE_REPLAY_OUTPUT')
            return {'status': 'PASS_VERIFIED_MANUFACTURED_SAVED_OUTPUTS' if manufactured_control else
                    'PASS_VERIFIED_SAVED_RETAINED_ENDPOINT_OUTPUTS', 'manufactured_control': manufactured_control,
                    'new_array_decodes': 0, 'new_source_callbacks': 0, 'saved_custody': custody,
                    'saved_checks': checks, 'saved_comparison': comparison,
                    'original_dependency_files_authenticated': 0 if manufactured_control else len(dependencies['opaque_byte_pins'])}
        require(output_root is not None and os.getuid() > 0 and os.geteuid() > 0,
                'FRESH_REPLAY_REQUIRES_NONROOT_FRESH_EXTERNAL_OUTPUT')
        out = Path(os.path.abspath(output_root)); directory(out.parent)
        forbidden = [v['root']] + ([] if manufactured_control else [dependencies['root']])
        require(all(out != path and path not in out.parents and out not in path.parents for path in forbidden),
                'FRESH_REPLAY_OUTPUT_MUST_BE_OUTSIDE_SOURCE_AND_INPUT_ROOTS')
        os.mkdir(out, 0o700)
        fresh_checks = {}; fresh_custody = {}
        for mode in ('normal', 'optimized'):
            reauthenticate_package(v); guard.reauthenticate(verified)
            if not manufactured_control:
                authenticate_dependencies(v, repository_root)
            run = out / mode
            fresh_custody[mode] = launch_registered(v, run, mode, out, repository_root, manufactured_control)
            fresh_checks[mode] = reader.validate(run)
            require(equal_json(fresh_checks[mode], load((run / 'OUTPUT_READBACK.json').read_bytes())),
                    'RECOMPUTED_FRESH_READBACK_DIFFERS:' + mode)
            compare_science(saved['normal'], run)
            guard.reauthenticate(verified); reauthenticate_package(v)
        fresh_comparison = compare_science(out / 'normal', out / 'optimized')
        if not manufactured_control:
            authenticate_dependencies(v, repository_root)
        result = {'status': 'PASS_COMPLETE_MANUFACTURED_SOURCE_REPLAY_CONTROL' if manufactured_control else
                  'PASS_COMPLETE_RETAINED_ENDPOINT_SOURCE_REPLAY', 'manufactured_control': manufactured_control,
                  'saved_custody': custody, 'saved_checks': checks, 'saved_comparison': comparison,
                  'fresh_custody': fresh_custody, 'fresh_checks': fresh_checks, 'fresh_comparison': fresh_comparison,
                  'exact_stable_files': list(STABLE), 'resource_entry_execution_are_separate_custody_evidence': True,
                  'new_physical_source_callbacks': 0 if manufactured_control else 256,
                  'new_retained_array_decodes': 0 if manufactured_control else 72}
        with (out / 'PORTABLE_REPLAY_RECEIPT.json').open('x') as stream:
            stream.write(json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + '\n')
        return result
    finally:
        sys.path[:] = original_path
        for name in module_names:
            sys.modules.pop(name, None)

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package', required=True)
    parser.add_argument('--payload-manifest-sha256', required=True)
    parser.add_argument('--helper-sha256', required=True)
    parser.add_argument('--replay-pins-sha256', required=True)
    parser.add_argument('--repository-root')
    parser.add_argument('--manufactured-control', action='store_true')
    parser.add_argument('--verify-only', action='store_true')
    parser.add_argument('--output-root')
    args = parser.parse_args(argv)
    require(args.verify_only == (args.output_root is None), 'VERIFY_ONLY_OR_FRESH_OUTPUT_REQUIRED')
    v = authenticate_package(args.package, args.payload_manifest_sha256, args.helper_sha256, args.replay_pins_sha256)
    result = replay(v, args.verify_only, args.repository_root, args.output_root, args.manufactured_control)
    print(json.dumps(result, sort_keys=True, separators=(',', ':'), allow_nan=False))
    return 0

if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (Stop, ValueError, OSError, KeyError, IndexError, importlib.metadata.PackageNotFoundError) as exc:
        print(json.dumps({'status': 'STOPPED', 'reason': str(exc)}), file=sys.stderr)
        raise SystemExit(2)

