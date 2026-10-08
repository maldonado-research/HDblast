"""Stdlib-only authentication for the new all-capsule comparison contract.

Only root-supplied external pins authorize actual work. This module does no
network I/O and cannot establish public publication by creating a receipt.
"""
from fractions import Fraction
from pathlib import Path, PurePosixPath
import hashlib
import json
import os
import platform
import re
import stat
import sys
from importlib import metadata

PREFIX = 'research/HDBLAST_CHECKPOINT_20261007_STORED_STATE_COMPARISON'
REMOTE_GO = 'PASS_REMOTE_REGISTERED_STORED_STATE_COMPARISON_GO'
INPUT_SCOPE = 'ALL_FOUR_RETAINED_INCOMING_CAPSULES_AT_FIXED_ANCHOR'
SOURCES = ('positive_B', 'signed_uB')
CAPSULES = tuple(s+'/'+g for s in SOURCES for g in ('coarse', 'fine'))
ZERO_KEYS = {'registered_source_callbacks', 'retained_array_decodes',
             'stored_state_comparisons', 'physical_trajectories', 'likelihood_evaluations'}
REQUIRED = {'REGISTRATION_CONTRACT.json', 'execution/registration_guard.py',
            'execution/run_registered.py', 'execution/comparison_worker.py',
            'execution/source_certificate.py', 'execution/readback.py',
            'execution/bounded_launcher.py', 'execution/kernel_guard.py',
            'decoder/exact_binary80.py', 'decoder/fabricated_fixtures.py',
            'theory/entire_target_engine.py', 'transport/discrete_transport.py',
            'source/route.py', 'source/source_algebra_baseline.py',
            'source/generic_operator_baseline.py', 'SOURCE_PINS.json',
            'ENVIRONMENT_PINS.json',
            'prior/BUDGET.json', 'prior/primary_DATA.json', 'prior/independent_DATA.json',
            'provenance/INPUT_SPEC.json', 'provenance/LAYOUT_PROVENANCE_BASIS.json',
            'provenance/ORIGINAL_ZIP_INVENTORIES.json'}

def require(ok, message):
    if not ok:
        raise ValueError(message)

def load(raw):
    def unique(items):
        out = {}
        for key, value in items:
            require(key not in out, 'Duplicate JSON key')
            out[key] = value
        return out
    def bad(value):
        raise ValueError('Nonfinite JSON constant: '+value)
    return json.loads(raw, object_pairs_hook=unique, parse_constant=bad)

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def digest(value):
    require(type(value) is str and re.fullmatch('[0-9a-f]{64}', value), 'External SHA256 required')
    return value

def rational(value):
    require(type(value) is str and re.fullmatch(r'-?(?:0|[1-9][0-9]*)/[1-9][0-9]*', value), 'Canonical rational required')
    result = Fraction(value)
    require(value == f'{result.numerator}/{result.denominator}', 'Noncanonical rational')
    return result

def real_directory(path):
    path = Path(os.path.abspath(path))
    require(path.is_dir() and path.resolve() == path, 'Real directory required')
    require(all(not p.is_symlink() for p in (path, *path.parents)), 'Symlink ancestor forbidden')
    return path

def safe_file(root, name):
    require(type(name) is str and name and '\\' not in name, 'Relative filename required')
    parts = PurePosixPath(name)
    require(not parts.is_absolute() and parts.as_posix() == name and
            all(p not in ('', '.', '..') for p in name.split('/')), 'Unsafe path')
    path = root
    for part in parts.parts:
        path = path/part
        require(not path.is_symlink(), 'Payload symlink forbidden')
    require(path.is_file() and stat.S_ISREG(path.stat().st_mode), 'Regular file required: '+name)
    return path

def zero_counts(value):
    return type(value) is dict and set(value) == ZERO_KEYS and all(type(x) is int and x == 0 for x in value.values())

FIXED_CONTRACT = {
    'schema_version': 1, 'input_scope': INPUT_SCOPE,
    'output_scope': 'ALL_NODE_EXACT_STORED_MINUS_CERTIFIED_BD_INCOMING_TARGET',
    'target_universe': 'PINNED_KNOWN_INPUT_49152_CAPSULE_NODES_NOT_ADAPTIVE',
    'incoming_anchor_eta': '-9/2', 'sources': list(SOURCES), 'capsules': list(CAPSULES),
    'distinct_capsule_nodes': 49152, 'prefix_cases': 12, 'prefix_node_appearances': 86016,
    'selected_arrays': 20, 'selected_scalar_binary80_slots': 294936,
    'physical_source_callbacks': 22, 'diagnostic_target_evaluations': 18,
    'target_node_evaluations': 49152, 'source_degree': 112, 'geometric_parents': 11,
    'coefficient_bits': 512, 'momentum_degree': 832, 'export_bits': 96,
    'max_complete_exported_L1_radius': '1/1000000000000000000',
    'represented_epsilon': '3777893186295716171/37778931862957161709568',
    'represented_Pi': '14488038916154245685/4611686018427387904',
    'source_error_rule': 'max(rebuilt_error,prior_uniform_error,prior_analytic_tail+prior_coefficient_error)',
    'source_builder': 'EXACT_THREE_IMMUTABLE_PRIOR_INDEPENDENT_FILES',
    'source_coefficient_pin': 'prior/BUDGET.json chosen_coefficients_sha256',
    'array_decode_policy': 'VERIFY_ALL_FOUR_CAPSULES_AND_ALL_19_MEMBERS_BEFORE_FIRST_SELECTED_DECODE',
    'other_members_decoded': 0, 'state_projection': False,
    'source_coordinate_equality': 'COMPARE_EXACT_DECODED_VALUES_AFTER_GO_ONLY',
    'diagnostic_crosscheck': '18_ROWS_36_RECTANGLE_PAIRS_72_COMPONENTS_BOTH_PRIOR_METHODS',
    'diagnostic_interpolation_certifies_all_nodes': False,
    'cap_delta': '1/128', 'cap_method': 'UNCHANGED_INDEPENDENT_FLAT_CAP_IBP_ENDPOINT_BOUND',
    'physical_trajectories': 0, 'likelihood_evaluations': 0,
    'higher_dimensional_Big_Bang_origin': 'NOT_ESTABLISHED',
    'full_twelve_case_continuous_pressure_contacts': 'UNRESOLVED',
    'metric_calibration': 'FAIL', 'external_novelty': 'NOT_ASSESSED',
    'resources': {'wall_seconds': 900, 'cpu_seconds': 900, 'address_space_bytes': 536870912,
                  'aggregate_output_bytes': 134217728, 'file_size_bytes': 134217728,
                  'nonroot': True, 'single_threads': True, 'nproc': 0, 'core': 0,
                  'network_allowed': False, 'child_processes_allowed': False},
}

def validate_contract(contract):
    require(type(contract) is dict, 'Contract object required')
    for key, expected in FIXED_CONTRACT.items():
        require(type(contract.get(key)) is type(expected) and contract[key] == expected,
                'Comparison contract differs: '+key)
    require(zero_counts(contract.get('new_evaluations_before_public_GO')), 'Nonprospective counters')
    return contract

def runtime_identity():
    require(sys.version_info[:3] == (3, 12, 14), 'Python 3.12.14 required')
    versions = {name: metadata.version(name) for name in ('python-flint', 'sympy', 'mpmath')}
    require(versions == {'python-flint': '0.9.0', 'sympy': '1.14.0', 'mpmath': '1.3.0'}, 'Pinned runtime differs')
    executable = Path(sys.executable).resolve()
    native = {}
    for name in ('python-flint','sympy','mpmath'):
        distribution = metadata.distribution(name)
        hashes = {}
        for member in distribution.files or ():
            if str(member).endswith(('.pyc','.pyo')) or str(member).endswith('/RECORD'):
                continue
            path = Path(distribution.locate_file(member))
            if path.is_file():
                hashes[str(member)] = sha(path.read_bytes())
        native[name] = sha(json.dumps(hashes,sort_keys=True,separators=(',',':')).encode())
    return {'python': sys.version, 'machine': platform.machine(), 'platform': sys.platform, 'packages': versions,
            'observed_python_executable':str(executable),'observed_python_executable_sha256':sha(executable.read_bytes()),
            'observed_distribution_content_sha256':native,
            'observed_hash_scope':'CURRENT_RUNTIME_PROVENANCE_NOT_CROSS_HOST_BYTE_EQUALITY_REQUIREMENT'}

def authenticate_runtime(root):
    pins = load(safe_file(root, 'ENVIRONMENT_PINS.json').read_bytes())
    require(type(pins) is dict and pins.get('schema_version') == 1 and
            pins.get('python_version') == '3.12.14' and
            pins.get('platform') == 'linux' and
            pins.get('native_guard') == 'libseccomp.so.2 fail-closed syscall filter',
            'Registered runtime/ABI requirements differ')
    require(pins.get('packages') == {'python-flint':'0.9.0','sympy':'1.14.0','mpmath':'1.3.0'},
            'Dependency pins differ')
    require(sys.version_info[:3] == (3,12,14) and sys.platform == 'linux', 'Registered Python/Linux required')
    versions={name:metadata.version(name) for name in pins['packages']}
    require(versions==pins['packages'],'Registered runtime package versions differ')
    require(pins.get('byte_hash_scope') == 'OBSERVED_PREPARATION_RUNTIME_PROVENANCE_ONLY',
            'Runtime hashes must retain their observation scope')
    return pins

def _local(candidate, registration_sha256):
    root = real_directory(candidate)
    digest(registration_sha256)
    raw = safe_file(root, 'FULL_REGISTRATION.json').read_bytes()
    require(sha(raw) == registration_sha256, 'Registration external pin differs')
    registration = load(raw)
    require(type(registration) is dict and set(registration) == {'schema_version', 'files'} and
            type(registration['schema_version']) is int and registration['schema_version'] == 1,
            'Exact registration schema required')
    files = registration['files']
    require(type(files) is dict and REQUIRED <= set(files) and 'FULL_REGISTRATION.json' not in files,
            'Incomplete or circular registration')
    for name, item in files.items():
        require(type(item) is dict and set(item) == {'bytes', 'sha256'} and type(item['bytes']) is int and
                0 <= item['bytes'] <= 20*1024*1024, 'Invalid registration member')
        digest(item['sha256'])
        path = safe_file(root, name)
        require(path.stat().st_size == item['bytes'] and sha(path.read_bytes()) == item['sha256'],
                'Registered file changed: '+name)
    for path in root.rglob('*'):
        require(not path.is_symlink(), 'Candidate symlink forbidden')
        if path.is_file():
            require(path.suffix not in ('.pyc', '.pyo', '.so', '.pyd', '.dll') and '__pycache__' not in path.parts,
                    'Cached/native candidate code forbidden')
            if path.suffix == '.py':
                require(path.relative_to(root).as_posix() in files, 'Unregistered Python source')
    contract = validate_contract(load(safe_file(root, 'REGISTRATION_CONTRACT.json').read_bytes()))
    authenticate_runtime(root)
    return {'root': root, 'registration': registration, 'contract': contract,
            'registration_sha256': registration_sha256, 'registration_raw': raw}

def authenticate_local(candidate, registration_sha256):
    verified = _local(candidate, registration_sha256)
    verified['fabricated_only'] = True
    return verified

def _verify_remote(item, raw, commit, name):
    require(type(item) is dict and type(item.get('bytes')) is int and item['bytes'] == len(raw) and
            item.get('sha256') == sha(raw), 'Remote size/SHA differs: '+name)
    git_blob = hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
    require(item.get('git_blob') == git_blob and item.get('independent_public_blob_download_sha256_verified') is True,
            'Independent public Git identity missing: '+name)
    require(item.get('immutable_contents_url') ==
            f'https://api.github.com/repos/maldonado-research/HDblast/contents/{PREFIX}/{name}?ref={commit}',
            'Immutable public path binding differs: '+name)

def authenticate(candidate, registration_sha256, public_go, public_go_sha256):
    digest(public_go_sha256)
    receipt_path = Path(os.path.abspath(public_go))
    require(receipt_path.is_file() and not receipt_path.is_symlink() and
            all(not p.is_symlink() for p in receipt_path.parents), 'Real external PUBLIC_GO receipt required')
    require(receipt_path.stat().st_size <= 20*1024*1024, 'Oversized PUBLIC_GO receipt')
    raw = receipt_path.read_bytes()
    require(sha(raw) == public_go_sha256, 'PUBLIC_GO external pin differs')
    receipt = load(raw)
    require(type(receipt) is dict and receipt.get('status') == REMOTE_GO and
            receipt.get('input_scope') == INPUT_SCOPE and
            receipt.get('public_repository') == 'https://github.com/maldonado-research/HDblast' and
            receipt.get('repository_checkpoint_path') == PREFIX, 'PUBLIC_GO status/scope/path differs')
    require(receipt.get('all_registered_public_blobs_independently_downloaded') is True and
            zero_counts(receipt.get('new_target_evaluations_before_public_verification')),
            'Premature work or incomplete public verification')
    commit = receipt.get('freeze_commit')
    require(type(commit) is str and re.fullmatch('[0-9a-f]{40}', commit), 'Immutable public commit required')
    require(receipt.get('registration_sha256') == registration_sha256, 'PUBLIC_GO registration differs')
    v = _local(candidate, registration_sha256)
    require(v['root'] not in receipt_path.parents, 'PUBLIC_GO must be outside candidate')
    remote = receipt.get('remote_registered_files')
    files = v['registration']['files']
    require('evidence/preparation/FINAL_PRE_FREEZE_REVIEW.json' in files,
            'Independent prefreeze review must be registered and publicly authenticated')
    require(type(remote) is dict and set(remote) == set(files)|{'FULL_REGISTRATION.json'} and
            type(receipt.get('remote_files_verified')) is int and receipt['remote_files_verified'] == len(remote),
            'Complete immutable public inventory required')
    for name in files:
        _verify_remote(remote[name], safe_file(v['root'], name).read_bytes(), commit, name)
    _verify_remote(remote['FULL_REGISTRATION.json'], v['registration_raw'], commit, 'FULL_REGISTRATION.json')
    review = load(safe_file(v['root'], 'evidence/preparation/FINAL_PRE_FREEZE_REVIEW.json').read_bytes())
    require(type(review) is dict and review.get('status') == 'GO_FOR_PROSPECTIVE_PUBLIC_FREEZE',
            'Independent prefreeze review GO required')
    reviewed = review.get('reviewed_source_sha256')
    require(type(reviewed) is dict and {n for n in files if n.endswith('.py')} <= set(reviewed),
            'Independent review must pin every Python file')
    for name, pin in reviewed.items():
        require(name in files and files[name]['sha256'] == pin, 'Reviewed source changed: '+name)
    v.update({'receipt': receipt, 'receipt_path': receipt_path, 'public_go_sha256': public_go_sha256,
              'fabricated_only': False})
    return v

def reauthenticate(v):
    if v['fabricated_only']:
        return authenticate_local(v['root'], v['registration_sha256'])
    return authenticate(v['root'], v['registration_sha256'], v['receipt_path'], v['public_go_sha256'])

def journal_event(path, event):
    with Path(path).open('ab') as handle:
        handle.write((json.dumps(event, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode())
        handle.flush()
        os.fsync(handle.fileno())

class SourceAuthorization:
    def __init__(self, verified, journal, centers):
        require(not verified['fabricated_only'], 'Fabricated guard cannot authorize real sources')
        self.verified, self.journal = verified, Path(journal)
        self.journal.open('xb').close()
        self.centers = tuple(centers)
        require(len(self.centers) == 11 and len(set(self.centers)) == 11, 'Eleven unique source parents required')
        self.seen = set()

    def __call__(self, event, detail):
        require(event == 'real_source_taylor' and type(detail) is dict and
                set(detail) == {'source', 'center', 'scope'}, 'Unexpected physical callback')
        require(detail['scope'] == 'UNIQUE_ANALYTIC_BD_PREHISTORY_ONLY' and detail['source'] in SOURCES,
                'Changed physical source identity')
        center = rational(detail['center'])
        key = (detail['source'], center)
        require(center in self.centers and key not in self.seen, 'Repeated/unregistered source callback')
        reauthenticate(self.verified)
        journal_event(self.journal, {'index': len(self.seen)+1, 'event': event, **detail})
        self.seen.add(key)

    def complete(self):
        require(self.seen == {(s,c) for s in SOURCES for c in self.centers}, 'Incomplete physical source construction')
        return len(self.seen)

class DecodeAuthorization:
    def __init__(self, verified, journal, snapshots):
        require(not verified['fabricated_only'], 'Fabricated guard cannot authorize retained decode')
        require(tuple(snapshots) == CAPSULES, 'All four opaque snapshots required before decode')
        self.verified, self.journal = verified, Path(journal)
        self.journal.open('xb').close()
        self.seen = set()

    def before(self, capsule, member):
        require(capsule in CAPSULES and member in ('k.npy','momentum_weights.npy','observation_eta.npy','u_1.npy','w_1.npy'),
                'Unregistered retained decode')
        key = (capsule, member)
        require(key not in self.seen, 'Repeated retained-array decode')
        reauthenticate(self.verified)
        journal_event(self.journal, {'index': len(self.seen)+1, 'capsule': capsule, 'member': member,
                                    'event': 'selected_retained_array_decode_attempt'})
        self.seen.add(key)

    def complete(self):
        require(len(self.seen) == 20, 'Incomplete selected-array decode universe')
        return len(self.seen)
