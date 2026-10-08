"""Stdlib guard for an externally pinned, root-authenticated public freeze.

No network operation occurs here. A supplied hash authenticates supplied bytes;
it does not turn an arbitrary self-produced receipt into independent public
authentication. The custodian must obtain expected receipt/registration hashes
from the root's authenticated immutable GitHub readback and independent review.
"""
from fractions import Fraction
from pathlib import Path, PurePosixPath
import hashlib
import json
import os
import re
import stat
import sys
import platform
from importlib import metadata

SOURCES = ('positive_B', 'signed_uB')
MOMENTA = ('0/1', '1/1099511627776', '1/4096', '1/4', '1/1', '16/1', '64/1', '128/1', '256/1')
INPUT_SCOPE = 'NO_RETAINED_ARRAYS_BD_PREHISTORY_ONLY'
OUTPUT_SCOPE = 'BD_PREHISTORY_TARGET_AT_FIXED_RATIONAL_PROBES'
REMOTE_GO = 'PASS_REMOTE_REGISTERED_BD_PREHISTORY_GO'
ZERO_COUNTER_KEYS = {'registered_source_callbacks', 'retained_array_decodes', 'stored_state_comparisons', 'physical_trajectories', 'likelihood_evaluations'}
CAP_FORMULAE = {
    'primary': 'Direct L1: D=delta*(2-delta); ell=1/(5-delta); H=(3/8)^63. G_positive=delta*H*(4*ell^2+4*ell/D^2+4/D^4); G_signed=delta*H*(4*ell^2+2*ell+(4*ell+4)/D^2+4/D^4). Ucap=G_source/2, Wcap=G_source.',
    'independent': 'IBP: D=delta*(2-delta); ell=1/(5-delta); H=(3/8)^63; p=2*(1-delta)/D^2; rho=0 for positive_B and 1 for signed_uB. Ucap=H*((p+rho+2*ell)/2+1+delta*(3*ell^2+2*k+2*ell)); Wcap=H*(p+rho+2*ell+2*k+delta*(4*k^2+4*k*ell+6*ell^2)).',
}
REQUIRED = {
    'REGISTRATION_CONTRACT.json', 'primary/route.py', 'independent/route.py',
    'execution/registration_guard.py', 'execution/run_registered.py',
    'execution/execute_bounded.py', 'execution/validate_outputs.py',
    'execution/test_controls.py', 'evidence/preparation/FINAL_PRE_FREEZE_REVIEW.json',
}

def require(ok, message):
    if not ok:
        raise ValueError(message)

def load(raw):
    def unique(items):
        value = {}
        for key, item in items:
            require(key not in value, 'Duplicate JSON key')
            value[key] = item
        return value
    def invalid(value):
        raise ValueError('Nonfinite JSON constant: ' + value)
    return json.loads(raw, object_pairs_hook=unique, parse_constant=invalid)

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def rational(value):
    require(type(value) is str and re.fullmatch(r'-?(?:0|[1-9][0-9]*)/[1-9][0-9]*', value), 'Canonical finite rational required')
    result = Fraction(value)
    require(value == f'{result.numerator}/{result.denominator}', 'Noncanonical rational')
    return result

def real_directory(path):
    path = Path(os.path.abspath(path))
    require(path.is_dir() and path.resolve() == path, 'Real directory required')
    require(all(not parent.is_symlink() for parent in (path, *path.parents)), 'Symlink ancestor forbidden')
    return path

def safe_file(root, name):
    require(type(name) is str and name and '\\' not in name, 'Canonical relative filename required')
    parsed = PurePosixPath(name)
    require(not parsed.is_absolute() and parsed.as_posix() == name and all(part not in ('', '.', '..') for part in name.split('/')), 'Unsafe relative filename')
    path = root
    for part in parsed.parts:
        path = path / part
        require(not path.is_symlink(), 'Payload symlink forbidden')
    require(path.is_file() and stat.S_ISREG(path.stat().st_mode), 'Regular frozen file required: ' + name)
    return path

def zero_counts(value):
    return type(value) is dict and set(value) == ZERO_COUNTER_KEYS and all(type(item) is int and item == 0 for item in value.values())

def validate_contract(contract):
    require(type(contract) is dict, 'Contract object required')
    fixed = {
        'schema_version': 1, 'input_scope': INPUT_SCOPE, 'output_scope': OUTPUT_SCOPE,
        'freeze_authorization': 'NOT_AUTHORIZED_BEFORE_REMOTE_BYTE_VERIFICATION',
        'scientific_outcome_before_freeze': 'UNCOMPUTED_REGISTERED_BD_PREHISTORY',
        'archive_array_decoding_allowed': False, 'original_binary80_state_error': 'NOT_ENCLOSED',
        'full_twelve_case_certificate': 'UNRESOLVED', 'metric_calibration': 'FAIL',
        'higher_dimensional_Big_Bang_origin': 'NOT_ESTABLISHED', 'external_novelty': 'NOT_ASSESSED',
        'domain': ['-5/1', '-9/2'], 'sources': list(SOURCES), 'probe_momenta': list(MOMENTA),
        'target': 'UNIQUE_CONTINUUM_BD_PREHISTORY_AT_FIXED_RATIONAL_PROBES',
        'zero_anchor_eta': '-6/1', 'zero_through_eta': '-5/1',
        'forcing_operator': 'g=4*L^2*h-2*L*h_prime-h_second; L=-1/eta',
        'represented_epsilon': '3777893186295716171/37778931862957161709568',
        'represented_Pi': '14488038916154245685/4611686018427387904',
        'whole_rows': 18, 'complex_rectangles': 36, 'max_complete_exported_L1_radius': '1/100000000000000000000',
        'cap_delta': '1/128', 'geometric_panels': 11, 'max_child_width': '1/64',
        'source_degree': 112, 'coefficient_bits': 512, 'primary_phase_degree': 128,
        'independent_mode_degree': 160, 'independent_exp_degree': 200,
    }
    for key, expected in fixed.items():
        require(type(contract.get(key)) is type(expected) and contract[key] == expected, 'Contract differs: ' + key)
    require(contract.get('cap_formulae') == CAP_FORMULAE, 'Exact per-route cap formulae differ')
    require(contract.get('cap_method_names') == {'primary': 'PRIMARY_FLAT_CAP_L1_DERIVATIVE_BOUND', 'independent': 'INDEPENDENT_FLAT_CAP_IBP_ENDPOINT_BOUND'}, 'Exact per-route cap method names differ')
    require(contract.get('cap_parameters') == {'delta': '1/128', 'D': '255/16384', 'H_base': '3/8', 'H_exponent': 63, 'lambda': '128/639'}, 'Exact cap parameters differ')
    require(rational(contract['represented_Pi']) >= 3, 'Represented Pi lower premise fails')
    require(contract.get('runtime') == {'python': '3.12.14', 'python-flint': '0.9.0', 'sympy': '1.14.0', 'mpmath': '1.3.0'}, 'Exact registered runtime versions differ')
    require(zero_counts(contract.get('new_target_evaluations_before_freeze')), 'Nonprospective contract counters')
    require(contract.get('resources') == {
        'each_route_wall_seconds': 900, 'each_route_peak_rss_kib': 524288,
        'each_route_output_bytes': 20971520, 'nonroot': True, 'single_threads': True,
        'child_processes_allowed': False, 'network_allowed': False,
    }, 'Resource contract differs')
    return contract

def runtime_identity():
    require(sys.version_info[:3] == (3, 12, 14), 'Registered Python 3.12.14 required')
    versions = {name: metadata.version(name) for name in ('python-flint', 'sympy', 'mpmath')}
    require(versions == {'python-flint': '0.9.0', 'sympy': '1.14.0', 'mpmath': '1.3.0'}, 'Registered dependency versions differ')
    return {'sys_version': sys.version, 'machine': platform.machine(), 'platform': sys.platform, 'package_versions': versions}

def verify_remote(item, raw, commit, name, prefix):
    require(type(item) is dict, 'Remote identity object required')
    require(type(item.get('bytes')) is int and item['bytes'] == len(raw) and item.get('sha256') == sha(raw), 'Remote size/SHA256 differs: ' + name)
    expected_git = hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest()
    require(item.get('git_blob') == expected_git and item.get('independent_public_blob_download_sha256_verified') is True, 'Remote Git identity differs: ' + name)
    expected_url = f'https://api.github.com/repos/maldonado-research/HDblast/contents/{prefix}/{name}?ref={commit}'
    require(item.get('immutable_contents_url') == expected_url, 'Immutable commit/path binding differs: ' + name)

def authenticate(root, receipt_path, receipt_sha256, registration_sha256):
    root = real_directory(root)
    require(re.fullmatch('[0-9a-f]{64}', receipt_sha256 or '') and re.fullmatch('[0-9a-f]{64}', registration_sha256 or ''), 'Explicit external SHA256 pins required')
    receipt_path = Path(os.path.abspath(receipt_path))
    require(receipt_path.is_file() and not receipt_path.is_symlink() and all(not p.is_symlink() for p in receipt_path.parents), 'Real external receipt required')
    require(receipt_path.stat().st_size <= 20 * 1024 * 1024, 'Oversized freeze receipt')
    receipt_raw = receipt_path.read_bytes()
    require(sha(receipt_raw) == receipt_sha256, 'Root public receipt hash differs')
    receipt = load(receipt_raw)
    require(type(receipt) is dict and receipt.get('status') == REMOTE_GO and receipt.get('public_repository') == 'https://github.com/maldonado-research/HDblast' and receipt.get('input_scope') == INPUT_SCOPE, 'Unverified freeze scope/status')
    require(receipt.get('all_registered_public_blobs_independently_downloaded') is True and zero_counts(receipt.get('new_target_evaluations_before_public_verification')), 'Preverification evaluation or incomplete public verification')
    commit = receipt.get('freeze_commit')
    require(type(commit) is str and re.fullmatch('[0-9a-f]{40}', commit), 'Exact immutable public commit required')
    require(receipt.get('registration_sha256') == registration_sha256, 'Receipt registration binding differs')
    prefix = receipt.get('repository_checkpoint_path')
    require(prefix == 'research/HDBLAST_CHECKPOINT_20261007_BD_PREHISTORY', 'Registered public path differs')
    registration_raw = safe_file(root, 'FULL_REGISTRATION.json').read_bytes()
    require(sha(registration_raw) == registration_sha256, 'Full registration hash differs')
    registration = load(registration_raw)
    require(type(registration) is dict and set(registration) == {'schema_version', 'files'} and type(registration['schema_version']) is int and registration['schema_version'] == 1, 'Exact registration schema required')
    files = registration['files']
    require(type(files) is dict and REQUIRED <= set(files) and 'FULL_REGISTRATION.json' not in files, 'Incomplete or circular source inventory')
    remote = receipt.get('remote_registered_files')
    require(type(remote) is dict and set(remote) == set(files) | {'FULL_REGISTRATION.json'} and type(receipt.get('remote_files_verified')) is int and receipt['remote_files_verified'] == len(remote), 'Complete immutable remote inventory required')
    for name, item in files.items():
        require(type(item) is dict and set(item) == {'bytes', 'sha256'} and type(item['bytes']) is int and 0 <= item['bytes'] <= 20 * 1024 * 1024 and re.fullmatch('[0-9a-f]{64}', item['sha256'] or ''), 'Invalid registration member')
        path = safe_file(root, name)
        require(path.stat().st_size == item['bytes'], 'Frozen size differs: ' + name)
        raw = path.read_bytes()
        require(sha(raw) == item['sha256'], 'Frozen SHA256 differs: ' + name)
        verify_remote(remote[name], raw, commit, name, prefix)
    verify_remote(remote['FULL_REGISTRATION.json'], registration_raw, commit, 'FULL_REGISTRATION.json', prefix)
    for path in root.rglob('*'):
        require(not path.is_symlink(), 'Checkpoint symlink forbidden')
        if path.is_file():
            require(path.suffix not in ('.pyc', '.pyo', '.so', '.pyd', '.dll') and '__pycache__' not in path.parts, 'Cached/native checkpoint code forbidden')
            if path.suffix == '.py':
                require(path.relative_to(root).as_posix() in files, 'Unregistered Python source')
    contract = validate_contract(load(safe_file(root, 'REGISTRATION_CONTRACT.json').read_bytes()))
    review = load(safe_file(root, 'evidence/preparation/FINAL_PRE_FREEZE_REVIEW.json').read_bytes())
    require(type(review) is dict and review.get('status') == 'GO_FOR_PROSPECTIVE_PUBLIC_FREEZE', 'Independent prefreeze GO required')
    reviewed = review.get('reviewed_source_sha256')
    require(type(reviewed) is dict and {name for name in files if name.endswith('.py')} <= set(reviewed), 'All prospective Python sources must be independently reviewed')
    for name, pin in reviewed.items():
        require(name in files and files[name]['sha256'] == pin, 'Reviewed source changed: ' + name)
    return {'root': root, 'receipt_path': receipt_path, 'receipt': receipt, 'registration': registration, 'contract': contract, 'receipt_sha256': receipt_sha256, 'registration_sha256': registration_sha256}

def authenticate_local_fabricated(root, registration_sha256):
    """Local preparation only: no remote GO, physical authorization or review claim."""
    root = real_directory(root)
    require(re.fullmatch('[0-9a-f]{64}', registration_sha256 or ''), 'External local registration pin required')
    raw = safe_file(root, 'FULL_REGISTRATION.json').read_bytes()
    require(sha(raw) == registration_sha256, 'Local fabricated registration hash differs')
    registration = load(raw)
    require(type(registration) is dict and set(registration) == {'schema_version', 'files'} and type(registration['schema_version']) is int and registration['schema_version'] == 1, 'Exact local registration schema required')
    files = registration['files']
    require(type(files) is dict and (REQUIRED - {'evidence/preparation/FINAL_PRE_FREEZE_REVIEW.json'}) <= set(files) and 'FULL_REGISTRATION.json' not in files, 'Incomplete local source inventory')
    for name, item in files.items():
        require(type(item) is dict and set(item) == {'bytes', 'sha256'} and type(item['bytes']) is int and 0 <= item['bytes'] <= 20 * 1024 * 1024 and re.fullmatch('[0-9a-f]{64}', item['sha256'] or ''), 'Invalid local registration member')
        path = safe_file(root, name)
        require(path.stat().st_size == item['bytes'] and sha(path.read_bytes()) == item['sha256'], 'Local frozen source differs: ' + name)
    for path in root.rglob('*'):
        require(not path.is_symlink(), 'Local checkpoint symlink forbidden')
        if path.is_file():
            require(path.suffix not in ('.pyc', '.pyo', '.so', '.pyd', '.dll') and '__pycache__' not in path.parts, 'Cached/native local checkpoint code forbidden')
            if path.suffix == '.py':
                require(path.relative_to(root).as_posix() in files, 'Unregistered local Python source')
    contract = validate_contract(load(safe_file(root, 'REGISTRATION_CONTRACT.json').read_bytes()))
    return {'root': root, 'registration': registration, 'contract': contract, 'registration_sha256': registration_sha256,
            'scope': 'LOCAL_FABRICATED_ONLY_NO_PUBLIC_FREEZE'}

def expected_centers():
    left = Fraction(1, 128)
    centers = []
    while left < Fraction(1, 2):
        right = min(3 * left / 2, Fraction(1, 2))
        centers.append((left + right) / 2 - 5)
        left = right
    require(len(centers) == 11, 'Wrong fixed parent geometry')
    return set(centers)

class SourceAuthorization:
    def __init__(self, verified, route, journal):
        require(route in ('primary', 'independent'), 'Known route required')
        self.verified = verified
        self.route = route
        self.journal = Path(journal)
        self.journal.open('xb').close()
        self.events = []
        self.seen = set()

    def __call__(self, event, detail):
        require(event == 'real_source_taylor' and type(detail) is dict and set(detail) == {'source', 'center', 'scope'}, 'Unexpected physical callback')
        require(detail['scope'] == 'UNIQUE_ANALYTIC_BD_PREHISTORY_ONLY' and detail['source'] in SOURCES, 'Callback target differs')
        center = rational(detail['center'])
        require(center in expected_centers(), 'Callback outside fixed geometric parents')
        key = (detail['source'], center)
        require(key not in self.seen, 'Repeated physical source construction')
        # Recheck frozen source bytes before every permitted real callback.
        v = self.verified
        authenticate(v['root'], v['receipt_path'], v['receipt_sha256'], v['registration_sha256'])
        record = {'index': len(self.events) + 1, 'route': self.route, **detail}
        with self.journal.open('ab') as handle:
            handle.write((json.dumps(record, sort_keys=True) + '\n').encode())
            handle.flush()
            os.fsync(handle.fileno())
        self.seen.add(key)
        self.events.append(record)

    def complete(self):
        require(self.seen == {(source, center) for source in SOURCES for center in expected_centers()}, 'Incomplete registered source construction')
        return {'route': self.route, 'physical_source_callbacks': 22, 'archive_arrays_decoded': 0, 'events': self.events}
