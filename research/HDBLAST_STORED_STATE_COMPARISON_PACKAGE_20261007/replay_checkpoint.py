#!/usr/bin/env python3
"""Authenticate the whole portable package before importing any study code.

External SHA256 values for this helper, PAYLOAD_MANIFEST and REPLAY_PINS must
come from the independently sealed release, not from the extracted package.
No network, extraction, native-array loading or public-GO creation occurs here.
"""
from pathlib import Path, PurePosixPath
import argparse
import hashlib
import importlib.util
import json
import math
import os
import re
import stat
import subprocess
import signal
import sys

sys.dont_write_bytecode = True
PREFIX = 'research/HDBLAST_CHECKPOINT_20261007_STORED_STATE_COMPARISON'
RESERVED = {'replay_checkpoint.py', 'PAYLOAD_MANIFEST.json', 'REPLAY_PINS.json'}
SCIENCE = ('DATA.json', 'SCIENCE_SUMMARY.json', 'SOURCE_CERTIFICATE.json',
           'DIAGNOSTIC_CROSSCHECK.json', 'NODE_TARGETS.jsonl.gz',
           'SOURCE_ATTEMPTS.jsonl', 'DECODE_ATTEMPTS.jsonl')
STABLE = SCIENCE + ('OUTPUT_READBACK.json',)
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
                            Path(name).suffix not in ('.pyc', '.pyo', '.so', '.pyd', '.dll'),
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

def validate_pins(root, rows, pins, manifest_sha256, helper_sha256):
    keys = {'schema_version', 'sealed', 'checkpoint_path', 'public_go_path',
            'registration_sha256', 'public_go_sha256', 'payload_manifest_sha256',
            'helper_sha256', 'runtime', 'saved_runs', 'stable_files'}
    require(type(pins) is dict and set(pins) == keys and type(pins['schema_version']) is int and
            pins['schema_version'] == 1 and pins['sealed'] is True, 'FINAL_SEALED_REPLAY_PINS_REQUIRED')
    require(pins['checkpoint_path'] == PREFIX and pins['payload_manifest_sha256'] == manifest_sha256 and
            pins['helper_sha256'] == helper_sha256 and pins['runtime'] == RUNTIME and
            pins['stable_files'] == list(STABLE), 'REPLAY_PIN_BINDING_OR_POLICY_DIFFERS')
    for key in ('registration_sha256', 'public_go_sha256'):
        digest(pins[key])
    go = relative(pins['public_go_path']).as_posix()
    require(go in rows and go != PREFIX and not go.startswith(PREFIX + '/'),
            'PUBLIC_GO_MUST_BE_REGISTERED_OUTSIDE_CHECKPOINT')
    require(rows[go]['sha256'] == pins['public_go_sha256'] and
            rows[PREFIX + '/FULL_REGISTRATION.json']['sha256'] == pins['registration_sha256'],
            'REGISTRATION_OR_PUBLIC_GO_MANIFEST_BINDING_DIFFERS')
    require(type(pins['saved_runs']) is dict and set(pins['saved_runs']) == {'normal', 'optimized'},
            'BOTH_SAVED_MODES_REQUIRED')
    paths = [relative(pins['saved_runs'][mode]).as_posix() for mode in ('normal', 'optimized')]
    require(paths[0] != paths[1] and all(not path.startswith(PREFIX + '/') for path in paths),
            'DISTINCT_SAVED_OUTPUTS_OUTSIDE_CHECKPOINT_REQUIRED')
    for path in paths:
        for name in STABLE + ('RESOURCE.json', 'ENTRY_RECEIPT.json', 'EXECUTION.json', 'child.log'):
            require(path + '/' + name in rows, 'COMPLETE_SAVED_EXECUTION_REQUIRED:' + name)
    for name in STABLE:
        require(rows[paths[0] + '/' + name] == rows[paths[1] + '/' + name],
                'SAVED_NORMAL_OPTIMIZED_SCIENCE_PINS_DIFFER:' + name)
    contract = load((root / PREFIX / 'REGISTRATION_CONTRACT.json').read_bytes())
    require(json.dumps(contract.get('resources'), sort_keys=True, allow_nan=False) ==
            json.dumps(LIMITS, sort_keys=True, allow_nan=False), 'REGISTERED_FIXED_RESOURCE_POLICY_DIFFERS')
    spec = load((root / PREFIX / 'provenance/INPUT_SPEC.json').read_bytes())
    required = [(spec['manifest_repository_path'], spec['manifest_sha256']),
                (spec['producer_repository_path'], spec['producer_sha256']),
                (spec['prior_static_audit']['repository_path'], spec['prior_static_audit']['sha256'])]
    require(type(spec['capsules']) is list and len(spec['capsules']) == 4, 'FOUR_CAPSULES_REQUIRED')
    for capsule in spec['capsules']:
        required.extend([(capsule['repository_path'], capsule['reader_inputspec']['file_sha256']),
                         (capsule['original']['repository_path'], capsule['original']['sha256'])])
    require(len({name for name, pin in required}) == 11, 'ALL_EIGHT_INPUTS_AND_THREE_PROVENANCE_FILES_REQUIRED')
    for name, pin in required:
        relative(name); digest(pin)
        require(name in rows and rows[name]['sha256'] == pin, 'INPUT_OR_PROVENANCE_PAYLOAD_BINDING_DIFFERS:' + name)

def reauthenticate_package(v):
    return authenticate_package(v['root'], *v['external'])

def study_module(name, path):
    existing = sys.modules.get(name)
    require(existing is None or (isinstance(getattr(existing, '__file__', None), str) and
                                  Path(existing.__file__).absolute() == path.absolute()),
            'FOREIGN_STUDY_MODULE_ALREADY_IMPORTED:' + name)
    if existing is not None:
        return existing
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec); sys.modules[name] = module
    spec.loader.exec_module(module)
    return module

def validate_custodian(output, mode):
    require(mode in ('normal','optimized'), 'REGISTERED_REPLAY_MODE_REQUIRED')
    output=directory(output)
    fd=os.open(output,os.O_RDONLY|os.O_DIRECTORY|os.O_NOFOLLOW|os.O_CLOEXEC)
    try:
        execution_pin,raw=read_member(fd,'EXECUTION.json',capture=True)
        require(execution_pin['bytes']<=65536,'CUSTODIAN_RECEIPT_LIMIT')
        receipt=load(raw)
        entry_pin,unused=read_member(fd,'ENTRY_RECEIPT.json')
        log_pin,unused=read_member(fd,'child.log')
    finally:os.close(fd)
    require(type(receipt) is dict and receipt.get('status')=='PASS_BOUNDED_REGISTERED_EXECUTION' and
            receipt.get('fabricated_only') is False and receipt.get('optimized') is (mode=='optimized'),
            'SAVED_OR_FRESH_CUSTODY_SCOPE_OR_MODE_DIFFERS:'+mode)
    require(all(type(receipt.get(key)) is int and receipt[key]==0
                for key in ('exit_code','wrapper_exit_code','wait_status')) and
            receipt.get('stop_reason') is None and receipt.get('monitor_error') is None,
            'CUSTODIAN_WAIT4_OUTCOME_NOT_SUCCESS:'+mode)
    limits={'wall_seconds':900,'cpu_seconds':900,'address_space_bytes':536870912,
            'aggregate_output_bytes':134217728,'file_size_bytes':134217728,
            'process_limit':0,'core_bytes':0,'custodian_receipt_reserve_bytes':65536}
    require(json.dumps(receipt.get('limits'),sort_keys=True,allow_nan=False)==
            json.dumps(limits,sort_keys=True,allow_nan=False),'CUSTODIAN_LIMITS_DIFFER:'+mode)
    require(receipt.get('entry_receipt')==entry_pin and
            receipt.get('captured_child_log_sha256')==log_pin['sha256'] and
            receipt.get('observed_child_log_sha256')==log_pin['sha256'] and
            receipt.get('preserved_child_log') is None,'CUSTODIAN_ENTRY_OR_LOG_LINKAGE_DIFFERS:'+mode)
    require(type(receipt.get('uid')) is int and receipt['uid']>0 and
            type(receipt.get('peak_rss_kib_wait4')) is int and
            0<receipt['peak_rss_kib_wait4']*1024<=536870912,'CUSTODIAN_UID_OR_RSS_DIFFERS:'+mode)
    for key in ('wall_seconds','cpu_seconds_wait4'):
        value=receipt.get(key)
        require(type(value) in (int,float) and math.isfinite(value) and 0<=value<=900,
                'CUSTODIAN_MEASURED_RESOURCE_DIFFERS:'+mode)
    summary=receipt.get('worker_output')
    require(type(summary) is dict and type(summary.get('bytes')) is int and
            0<summary['bytes']<=134217728-65536 and type(summary.get('file_count')) is int and
            summary['file_count']>=11,'CUSTODIAN_AGGREGATE_ACCEPTANCE_DIFFERS:'+mode)
    return {'status':'PASS_AUTHENTICATED_CUSTODIAN_OUTCOME','mode':mode,
            'execution_receipt':execution_pin,'entry_receipt':entry_pin,'child_log':log_pin}

def registered_command(v, output, mode):
    require(mode in ('normal', 'optimized'), 'REGISTERED_REPLAY_MODE_REQUIRED')
    pins=v['pins'];candidate=v['root']/PREFIX
    command=[sys.executable, '-I', '-B', str(candidate/'execution/bounded_launcher.py'),
             '--candidate', str(candidate), '--output', str(output),
             '--registration-sha256', pins['registration_sha256'],
             '--public-go', str(v['root']/pins['public_go_path']),
             '--public-go-sha256', pins['public_go_sha256'],
             '--repository-root', str(v['root'])]
    if mode=='optimized':command.append('--optimized')
    return command

def launch_registered(v, output, mode, logs):
    # Use a fresh stdlib-only custodian process. This helper's preceding saved
    # readback may have imported numerical libraries, so it must not fork the
    # bounded worker itself through an in-process launcher import.
    command=registered_command(v, output, mode)
    with (logs/(mode+'-custodian.stdout.log')).open('xb') as stdout, \
         (logs/(mode+'-custodian.stderr.log')).open('xb') as stderr:
        process=subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=stdout, stderr=stderr,
                                 env={'PATH':'/usr/bin:/bin','LANG':'C.UTF-8','LC_ALL':'C.UTF-8'})
        try:
            code=process.wait()
        except BaseException:
            # SIGINT lets the Python custodian run its child termination/reap
            # cleanup, including the worker's separate process group.
            process.send_signal(signal.SIGINT)
            try:process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                process.kill();process.wait()
            raise
    require(code==0, 'BOUNDED_REPLAY_PROCESS_FAILED:'+mode)
    path=output/'EXECUTION.json'
    require(path.is_file() and not path.is_symlink() and path.stat().st_size<=65536,
            'BOUNDED_EXECUTION_RECEIPT_REQUIRED')
    validate_custodian(output,mode)
    receipt=load(path.read_bytes())
    require(receipt.get('status')=='PASS_BOUNDED_REGISTERED_EXECUTION' and
            receipt.get('fabricated_only') is False and receipt.get('optimized') is (mode=='optimized') and
            type(receipt.get('wrapper_exit_code')) is int and receipt['wrapper_exit_code']==0,
            'BOUNDED_REPLAY_FAILED:'+mode)
    return receipt

def replay(v, verify_only, output_root=None):
    """Only called after complete package byte authentication."""
    require(sys.flags.isolated, 'USE_PYTHON_ISOLATED_MODE_-I')
    reauthenticate_package(v)
    pins = v['pins']; candidate = v['root'] / PREFIX; execution = candidate / 'execution'
    sys.path.insert(0, str(execution))
    try:
        guard = study_module('registration_guard', execution / 'registration_guard.py')
        verified = guard.authenticate(candidate, pins['registration_sha256'],
                                      v['root'] / pins['public_go_path'], pins['public_go_sha256'])
        saved = {mode: v['root'] / pins['saved_runs'][mode] for mode in ('normal', 'optimized')}
        custody = {mode: validate_custodian(path,mode) for mode,path in saved.items()}
        reader = study_module('readback', execution / 'readback.py')
        checks = {mode: reader.validate_output(path, candidate, fabricated=False, verified=verified)
                  for mode, path in saved.items()}
        comparison = reader.compare_replays(saved['normal'], saved['optimized'])
        guard.reauthenticate(verified); reauthenticate_package(v)
        if verify_only:
            return {'status': 'PASS_VERIFIED_SAVED_NORMAL_OPTIMIZED_PACKAGE', 'new_array_decodes': 0,
                    'new_source_callbacks': 0, 'saved_custody': custody, 'saved_checks': checks, 'saved_comparison': comparison}
        require(output_root is not None and os.getuid() > 0 and os.geteuid() > 0,
                'FRESH_REPLAY_REQUIRES_NONROOT_AND_EXTERNAL_OUTPUT')
        out = Path(os.path.abspath(output_root)); directory(out.parent)
        require(out != v['root'] and v['root'] not in out.parents and out not in v['root'].parents,
                'REPLAY_OUTPUT_MUST_BE_OUTSIDE_PACKAGE')
        os.mkdir(out, 0o700)
        receipts = {}; fresh_checks = {}
        for mode in ('normal', 'optimized'):
            reauthenticate_package(v); guard.reauthenticate(verified)
            run = out / mode
            receipt = launch_registered(v, run, mode, out)
            require(receipt['status'] == 'PASS_BOUNDED_REGISTERED_EXECUTION' and receipt['wrapper_exit_code'] == 0,
                    'BOUNDED_REPLAY_FAILED:' + mode)
            receipts[mode] = receipt
            fresh_checks[mode] = reader.validate_output(run, candidate, fabricated=False, verified=verified)
            reader.compare_replays(saved['normal'], run)
            guard.reauthenticate(verified); reauthenticate_package(v)
        fresh_comparison = reader.compare_replays(out / 'normal', out / 'optimized')
        result = {'status': 'PASS_COMPLETE_PORTABLE_STORED_STATE_COMPARISON_REPLAY',
                  'saved_custody': custody, 'saved_checks': checks, 'saved_comparison': comparison, 'fresh_checks': fresh_checks,
                  'fresh_comparison': fresh_comparison, 'exact_scientific_files': list(SCIENCE),
                  'other_exact_file': 'OUTPUT_READBACK.json', 'resource_time_RSS_are_nonidentity': True,
                  'execution_receipts': {mode: str(out / mode / 'EXECUTION.json') for mode in receipts}}
        (out / 'PORTABLE_REPLAY_RECEIPT.json').write_text(json.dumps(result, sort_keys=True, indent=2) + '\n')
        return result
    finally:
        sys.path.remove(str(execution))

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package', required=True)
    parser.add_argument('--payload-manifest-sha256', required=True)
    parser.add_argument('--helper-sha256', required=True)
    parser.add_argument('--replay-pins-sha256', required=True)
    parser.add_argument('--verify-only', action='store_true')
    parser.add_argument('--output-root')
    args = parser.parse_args(argv)
    require(args.verify_only == (args.output_root is None), 'VERIFY_ONLY_OR_FRESH_OUTPUT_REQUIRED')
    v = authenticate_package(args.package, args.payload_manifest_sha256,
                             args.helper_sha256, args.replay_pins_sha256)
    result = replay(v, args.verify_only, args.output_root)
    print(json.dumps(result, sort_keys=True, separators=(',', ':'), allow_nan=False))
    return 0

if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (Stop, ValueError, OSError, KeyError) as exc:
        print(json.dumps({'status': 'STOPPED', 'reason': str(exc)}), file=sys.stderr)
        raise SystemExit(2)
