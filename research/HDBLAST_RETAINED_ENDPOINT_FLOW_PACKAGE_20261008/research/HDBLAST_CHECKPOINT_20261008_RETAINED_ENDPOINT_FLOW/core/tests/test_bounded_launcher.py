"""Manufactured-only kernel/launcher checks; no archive or physical source access.

The isolated test harness replaces candidate authentication with a stub for
its own temporary, self-written probe worker. Production authentication is
unchanged. A lowered wall/CPU threshold exercises timeout handling quickly;
the normal probe separately verifies the exact registered hard limits.
"""
from pathlib import Path
from contextlib import nullcontext
import errno
import json
import os
import shutil
import signal
import stat
import sys
import tempfile
import time
import unittest
from unittest import mock

sys.dont_write_bytecode = True
import types
_launcher_path = Path(__file__).resolve().parents[1] / 'execution' / 'bounded_launcher.py'
bounded_launcher = types.ModuleType('bounded_launcher')
bounded_launcher.__file__ = str(_launcher_path)
sys.modules['bounded_launcher'] = bounded_launcher
exec(compile(_launcher_path.read_bytes(), str(_launcher_path), 'exec'), bounded_launcher.__dict__)
kernel_guard = sys.modules['kernel_guard']

PROBE = r'''
from pathlib import Path
import argparse, ctypes, errno, hashlib, json, os, resource, signal, socket, sys, time
from decimal import Decimal
from fractions import Fraction
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kernel_guard
parser = argparse.ArgumentParser()
parser.add_argument('--root', required=True)
parser.add_argument('--registration', required=True)
parser.add_argument('--output', required=True)
parser.add_argument('--registration-sha256', required=True)
parser.add_argument('--manufactured-only', action='store_true', required=True)
parser.add_argument('--probe-mode', required=True)
args = parser.parse_args()
output = Path(args.output)
mode = args.probe_mode
if mode != 'cpu':
    kernel_guard.require_worker_limits()
limits = {name: list(resource.getrlimit(getattr(resource, 'RLIMIT_' + name)))
          for name in ('AS', 'CPU', 'FSIZE', 'NPROC', 'CORE')}
libc = ctypes.CDLL(None, use_errno=True)
libc.syscall.restype = ctypes.c_long
seccomp = ctypes.CDLL('libseccomp.so.2', use_errno=True)
seccomp.seccomp_syscall_resolve_name.argtypes = [ctypes.c_char_p]
seccomp.seccomp_syscall_resolve_name.restype = ctypes.c_int
numbers = {name: seccomp.seccomp_syscall_resolve_name(name.encode())
           for name in kernel_guard.DENIED_SYSCALLS}
preexisting_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
denied = kernel_guard.install()
print('MANUFACTURED_PROBE_DIAGNOSTIC ' + mode, flush=True)
data = {'mode': mode, 'limits': limits, 'uid': os.getuid(), 'denied_syscalls': denied,
        'optimized': sys.flags.optimize, 'stdin_eof': os.read(0, 1) == b'',
        'credentials_absent': all(name not in os.environ for name in
                                  ('GITHUB_TOKEN', 'AWS_SECRET_ACCESS_KEY', 'PROBE_SECRET')),
        'single_thread_env': all(os.environ.get(name) == '1'
                                 for name in kernel_guard.THREAD_VARIABLES)}
if mode == 'security':
    observations = {}
    # Invoke every resolvable denied call directly with inert/invalid arguments.
    # EPERM therefore proves kernel denial, rather than an API wrapper or NPROC.
    for name, number in numbers.items():
        if number >= 0:
            ctypes.set_errno(0)
            result = libc.syscall(ctypes.c_long(number), ctypes.c_long(-1),
                                  ctypes.c_long(0), ctypes.c_long(0),
                                  ctypes.c_long(0), ctypes.c_long(0), ctypes.c_long(0))
            observations[name] = {'return': result, 'errno': ctypes.get_errno()}
    try:
        preexisting_socket.connect(('127.0.0.1', 9))
        observations['python_connect'] = 'UNEXPECTED_SUCCESS'
    except OSError as exc:
        observations['python_connect'] = exc.errno
    data['kernel_observations'] = observations
elif mode == 'failure':
    raise SystemExit(7)
elif mode == 'signal':
    os.kill(os.getpid(), signal.SIGTERM)
elif mode == 'wall':
    time.sleep(10)
elif mode == 'cpu':
    while True:
        pass
elif mode == 'symlink':
    (output / 'unsafe-link').symlink_to('/etc/passwd')
elif mode == 'fifo':
    os.mkfifo(output / 'unsafe-pipe')
elif mode == 'hardlink':
    (output / 'original').write_text('manufactured')
    os.link(output / 'original', output / 'alias')
elif mode == 'aggregate':
    for name in ('part-one', 'part-two'):
        with (output / name).open('xb') as handle:
            handle.truncate(70 * 1024 * 1024)
elif mode == 'fsize':
    signal.signal(signal.SIGXFSZ, signal.SIG_DFL)
    with (output / 'oversized').open('xb') as handle:
        handle.truncate(kernel_guard.OUTPUT_BYTES + 1)
elif mode == 'as':
    try:
        allocation = bytearray(kernel_guard.ADDRESS_SPACE_BYTES * 2)
        data['allocation'] = 'UNEXPECTED_SUCCESS'
    except MemoryError:
        data['allocation'] = 'DENIED_MEMORY_ERROR'
elif mode == 'log-removal':
    (output / 'child.log').unlink()
elif mode == 'log-overwrite':
    (output / 'child.log').write_text('manufactured log contamination\n')
elif mode == 'receipt-tamper':
    (output / 'EXECUTION.json').symlink_to('/etc/passwd')
elif mode == 'directory-replacement':
    renamed = output.with_name(output.name + '-renamed')
    output.rename(renamed)
    output.mkdir()
    output = renamed
(output / 'PROBE.json').write_text(json.dumps(data, sort_keys=True) + '\n')

if mode != 'false-zero':
    core = ('SOURCE_CERTIFICATE.json', 'TARGET_BUDGETS.json', 'NODE_CERTIFICATES.jsonl.gz',
            'DATA.json', 'SCIENCE_SUMMARY.json', 'RESOURCE.json', 'SOURCE_ATTEMPTS.jsonl',
            'DECODE_ATTEMPTS.jsonl', 'INPUT_AUTHENTICATION.json', 'OUTPUT_READBACK.json')
    manifest = {}
    for name in core:
        raw = b'manufactured custodian-contract fixture only\n'
        if name == 'OUTPUT_READBACK.json':
            readback={'status':'PASS_STANDALONE_ENDPOINT_READBACK','nodes':49152,
                      'endpoint_comparisons':98304,'prefix_cases':24}
            if mode == 'readback-status-tamper':readback['status']='PASS_UNRELATED'
            if mode == 'readback-count-tamper':readback['prefix_cases']=True
            raw=(json.dumps(readback,sort_keys=True)+'\n').encode()
        (output / name).write_bytes(raw)
        manifest[name] = {'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}
    entry = {'status':'PASS_MANUFACTURED_ENDPOINT_ENTRY', 'manufactured':True,
             'registration_sha256':args.registration_sha256, 'go_sha256':None,
             'scope':'LATER_RETAINED_ENDPOINT_FLOW_CERTIFICATE', 'nodes':49152,
             'endpoint_comparisons':98304, 'outputs':manifest,
             'kernel_denied_syscalls':denied}
    if mode == 'entry-hash-tamper':
        entry['outputs']['DATA.json']['sha256'] = 'f' * 64
    elif mode == 'entry-status-tamper':
        entry['status'] = 'PASS_REGISTERED_ENDPOINT_ENTRY'
    elif mode == 'entry-no-kernel':
        entry['kernel_denied_syscalls'] = []
    elif mode == 'entry-bool-count':
        entry['nodes'] = True
    (output / 'ENTRY_RECEIPT.json').write_text(json.dumps(entry,sort_keys=True)+'\n')
'''


class ManufacturedLauncherTests(unittest.TestCase):
    def setUp(self):
        self.artifacts=Path(os.environ.get('HDBLAST_ENDPOINT_TEST_OUTPUT',str(Path.cwd())))
        self.artifacts.mkdir(parents=True,exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(prefix='manufactured-launcher-', dir=self.artifacts)
        self.base = Path(self.temporary.name)
        self.candidate = self.base / 'candidate'
        execution = self.candidate / 'execution'
        execution.mkdir(parents=True)
        (execution / 'run_endpoint.py').write_text(PROBE)
        shutil.copyfile(Path(kernel_guard.__file__), execution / 'kernel_guard.py')
        self.counter = 0

    def tearDown(self):
        # Preserve the custodian's small receipt/log/probe evidence before
        # removing self-manufactured hazards (sparse quota probes, links/FIFOs).
        evidence=self.artifacts/('preserved-'+self.base.name)
        evidence.mkdir()
        for directory,unused,filenames in os.walk(self.base,followlinks=False):
            for name in filenames:
                path=Path(directory)/name
                metadata=path.lstat()
                if (stat.S_ISREG(metadata.st_mode) and metadata.st_size<=1024*1024 and
                        (name.endswith('.log') or name in ('EXECUTION.json','CUSTODIAN_FAILURE.json','PROBE.json'))):
                    target=evidence/path.relative_to(self.base)
                    target.parent.mkdir(parents=True,exist_ok=True)
                    shutil.copyfile(path,target)
        self.temporary.cleanup()

    def run_probe(self, mode, optimized=False, **overrides):
        self.counter += 1
        output = self.base / ('output-' + str(self.counter))
        with mock.patch.object(bounded_launcher, '_authenticate_candidate',
                               return_value={'scope': 'MANUFACTURED_TEST_HARNESS_ONLY'}):
            with mock.patch.dict(os.environ, {'PROBE_SECRET': 'manufactured-not-a-secret',
                                         'GITHUB_TOKEN': 'manufactured-not-a-token',
                                         'AWS_SECRET_ACCESS_KEY': 'manufactured-not-a-key'}):
                with mock.patch.multiple(bounded_launcher, **overrides) if overrides else nullcontext():
                    receipt = bounded_launcher.launch(
                        self.candidate, output,
                        ['--registration', str(self.base / 'manufactured-registration.json'),
                         '--registration-sha256', '0' * 64, '--probe-mode', mode],
                        optimized=optimized, manufactured=True)
        return receipt, output

    def test_exact_limits_environment_and_optimization(self):
        for optimized in (False, True):
            receipt, output = self.run_probe('success', optimized)
            self.assertEqual(receipt['status'], 'PASS_BOUNDED_MANUFACTURED_EXECUTION')
            data = json.loads((output / 'PROBE.json').read_text())
            self.assertEqual(data['limits'], {'AS': [536870912] * 2, 'CPU': [900] * 2,
                             'FSIZE': [134217728] * 2, 'NPROC': [0, 0], 'CORE': [0, 0]})
            self.assertTrue(data['credentials_absent'])
            self.assertTrue(data['stdin_eof'])
            self.assertTrue(data['single_thread_env'])
            self.assertEqual(data['optimized'], int(optimized))
            self.assertLessEqual(kernel_guard.check_output(output)['bytes'], kernel_guard.OUTPUT_BYTES)

    def test_kernel_network_process_exec_and_io_uring_denial(self):
        receipt, output = self.run_probe('security')
        self.assertEqual(receipt['wrapper_exit_code'], 0)
        data = json.loads((output / 'PROBE.json').read_text())
        observations = data['kernel_observations']
        for name in data['denied_syscalls']:
            self.assertEqual(observations[name], {'return': -1, 'errno': errno.EPERM}, name)
        self.assertEqual(observations['python_connect'], errno.EPERM)

    def test_authoritative_exit_and_failure_log(self):
        for mode, code, wrapper in (('failure', 7, 7), ('signal', -signal.SIGTERM, 143)):
            receipt, output = self.run_probe(mode)
            self.assertEqual(receipt['exit_code'], code)
            self.assertEqual(receipt['wrapper_exit_code'], wrapper)
            self.assertEqual(os.waitstatus_to_exitcode(receipt['wait_status']), code)
            self.assertIn('MANUFACTURED_PROBE_DIAGNOSTIC', (output / 'child.log').read_text())

    def test_wall_timeout_and_cpu_limit_are_reaped(self):
        receipt, output = self.run_probe('wall', WALL_SECONDS=0.25)
        self.assertEqual(receipt['stop_reason'], 'WALL_TIMEOUT')
        self.assertEqual(receipt['exit_code'], -signal.SIGKILL)
        self.assertLess(receipt['wall_seconds'], 2)
        receipt, output = self.run_probe('cpu', CPU_SECONDS=1)
        self.assertEqual(receipt['exit_code'], -signal.SIGKILL)
        self.assertGreaterEqual(receipt['cpu_seconds_wait4'], 0.8)

    def test_address_space_and_file_size_limits(self):
        receipt, output = self.run_probe('as')
        self.assertEqual(receipt['wrapper_exit_code'], 0)
        self.assertEqual(json.loads((output / 'PROBE.json').read_text())['allocation'],
                         'DENIED_MEMORY_ERROR')
        receipt, output = self.run_probe('fsize')
        self.assertEqual(receipt['exit_code'], -signal.SIGXFSZ)
        self.assertIn('MANUFACTURED_PROBE_DIAGNOSTIC', (output / 'child.log').read_text())

    def test_aggregate_cap_rejects_individually_legal_files(self):
        receipt, output = self.run_probe('aggregate')
        self.assertEqual(receipt['status'], 'EXECUTION_FAILED')
        self.assertEqual(receipt['stop_reason'], 'OUTPUT_LIMIT')
        self.assertLessEqual((output / 'part-one').stat().st_size, kernel_guard.OUTPUT_BYTES)
        self.assertLessEqual((output / 'part-two').stat().st_size, kernel_guard.OUTPUT_BYTES)
        with self.assertRaises(ValueError):
            kernel_guard.check_output(output)

    def test_output_links_special_files_and_custodian_tampering(self):
        for mode in ('symlink', 'fifo', 'hardlink', 'log-removal', 'log-overwrite',
                     'receipt-tamper', 'directory-replacement'):
            receipt, output = self.run_probe(mode)
            self.assertEqual(receipt['status'], 'EXECUTION_FAILED', mode)
            self.assertEqual(receipt['stop_reason'], 'UNSAFE_OUTPUT', mode)
            if mode == 'log-removal':
                self.assertIn('MANUFACTURED_PROBE_DIAGNOSTIC',
                              (output / receipt['preserved_child_log']).read_text())
            if mode == 'receipt-tamper':
                self.assertTrue((output / 'CUSTODIAN_FAILURE.json').is_file())
                self.assertTrue((output / 'EXECUTION.json').is_symlink())

    def test_timeout_before_process_group_setup_kills_direct_child(self):
        with mock.patch.object(bounded_launcher.os, 'setsid',
                               side_effect=lambda: time.sleep(5)):
            receipt, output = self.run_probe('success', WALL_SECONDS=0.15)
        self.assertEqual(receipt['stop_reason'], 'WALL_TIMEOUT')
        self.assertEqual(receipt['exit_code'], -signal.SIGKILL)
        self.assertLess(receipt['wall_seconds'], 1)

    def test_unexpected_monitor_exception_terminates_and_reaps(self):
        original = bounded_launcher.check_output
        calls = []
        def fail_once(*args, **kwargs):
            if not calls:
                calls.append(True)
                raise RuntimeError('manufactured monitor fault')
            return original(*args, **kwargs)
        with mock.patch.object(bounded_launcher, 'check_output', side_effect=fail_once):
            with mock.patch.object(bounded_launcher, '_terminate',
                                   wraps=bounded_launcher._terminate) as terminate:
                receipt, output = self.run_probe('wall')
                pid = terminate.call_args.args[0]
        self.assertEqual(receipt['stop_reason'], 'CUSTODIAN_ERROR')
        self.assertEqual(receipt['exit_code'], -signal.SIGKILL)
        with self.assertRaises(ChildProcessError):
            os.waitpid(pid, os.WNOHANG)
        self.assertTrue((output / 'child.log').is_file())

    def test_pinned_output_descriptor_reusable_and_wrong_root_rejected(self):
        root = self.base / 'held-output'
        replacement = self.base / 'replacement-output'
        root.mkdir()
        replacement.mkdir()
        (root / 'original.bin').write_bytes(b'held original output bytes')
        (replacement / 'replacement.bin').write_bytes(b'replacement')
        unused, held_fd = kernel_guard.open_real_directory(root)
        unused, wrong_fd = kernel_guard.open_real_directory(replacement)
        try:
            first = kernel_guard.check_output(root, hash_files=True, directory_fd=held_fd)
            second = kernel_guard.check_output(root, hash_files=True, directory_fd=held_fd)
            self.assertEqual(first, second)
            self.assertEqual(first['bytes'], len(b'held original output bytes'))
            self.assertEqual(first['file_count'], 1)
            self.assertEqual(os.fstat(held_fd).st_ino, root.stat().st_ino)
            with self.assertRaisesRegex(ValueError, 'does not name held root'):
                kernel_guard.check_output(root, directory_fd=wrong_fd)
            self.assertEqual(os.fstat(wrong_fd).st_ino, replacement.stat().st_ino)
        finally:
            os.close(wrong_fd)
            os.close(held_fd)

    def test_transient_root_swap_cannot_substitute_scanned_directory(self):
        root = self.base / 'held-output'
        replacement = self.base / 'replacement-output'
        renamed = self.base / 'held-output-renamed'
        root.mkdir()
        replacement.mkdir()
        (root / 'original.bin').write_bytes(b'original output')
        (replacement / 'replacement.bin').write_bytes(b'x')
        original_open = kernel_guard.open_real_directory
        unused, held_fd = original_open(root)

        def swap_open_and_restore(path):
            # Both bracketing pathname checks see the original directory, but
            # this root-open boundary briefly returns the replacement inode.
            # A scanner that chooses its root from this fresh pathname would
            # inspect the wrong tree and could accept a smaller replacement.
            self.assertEqual(Path(path), root)
            root.rename(renamed)
            replacement.rename(root)
            try:
                return original_open(root)
            finally:
                root.rename(replacement)
                renamed.rename(root)

        try:
            bounded_launcher._same_output_directory(held_fd, root)
            with mock.patch.object(kernel_guard, 'open_real_directory',
                                   side_effect=swap_open_and_restore):
                with self.assertRaisesRegex(ValueError, 'does not name held root'):
                    kernel_guard.check_output(root, hash_files=True, directory_fd=held_fd)
            bounded_launcher._same_output_directory(held_fd, root)
            self.assertEqual(os.fstat(held_fd).st_ino, root.stat().st_ino)
            self.assertEqual(kernel_guard.check_output(root, directory_fd=held_fd)['bytes'],
                             len(b'original output'))
        finally:
            os.close(held_fd)

    def test_launcher_pins_live_and_final_inspections_to_held_descriptor(self):
        with mock.patch.object(bounded_launcher, 'check_output',
                               wraps=bounded_launcher.check_output) as inspected:
            receipt, output = self.run_probe('success')
        self.assertEqual(receipt['wrapper_exit_code'], 0)
        self.assertGreaterEqual(inspected.call_count, 2)
        held = {call.kwargs.get('directory_fd') for call in inspected.call_args_list}
        self.assertEqual(len(held), 1)
        self.assertTrue(all(type(fd) is int for fd in held))
        self.assertTrue(any(call.kwargs.get('hash_files') for call in inspected.call_args_list))

    def test_output_must_be_fresh_external_real_directory(self):
        with mock.patch.object(bounded_launcher, '_authenticate_candidate', return_value={}):
            for output in (self.candidate / 'inside', self.base):
                with self.assertRaises((ValueError, FileExistsError)):
                    bounded_launcher.launch(self.candidate, output, [], manufactured=True)
            link = self.base / 'parent-link'
            link.symlink_to(self.base, target_is_directory=True)
            with self.assertRaises(ValueError):
                bounded_launcher.launch(self.candidate, link / 'new-output', [], manufactured=True)

    def test_authentication_rejection_precedes_output_and_fork(self):
        output = self.base / 'not-created'
        with mock.patch.object(bounded_launcher, '_authenticate_candidate',
                               side_effect=ValueError('manufactured denied pin')):
            with mock.patch.object(bounded_launcher.os, 'fork') as fork:
                with self.assertRaisesRegex(ValueError, 'denied pin'):
                    bounded_launcher.launch(self.candidate, output, [], manufactured=True)
                fork.assert_not_called()
        self.assertFalse(output.exists())


    def test_zero_exit_cannot_replace_bound_worker_entry(self):
        for mode in ('false-zero', 'entry-hash-tamper', 'entry-status-tamper',
                     'entry-no-kernel', 'entry-bool-count', 'readback-status-tamper',
                     'readback-count-tamper'):
            receipt, output = self.run_probe(mode)
            self.assertEqual(receipt['exit_code'], 0, mode)
            self.assertEqual(receipt['status'], 'EXECUTION_FAILED', mode)
            self.assertEqual(receipt['stop_reason'], 'INVALID_ENTRY_OUTPUT', mode)
            self.assertEqual(receipt['wrapper_exit_code'], 1, mode)

    def test_manufactured_authentication_cannot_claim_physical_scope(self):
        arguments = ['--registration', str(self.base/'registration.json'),
                     '--registration-sha256', '0'*64]
        for flag in ('--go', '--go-sha256', '--repository-root'):
            with self.assertRaisesRegex(ValueError, 'must not claim'):
                bounded_launcher._authenticate_candidate(self.candidate,
                    arguments+[flag, str(self.base/'fake') if flag!='--go-sha256' else '0'*64], True)
        with self.assertRaisesRegex(ValueError, 'requires public GO'):
            bounded_launcher._authenticate_candidate(self.candidate, arguments, False)

    def test_real_source_authentication_rejects_mutation_before_fork(self):
        import hashlib
        owned = self.base/'owned-source'
        source = Path(__file__).resolve().parents[1]
        shutil.copytree(source, owned, ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
        contract = {'arithmetic_bits':1024,'momentum_degree':2048,'source_degree':24,
            'source_coefficient_bits':512,'source_cells':64,'export_bits':96,
            'export_gate':'1/1000000000000000000','anchor':'-9/2','endpoints':['-4','-7/2'],
            'selected_decode_members':['k.npy','momentum_weights.npy','observation_eta.npy',
                'u_1.npy','w_1.npy','u_2.npy','w_2.npy','u_3.npy','w_3.npy'],
            'wall_seconds':900,'rss_kib':524288,'output_bytes':134217728}
        files={}
        for path in sorted(owned.rglob('*')):
            if path.is_file():
                raw=path.read_bytes()
                files[str(path.relative_to(owned))]={'bytes':len(raw),'sha256':hashlib.sha256(raw).hexdigest()}
        registration=self.base/'closure.json'
        raw=(json.dumps({'schema_version':1,'scope':'LATER_RETAINED_ENDPOINT_FLOW_CERTIFICATE',
                         'files':files,'contract':contract},sort_keys=True)+'\n').encode()
        registration.write_bytes(raw)
        pin=hashlib.sha256(raw).hexdigest()
        result=bounded_launcher._authenticate_candidate(owned,
            ['--registration',str(registration),'--registration-sha256',pin],True)
        self.assertTrue(result['manufactured'])
        target=owned/'execution'/'run_endpoint.py'
        target.write_bytes(target.read_bytes()+b'\n# manufactured byte mutation\n')
        output=self.base/'never-created-after-source-mismatch'
        with mock.patch.object(bounded_launcher.os,'fork') as fork:
            with self.assertRaisesRegex(ValueError,'frozen source byte mismatch'):
                bounded_launcher.launch(owned,output,
                    ['--registration',str(registration),'--registration-sha256',pin], manufactured=True)
            fork.assert_not_called()
        self.assertFalse(output.exists())


    def test_bootstrap_rejects_executable_cache_before_local_import(self):
        import importlib._bootstrap_external
        import subprocess
        root=self.base/'cache-bootstrap'
        execution=root/'execution'
        execution.mkdir(parents=True)
        shutil.copyfile(Path(bounded_launcher.__file__),execution/'bounded_launcher.py')
        guard=execution/'kernel_guard.py'
        shutil.copyfile(Path(kernel_guard.__file__),guard)
        cachedir=execution/'__pycache__'
        cachedir.mkdir()
        marker='MANUFACTURED_CACHED_KERNEL_EXECUTION_MUST_NEVER_HAPPEN'
        code=compile('raise RuntimeError('+repr(marker)+')',str(guard),'exec')
        meta=guard.stat()
        cached=cachedir/('kernel_guard.'+sys.implementation.cache_tag+'.pyc')
        cached.write_bytes(importlib._bootstrap_external._code_to_timestamp_pyc(
            code,int(meta.st_mtime),meta.st_size))
        result=subprocess.run([sys.executable,'-I','-B',str(execution/'bounded_launcher.py'),'--help'],
            capture_output=True,text=True,timeout=10,env={'PATH':'/usr/bin:/bin','LANG':'C.UTF-8',
                                                        'PYTHONDONTWRITEBYTECODE':'1'})
        self.assertNotEqual(result.returncode,0)
        self.assertIn('bytecode caches are forbidden before local imports',result.stderr)
        self.assertNotIn(marker,result.stderr)
        # Even cache entries outside __pycache__ must fail before importing a guard.
        other=self.base/'standalone-cache'
        other.mkdir()
        (other/'standalone.pyc').write_bytes(cached.read_bytes())
        with self.assertRaisesRegex(ValueError,'bytecode caches'):
            bounded_launcher._reject_local_bytecode(other)

        native=self.base/'native-bootstrap'/'execution'
        native.mkdir(parents=True)
        shutil.copyfile(Path(bounded_launcher.__file__),native/'bounded_launcher.py')
        shutil.copyfile(Path(kernel_guard.__file__),native/'kernel_guard.py')
        (native/'kernel_guard.so').write_bytes(b'nonexecuted manufactured native alias')
        result=subprocess.run([sys.executable,'-I','-B',str(native/'bounded_launcher.py'),'--help'],
            capture_output=True,text=True,timeout=10,env={'PATH':'/usr/bin:/bin','LANG':'C.UTF-8',
                                                        'PYTHONDONTWRITEBYTECODE':'1'})
        self.assertNotEqual(result.returncode,0)
        self.assertIn('compiled local modules are forbidden before local imports',result.stderr)


    def test_bootstrap_rejects_python_packages_and_stdlib_aliases(self):
        import subprocess
        marker='MANUFACTURED_PYTHON_ALIAS_EXECUTION_MUST_NEVER_HAPPEN'
        for number,relative in enumerate(('kernel_guard/__init__.py',
                'registration_guard/__init__.py','ctypes.py','resource.py')):
            execution=self.base/('python-alias-'+str(number))/'execution'
            execution.mkdir(parents=True)
            shutil.copyfile(Path(bounded_launcher.__file__),execution/'bounded_launcher.py')
            shutil.copyfile(Path(kernel_guard.__file__),execution/'kernel_guard.py')
            alias=execution/relative
            alias.parent.mkdir(parents=True,exist_ok=True)
            alias.write_text('raise RuntimeError('+repr(marker)+')\n')
            result=subprocess.run([sys.executable,'-I','-B',str(execution/'bounded_launcher.py'),'--help'],
                capture_output=True,text=True,timeout=10,env={'PATH':'/usr/bin:/bin','LANG':'C.UTF-8',
                                                            'PYTHONDONTWRITEBYTECODE':'1'})
            self.assertNotEqual(result.returncode,0,relative)
            self.assertIn('aliases are forbidden before local imports',result.stderr,relative)
            self.assertNotIn(marker,result.stderr,relative)

    def test_bootstrap_unreadable_package_stops_before_known_child_import(self):
        import subprocess
        root=self.base/'unreadable-bootstrap'
        execution=root/'execution'
        execution.mkdir(parents=True)
        shutil.copyfile(Path(bounded_launcher.__file__),execution/'bounded_launcher.py')
        shutil.copyfile(Path(kernel_guard.__file__),execution/'kernel_guard.py')
        package=root/'source'/'later_source'
        package.mkdir(parents=True)
        marker='MANUFACTURED_EXECUTE_ONLY_PACKAGE_MUST_NEVER_EXECUTE'
        (package/'__init__.py').write_text('raise RuntimeError('+repr(marker)+')\n')
        # A known child is readable through an execute-only directory while
        # enumeration fails. This reproduced the silent-rglob omission class.
        package.chmod(0o111)
        try:
            self.assertGreater(os.geteuid(),0)
            self.assertIn(marker,(package/'__init__.py').read_text())
            with self.assertRaisesRegex(ValueError,'Candidate tree inspection failed'):
                bounded_launcher._reject_local_bytecode(root)
            result=subprocess.run([sys.executable,'-I','-B',str(execution/'bounded_launcher.py'),'--help'],
                capture_output=True,text=True,timeout=10,env={'PATH':'/usr/bin:/bin','LANG':'C.UTF-8'})
            self.assertNotEqual(result.returncode,0)
            self.assertIn('Candidate tree inspection failed',result.stderr)
            self.assertNotIn(marker,result.stderr)
        finally:
            package.chmod(0o700)
        with self.assertRaisesRegex(ValueError,'package aliases are forbidden'):
            bounded_launcher._reject_local_bytecode(root)

    def test_bootstrap_descriptor_scan_rejects_links_replacement_and_scan_errors(self):
        root=self.base/'descriptor-bootstrap'
        root.mkdir()
        plain=root/'plain'
        plain.mkdir()
        (plain/'safe.txt').write_bytes(b'manufactured tree bytes')
        bounded_launcher._reject_local_bytecode(root)
        link=root/'linked'
        link.symlink_to(plain,target_is_directory=True)
        try:
            with self.assertRaisesRegex(ValueError,'regular single-link'):
                bounded_launcher._reject_local_bytecode(root)
        finally:
            link.unlink()
        os.link(plain/'safe.txt',root/'hardlink.txt')
        try:
            with self.assertRaisesRegex(ValueError,'regular single-link'):
                bounded_launcher._reject_local_bytecode(root)
        finally:
            (root/'hardlink.txt').unlink()
        with mock.patch.object(bounded_launcher.os,'scandir',side_effect=PermissionError('manufactured scan failure')):
            with self.assertRaisesRegex(ValueError,'Candidate tree inspection failed'):
                bounded_launcher._reject_local_bytecode(root)
        actual_open=os.open
        swapped=False
        def replace_child_after_open(path,flags,*args,**kwargs):
            nonlocal swapped
            fd=actual_open(path,flags,*args,**kwargs)
            if path=='plain' and kwargs.get('dir_fd') is not None and not swapped:
                swapped=True
                plain.rename(root/'displaced')
                plain.mkdir()
            return fd
        with mock.patch.object(bounded_launcher.os,'open',side_effect=replace_child_after_open):
            with self.assertRaisesRegex(ValueError,'Candidate directory (changed before|replaced during) inspection'):
                bounded_launcher._reject_local_bytecode(root)
        self.assertTrue(swapped)

    def test_nonisolated_cli_refuses_before_shadowed_stdlib_import(self):
        import subprocess
        execution=self.base/'nonisolated-bootstrap'/'execution'
        execution.mkdir(parents=True)
        shutil.copyfile(Path(bounded_launcher.__file__),execution/'bounded_launcher.py')
        marker='MANUFACTURED_RESOURCE_SHADOW_MUST_NEVER_EXECUTE'
        (execution/'resource.py').write_text('raise RuntimeError('+repr(marker)+')\n')
        result=subprocess.run([sys.executable,'-B',str(execution/'bounded_launcher.py'),'--help'],
            capture_output=True,text=True,timeout=10,env={'PATH':'/usr/bin:/bin','LANG':'C.UTF-8',
                                                        'PYTHONDONTWRITEBYTECODE':'1'})
        self.assertNotEqual(result.returncode,0)
        self.assertIn('Isolated Python -I invocation is required before bootstrap imports',result.stderr)
        self.assertNotIn(marker,result.stderr)

    def test_registration_guard_snapshot_authentication_precedes_its_execution(self):
        import hashlib
        root=self.base/'guard-snapshot'
        execution=root/'execution'
        execution.mkdir(parents=True)
        marker='MANUFACTURED_UNAUTHENTICATED_GUARD_MUST_NEVER_EXECUTE'
        raw=('raise RuntimeError('+repr(marker)+')\n').encode()
        (execution/'registration_guard.py').write_bytes(raw)
        document={'schema_version':1,'scope':'LATER_RETAINED_ENDPOINT_FLOW_CERTIFICATE',
                  'files':{'execution/registration_guard.py':{'bytes':len(raw),'sha256':'0'*64}}}
        registration=self.base/'guard-snapshot-registration.json'
        serialized=(json.dumps(document,sort_keys=True)+'\n').encode()
        registration.write_bytes(serialized)
        before=list(sys.path)
        with self.assertRaisesRegex(ValueError,'Bootstrap registration guard byte pin differs'):
            bounded_launcher._authenticate_candidate(root,
                ['--registration',str(registration),'--registration-sha256',
                 hashlib.sha256(serialized).hexdigest()],True)
        self.assertEqual(sys.path,before)
        with self.assertRaisesRegex(ValueError,'Bootstrap registration SHA mismatch'):
            bounded_launcher._authenticate_candidate(root,
                ['--registration',str(registration),'--registration-sha256','f'*64],True)


if __name__ == '__main__':
    unittest.main(verbosity=2)
