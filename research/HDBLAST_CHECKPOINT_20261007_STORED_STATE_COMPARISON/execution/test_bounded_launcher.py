"""Fabricated-only kernel/launcher checks; no archive or physical source access.

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
sys.path.insert(0, str(Path(__file__).resolve().parent))
import bounded_launcher
import kernel_guard

PROBE = r'''
from pathlib import Path
import argparse, ctypes, errno, json, os, resource, signal, socket, sys, time
from decimal import Decimal
from fractions import Fraction
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kernel_guard
parser = argparse.ArgumentParser()
parser.add_argument('--candidate', required=True)
parser.add_argument('--output', required=True)
parser.add_argument('--registration-sha256', required=True)
parser.add_argument('--fabricated-only', action='store_true', required=True)
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
print('FABRICATED_PROBE_DIAGNOSTIC ' + mode, flush=True)
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
    (output / 'original').write_text('fabricated')
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
    (output / 'child.log').write_text('fabricated log contamination\n')
elif mode == 'receipt-tamper':
    (output / 'EXECUTION.json').symlink_to('/etc/passwd')
elif mode == 'directory-replacement':
    renamed = output.with_name(output.name + '-renamed')
    output.rename(renamed)
    output.mkdir()
    output = renamed
(output / 'PROBE.json').write_text(json.dumps(data, sort_keys=True) + '\n')
'''


class FabricatedLauncherTests(unittest.TestCase):
    def setUp(self):
        self.artifacts=Path(os.environ.get('HDBLAST_COMPARISON_TEST_OUTPUT',str(Path.cwd())))
        self.artifacts.mkdir(parents=True,exist_ok=True)
        self.temporary = tempfile.TemporaryDirectory(prefix='fabricated-launcher-', dir=self.artifacts)
        self.base = Path(self.temporary.name)
        self.candidate = self.base / 'candidate'
        execution = self.candidate / 'execution'
        execution.mkdir(parents=True)
        (execution / 'run_registered.py').write_text(PROBE)
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
                               return_value={'scope': 'FABRICATED_TEST_HARNESS_ONLY'}):
            with mock.patch.object(bounded_launcher, '_validate_successful_worker', return_value=None):
                with mock.patch.dict(os.environ, {'PROBE_SECRET': 'fabricated-not-a-secret',
                                             'GITHUB_TOKEN': 'fabricated-not-a-token',
                                             'AWS_SECRET_ACCESS_KEY': 'fabricated-not-a-key'}):
                    with mock.patch.multiple(bounded_launcher, **overrides) if overrides else nullcontext():
                        receipt = bounded_launcher.launch(
                            self.candidate, output,
                            ['--registration-sha256', '0' * 64, '--probe-mode', mode],
                            optimized=optimized, fabricated=True)
        return receipt, output

    def test_exact_limits_environment_and_optimization(self):
        for optimized in (False, True):
            receipt, output = self.run_probe('success', optimized)
            self.assertEqual(receipt['status'], 'PASS_BOUNDED_FABRICATED_EXECUTION')
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
            self.assertIn('FABRICATED_PROBE_DIAGNOSTIC', (output / 'child.log').read_text())

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
        self.assertIn('FABRICATED_PROBE_DIAGNOSTIC', (output / 'child.log').read_text())

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
                self.assertIn('FABRICATED_PROBE_DIAGNOSTIC',
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
                raise RuntimeError('fabricated monitor fault')
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
                    bounded_launcher.launch(self.candidate, output, [], fabricated=True)
            link = self.base / 'parent-link'
            link.symlink_to(self.base, target_is_directory=True)
            with self.assertRaises(ValueError):
                bounded_launcher.launch(self.candidate, link / 'new-output', [], fabricated=True)

    def test_authentication_rejection_precedes_output_and_fork(self):
        output = self.base / 'not-created'
        with mock.patch.object(bounded_launcher, '_authenticate_candidate',
                               side_effect=ValueError('fabricated denied pin')):
            with mock.patch.object(bounded_launcher.os, 'fork') as fork:
                with self.assertRaisesRegex(ValueError, 'denied pin'):
                    bounded_launcher.launch(self.candidate, output, [], fabricated=True)
                fork.assert_not_called()
        self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main(verbosity=2)
