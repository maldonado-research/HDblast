"""Nonroot custodian for the registered, offline comparison worker.

The worker is exec'd before it installs seccomp after numeric imports. The
custodian imports no numerical dependency, authenticates first, applies exact
hard rlimits before exec, and records the authoritative Linux wait4 outcome.
Aggregate output is polled and checked again after exit; it is an acceptance
limit, not a filesystem quota. No physical authorization is created here.
"""
from pathlib import Path
import argparse
import errno
import hashlib
import json
import os
import resource
import selectors
import signal
import stat
import sys
import time

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from kernel_guard import (ADDRESS_SPACE_BYTES, CPU_SECONDS, OUTPUT_BYTES,
                          THREAD_VARIABLES, WALL_SECONDS, check_output,
                          no_new_privileges, nonroot, open_real_directory,
                          real_directory, require)

RECEIPT_RESERVE_BYTES = 64 * 1024
WORKER_OUTPUT_BYTES = OUTPUT_BYTES - RECEIPT_RESERVE_BYTES
POLL_SECONDS = 0.05
REQUIRED_WORKER_OUTPUTS = (
    'SOURCE_CERTIFICATE.json', 'DIAGNOSTIC_CROSSCHECK.json', 'NODE_TARGETS.jsonl.gz',
    'DATA.json', 'SCIENCE_SUMMARY.json', 'RESOURCE.json', 'SOURCE_ATTEMPTS.jsonl',
    'DECODE_ATTEMPTS.jsonl', 'OUTPUT_READBACK.json', 'ENTRY_RECEIPT.json',
)


def _argument(arguments, flag, required=False):
    positions = [index for index, value in enumerate(arguments) if value == flag]
    require(len(positions) <= 1, 'Repeated worker argument: ' + flag)
    if not positions:
        require(not required, 'Required worker argument: ' + flag)
        return None
    index = positions[0]
    require(index + 1 < len(arguments) and not arguments[index + 1].startswith('--'),
            'Missing worker argument value: ' + flag)
    return arguments[index + 1]


def _authenticate_candidate(candidate, arguments, fabricated):
    # No numeric imports or source callbacks occur in this authentication.
    from registration_guard import authenticate, authenticate_local
    registration = _argument(arguments, '--registration-sha256', True)
    public_go = _argument(arguments, '--public-go')
    public_go_sha = _argument(arguments, '--public-go-sha256')
    if fabricated:
        require(public_go is None and public_go_sha is None and
                _argument(arguments, '--repository-root') is None,
                'Fabricated launch must not claim physical/public authorization')
        return authenticate_local(candidate, registration)
    require(public_go is not None and public_go_sha is not None,
            'Physical launch requires an externally pinned public GO receipt')
    return authenticate(candidate, registration, public_go, public_go_sha)


def _terminate(pid):
    # The direct kill also covers the race before the child calls setsid().
    for kill in (lambda: os.killpg(pid, signal.SIGKILL),
                 lambda: os.kill(pid, signal.SIGKILL)):
        try:
            kill()
        except ProcessLookupError:
            pass


def _wait4(pid, options):
    while True:
        try:
            return os.wait4(pid, options)
        except InterruptedError:
            continue


def _fresh_output(candidate, output):
    output = Path(os.path.abspath(os.fspath(output)))
    require(output != candidate and candidate not in output.parents,
            'Fresh external worker output directory required')
    parent, parent_fd = open_real_directory(output.parent)
    try:
        # mkdir is exclusive: an existing empty directory is contamination too.
        os.mkdir(output.name, 0o700, dir_fd=parent_fd)
        fd = os.open(output.name, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW |
                     os.O_CLOEXEC, dir_fd=parent_fd)
        try:
            require(os.fstat(fd).st_uid == os.getuid(), 'Output custodian owner required')
        except BaseException:
            os.close(fd)
            raise
        return output, fd
    finally:
        os.close(parent_fd)


def _same_named_file(directory_fd, name, original):
    current = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
    require(stat.S_ISREG(current.st_mode) and current.st_nlink == 1 and
            (current.st_dev, current.st_ino) == (original.st_dev, original.st_ino),
            'Custodian log was removed or replaced')


def _same_output_directory(directory_fd, output):
    unused_path, named_fd = open_real_directory(output)
    try:
        named = os.fstat(named_fd)
        held = os.fstat(directory_fd)
        require((named.st_dev, named.st_ino) == (held.st_dev, held.st_ino),
                'Custodian output directory was renamed or replaced')
    finally:
        os.close(named_fd)


def _validate_successful_worker(directory_fd, arguments, fabricated):
    """Verify the concrete successful-entry artifacts without numeric imports."""
    from registration_guard import load
    raw_files = {}
    for name in REQUIRED_WORKER_OUTPUTS:
        fd = os.open(name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK |
                     os.O_CLOEXEC, dir_fd=directory_fd)
        try:
            metadata = os.fstat(fd)
            require(stat.S_ISREG(metadata.st_mode) and metadata.st_nlink == 1,
                    'Required worker output must be a regular unaliased file: ' + name)
            require(metadata.st_size > 0 or name in ('SOURCE_ATTEMPTS.jsonl', 'DECODE_ATTEMPTS.jsonl'),
                    'Required worker output is empty: ' + name)
            digest = hashlib.sha256()
            receipt_raw = bytearray()
            bytes_read = 0
            while True:
                raw = os.read(fd, 1024 * 1024)
                if not raw:
                    break
                bytes_read += len(raw)
                require(bytes_read <= WORKER_OUTPUT_BYTES, 'Required worker output exceeds bound')
                digest.update(raw)
                if name == 'ENTRY_RECEIPT.json':
                    require(bytes_read <= RECEIPT_RESERVE_BYTES, 'Worker entry receipt exceeds bound')
                    receipt_raw.extend(raw)
            require(bytes_read == metadata.st_size, 'Required worker output changed')
            raw_files[name] = {'bytes': bytes_read, 'sha256': digest.hexdigest()}
            if name == 'ENTRY_RECEIPT.json':
                entry = load(bytes(receipt_raw))
        finally:
            os.close(fd)
    expected = ('PASS_FABRICATED_STORED_COMPARISON' if fabricated else
                'PASS_AUTHENTICATED_STORED_COMPARISON')
    require(type(entry) is dict and entry.get('status') == expected and
            entry.get('fabricated_only') is bool(fabricated) and
            entry.get('registration_sha256') == _argument(arguments, '--registration-sha256', True) and
            entry.get('public_go_sha256') == _argument(arguments, '--public-go-sha256'),
            'Worker entry status or authorization binding differs')
    require(entry.get('files') == {name: item for name, item in raw_files.items()
                                 if name != 'ENTRY_RECEIPT.json'},
            'Worker entry hashes do not match output bytes')
    return raw_files['ENTRY_RECEIPT.json']


def _preserve_log(directory_fd, log_fd):
    """Recover the held log inode when a fabricated worker unlinks its name."""
    name = 'CUSTODIAN_CHILD.log'
    try:
        recovered_fd = os.open(name, os.O_WRONLY | os.O_CREAT | os.O_EXCL |
                               os.O_NOFOLLOW | os.O_CLOEXEC, 0o600, dir_fd=directory_fd)
    except FileExistsError:
        name = 'CUSTODIAN_CHILD_' + str(time.monotonic_ns()) + '.log'
        recovered_fd = os.open(name, os.O_WRONLY | os.O_CREAT | os.O_EXCL |
                               os.O_NOFOLLOW | os.O_CLOEXEC, 0o600, dir_fd=directory_fd)
    try:
        offset = 0
        while offset < WORKER_OUTPUT_BYTES:
            raw = os.pread(log_fd, min(1024 * 1024, WORKER_OUTPUT_BYTES - offset), offset)
            if not raw:
                break
            written = 0
            while written < len(raw):
                written += os.write(recovered_fd, raw[written:])
            offset += len(raw)
        os.fsync(recovered_fd)
    finally:
        os.close(recovered_fd)
    return name


def launch(candidate, output, worker_args, optimized=False, fabricated=False):
    """Return an execution receipt; never turn local preparation into public GO.

    ``worker_args`` supplies the registration pin and, only for actual work,
    public GO pins and repository root. Candidate/output/scope flags belong to
    this custodian. Launch failures leave the fresh directory and child log.
    """
    owned_fds = set()
    owned_selectors = []
    try:
        return _launch(candidate, output, worker_args, optimized, fabricated,
                       owned_fds, owned_selectors)
    finally:
        for selector in owned_selectors:
            selector.close()
        for fd in owned_fds:
            try:
                os.close(fd)
            except OSError as exc:
                if exc.errno != errno.EBADF:
                    raise


def _launch(candidate, output, worker_args, optimized, fabricated,
            owned_fds, owned_selectors):
    nonroot()
    require(sys.platform == 'linux', 'Linux wait4/seccomp execution required')
    candidate = real_directory(candidate)
    require(type(worker_args) in (list, tuple) and
            all(type(value) is str and value and '\0' not in value for value in worker_args),
            'Worker arguments must be finite strings')
    arguments = list(worker_args)
    require(sum(len(value.encode()) for value in arguments) <= 16 * 1024,
            'Worker arguments exceed custodian receipt bound')
    require(not any(value in ('--candidate', '--output', '--fabricated-only')
                    or value.startswith(('--candidate=', '--output=', '--fabricated-only='))
                    for value in arguments), 'Custodian-owned worker argument supplied')
    _authenticate_candidate(candidate, arguments, fabricated)
    worker = candidate / 'execution' / 'run_registered.py'
    require(worker.is_file() and not worker.is_symlink(), 'Authenticated worker required')
    output, directory_fd = _fresh_output(candidate, output)
    owned_fds.add(directory_fd)
    command = [sys.executable, '-I', '-B'] + (['-O'] if optimized else []) + [
        str(worker), '--candidate', str(candidate), '--output', str(output), *arguments]
    if fabricated:
        command.append('--fabricated-only')
    environment = {'PATH': '/usr/local/bin:/usr/bin:/bin', 'LANG': 'C.UTF-8',
                   'PYTHONDONTWRITEBYTECODE': '1'}
    environment.update({name: '1' for name in THREAD_VARIABLES})
    log_fd = os.open('child.log', os.O_RDWR | os.O_CREAT | os.O_EXCL |
                     os.O_NOFOLLOW | os.O_CLOEXEC, 0o600, dir_fd=directory_fd)
    owned_fds.add(log_fd)
    log_identity = os.fstat(log_fd)
    read_fd, write_fd = os.pipe2(os.O_CLOEXEC)
    owned_fds.update((read_fd, write_fd))
    selector = selectors.DefaultSelector()
    owned_selectors.append(selector)
    os.set_blocking(read_fd, False)
    selector.register(read_fd, selectors.EVENT_READ)
    started = time.monotonic()
    pid = None
    status = usage = None
    stop_reason = None
    monitor_error = None
    summary = {'bytes': 0, 'file_count': 1, 'directory_count': 0}
    log_bytes = 0
    log_digest = hashlib.sha256()
    try:
        pid = os.fork()
        if pid == 0:
            try:
                os.setsid()
                null_fd = os.open('/dev/null', os.O_RDONLY)
                os.dup2(null_fd, 0)
                os.dup2(write_fd, 1)
                os.dup2(write_fd, 2)
                # Enumerate actual descriptors: an inherited fd can exceed a
                # subsequently lowered RLIMIT_NOFILE value.
                inherited_fds = [int(name) for name in os.listdir('/proc/self/fd')]
                for inherited_fd in inherited_fds:
                    if inherited_fd > 2:
                        try:
                            os.close(inherited_fd)
                        except OSError as exc:
                            if exc.errno != errno.EBADF:
                                raise
                for limit, expected in ((resource.RLIMIT_CORE, 0),
                                        (resource.RLIMIT_AS, ADDRESS_SPACE_BYTES),
                                        (resource.RLIMIT_CPU, CPU_SECONDS),
                                        (resource.RLIMIT_FSIZE, OUTPUT_BYTES),
                                        (resource.RLIMIT_NPROC, 0)):
                    resource.setrlimit(limit, (expected, expected))
                no_new_privileges()
                os.chdir(output)
                os.execve(sys.executable, command, environment)
            except BaseException as exc:
                print('CHILD_LAUNCH_FAILED ' + type(exc).__name__ + ': ' + str(exc),
                      file=sys.stderr, flush=True)
                os._exit(126)
        os.close(write_fd)
        owned_fds.discard(write_fd)
        write_fd = None
        while status is None:
            for key, unused in selector.select(POLL_SECONDS):
                try:
                    raw = os.read(key.fd, 64 * 1024)
                except BlockingIOError:
                    continue
                if not raw:
                    selector.unregister(key.fd)
                    continue
                remaining = max(0, WORKER_OUTPUT_BYTES - log_bytes)
                kept = raw[:remaining]
                offset = 0
                while offset < len(kept):
                    offset += os.write(log_fd, kept[offset:])
                log_bytes += len(kept)
                log_digest.update(kept)
                if len(kept) != len(raw):
                    stop_reason = 'OUTPUT_LIMIT'
            if time.monotonic() - started > WALL_SECONDS:
                stop_reason = 'WALL_TIMEOUT'
            if stop_reason is None:
                try:
                    _same_output_directory(directory_fd, output)
                    _same_named_file(directory_fd, 'child.log', log_identity)
                    summary = check_output(output, WORKER_OUTPUT_BYTES,
                                           forbidden_roots=(candidate,),
                                           deadline=started + WALL_SECONDS,
                                           directory_fd=directory_fd)
                    _same_output_directory(directory_fd, output)
                    require(not os.path.lexists(output / 'EXECUTION.json'),
                            'Worker contaminated custodian receipt filename')
                except (ValueError, OSError) as exc:
                    stop_reason = ('OUTPUT_LIMIT' if 'exceeds bound' in str(exc)
                                   else 'UNSAFE_OUTPUT')
                    monitor_error = str(exc)
            if stop_reason:
                _terminate(pid)
                unused_pid, status, usage = _wait4(pid, 0)
            else:
                done, candidate_status, candidate_usage = _wait4(pid, os.WNOHANG)
                if done:
                    status, usage = candidate_status, candidate_usage
        # Drain already-produced diagnostics after the authoritative child exit.
        while selector.get_map():
            try:
                raw = os.read(read_fd, 64 * 1024)
            except BlockingIOError:
                break
            if not raw:
                break
            remaining = max(0, WORKER_OUTPUT_BYTES - log_bytes)
            kept = raw[:remaining]
            offset = 0
            while offset < len(kept):
                offset += os.write(log_fd, kept[offset:])
            log_bytes += len(kept)
            log_digest.update(kept)
            if len(kept) != len(raw):
                stop_reason = 'OUTPUT_LIMIT'
    except Exception as exc:
        if pid:
            if status is None:
                _terminate(pid)
                unused_pid, status, usage = _wait4(pid, 0)
            stop_reason = 'CUSTODIAN_ERROR'
            monitor_error = type(exc).__name__ + ': ' + str(exc)
        else:
            os.close(log_fd)
            owned_fds.discard(log_fd)
            os.close(directory_fd)
            owned_fds.discard(directory_fd)
            raise
    except BaseException:
        if pid and status is None:
            _terminate(pid)
            _wait4(pid, 0)
        os.close(log_fd)
        owned_fds.discard(log_fd)
        os.close(directory_fd)
        owned_fds.discard(directory_fd)
        raise
    finally:
        selector.close()
        os.close(read_fd)
        owned_fds.discard(read_fd)
        if write_fd is not None:
            os.close(write_fd)
            owned_fds.discard(write_fd)
        if status is not None:
            try:
                os.fsync(log_fd)
            except OSError as exc:
                stop_reason = 'CUSTODIAN_IO_ERROR'
                monitor_error = str(exc)
    elapsed = time.monotonic() - started
    code = os.waitstatus_to_exitcode(status)
    preserved_log = None
    try:
        _same_output_directory(directory_fd, output)
    except (ValueError, OSError) as exc:
        stop_reason = 'UNSAFE_OUTPUT'
        monitor_error = str(exc)
    try:
        _same_named_file(directory_fd, 'child.log', log_identity)
    except (ValueError, OSError) as exc:
        stop_reason = 'UNSAFE_OUTPUT'
        monitor_error = str(exc)
        preserved_log = _preserve_log(directory_fd, log_fd)
    observed_log_digest = hashlib.sha256()
    observed_log_bytes = 0
    while observed_log_bytes <= WORKER_OUTPUT_BYTES:
        raw = os.pread(log_fd, 1024 * 1024, observed_log_bytes)
        if not raw:
            break
        observed_log_bytes += len(raw)
        observed_log_digest.update(raw)
    if observed_log_bytes != log_bytes or observed_log_digest.digest() != log_digest.digest():
        stop_reason = 'UNSAFE_OUTPUT'
        monitor_error = 'Custodian log bytes were modified by worker'
    try:
        summary = check_output(output, WORKER_OUTPUT_BYTES, hash_files=True,
                               forbidden_roots=(candidate,),
                               deadline=time.monotonic() + 10,
                               directory_fd=directory_fd)
        require(not os.path.lexists(output / 'EXECUTION.json'),
                'Worker contaminated custodian receipt filename')
    except (ValueError, OSError) as exc:
        stop_reason = ('OUTPUT_LIMIT' if 'exceeds bound' in str(exc) else 'UNSAFE_OUTPUT')
        monitor_error = str(exc)
    entry_identity = None
    if code == 0 and stop_reason is None:
        try:
            entry_identity = _validate_successful_worker(directory_fd, arguments, fabricated)
            _same_output_directory(directory_fd, output)
        except (ValueError, OSError) as exc:
            stop_reason = 'INVALID_ENTRY_OUTPUT'
            monitor_error = str(exc)
    passed = (code == 0 and stop_reason is None and elapsed <= WALL_SECONDS and
              usage.ru_maxrss * 1024 <= ADDRESS_SPACE_BYTES and
              usage.ru_utime + usage.ru_stime <= CPU_SECONDS)
    receipt = {
        'status': 'PASS_BOUNDED_FABRICATED_EXECUTION' if passed and fabricated else
                  'PASS_BOUNDED_REGISTERED_EXECUTION' if passed else 'EXECUTION_FAILED',
        'fabricated_only': bool(fabricated), 'optimized': bool(optimized),
        'uid': os.getuid(), 'exit_code': code, 'wait_status': status,
        'wrapper_exit_code': code if code > 0 else 128 - code if code < 0 else 0 if passed else 1,
        'stop_reason': stop_reason, 'monitor_error': monitor_error,
        'preserved_child_log': preserved_log,
        'captured_child_log_sha256': log_digest.hexdigest(),
        'observed_child_log_sha256': observed_log_digest.hexdigest(),
        'entry_receipt': entry_identity,
        'custodian_output_directory': os.readlink('/proc/self/fd/' + str(directory_fd)),
        'wall_seconds': elapsed, 'peak_rss_kib_wait4': usage.ru_maxrss,
        'cpu_seconds_wait4': usage.ru_utime + usage.ru_stime,
        'worker_output': summary, 'command': command,
        'limits': {'wall_seconds': WALL_SECONDS, 'cpu_seconds': CPU_SECONDS,
                   'address_space_bytes': ADDRESS_SPACE_BYTES,
                   'aggregate_output_bytes': OUTPUT_BYTES, 'file_size_bytes': OUTPUT_BYTES,
                   'process_limit': 0, 'core_bytes': 0,
                   'custodian_receipt_reserve_bytes': RECEIPT_RESERVE_BYTES},
        'aggregate_limit_enforcement': 'POLL_AND_FINAL_ACCEPTANCE_CHECK_WITH_RECEIPT_RESERVE',
    }
    raw = (json.dumps(receipt, sort_keys=True, indent=2, allow_nan=False) + '\n').encode()
    require(len(raw) <= RECEIPT_RESERVE_BYTES, 'Custodian receipt exceeds reserved bound')
    try:
        # Never follow or replace worker-created receipt contamination.
        receipt_fd = os.open('EXECUTION.json', os.O_WRONLY | os.O_CREAT | os.O_EXCL |
                             os.O_NOFOLLOW | os.O_CLOEXEC, 0o600, dir_fd=directory_fd)
    except FileExistsError:
        try:
            receipt_fd = os.open('CUSTODIAN_FAILURE.json', os.O_WRONLY | os.O_CREAT | os.O_EXCL |
                                 os.O_NOFOLLOW | os.O_CLOEXEC, 0o600, dir_fd=directory_fd)
        except FileExistsError:
            receipt_fd = os.open('CUSTODIAN_FAILURE_' + str(time.monotonic_ns()) + '.json',
                                 os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW |
                                 os.O_CLOEXEC, 0o600, dir_fd=directory_fd)
    owned_fds.add(receipt_fd)
    try:
        offset = 0
        while offset < len(raw):
            offset += os.write(receipt_fd, raw[offset:])
        os.fsync(receipt_fd)
    finally:
        os.close(receipt_fd)
        owned_fds.discard(receipt_fd)
        os.close(log_fd)
        owned_fds.discard(log_fd)
        os.close(directory_fd)
        owned_fds.discard(directory_fd)
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--registration-sha256', required=True)
    parser.add_argument('--public-go')
    parser.add_argument('--public-go-sha256')
    parser.add_argument('--repository-root')
    parser.add_argument('--optimized', action='store_true')
    parser.add_argument('--fabricated-only', action='store_true')
    args = parser.parse_args()
    worker_args = ['--registration-sha256', args.registration_sha256]
    for name in ('public_go', 'public_go_sha256', 'repository_root'):
        value = getattr(args, name)
        if value is not None:
            worker_args.extend(['--' + name.replace('_', '-'), value])
    receipt = launch(args.candidate, args.output, worker_args, args.optimized, args.fabricated_only)
    print(json.dumps(receipt, sort_keys=True), flush=True)
    raise SystemExit(receipt['wrapper_exit_code'])


if __name__ == '__main__':
    main()
