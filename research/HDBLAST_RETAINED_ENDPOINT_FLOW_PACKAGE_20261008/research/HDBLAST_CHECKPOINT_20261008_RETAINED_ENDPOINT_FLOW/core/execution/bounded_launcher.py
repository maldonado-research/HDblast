"""Nonroot custodian for the registered, offline endpoint worker.

The worker is exec'd before it installs seccomp after numeric imports. The
custodian imports no numerical dependency, authenticates first, applies exact
hard rlimits before exec, and records the authoritative Linux wait4 outcome.
Aggregate output is polled and checked again after exit; it is an acceptance
limit, not a filesystem quota. No physical authorization is created here.
"""
import sys
if __name__ == '__main__' and not sys.flags.isolated:
    raise ValueError('Isolated Python -I invocation is required before bootstrap imports')

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
import types

sys.dont_write_bytecode = True


def _reject_local_bytecode(root):
    """Fail closed on every directory, alias and traversal error before imports.

    Enumeration uses held directory descriptors. An execute-only directory
    must fail inspection even when its known children remain importable.
    """
    root = Path(os.path.abspath(os.fspath(root)))
    directory_flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
    identity = lambda metadata: (metadata.st_dev, metadata.st_ino, metadata.st_mode,
                                 metadata.st_mtime_ns, metadata.st_ctime_ns)
    entries = 0

    def walk(fd, depth):
        nonlocal entries
        if depth > 24:
            raise ValueError('Candidate tree inspection depth exceeds bound')
        before = os.fstat(fd)
        with os.scandir(fd) as iterator:
            names = [entry.name for entry in iterator]
        for name in names:
            entries += 1
            if entries > 8192:
                raise ValueError('Candidate tree inspection entry count exceeds bound')
            if name == '__pycache__' or name.endswith(('.pyc', '.pyo')):
                raise ValueError('Candidate bytecode caches are forbidden before local imports')
            if name.endswith(('.so', '.pyd', '.dll', '.dylib')):
                raise ValueError('Candidate compiled local modules are forbidden before local imports')
            if name == '__init__.py':
                raise ValueError('Candidate package aliases are forbidden before local imports')
            metadata = os.stat(name, dir_fd=fd, follow_symlinks=False)
            if (name.endswith('.py') and name[:-3] in sys.stdlib_module_names) or (
                    stat.S_ISDIR(metadata.st_mode) and name in sys.stdlib_module_names):
                raise ValueError('Candidate standard-library aliases are forbidden before local imports')
            if stat.S_ISDIR(metadata.st_mode):
                child = os.open(name, directory_flags, dir_fd=fd)
                try:
                    if identity(metadata) != identity(os.fstat(child)):
                        raise ValueError('Candidate directory changed before inspection')
                    walk(child, depth + 1)
                    if identity(os.fstat(child)) != identity(os.stat(name, dir_fd=fd, follow_symlinks=False)):
                        raise ValueError('Candidate directory replaced during inspection')
                finally:
                    os.close(child)
            elif not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1:
                raise ValueError('Candidate tree requires regular single-link source files')
        if identity(before) != identity(os.fstat(fd)):
            raise ValueError('Candidate directory changed during inspection')

    fd = None
    try:
        for ancestor in (root, *root.parents):
            if not stat.S_ISDIR(os.lstat(ancestor).st_mode):
                raise ValueError('Real candidate ancestors without symlinks required')
        fd = os.open(root, directory_flags)
        walk(fd, 0)
        named = os.open(root, directory_flags)
        try:
            if identity(os.fstat(fd)) != identity(os.fstat(named)):
                raise ValueError('Candidate root replaced during inspection')
        finally:
            os.close(named)
    except OSError as error:
        raise ValueError('Candidate tree inspection failed: ' + str(error)) from error
    finally:
        if fd is not None:
            os.close(fd)


def _strict_json(raw):
    def unique_pairs(items):
        result = {}
        for name, value in items:
            if name in result:
                raise ValueError('Duplicate JSON key in bootstrap document')
            result[name] = value
        return result
    def reject_nonfinite(value):
        raise ValueError('Nonfinite JSON value in bootstrap document')
    return json.loads(raw, object_pairs_hook=unique_pairs, parse_constant=reject_nonfinite)


def _source_bytes(path):
    """Read a bounded, stable, unaliased regular source snapshot without imports."""
    path = Path(os.path.abspath(os.fspath(path)))
    for ancestor in reversed(path.parents):
        if not stat.S_ISDIR(os.lstat(ancestor).st_mode):
            raise ValueError('Real source ancestors without symlinks required')
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC)
    try:
        before = os.fstat(fd)
        if not stat.S_ISREG(before.st_mode) or before.st_nlink != 1 or before.st_size > 2*1024*1024:
            raise ValueError('Bounded unaliased regular bootstrap source required')
        data = bytearray()
        while True:
            block = os.read(fd, 65536)
            if not block:
                break
            data.extend(block)
            if len(data) > 2*1024*1024:
                raise ValueError('Bootstrap source snapshot exceeds bound')
        after = os.fstat(fd)
        if len(data) != before.st_size or (before.st_size, before.st_mtime_ns, before.st_ctime_ns) != (
                after.st_size, after.st_mtime_ns, after.st_ctime_ns):
            raise ValueError('Bootstrap source changed during snapshot')
        return bytes(data)
    finally:
        os.close(fd)


def _captured_module(name, path, raw):
    """Execute already authenticated bytes; never ask an import finder for them."""
    module = types.ModuleType(name)
    module.__file__ = str(path)
    module.__package__ = ''
    previous = sys.modules.get(name)
    sys.modules[name] = module
    try:
        exec(compile(raw, str(path), 'exec'), module.__dict__)
    except BaseException:
        if previous is None:
            sys.modules.pop(name, None)
        else:
            sys.modules[name] = previous
        raise
    return module


BOOTSTRAP_ROOT = Path(__file__).resolve().parents[1]
_reject_local_bytecode(BOOTSTRAP_ROOT)
_kernel_path = BOOTSTRAP_ROOT / 'execution' / 'kernel_guard.py'
_kernel_raw = _source_bytes(_kernel_path)
if hashlib.sha256(_kernel_raw).hexdigest() != '886a0bc7148e7c49d57b68303f8113c657a455601713c974513a90abd0bf2b43':
    raise ValueError('Immutable kernel bootstrap source hash differs')
_kernel = _captured_module('kernel_guard', _kernel_path, _kernel_raw)
for _name in ('ADDRESS_SPACE_BYTES','CPU_SECONDS','OUTPUT_BYTES','THREAD_VARIABLES',
              'WALL_SECONDS','DENIED_SYSCALLS','check_output','no_new_privileges',
              'nonroot','open_real_directory','real_directory','require'):
    globals()[_name] = getattr(_kernel, _name)

RECEIPT_RESERVE_BYTES = 64 * 1024
WORKER_OUTPUT_BYTES = OUTPUT_BYTES - RECEIPT_RESERVE_BYTES
POLL_SECONDS = 0.05
REQUIRED_WORKER_OUTPUTS = (
    'SOURCE_CERTIFICATE.json', 'TARGET_BUDGETS.json', 'NODE_CERTIFICATES.jsonl.gz',
    'DATA.json', 'SCIENCE_SUMMARY.json', 'RESOURCE.json', 'SOURCE_ATTEMPTS.jsonl',
    'DECODE_ATTEMPTS.jsonl', 'INPUT_AUTHENTICATION.json', 'OUTPUT_READBACK.json', 'ENTRY_RECEIPT.json',
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


def _authenticate_candidate(candidate, arguments, manufactured):
    # Authenticated source snapshots precede all candidate code and numeric imports.
    _reject_local_bytecode(candidate)
    registration_path = _argument(arguments, '--registration', True)
    registration_sha = _argument(arguments, '--registration-sha256', True)
    public_go = _argument(arguments, '--go')
    public_go_sha = _argument(arguments, '--go-sha256')
    repository_root = _argument(arguments, '--repository-root')
    require(Path(registration_path).is_absolute(), 'Absolute registration path required')
    if manufactured:
        require(public_go is None and public_go_sha is None and repository_root is None,
                'Manufactured launch must not claim physical/public authorization')
    else:
        require(public_go is not None and public_go_sha is not None and repository_root is not None,
                'Physical launch requires public GO pins and retained repository root')
        require(Path(public_go).is_absolute() and Path(repository_root).is_absolute(),
                'Absolute GO and repository paths required')
    registration_raw = _source_bytes(registration_path)
    require(type(registration_sha) is str and len(registration_sha) == 64 and
            all(c in '0123456789abcdef' for c in registration_sha) and
            hashlib.sha256(registration_raw).hexdigest() == registration_sha,
            'Bootstrap registration SHA mismatch')
    registration = _strict_json(registration_raw)
    require(type(registration) is dict and registration.get('schema_version') == 1 and
            registration.get('scope') == 'LATER_RETAINED_ENDPOINT_FLOW_CERTIFICATE' and
            type(registration.get('files')) is dict, 'Bootstrap registration scope differs')
    relative = 'execution/registration_guard.py'
    guard_path = candidate / relative
    guard_raw = _source_bytes(guard_path)
    pin = registration['files'].get(relative)
    require(type(pin) is dict and type(pin.get('bytes')) is int and
            pin.get('bytes') == len(guard_raw) and
            pin.get('sha256') == hashlib.sha256(guard_raw).hexdigest(),
            'Bootstrap registration guard byte pin differs')
    guard = _captured_module('_hdblast_endpoint_registration_guard', guard_path, guard_raw)
    return guard.authenticate(candidate, registration_path, registration_sha, manufactured,
                              public_go, public_go_sha)


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


def _validate_successful_worker(directory_fd, arguments, manufactured):
    """Verify the concrete successful-entry artifacts without numeric imports."""
    load = _strict_json
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
                if name in ('ENTRY_RECEIPT.json', 'OUTPUT_READBACK.json'):
                    require(bytes_read <= RECEIPT_RESERVE_BYTES, 'Worker entry receipt exceeds bound')
                    receipt_raw.extend(raw)
            require(bytes_read == metadata.st_size, 'Required worker output changed')
            raw_files[name] = {'bytes': bytes_read, 'sha256': digest.hexdigest()}
            if name == 'ENTRY_RECEIPT.json':
                entry = load(bytes(receipt_raw))
            elif name == 'OUTPUT_READBACK.json':
                readback = load(bytes(receipt_raw))
        finally:
            os.close(fd)
    expected = ('PASS_MANUFACTURED_ENDPOINT_ENTRY' if manufactured else
                'PASS_REGISTERED_ENDPOINT_ENTRY')
    require(type(entry) is dict and entry.get('status') == expected and
            entry.get('manufactured') is bool(manufactured) and
            entry.get('registration_sha256') == _argument(arguments, '--registration-sha256', True) and
            entry.get('go_sha256') == _argument(arguments, '--go-sha256') and
            entry.get('scope') == 'LATER_RETAINED_ENDPOINT_FLOW_CERTIFICATE' and
            type(entry.get('nodes')) is int and entry.get('nodes') == 49152 and
            type(entry.get('endpoint_comparisons')) is int and entry.get('endpoint_comparisons') == 98304,
            'Worker entry status or authorization binding differs')
    require(entry.get('outputs') == {name: item for name, item in raw_files.items()
                                 if name != 'ENTRY_RECEIPT.json'},
            'Worker entry hashes do not match output bytes')
    require(type(readback) is dict and
            readback.get('status') == 'PASS_STANDALONE_ENDPOINT_READBACK' and
            type(readback.get('nodes')) is int and readback.get('nodes') == 49152 and
            type(readback.get('endpoint_comparisons')) is int and
            readback.get('endpoint_comparisons') == 98304 and
            type(readback.get('prefix_cases')) is int and readback.get('prefix_cases') == 24,
            'Worker standalone readback status or coverage differs')
    denied = entry.get('kernel_denied_syscalls')
    require(type(denied) is list and all(type(name) is str for name in denied) and
            len(denied) == len(set(denied)) and set(denied) <= set(DENIED_SYSCALLS) and
            {'socket', 'socketpair', 'connect', 'clone', 'execve', 'io_uring_setup',
             'io_uring_enter', 'io_uring_register'} <= set(denied),
            'Worker entry omits mandatory kernel restriction evidence')
    return raw_files['ENTRY_RECEIPT.json']


def _preserve_log(directory_fd, log_fd):
    """Recover the held log inode when a manufactured worker unlinks its name."""
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


def launch(candidate, output, worker_args, optimized=False, manufactured=False):
    """Return an execution receipt; never turn local preparation into public GO.

    ``worker_args`` supplies the registration pin and, only for actual work,
    public GO pins and repository root. Candidate/output/scope flags belong to
    this custodian. Launch failures leave the fresh directory and child log.
    """
    owned_fds = set()
    owned_selectors = []
    try:
        return _launch(candidate, output, worker_args, optimized, manufactured,
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


def _launch(candidate, output, worker_args, optimized, manufactured,
            owned_fds, owned_selectors):
    require(sys.flags.isolated, 'Isolated Python -I execution required')
    nonroot()
    require(sys.platform == 'linux', 'Linux wait4/seccomp execution required')
    candidate = real_directory(candidate)
    require(type(worker_args) in (list, tuple) and
            all(type(value) is str and value and '\0' not in value for value in worker_args),
            'Worker arguments must be finite strings')
    arguments = list(worker_args)
    require(sum(len(value.encode()) for value in arguments) <= 16 * 1024,
            'Worker arguments exceed custodian receipt bound')
    require(not any(value in ('--candidate', '--root', '--output', '--manufactured-only', '--fabricated-only')
                    or value.startswith(('--candidate=', '--root=', '--output=', '--manufactured-only=', '--fabricated-only='))
                    for value in arguments), 'Custodian-owned worker argument supplied')
    _authenticate_candidate(candidate, arguments, manufactured)
    worker = candidate / 'execution' / 'run_endpoint.py'
    require(worker.is_file() and not worker.is_symlink(), 'Authenticated worker required')
    output, directory_fd = _fresh_output(candidate, output)
    owned_fds.add(directory_fd)
    command = [sys.executable, '-I', '-B'] + (['-O'] if optimized else []) + [
        str(worker), '--root', str(candidate), '--output', str(output), *arguments]
    if manufactured:
        command.append('--manufactured-only')
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
            entry_identity = _validate_successful_worker(directory_fd, arguments, manufactured)
            _same_output_directory(directory_fd, output)
        except (ValueError, OSError) as exc:
            stop_reason = 'INVALID_ENTRY_OUTPUT'
            monitor_error = str(exc)
    passed = (code == 0 and stop_reason is None and elapsed <= WALL_SECONDS and
              usage.ru_maxrss * 1024 <= ADDRESS_SPACE_BYTES and
              usage.ru_utime + usage.ru_stime <= CPU_SECONDS)
    receipt = {
        'status': 'PASS_BOUNDED_MANUFACTURED_EXECUTION' if passed and manufactured else
                  'PASS_BOUNDED_REGISTERED_EXECUTION' if passed else 'EXECUTION_FAILED',
        'manufactured_only': bool(manufactured), 'optimized': bool(optimized),
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
    parser.add_argument('--optimized', action='store_true')
    parser.add_argument('--manufactured-only', action='store_true')
    parser.add_argument('worker_arguments', nargs=argparse.REMAINDER,
                        help='After --, supply registration and optional physical GO/repository arguments')
    args = parser.parse_args()
    arguments = args.worker_arguments
    require(arguments and arguments[0] == '--', 'Worker arguments must follow explicit -- separator')
    receipt = launch(args.candidate, args.output, arguments[1:],
                     args.optimized, args.manufactured_only)
    print(json.dumps(receipt, sort_keys=True), flush=True)
    raise SystemExit(receipt['wrapper_exit_code'])


if __name__ == '__main__':
    main()
