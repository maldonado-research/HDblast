"""Linux kernel restrictions and bounded, symlink-free output inspection.

This module uses only the standard library. ``install`` must run after all
numeric imports and before either fabricated work or an authorized source
callback. Failure to load the kernel policy is fatal; there is no fallback.
"""
from pathlib import Path
import ctypes
import errno
import hashlib
import os
import resource
import stat
import sys
import time

WALL_SECONDS = 900
CPU_SECONDS = 900
ADDRESS_SPACE_BYTES = 512 * 1024 * 1024
OUTPUT_BYTES = 128 * 1024 * 1024
MAX_OUTPUT_ENTRIES = 4096
MAX_OUTPUT_DEPTH = 32
THREAD_VARIABLES = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                    'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'BLIS_NUM_THREADS')
DENIED_SYSCALLS = (
    'socket', 'socketpair', 'connect', 'bind', 'listen', 'accept', 'accept4',
    'fork', 'vfork', 'clone', 'clone3', 'execve', 'execveat',
    # io_uring can otherwise issue socket/connect operations internally.
    'io_uring_setup', 'io_uring_enter', 'io_uring_register',
)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def nonroot():
    require(os.getuid() != 0 and os.geteuid() != 0, 'Nonroot execution required')
    require(os.getuid() == os.geteuid() and os.getgid() == os.getegid(),
            'Set-ID execution forbidden')
    # NPROC is not enforced for processes with either of these capabilities.
    if sys.platform == 'linux':
        for line in Path('/proc/self/status').read_text().splitlines():
            if line.startswith('CapEff:'):
                capabilities = int(line.split()[1], 16)
                require(not capabilities & ((1 << 21) | (1 << 24)),
                        'Resource-exempt capabilities forbidden')
                break
        else:
            raise ValueError('Linux effective capability status unavailable')


def no_new_privileges():
    require(sys.platform == 'linux', 'Linux kernel restrictions required')
    libc = ctypes.CDLL(None, use_errno=True)
    libc.prctl.argtypes = [ctypes.c_int, ctypes.c_ulong, ctypes.c_ulong,
                          ctypes.c_ulong, ctypes.c_ulong]
    libc.prctl.restype = ctypes.c_int
    require(libc.prctl(38, 1, 0, 0, 0) == 0, 'PR_SET_NO_NEW_PRIVS failed')


def require_worker_limits():
    nonroot()
    for limit, expected in ((resource.RLIMIT_AS, ADDRESS_SPACE_BYTES),
                            (resource.RLIMIT_CPU, CPU_SECONDS),
                            (resource.RLIMIT_FSIZE, OUTPUT_BYTES),
                            (resource.RLIMIT_NPROC, 0), (resource.RLIMIT_CORE, 0)):
        require(resource.getrlimit(limit) == (expected, expected),
                'Exact inherited resource limits required')
    require(all(os.environ.get(name) == '1' for name in THREAD_VARIABLES),
            'Single-thread environment required')


def install():
    """Deny networking and process/thread creation in the kernel, or fail."""
    nonroot()
    require(sys.platform == 'linux', 'Linux seccomp required')
    require(len(os.listdir('/proc/self/task')) == 1,
            'Numeric imports must leave exactly one thread before seccomp')
    no_new_privileges()
    # Loading by soname avoids find_library's possible child process.
    lib = ctypes.CDLL('libseccomp.so.2', use_errno=True)
    lib.seccomp_init.argtypes = [ctypes.c_uint32]
    lib.seccomp_init.restype = ctypes.c_void_p
    lib.seccomp_syscall_resolve_name.argtypes = [ctypes.c_char_p]
    lib.seccomp_syscall_resolve_name.restype = ctypes.c_int
    lib.seccomp_rule_add.argtypes = [ctypes.c_void_p, ctypes.c_uint32,
                                   ctypes.c_int, ctypes.c_uint]
    lib.seccomp_rule_add.restype = ctypes.c_int
    lib.seccomp_load.argtypes = [ctypes.c_void_p]
    lib.seccomp_load.restype = ctypes.c_int
    lib.seccomp_release.argtypes = [ctypes.c_void_p]
    context = lib.seccomp_init(0x7fff0000)  # default ALLOW
    require(bool(context), 'seccomp allocation failed')
    denied = []
    try:
        action = 0x00050000 | errno.EPERM
        for name in DENIED_SYSCALLS:
            number = lib.seccomp_syscall_resolve_name(name.encode('ascii'))
            # Nonexistent architecture-specific calls cannot be invoked. All
            # essential current-Linux calls must nevertheless resolve.
            if number >= 0:
                require(lib.seccomp_rule_add(context, action, number, 0) == 0,
                        'seccomp rule failed: ' + name)
                denied.append(name)
        require({'socket', 'socketpair', 'connect', 'clone', 'execve',
                 'io_uring_setup', 'io_uring_enter', 'io_uring_register'} <= set(denied),
                'Essential offline/no-child syscalls unresolved')
        require(lib.seccomp_load(context) == 0, 'seccomp policy load failed')
    finally:
        lib.seccomp_release(context)
    return denied


def real_directory(directory):
    path = Path(os.path.abspath(os.fspath(directory)))
    for ancestor in reversed((path, *path.parents)):
        require(stat.S_ISDIR(os.lstat(ancestor).st_mode),
                'Real directory without symlink ancestors required')
    return path


def open_real_directory(directory):
    """Open every absolute path component without following any symlink."""
    path = Path(os.path.abspath(os.fspath(directory)))
    flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
    fd = os.open('/', flags)
    try:
        for part in path.parts[1:]:
            next_fd = os.open(part, flags, dir_fd=fd)
            os.close(fd)
            fd = next_fd
        return path, fd
    except OSError as exc:
        os.close(fd)
        raise ValueError('Real directory without symlink ancestors required: ' + str(exc)) from exc
    except BaseException:
        os.close(fd)
        raise


def check_output(directory, max_bytes=OUTPUT_BYTES, hash_files=False,
                 forbidden_roots=(), deadline=None, directory_fd=None):
    """Inspect regular output files using dirfds; reject links/special files.

    During a live run this returns a bounded size snapshot. After wait4, use
    ``hash_files=True`` to require stable bytes and obtain a tree SHA256.
    ``deadline`` is an absolute monotonic deadline for potentially large scans.
    ``directory_fd`` pins the scan to an already held output root. This function
    owns only a duplicate and also verifies that the pathname names that inode.
    """
    require(type(max_bytes) is int and max_bytes >= 0, 'Finite output cap required')
    root = real_directory(directory)
    for forbidden in forbidden_roots:
        forbidden = real_directory(forbidden)
        require(root != forbidden and forbidden not in root.parents,
                'Worker output must be external to candidate')
    flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
    summary = {'bytes': 0, 'file_count': 0, 'directory_count': 0}
    entries = 0
    tree = hashlib.sha256()
    seen = set()
    root_fd = None
    try:
        if directory_fd is None:
            root, root_fd = open_real_directory(root)
        else:
            require(type(directory_fd) is int and directory_fd >= 0,
                    'Held output directory descriptor required')
            root_fd = os.dup(directory_fd)
        root_stat = os.fstat(root_fd)
        require(stat.S_ISDIR(root_stat.st_mode), 'Output descriptor must be a directory')
    except BaseException:
        if root_fd is not None:
            os.close(root_fd)
        raise

    def require_named_root():
        unused_root, named_fd = open_real_directory(root)
        try:
            named = os.fstat(named_fd)
            require((named.st_dev, named.st_ino) == (root_stat.st_dev, root_stat.st_ino),
                    'Output pathname does not name held root directory')
        finally:
            os.close(named_fd)

    def scan(fd, prefix, depth):
        nonlocal entries
        require(depth <= MAX_OUTPUT_DEPTH, 'Output nesting exceeds bound')
        names = []
        with os.scandir(fd) as iterator:
            for entry in iterator:
                require(deadline is None or time.monotonic() <= deadline,
                        'Output inspection exceeded deadline')
                entries += 1
                require(entries <= MAX_OUTPUT_ENTRIES, 'Output entry count exceeds bound')
                names.append(entry.name)
        for name in sorted(names):
            metadata = os.stat(name, dir_fd=fd, follow_symlinks=False)
            require(metadata.st_dev == root_stat.st_dev, 'Output mount crossing forbidden')
            require(metadata.st_uid == os.getuid(), 'Output owner differs from custodian')
            identity = (metadata.st_dev, metadata.st_ino)
            require(identity not in seen, 'Aliased output inode forbidden')
            seen.add(identity)
            relative = prefix + name
            if stat.S_ISDIR(metadata.st_mode):
                child_fd = os.open(name, flags, dir_fd=fd)
                try:
                    current = os.fstat(child_fd)
                    require((current.st_dev, current.st_ino) == identity,
                            'Output directory changed during inspection')
                    summary['directory_count'] += 1
                    scan(child_fd, relative + '/', depth + 1)
                finally:
                    os.close(child_fd)
            else:
                require(stat.S_ISREG(metadata.st_mode),
                        'Only regular output files and directories allowed')
                require(metadata.st_nlink == 1, 'Output hardlink forbidden')
                file_fd = os.open(name, os.O_RDONLY | os.O_NOFOLLOW |
                                  os.O_NONBLOCK | os.O_CLOEXEC, dir_fd=fd)
                try:
                    before = os.fstat(file_fd)
                    require(stat.S_ISREG(before.st_mode) and before.st_nlink == 1 and
                            (before.st_dev, before.st_ino) == identity,
                            'Output file changed during inspection')
                    summary['bytes'] += before.st_size
                    summary['file_count'] += 1
                    require(summary['bytes'] <= max_bytes, 'Aggregate output exceeds bound')
                    if hash_files:
                        digest = hashlib.sha256()
                        read_bytes = 0
                        while True:
                            require(deadline is None or time.monotonic() <= deadline,
                                    'Output hashing exceeded deadline')
                            raw = os.read(file_fd, 1024 * 1024)
                            if not raw:
                                break
                            read_bytes += len(raw)
                            require(read_bytes <= before.st_size,
                                    'Output grew during final inspection')
                            digest.update(raw)
                        after = os.fstat(file_fd)
                        require(read_bytes == before.st_size and
                                (before.st_size, before.st_mtime_ns, before.st_ctime_ns) ==
                                (after.st_size, after.st_mtime_ns, after.st_ctime_ns),
                                'Output changed during final inspection')
                        tree.update(relative.encode('utf-8', 'surrogateescape') + b'\0' +
                                    str(read_bytes).encode() + b'\0' + digest.digest())
                finally:
                    os.close(file_fd)
    try:
        require_named_root()
        scan(root_fd, '', 0)
        require_named_root()
    except OSError as exc:
        raise ValueError('Unsafe or inaccessible output: ' + str(exc)) from exc
    finally:
        os.close(root_fd)
    if hash_files:
        summary['tree_sha256'] = tree.hexdigest()
    return summary
