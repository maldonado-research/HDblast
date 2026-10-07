"""Offline numeric entry: verify the public freeze before loading either route."""
from pathlib import Path
import argparse
import ctypes
import errno
import importlib.util
import json
import os
import resource
import sys

# The outer launcher has checked this complete directory before interpreter
# startup. The worker repeats authentication before any numeric route import.
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from registration_guard import SourceAuthorization, authenticate, authenticate_local_fabricated, require, runtime_identity, safe_file, sha
from validate_outputs import validate_output

def deny_network_and_children():
    """Fail-closed Linux seccomp policy denies sockets, process and thread spawn."""
    require(sys.platform == 'linux', 'Registered worker requires Linux seccomp')
    # Use the runtime's system dynamic loader directly; find_library can spawn
    # a subprocess and is incompatible with the no-child resource contract.
    lib = ctypes.CDLL('libseccomp.so.2', use_errno=True)
    lib.seccomp_init.argtypes = [ctypes.c_uint32]
    lib.seccomp_init.restype = ctypes.c_void_p
    lib.seccomp_syscall_resolve_name.argtypes = [ctypes.c_char_p]
    lib.seccomp_syscall_resolve_name.restype = ctypes.c_int
    lib.seccomp_rule_add.argtypes = [ctypes.c_void_p, ctypes.c_uint32, ctypes.c_int, ctypes.c_uint]
    lib.seccomp_rule_add.restype = ctypes.c_int
    lib.seccomp_load.argtypes = [ctypes.c_void_p]
    lib.seccomp_load.restype = ctypes.c_int
    lib.seccomp_release.argtypes = [ctypes.c_void_p]
    context = lib.seccomp_init(0x7fff0000)
    require(context, 'seccomp allocation failed')
    denied = []
    try:
        action = 0x00050000 | errno.EPERM
        for name in ('socket', 'socketpair', 'connect', 'bind', 'listen', 'accept', 'accept4', 'fork', 'vfork', 'clone', 'clone3', 'execve', 'execveat'):
            number = lib.seccomp_syscall_resolve_name(name.encode())
            if number >= 0:
                require(lib.seccomp_rule_add(context, action, number, 0) == 0, 'seccomp rule failed: ' + name)
                denied.append(name)
        require({'socket', 'connect', 'clone', 'execve'} <= set(denied), 'Essential offline/no-child syscalls unresolved')
        require(lib.seccomp_load(context) == 0, 'seccomp policy load failed')
    finally:
        lib.seccomp_release(context)
    return denied

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', required=True)
    parser.add_argument('--receipt')
    parser.add_argument('--receipt-sha256')
    parser.add_argument('--registration-sha256', required=True)
    parser.add_argument('--route', choices=('primary', 'independent'), required=True)
    parser.add_argument('--output-directory', required=True)
    parser.add_argument('--fabricated-only', action='store_true')
    args = parser.parse_args()
    require(os.getuid() != 0, 'Numeric worker must be nonroot')
    runtime = runtime_identity()
    for limit, expected in [(resource.RLIMIT_AS, 524288 * 1024), (resource.RLIMIT_CPU, 900),
                            (resource.RLIMIT_FSIZE, 20971520), (resource.RLIMIT_NPROC, 0), (resource.RLIMIT_CORE, 0)]:
        require(resource.getrlimit(limit) == (expected, expected), 'Worker must inherit the exact custodian resource limits')
    require(all(os.environ.get(name) == '1' for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')), 'Single-thread environment required')
    if args.fabricated_only:
        require(args.receipt is None and args.receipt_sha256 is None, 'Local fabricated entry must not claim remote authorization')
        verified = authenticate_local_fabricated(args.root, args.registration_sha256)
    else:
        require(args.receipt is not None and args.receipt_sha256 is not None, 'Physical entry requires externally pinned root public receipt')
        verified = authenticate(args.root, args.receipt, args.receipt_sha256, args.registration_sha256)
    root = verified['root']
    output = Path(os.path.abspath(args.output_directory))
    require(output.is_dir() and output.resolve() == output and output != root and root not in output.parents, 'Existing real external worker output directory required')
    require(all(not path.is_symlink() for path in (output, *output.parents)), 'Output symlink forbidden')
    denied = deny_network_and_children()
    relative = args.route + '/route.py'
    path = safe_file(root, relative)
    # Every Python helper below root was already pinned and independently
    # reviewed. Only this verified route directory is added for sibling imports.
    sys.path.insert(0, str(path.parent))
    spec = importlib.util.spec_from_file_location('registered_' + args.route, path)
    require(spec is not None and spec.loader is not None, 'Frozen route load failed')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    journal = output / 'SOURCE_ATTEMPTS.jsonl'
    if args.fabricated_only:
        journal.open('xb').close()
        data, budget = module.run_fabricated()
        callbacks = {'physical_source_callbacks': 0, 'archive_arrays_decoded': 0, 'events': []}
    else:
        authorization = SourceAuthorization(verified, args.route, journal)
        data, budget = module.run(authorization)
        callbacks = authorization.complete()
    review, unused_boxes = validate_output(data, budget, args.route, args.fabricated_only, enforce_gate=True)
    if args.fabricated_only:
        authenticate_local_fabricated(root, args.registration_sha256)
    else:
        authenticate(root, args.receipt, args.receipt_sha256, args.registration_sha256)
    bodies = {name: (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False) + '\n').encode() for name, value in [('DATA.json', data), ('BUDGET.json', budget), ('OUTPUT_REVIEW.json', review)]}
    require(sum(len(raw) for raw in bodies.values()) <= 20 * 1024 * 1024, 'Combined worker JSON exceeds 20 MiB')
    for name, raw in bodies.items():
        with (output / name).open('xb') as handle:
            handle.write(raw)
    audit = {'status': 'PASS_FABRICATED_ENTRY' if args.fabricated_only else 'PASS_AUTHENTICATED_OFFLINE_REGISTERED_ENTRY',
             'route': args.route, 'fabricated_only': args.fabricated_only,
             'freeze_commit': None if args.fabricated_only else verified['receipt']['freeze_commit'],
             'authorization_scope': 'LOCAL_FABRICATED_ONLY_NO_PUBLIC_FREEZE' if args.fabricated_only else verified['receipt']['input_scope'],
             'registration_sha256': args.registration_sha256,
             'receipt_sha256': args.receipt_sha256, 'uid': os.getuid(), 'denied_syscalls': denied,
             'runtime': runtime,
             'source_callbacks': callbacks, 'files': {name: {'bytes': len(raw), 'sha256': sha(raw)} for name, raw in bodies.items()},
             'original_binary80_state_error': 'NOT_ENCLOSED', 'full_twelve_case_certificate': 'UNRESOLVED'}
    with (output / 'ENTRY_RECEIPT.json').open('x') as handle:
        handle.write(json.dumps(audit, sort_keys=True, indent=2) + '\n')
    print(audit['status'], flush=True)

if __name__ == '__main__':
    main()
