"""Nonroot offline custodian with authoritative wait4 outcome and hard limits."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import resource
import signal
import sys
import time

sys.dont_write_bytecode = True
from registration_guard import authenticate, authenticate_local_fabricated, require, runtime_identity

WALL_SECONDS = 900
RSS_KIB = 524288
OUTPUT_BYTES = 20971520

def output_size(directory):
    total = 0
    for path in directory.rglob('*'):
        require(not path.is_symlink(), 'Worker output symlink forbidden')
        if path.is_file():
            total += path.stat().st_size
    return total

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', required=True)
    parser.add_argument('--receipt')
    parser.add_argument('--receipt-sha256')
    parser.add_argument('--registration-sha256', required=True)
    parser.add_argument('--python', required=True)
    parser.add_argument('--route', choices=('primary', 'independent'), required=True)
    parser.add_argument('--output-directory', required=True)
    parser.add_argument('--optimized', action='store_true')
    parser.add_argument('--fabricated-only', action='store_true')
    args = parser.parse_args()
    require(os.getuid() != 0, 'Launch as the existing nonroot user')
    runtime = runtime_identity()
    if args.fabricated_only:
        require(args.receipt is None and args.receipt_sha256 is None, 'Local fabricated launch must not claim remote authorization')
        verified = authenticate_local_fabricated(args.root, args.registration_sha256)
    else:
        require(args.receipt is not None and args.receipt_sha256 is not None, 'Physical launch requires externally pinned root public receipt')
        verified = authenticate(args.root, args.receipt, args.receipt_sha256, args.registration_sha256)
    root = verified['root']
    python = Path(os.path.abspath(args.python))
    require(python.is_file() and os.access(python, os.X_OK), 'Explicit executable Python required')
    out = Path(os.path.abspath(args.output_directory))
    require(out != root and root not in out.parents and not out.exists(), 'Fresh external output directory required')
    require(all(not p.is_symlink() for p in (out, *out.parents)), 'Output symlink ancestor forbidden')
    out.mkdir(parents=True)
    cmd = [str(python), '-I', '-B'] + (['-O'] if args.optimized else []) + [
        str(root / 'execution/run_registered.py'), '--root', str(root),
        '--registration-sha256', args.registration_sha256, '--route', args.route,
        '--output-directory', str(out)]
    if args.fabricated_only:
        cmd.append('--fabricated-only')
    else:
        cmd.extend(['--receipt', str(Path(args.receipt).absolute()), '--receipt-sha256', args.receipt_sha256])
    environment = {'PATH': '/usr/local/bin:/usr/bin:/bin', 'LANG': 'C.UTF-8', 'PYTHONHASHSEED': '0',
                   'PYTHONDONTWRITEBYTECODE': '1', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1',
                   'MKL_NUM_THREADS': '1', 'NUMEXPR_NUM_THREADS': '1'}
    started = time.monotonic()
    stop_reason = None
    with (out / 'child.log').open('xb') as log:
        pid = os.fork()
        if pid == 0:
            try:
                os.setsid()
                os.dup2(log.fileno(), 1)
                os.dup2(log.fileno(), 2)
                resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
                resource.setrlimit(resource.RLIMIT_AS, (RSS_KIB * 1024, RSS_KIB * 1024))
                resource.setrlimit(resource.RLIMIT_CPU, (WALL_SECONDS, WALL_SECONDS))
                resource.setrlimit(resource.RLIMIT_FSIZE, (OUTPUT_BYTES, OUTPUT_BYTES))
                resource.setrlimit(resource.RLIMIT_NPROC, (0, 0))
                os.chdir(out)
                os.execve(str(python), cmd, environment)
            except BaseException as exc:
                print('CHILD_LAUNCH_FAILED ' + type(exc).__name__ + ': ' + str(exc), file=sys.stderr, flush=True)
                os._exit(126)
        while True:
            done, status, usage = os.wait4(pid, os.WNOHANG)
            if done:
                break
            if time.monotonic() - started > WALL_SECONDS:
                stop_reason = 'WALL_TIMEOUT'
            try:
                if output_size(out) > OUTPUT_BYTES:
                    stop_reason = 'OUTPUT_LIMIT'
            except ValueError:
                stop_reason = 'UNSAFE_OUTPUT'
            if stop_reason:
                try:
                    os.killpg(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                _, status, usage = os.wait4(pid, 0)
                break
            time.sleep(0.1)
    elapsed = time.monotonic() - started
    code = os.waitstatus_to_exitcode(status)
    try:
        data_bytes = output_size(out)
    except ValueError:
        data_bytes = OUTPUT_BYTES + 1
        stop_reason = 'UNSAFE_OUTPUT'
    required_outputs = ['DATA.json', 'BUDGET.json', 'OUTPUT_REVIEW.json', 'ENTRY_RECEIPT.json', 'SOURCE_ATTEMPTS.jsonl']
    completed = code == 0 and stop_reason is None and elapsed <= WALL_SECONDS and usage.ru_maxrss <= RSS_KIB and data_bytes <= OUTPUT_BYTES and all((out / name).is_file() for name in required_outputs)
    receipt = {'status': 'PASS_BOUNDED_REGISTERED_ENTRY' if completed else 'EXECUTION_FAILED',
               'route': args.route, 'fabricated_only': args.fabricated_only, 'optimized': args.optimized,
               'exit_code': code, 'wait_status': status, 'stop_reason': stop_reason,
               'wall_seconds': elapsed, 'peak_rss_kib_wait4': usage.ru_maxrss,
               'cpu_seconds_wait4': usage.ru_utime + usage.ru_stime, 'worker_output_bytes': data_bytes,
               'outer_limits': {'wall_seconds': WALL_SECONDS, 'peak_rss_kib': RSS_KIB, 'output_bytes': OUTPUT_BYTES},
               'uid': os.getuid(), 'command': cmd, 'registration_sha256': args.registration_sha256,
               'runtime': runtime,
               'remote_receipt_sha256': args.receipt_sha256,
               'files': {name: {'bytes': (out / name).stat().st_size, 'sha256': hashlib.sha256((out / name).read_bytes()).hexdigest()} for name in ['child.log', *required_outputs] if (out / name).is_file() and not (out / name).is_symlink()}}
    with (out / 'EXECUTION.json').open('x') as handle:
        handle.write(json.dumps(receipt, sort_keys=True, indent=2) + '\n')
    print(json.dumps(receipt, sort_keys=True), flush=True)
    # Preserve the authoritative child outcome, including signal-style 128+n.
    if code:
        sys.exit(code if code > 0 else 128 - code)
    if not completed:
        sys.exit(1)

if __name__ == '__main__':
    main()
