"""Offline custodian launcher with authoritative wait4 resources and status."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import sys
import time

from registration_guard import authenticate, require


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', required=True)
    parser.add_argument('--receipt', required=True)
    parser.add_argument('--receipt-sha256', required=True)
    parser.add_argument('--registration-sha256', required=True)
    parser.add_argument('--python', required=True)
    parser.add_argument('--route', choices=('primary', 'independent'), required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--optimized', action='store_true')
    parser.add_argument('--fabricated-only', action='store_true')
    args = parser.parse_args()
    require(os.getuid() != 0, 'Launch as the existing nonroot user')
    root = Path(args.root).absolute()
    verified = authenticate(root, args.receipt, args.receipt_sha256,
                            args.registration_sha256)
    out = Path(args.output_dir).absolute()
    require(root not in out.parents and not out.exists(), 'Exclusive external run directory required')
    out.mkdir(parents=True)
    limits = verified['contract']['resources']
    wall = limits['each_route_wall_seconds']
    rss = limits['each_route_peak_rss_kib']
    require(type(wall) is int and wall > 0 and type(rss) is int and rss > 0,
            'Positive exact outer limits required')
    rawlog = out/'child.log'
    output = out/'OUTPUT.json'
    cmd = [args.python, '-B'] + (['-O'] if args.optimized else []) + [
        str(root/'execution/run_registered.py'), '--root', str(root),
        '--receipt', str(Path(args.receipt).absolute()),
        '--receipt-sha256', args.receipt_sha256,
        '--registration-sha256', args.registration_sha256,
        '--route', args.route, '--output', str(output)]
    if args.fabricated_only:
        cmd.append('--fabricated-only')
    environment = {'PATH': '/usr/local/bin:/usr/bin:/bin', 'LANG': 'C.UTF-8',
        'PYTHONHASHSEED': '0', 'PYTHONDONTWRITEBYTECODE': '1',
        'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
        'NUMEXPR_NUM_THREADS': '1'}
    started = time.monotonic()
    with rawlog.open('xb') as log:
        pid = os.fork()
        if pid == 0:
            try:
                os.setsid()
                os.dup2(log.fileno(), 1)
                os.dup2(log.fileno(), 2)
                resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
                resource.setrlimit(resource.RLIMIT_AS, (rss*1024, rss*1024))
                resource.setrlimit(resource.RLIMIT_CPU, (wall, wall))
                os.chdir(out)
                os.execve(args.python, cmd, environment)
            except BaseException as exc:
                print('CHILD_LAUNCH_FAILED '+type(exc).__name__+': '+str(exc),
                      file=sys.stderr, flush=True)
                os._exit(126)
        timed_out = False
        while True:
            done, status, usage = os.wait4(pid, os.WNOHANG)
            if done:
                break
            if time.monotonic()-started > wall:
                timed_out = True
                os.killpg(pid, signal.SIGKILL)
                _, status, usage = os.wait4(pid, 0)
                break
            time.sleep(0.1)
    elapsed = time.monotonic()-started
    code = os.waitstatus_to_exitcode(status)
    completed = (code == 0 and not timed_out and output.is_file() and
                 usage.ru_maxrss <= rss and elapsed <= wall)
    receipt = {'status': 'PASS_BOUNDED_ROUTE_EXECUTION' if completed else 'EXECUTION_FAILED',
        'route': args.route, 'command': cmd, 'uid': os.getuid(),
        'fabricated_only': args.fabricated_only,
        'exit_code': code, 'wait_status': status, 'timed_out': timed_out,
        'wall_seconds': elapsed, 'peak_rss_kib_wait4': usage.ru_maxrss,
        'cpu_seconds_wait4': usage.ru_utime+usage.ru_stime,
        'outer_limits': {'wall_seconds': wall, 'peak_rss_kib': rss},
        'registration_sha256': args.registration_sha256,
        'remote_receipt_sha256': args.receipt_sha256,
        'child_log_sha256': hashlib.sha256(rawlog.read_bytes()).hexdigest(),
        'source_attempt_journal_sha256': hashlib.sha256(
            (out/'SOURCE_ATTEMPTS.jsonl').read_bytes()).hexdigest()
            if (out/'SOURCE_ATTEMPTS.jsonl').is_file() else None,
        'source_attempt_journal_bytes': (out/'SOURCE_ATTEMPTS.jsonl').stat().st_size
            if (out/'SOURCE_ATTEMPTS.jsonl').is_file() else 0,
        'output_sha256': hashlib.sha256(output.read_bytes()).hexdigest()
            if output.is_file() else None,
        'output_bytes': output.stat().st_size if output.is_file() else 0}
    (out/'EXECUTION.json').write_text(json.dumps(receipt, indent=2, sort_keys=True)+'\n')
    print(json.dumps(receipt), flush=True)
    if not completed:
        sys.exit(1)


if __name__ == '__main__':
    main()
