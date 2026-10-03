#!/usr/bin/env python3
"""Authenticate and reproduce the complete registered active-source diagnostic.

Execution receipts describe actual child exits/resources. A completed replay
may report any registered scientific classification; earlier metric FAILs stay.
--plan-only authenticates the exact plan and launches no child processes.
"""
from __future__ import annotations

import argparse
import ctypes
from datetime import datetime, timezone
import importlib.metadata
import math
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import sys
import time
import traceback

from active_integrity import (MEMORY_LIMIT_KIB, TIME_LIMIT_SECONDS, VERSIONS,
                              file_sha, fresh_external, read_json, require, safe_file,
                              verify_frozen, write_json)

CLASSIFICATIONS = {'LEDGER_ERROR_DEMONSTRATED', 'NO_GATE_SCALE_ATTRIBUTION', 'CONSISTENCY_FAILURE'}


def stamp():
    return datetime.now(timezone.utc).isoformat()


def runtime_check():
    require(sys.flags.optimize == 0, 'Physical replay driver must run normal Python')
    require(sys.version_info[:3] == (3, 12, 14), 'Python 3.12.14 required')
    for name, version in VERSIONS.items():
        require(importlib.metadata.version(name) == version, 'Pinned dependency mismatch: ' + name)


def command_plan(work, fresh, inputs):
    experiment = read_json(safe_file(work, 'EXPERIMENT.json'))
    reproduction = experiment['reproduction']
    require(set(reproduction['routes']) == {'primary', 'independent'}, 'Exactly two routes required')
    require(reproduction['result_filename'] == 'diagnostic.json', 'Diagnostic filename differs')
    require(reproduction['validation_filename'] == 'VALIDATION.json', 'Validation filename differs')
    registered = inputs['frozen_files']
    pin_args = ['--checkpoint-root', str(work), '--registration-sha256', inputs['registration_sha256'],
                '--freeze-commit', inputs['public_freeze_commit']]
    commands = []
    for name in ('primary', 'independent'):
        entrypoint = reproduction['routes'][name]
        require(entrypoint in registered, 'Route omitted from registration: ' + entrypoint)
        commands.append({'name': name, 'physical_route': True,
                         'command': [sys.executable, str(safe_file(work, entrypoint)), *pin_args,
                                     '--output-dir', str(fresh / name)],
                         'timeout_seconds': TIME_LIMIT_SECONDS, 'memory_limit_kib': MEMORY_LIMIT_KIB})
    validator = reproduction['validator']
    require(validator in registered, 'Validator omitted from registration')
    for optimized in (False, True):
        name = 'validation_optimized' if optimized else 'validation'
        commands.append({'name': name, 'physical_route': False,
                         'command': [sys.executable, *(['-O'] if optimized else []),
                                     str(safe_file(work, validator)), *pin_args,
                                     '--primary', str(fresh / 'primary' / 'diagnostic.json'),
                                     '--independent', str(fresh / 'independent' / 'diagnostic.json'),
                                     '--execution-receipts', str(fresh.parent / 'EXECUTION.json'),
                                     '--output-dir', str(fresh / name)],
                         'timeout_seconds': TIME_LIMIT_SECONDS, 'memory_limit_kib': MEMORY_LIMIT_KIB})
    require(all('-O' not in c['command'] for c in commands if c['physical_route']),
            'Physical routes must run normal Python')
    return commands


def group_rss_kib(process_group):
    """Best-effort live Linux group RSS; wait4 peak RSS remains authoritative."""
    total = 0
    for entry in Path('/proc').iterdir():
        if not entry.name.isdigit():
            continue
        try:
            fields = (entry / 'stat').read_text().rsplit(')', 1)[1].split()
            if int(fields[2]) != process_group:
                continue
            for line in (entry / 'status').read_text().splitlines():
                if line.startswith('VmRSS:'):
                    total += int(line.split()[1])
                    break
        except (FileNotFoundError, ProcessLookupError, PermissionError):
            continue
    return total


def child_pids(parent):
    result = []
    for entry in Path('/proc').iterdir():
        if entry.name.isdigit():
            try:
                fields = (entry / 'stat').read_text().rsplit(')', 1)[1].split()
                if int(fields[1]) == parent:
                    result.append(int(entry.name))
            except (FileNotFoundError, ProcessLookupError, PermissionError):
                pass
    return result


def drain_descendants():
    """As a Linux subreaper, kill/reap orphans even if they escaped the group."""
    receipts = []
    deadline = time.monotonic() + 2
    while True:
        children = child_pids(os.getpid())
        for pid in children:
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        while True:
            try:
                pid, status, usage = os.wait4(-1, os.WNOHANG)
            except ChildProcessError:
                pid = 0
            if not pid:
                break
            receipts.append({'pid': pid, 'exit_code': os.waitstatus_to_exitcode(status),
                             'peak_rss_kib': int(usage.ru_maxrss),
                             'user_seconds': usage.ru_utime, 'system_seconds': usage.ru_stime})
        if not child_pids(os.getpid()):
            return receipts
        require(time.monotonic() < deadline, 'Could not reap all route descendants')
        time.sleep(0.01)


def run_command(command, logs, env):
    """Run once, preserve logs, kill the full process group on a budget breach."""
    record = {**command, 'started_utc': stamp(), 'status': 'RUNNING',
              'exit_code': None, 'timed_out': False, 'memory_exceeded': False,
              'process_creation_limit': 0, 'executor_uid': os.getuid()}
    started = time.monotonic()
    process = None
    usage = None
    peak_group_rss = 0
    def limits():
        # Frozen routes use a single process and one numeric thread. Prohibit
        # descendants/extra threads in-kernel rather than miss aggregate RSS.
        resource.setrlimit(resource.RLIMIT_NPROC, (0, 0))
        resource.setrlimit(resource.RLIMIT_CPU, (math.ceil(command['timeout_seconds']),
                                               math.ceil(command['timeout_seconds']) + 1))
    with (logs / (command['name'] + '.stdout.log')).open('xb') as stdout, \
         (logs / (command['name'] + '.stderr.log')).open('xb') as stderr:
        try:
            require(os.getuid() != 0, 'Nonroot Linux execution is required for RLIMIT_NPROC isolation')
            # PR_SET_CHILD_SUBREAPER: orphaned grandchildren belong to this driver.
            libc = ctypes.CDLL(None, use_errno=True)
            require(libc.prctl(36, 1, 0, 0, 0) == 0, 'Linux child-subreaper isolation unavailable')
            process = subprocess.Popen(command['command'], cwd=command['cwd'], env=env,
                                       stdout=stdout, stderr=stderr, start_new_session=True,
                                       preexec_fn=limits)
            record['pid'] = process.pid
            while True:
                waited, status, usage = os.wait4(process.pid, os.WNOHANG)
                if waited:
                    process.returncode = os.waitstatus_to_exitcode(status)
                    record['exit_code'] = process.returncode
                    survivors = child_pids(os.getpid())
                    if survivors:
                        record['surviving_descendant_pids'] = survivors
                        record['descendant_exit_records'] = drain_descendants()
                    break
                usage = None
                peak_group_rss = max(peak_group_rss, group_rss_kib(process.pid))
                elapsed = time.monotonic() - started
                record['timed_out'] = elapsed > command['timeout_seconds']
                record['memory_exceeded'] = peak_group_rss > command['memory_limit_kib']
                if record['timed_out'] or record['memory_exceeded']:
                    try:
                        os.killpg(process.pid, signal.SIGKILL)
                    except ProcessLookupError:
                        pass
                    _, status, usage = os.wait4(process.pid, 0)
                    process.returncode = os.waitstatus_to_exitcode(status)
                    record['exit_code'] = process.returncode
                    break
                time.sleep(0.02)
            record['peak_rss_kib'] = max(int(usage.ru_maxrss), peak_group_rss)
            record['wait4_peak_rss_kib'] = int(usage.ru_maxrss)
            record['observed_peak_group_rss_kib'] = peak_group_rss
            record['user_seconds'] = usage.ru_utime
            record['system_seconds'] = usage.ru_stime
            record['memory_exceeded'] |= record['peak_rss_kib'] > command['memory_limit_kib']
        except BaseException as exc:
            record['exception'] = repr(exc)
            if process is not None and process.returncode is None:
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                    _, status, usage = os.wait4(process.pid, 0)
                    process.returncode = os.waitstatus_to_exitcode(status)
                    record['exit_code'] = process.returncode
                except ProcessLookupError:
                    pass
        finally:
            if process is not None and child_pids(os.getpid()):
                record['surviving_descendant_pids'] = child_pids(os.getpid())
                try:
                    record['descendant_exit_records'] = drain_descendants()
                except RuntimeError as exc:
                    record['descendant_cleanup_error'] = repr(exc)
            record['elapsed_seconds'] = time.monotonic() - started
            record['completed_utc'] = stamp()
            record['timed_out'] |= record['elapsed_seconds'] > command['timeout_seconds']
            record['status'] = ('PASS_EXECUTION' if record['exit_code'] == 0
                                and not record['timed_out'] and not record['memory_exceeded']
                                and 'exception' not in record and 'surviving_descendant_pids' not in record
                                else 'FAIL_EXECUTION')
            write_json(logs / (command['name'] + '.exit.json'), record)
    return record


def compare_verdicts(fresh):
    reports = [read_json(fresh / name / 'VALIDATION.json')
               for name in ('validation', 'validation_optimized')]
    for report in reports:
        require(report.get('classification') in CLASSIFICATIONS, 'Missing authoritative diagnostic classification')
        require(report.get('old_metric_status') == 'FAIL', 'Earlier metric FAIL must remain explicit')
    normalized = [{k: v for k, v in report.items() if k != 'python_optimization'} for report in reports]
    require(normalized[0] == normalized[1], 'Normal/optimized validator reports differ')
    return reports[0]['classification']


def replay(root, output, registration_sha256, freeze_commit, plan_only=False):
    require(sys.flags.optimize == 0, 'Physical replay driver must run normal Python')
    root = Path(root).absolute()
    output = fresh_external(root, output)
    output.mkdir(parents=True, exist_ok=False)
    work, fresh, logs = output / 'checkpoint', output / 'fresh', output / 'logs'
    state = {'status': 'RUNNING', 'started_utc': stamp(), 'driver_sha256': file_sha(__file__),
             'commands': [], 'physical_routes_planned': 2, 'old_metric_status': 'FAIL',
             'classification': 'UNCOMPUTED', 'plan_only': plan_only}
    inputs = None
    write_json(output / 'REPLAY.json', state)
    try:
        runtime_check()
        inputs = verify_frozen(root, registration_sha256, freeze_commit)
        write_json(output / 'INPUT_MANIFEST.json', inputs)
        state['source_verification_before'] = 'PASS'
        shutil.copytree(root, work, ignore=shutil.ignore_patterns('__pycache__', '*.pyc', '*.pyo'))
        require(verify_frozen(work, registration_sha256, freeze_commit) == inputs,
                'Copied source authentication differs')
        fresh.mkdir(); logs.mkdir()
        commands = [{**command, 'cwd': str(work)} for command in command_plan(work, fresh, inputs)]
        write_json(output / 'COMMAND_PLAN.json', commands)
        state['planned_commands'] = len(commands)
        if plan_only:
            state['status'] = 'VERIFIED_PLAN_ONLY'
        else:
            env = os.environ.copy()
            env.update(PYTHONDONTWRITEBYTECODE='1', PYTHONOPTIMIZE='0', OPENBLAS_NUM_THREADS='1',
                       OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1',
                       VECLIB_MAXIMUM_THREADS='1', BLIS_NUM_THREADS='1', MPLBACKEND='Agg')
            for command in commands:
                if not command['physical_route'] and any(c['status'] != 'PASS_EXECUTION'
                                                         for c in state['commands'] if c['physical_route']):
                    state.setdefault('skipped_commands', []).append({'name': command['name'],
                                                                     'reason': 'Physical route execution failed'})
                    continue
                record = run_command(command, logs, env)
                if command['physical_route'] and record['status'] == 'PASS_EXECUTION':
                    try:
                        safe_file(fresh / command['name'], 'diagnostic.json')
                    except RuntimeError as exc:
                        record.update(status='FAIL_EXECUTION', result_error=repr(exc))
                        write_json(logs / (command['name'] + '.exit.json'), record)
                state['commands'].append(record)
                write_json(output / 'EXECUTION.json', state['commands'])
                write_json(output / 'REPLAY.json', state)
                require(verify_frozen(root, registration_sha256, freeze_commit) == inputs,
                        'Original frozen inputs changed during ' + command['name'])
                require(verify_frozen(work, registration_sha256, freeze_commit) == inputs,
                        'Copied frozen inputs changed during ' + command['name'])
                print(command['name'] + ' ' + record['status'], flush=True)
            require(len(state['commands']) == len(commands)
                    and all(c['status'] == 'PASS_EXECUTION' for c in state['commands']),
                    'Replay command failed; actual logs, receipts and progress retained')
            state['classification'] = compare_verdicts(fresh)
            state['status'] = 'COMPLETED_REPLAY'
    except BaseException as exc:
        state.update(status='FAIL_REPLAY', failure={'exception': repr(exc), 'traceback': traceback.format_exc()})
    finally:
        try:
            if inputs is not None:
                require(verify_frozen(root, registration_sha256, freeze_commit) == inputs,
                        'Original source changed during replay')
                if work.exists():
                    require(verify_frozen(work, registration_sha256, freeze_commit) == inputs,
                            'Copied source changed during replay')
                state['source_verification_after'] = 'PASS'
            else:
                state['source_verification_after'] = 'NOT_STARTED'
        except BaseException as exc:
            state.update(status='FAIL_REPLAY', source_verification_after='FAIL',
                         source_verification_failure=repr(exc))
        state.update(completed_utc=stamp(), attempted_commands=len(state['commands']),
                     successful_commands=sum(c['status'] == 'PASS_EXECUTION' for c in state['commands']),
                     physical_routes_attempted=sum(c['physical_route'] for c in state['commands']),
                     physical_routes_completed=sum(c['physical_route'] and c['status'] == 'PASS_EXECUTION'
                                                   for c in state['commands']))
        if fresh.exists():
            try:
                output_files = {}
                for path in sorted(fresh.rglob('*')):
                    require(not path.is_symlink(), 'Fresh output contains a symlink: ' + str(path))
                    if path.is_file():
                        name = path.relative_to(fresh).as_posix()
                        output_files[name] = file_sha(safe_file(fresh, name))
                state['fresh_output_sha256'] = output_files
            except RuntimeError as exc:
                state.update(status='FAIL_REPLAY', output_integrity_failure=repr(exc))
        if state['status'] == 'FAIL_REPLAY':
            if state['classification'] != 'UNCOMPUTED':
                state['unaccepted_validator_classification'] = state['classification']
            state['classification'] = 'UNCOMPUTED'
        write_json(output / 'REPLAY.json', state)
    return state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkpoint-root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--registration-sha256', required=True)
    parser.add_argument('--freeze-commit', required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--plan-only', action='store_true')
    args = parser.parse_args()
    state = replay(args.checkpoint_root, args.output_dir, args.registration_sha256,
                   args.freeze_commit, args.plan_only)
    print({k: state[k] for k in ('status', 'classification', 'attempted_commands', 'successful_commands')})
    return 0 if state['status'] in ('COMPLETED_REPLAY', 'VERIFIED_PLAN_ONLY') else 1


if __name__ == '__main__':
    raise SystemExit(main())
