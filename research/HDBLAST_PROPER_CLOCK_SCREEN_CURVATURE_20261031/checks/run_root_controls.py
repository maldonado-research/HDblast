#!/usr/bin/env python3
"""Strict root provenance-corruption subprocess controls, without proof replay."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
CAP = HERE.parent
SENTINEL = 'ROOT_CONTROL_SENTINEL_MUST_NOT_RUN'


def need(ok, witness):
    if not ok:
        raise ValueError(witness)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_bytes())


def relative(value):
    need(isinstance(value, str) and bool(value), 'NONEMPTY_RELATIVE_PATH')
    path = PurePosixPath(value)
    need(not path.is_absolute() and '..' not in path.parts and path.as_posix() == value
         and value != '.', 'SAFE_RELATIVE_PATH')
    return Path(*path.parts)


def parsed_failure(proc):
    lines = proc.stderr.strip().splitlines()
    line = lines[-1] if lines else ''
    prefix = 'ValueError: '
    return line[len(prefix):] if line.startswith(prefix) else None


def successful_stdout(proc):
    for line in proc.stdout.splitlines():
        try:
            value = json.loads(line)
        except (ValueError, TypeError):
            continue
        if isinstance(value, dict):
            status = value.get('status', '')
            if isinstance(status, str) and status.startswith(('PASS', 'CERTIFIED')):
                return True
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo-root', required=True, type=Path)
    ap.add_argument('--output', required=True, type=Path)
    args = ap.parse_args()
    need(not args.output.exists() and not args.output.is_symlink(), 'FRESH_CONTROL_OUTPUT')
    need(all(not p.is_symlink() for p in args.output.absolute().parents), 'OUTPUT_ANCESTOR_SYMLINK')
    own_pin = sha(__file__)
    manifest_pin = sha(CAP/'MANIFEST.json')
    plan = read(HERE/'REPLAY_PLAN.json')
    source_capsules = plan['source_capsules']
    need(bool(source_capsules), 'NONEMPTY_INHERITED_CAPSULE_ROSTER')
    old_name = sorted(source_capsules)[0]
    old_root = args.repo_root/'research'/relative(old_name)
    old_manifest_pin = sha(old_root/'MANIFEST.json')
    need(old_manifest_pin == source_capsules[old_name], 'BASELINE_INHERITED_MANIFEST_PIN')
    members = read(old_root/'MANIFEST.json')['files']
    if isinstance(members, dict):
        members = [{'path': key, **value} for key, value in members.items()]
    need(isinstance(members, list) and bool(members), 'NONEMPTY_INHERITED_MEMBER_ROSTER')
    old_member = sorted(row['path'] for row in members)[0]
    relative(old_member)
    results = []
    with tempfile.TemporaryDirectory(prefix='proper-clock-root-controls-') as scratch:
        scratch = Path(scratch)

        def reject(case, optimized, runner, repo, expected):
            mode = 'optimized' if optimized else 'normal'
            out = scratch/(case+'-'+mode+'-output')
            need(not out.exists(), 'FRESH_REJECTION_OUTPUT_DIRECTORY')
            py = [sys.executable, '-I', '-S', '-B']+(['-O'] if optimized else [])
            argv = py+[str(runner), '--repo-root', str(repo), '--output-dir', str(out),
                       '--manifest-sha256', manifest_pin, '--authenticate-only']
            proc = subprocess.run(argv, text=True, capture_output=True)
            observed = parsed_failure(proc)
            no_success_stdout = not successful_stdout(proc)
            need(proc.returncode != 0 and observed == expected and not out.exists()
                 and no_success_stdout and SENTINEL not in proc.stdout+proc.stderr,
                 'ROOT_CONTROL_NOT_REJECTED:'+case)
            results.append({'case': case, 'category': 'source_authentication', 'optimized': optimized,
                            'subprocess_exit': proc.returncode,
                            'expected_witness': expected, 'observed_witness': observed,
                            'success_artifact_absent': not out.exists(),
                            'success_status_absent_from_stdout': no_success_stdout,
                            'sentinel_not_executed': True,
                            'authentication_only': True})

        # Each new-root corruption gets its own full isolated capsule copy.
        for case in ('new_member_bytes', 'new_unlisted_member'):
            copied_cap = scratch/case/CAP.name
            shutil.copytree(CAP, copied_cap)
            if case == 'new_member_bytes':
                target = copied_cap/'numeric/check_proper_clock.py'
                with target.open('a') as handle:
                    handle.write('\nraise RuntimeError("'+SENTINEL+'")\n')
                expected = 'CLOSED_MEMBER:'+CAP.name+'/numeric/check_proper_clock.py'
            else:
                target = copied_cap/'ROOT_CONTROL_UNLISTED.txt'
                need(not target.exists(), 'UNIQUE_UNLISTED_CONTROL_MEMBER')
                target.write_text('Deliberate root inventory corruption.\n')
                expected = 'CLOSED_INVENTORY:'+CAP.name
            for optimized in (False, True):
                reject(case, optimized, copied_cap/'checks/run_checkpoint.py', args.repo_root, expected)

        # Only the first inherited capsule is needed: authentication must reject
        # it before the runner attempts any later capsule in the sorted roster.
        for case in ('old_manifest_bytes', 'old_member_bytes'):
            copied_repo = scratch/case/'repository'
            copied_old = copied_repo/'research'/relative(old_name)
            shutil.copytree(old_root, copied_old)
            if case == 'old_manifest_bytes':
                target = copied_old/'MANIFEST.json'
                target.write_bytes(target.read_bytes()+b' ')
                expected = 'EXTERNAL_MANIFEST_PIN:'+old_name
            else:
                target = copied_old/relative(old_member)
                target.write_bytes(target.read_bytes()+b' ')
                expected = 'CLOSED_MEMBER:'+old_name+'/'+old_member
            for optimized in (False, True):
                reject(case, optimized, HERE/'run_checkpoint.py', copied_repo, expected)
    need(sha(__file__) == own_pin and sha(CAP/'MANIFEST.json') == manifest_pin
         and sha(old_root/'MANIFEST.json') == old_manifest_pin, 'ROOT_CONTROL_INPUTS_CHANGED')
    report = {'schema': 'proper-clock-root-strict-provenance-controls-v1',
              'status': 'PASS_ROOT_PROVENANCE_CORRUPTION_CONTROLS',
              'controls_source_sha256': own_pin,
              'capsule_manifest_sha256': manifest_pin,
              'distinct_controls': 4, 'subprocess_rejections': len(results),
              'both_normal_and_optimized': True, 'results': results,
              'scope': 'Only source provenance is exercised here; authenticate-only performs no scientific imports or checkpoint proof replay.'}
    with args.output.open('x') as handle:
        handle.write(json.dumps(report, sort_keys=True, indent=2)+'\n')
    print(json.dumps({'status': report['status'], 'subprocess_rejections': len(results)}, sort_keys=True))


if __name__ == '__main__':
    main()
