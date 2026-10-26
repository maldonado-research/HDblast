#!/usr/bin/env python3
"""Verify evaluator I/O rejection before authentication or source evaluation."""
import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

CAPSULES = (
    'HDBLAST_COUPLED_MEAN_PILOT_20261025',
    'HDBLAST_INPUT_PRECISION_20261023',
    'HDBLAST_COMMON_GAUGE_FORCING_20261016',
    'HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016',
    'HDBLAST_GEOMETRIC_BULK_SOURCES_20261016',
    'HDBLAST_FORCED_RESPONSE_BOUNDS_20261023',
    'HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023',
    'HDBLAST_SHELL_FEEDBACK_PRECISION_20261025',
)


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rejected(evaluator, flags, repo, output, message):
    result = subprocess.run([sys.executable, '-B', *flags, str(evaluator),
                             '--repo-root', str(repo), '--output', str(output)],
                            capture_output=True, text=True, timeout=10)
    need(result.returncode != 0 and message in result.stderr, 'GUARD_NOT_REJECTED')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', required=True, type=Path)
    parser.add_argument('--evaluator', required=True, type=Path)
    parser.add_argument('--report', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    repo, target = args.repo_root.resolve(), args.output.resolve()
    protected = [(repo/'research'/name).resolve() for name in CAPSULES]
    need(not args.output.exists() and not args.output.is_symlink(), 'FRESH_OUTPUT')
    need(all(target != root and root not in target.parents for root in protected),
         'OUTPUT_INSIDE_PROTECTED_CAPSULE')
    old_hash, evaluator_hash = sha(args.report), sha(args.evaluator)
    controls = []
    for flags, mode in (([], 'normal'), (['-O'], 'optimized')):
        for capsule in CAPSULES:
            output = repo/'research'/capsule/'_TT_IO_GUARD_MUST_NOT_WRITE_.json'
            need(not output.exists(), 'GUARD_PROBE_COLLISION')
            rejected(args.evaluator, flags, repo, output, 'OUTPUT_INSIDE_FROZEN_CAPSULE')
            need(not output.exists(), 'PROTECTED_OUTPUT_WRITTEN')
            controls.append(mode + ':' + capsule)
        rejected(args.evaluator, flags, repo, args.report, 'FRESH_OUTPUT_REQUIRED')
        controls.append(mode + ':existing_output')
        with tempfile.TemporaryDirectory(prefix='tt-guard-') as temporary:
            fake_repo = Path(temporary)/'repo'
            (fake_repo/'research').mkdir(parents=True)
            linked = fake_repo/'research'/CAPSULES[-1]
            linked.symlink_to(protected[-1], target_is_directory=True)
            output = linked/'_TT_IO_GUARD_MUST_NOT_WRITE_.json'
            rejected(args.evaluator, flags, fake_repo, output, 'OUTPUT_INSIDE_FROZEN_CAPSULE')
            need(not output.exists(), 'SYMLINKED_PROTECTED_OUTPUT_WRITTEN')
        controls.append(mode + ':symlinked_helper_capsule')
    need(old_hash == sha(args.report) and evaluator_hash == sha(args.evaluator), 'INPUT_CHANGED')
    out = {'schema': 'HDBLAST_TT_EVALUATOR_IO_GUARD_REVIEW_V1',
           'status': 'PASS_ALL_EVALUATOR_IO_GUARDS',
           'rejected_cases': controls, 'rejected_case_count': len(controls),
           'source_oracle_executed': False,
           'files_sha256': {'tensor/tt_response.py': evaluator_hash,
                            'tensor/TT_TOTAL_RESPONSE_NORMAL.json': old_hash,
                            'tensor/review/check_tt_io_guards.py': sha(Path(__file__))}}
    with args.output.open('x') as handle:
        json.dump(out, handle, indent=2, sort_keys=True)
        handle.write('\n')
    print(json.dumps({'status': out['status'], 'rejected_cases': len(controls)}))


if __name__ == '__main__':
    main()
