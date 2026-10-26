#!/usr/bin/env python3
"""Check byte-identical normal and optimized full23 response outputs."""
import argparse
import importlib.util
import json
import sys
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
COMPUTE = HERE.parent
DEFAULT_REPO = Path('/workspace/hdblast-research-work/shell-observable-20261025/repo')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, default=DEFAULT_REPO)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location('full23_response_mode_review', COMPUTE/'full23_response.py')
    producer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(producer)
    response = producer.Response(args.repo_root)
    files = {}
    for name in ('RESPONSE_REPORT.json', 'RESPONSE_PROFILES.json', 'EVALUATED_CANDIDATE.json'):
        normal, optimized = COMPUTE/'response-normal'/name, COMPUTE/'response-optimized'/name
        producer.need(normal.read_bytes() == optimized.read_bytes(), 'NORMAL_OPTIMIZED_BYTES:'+name)
        files[name] = {'sha256': producer.sha(normal), 'bytes': normal.stat().st_size}
    report = json.loads((COMPUTE/'response-normal/RESPONSE_REPORT.json').read_text())
    producer.need(report['passes_target'] and F(report['ledger']['total_X_error_upper']) < F(report['target_X_error']),
                  'EXACT_TARGET_GATE')
    producer.need(report['regression']['exact_mean_feedback_coefficients_equal'] == 90,
                  'NINETY_EXACT_MEAN_FEEDBACK_REGRESSIONS')
    response.context.reauthenticate()
    producer.need(producer.authenticate_local() == response.local_before, 'LOCAL_DEPENDENCY_CHANGE')
    result = {
        'status': 'PASS_BYTE_IDENTICAL_NORMAL_OPTIMIZED_FULL23_RESPONSE_OUTPUTS',
        'producer_sha256': producer.sha(COMPUTE/'full23_response.py'),
        'files': files,
        'exact_target_gate': True,
        'complete_roots_and_local_dependencies_authenticated_before_and_after': True,
        'scaled_derivative_scope': 'The pointwise derivative error in the native report refers to coordinate scaled x/r derivatives. Angular covariant derivatives require separate basis and frame factors.',
    }
    (HERE/'MODE_PARITY_REVIEW.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
