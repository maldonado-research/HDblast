#!/usr/bin/env python3
"""Check byte-identical normal/optimized full23 initial production outputs."""
import argparse
import importlib.util
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_REPO = Path('/workspace/hdblast-research-work/shell-observable-20261025/repo')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, default=DEFAULT_REPO)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location('mode_parity_initial23', ROOT/'full23_initial.py')
    producer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(producer)
    context = producer.Context(args.repo_root)
    files = {}
    for name in ('INITIAL_REPORT.json', 'INITIAL_PROFILES.json'):
        normal, optimized = ROOT/'results'/name, ROOT/'optimized-results'/name
        producer.need(normal.read_bytes() == optimized.read_bytes(), 'NORMAL_OPTIMIZED_BYTES:'+name)
        files[name] = {'sha256': producer.sha(normal), 'bytes': normal.stat().st_size}
    exact_review = json.loads((HERE/'EXACT_ADAPTER_REVIEW.json').read_text())
    output_review = json.loads((HERE/'OUTPUT_REVIEW.json').read_text())
    producer.need(exact_review['status'] == 'PASS_EXACT_RATIONAL_JET_FULL23_INITIAL_ADAPTER_REVIEW', 'EXACT_REVIEW')
    producer.need(exact_review['producer_sha256'] == producer.sha(ROOT/'full23_initial.py'), 'EXACT_REVIEW_SOURCE_CHANGED')
    producer.need(output_review['report_sha256'] == files['INITIAL_REPORT.json']['sha256'], 'REPORT_CHANGED_AFTER_REVIEW')
    producer.need(output_review['profiles_sha256'] == files['INITIAL_PROFILES.json']['sha256'], 'PROFILES_CHANGED_AFTER_REVIEW')
    context.reauthenticate()
    result = {
        'status': 'PASS_BYTE_IDENTICAL_NORMAL_OPTIMIZED_FULL23_INITIAL_OUTPUTS',
        'files': files,
        'producer_sha256': producer.sha(ROOT/'full23_initial.py'),
        'all_scientific_decisions_use_explicit_need_not_assert': True,
        'independent_exact_adapter_review_sha256': producer.sha(HERE/'EXACT_ADAPTER_REVIEW.json'),
        'output_support_and_mean_review_sha256': producer.sha(HERE/'OUTPUT_REVIEW.json'),
        'complete_roots_authenticated_before_and_after': True,
    }
    (HERE/'MODE_PARITY_REVIEW.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
