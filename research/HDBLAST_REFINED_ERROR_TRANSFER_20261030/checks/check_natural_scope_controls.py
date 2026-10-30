#!/usr/bin/env python3
"""Execute corrupted inputs and formula variants through the natural checker."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def need(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--c2-report', type=Path, required=True)
    ap.add_argument('--c2-sha256', required=True)
    ap.add_argument('--checker', type=Path, required=True)
    ap.add_argument('--checker-sha256', required=True)
    a = ap.parse_args()
    original_raw = a.c2_report.read_bytes()
    checker_raw = a.checker.read_bytes()
    need(hashlib.sha256(original_raw).hexdigest() == a.c2_sha256, 'authenticate original C2 input')
    need(hashlib.sha256(checker_raw).hexdigest() == a.checker_sha256, 'authenticate natural checker')
    original = json.loads(original_raw)
    flags = [sys.executable, '-I', '-S', '-B']+(['-O'] if sys.flags.optimize else [])
    premise_flags = ['--assume-first-harmonic-neighbourhood', '--assume-formal-subsidiary-completion']
    results = []

    def reject(name, report, digest, expected, extra=(), premises=premise_flags):
        result = subprocess.run(flags+[str(a.checker), '--repo-root', str(a.repo_root),
            '--c2-report', str(report), '--c2-sha256', digest]+list(premises)+list(extra),
            capture_output=True, timeout=60)
        need(result.returncode != 0 and not result.stdout and expected in result.stderr.decode(errors='replace'),
             'corrupted natural input not rejected: '+name)
        results.append({'name': name, 'rejected': True, 'guard': expected})

    reject('missing first harmonic premise', a.c2_report, a.c2_sha256,
        'additional retained first-harmonic neighbourhood premise must be explicit',
        premises=['--assume-formal-subsidiary-completion'])
    reject('missing formal target subsidiary premise', a.c2_report, a.c2_sha256,
        'additional retained formal subsidiary completion premise must be explicit',
        premises=['--assume-first-harmonic-neighbourhood'])
    variants = [
        ('NOT_PASS input', lambda d: d.update(status='NOT_PASS'), 'complete original-target absolute C2 error certificate'),
        ('filtered candidate name', lambda d: d.update(candidate='Z_epsilon'), 'unchanged original-source W1_quadratic candidate'),
        ('changed candidate', lambda d: d.update(candidate_changed=True), 'same candidate target and actual source'),
        ('changed source', lambda d: d.update(actual_sources_changed=True), 'same candidate target and actual source'),
        ('changed target', lambda d: d.update(actual_target_unchanged=False), 'same candidate target and actual source'),
        ('wrong domain', lambda d: d.update(domain='outside physical collar'), 'same original target domain and coordinate norm'),
        ('wrong Gram', lambda d: d['Gram'].__setitem__(0, '2'), 'same complete full23 H norm'),
        ('missing second derivative', lambda d: d['original_target_second_coordinate_error_H_caps'].pop('rr'),
         'all three nonnegative absolute second-error caps'),
        ('negative second derivative cap', lambda d: d['original_target_second_coordinate_error_H_caps'].__setitem__('rr', '-1'),
         'all three nonnegative absolute second-error caps'),
        ('zero reduced defect', lambda d: d.update(complete_original_reduced_defect_H_cap='0'),
         'retained C1 and original reduced-defect caps'),
        ('changed C1 cap', lambda d: d.update(complete_original_C1_error_X_cap='1'),
         'retained C1 and original reduced-defect caps'),
    ]
    with tempfile.TemporaryDirectory(prefix='hdblast-natural-transfer-controls-') as tmp:
        for name, mutate, expected in variants:
            data = copy.deepcopy(original)
            mutate(data)
            raw = (json.dumps(data, sort_keys=True)+'\n').encode()
            path = Path(tmp)/('case_'+str(len(results))+'.json')
            path.write_bytes(raw)
            reject(name, path, hashlib.sha256(raw).hexdigest(), expected)
    for fault, expected in (
        ('symmetrization', 'exact shared-source one-half symmetric gradient'),
        ('metric-normalization', 'exact original metric forcing normalization'),
        ('radial-product-sign', 'retain radial physical product-rule sign'),
        ('filtered-source-charge', 'original W1q reduced defect without filtered-source charge'),
        ('value-bound', 'value bound uses retained X definition'),
    ):
        reject('formula '+fault, a.c2_report, a.c2_sha256, expected, ['--inject-error', fault])
    need(a.c2_report.read_bytes() == original_raw and a.checker.read_bytes() == checker_raw,
         'retained natural-transfer inputs unchanged on readback')
    report = {'status': 'PASS_NATURAL_W1Q_ERROR_TRANSFER_SCOPE_CONTROLS',
        'input_sha256': {'new_C2': a.c2_sha256, 'checker': a.checker_sha256},
        'corrupted_runs_rejected': results, 'count': len(results),
        'scope': 'Explicit premise, candidate/source/target and formula guard executions; flags record conditional premises and do not prove physical constraints'}
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
