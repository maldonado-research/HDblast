#!/usr/bin/env python3
"""Run corrupted refinement reports through the actual C2 checker."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def need(ok, label):
    if not ok:
        raise RuntimeError(label)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--initial-refinement', type=Path, required=True)
    ap.add_argument('--initial-refinement-sha256', required=True)
    ap.add_argument('--checker', type=Path, required=True)
    ap.add_argument('--checker-sha256', required=True)
    a = ap.parse_args()
    raw = a.initial_refinement.read_bytes()
    need(hashlib.sha256(raw).hexdigest() == a.initial_refinement_sha256,
         'authenticate original refinement report')
    code = a.checker.read_bytes()
    need(hashlib.sha256(code).hexdigest() == a.checker_sha256,
         'authenticate guarded C2 checker')
    original = json.loads(raw)
    variants = [
        ('NOT_PASS_status', lambda d: d.update(status='NOT_PASS'), 'successful refined initial report'),
        ('changed_candidate', lambda d: d.update(candidate_changed=True), 'refinement changes estimates only, not candidate'),
        ('other_candidate', lambda d: d.update(same_candidate_W1_quadratic=False), 'refinement changes estimates only, not candidate'),
        ('lower_bound_role', lambda d: d.update(refined_u_rr_H_cap_is_upper_bound=False), 'refined input is an upper bound'),
        ('point_value_role', lambda d: d.update(bound_quantity='point_signed_actual_minus_smoothed_initial_u_rr'), 'refined input is an upper bound'),
        ('missing_component', lambda d: d.update(labels=d['labels'][:-1]), 'refined input has the same complete positive H norm'),
        ('altered_Gram', lambda d: d['Gram'].__setitem__(0, '2'), 'refined input has the same complete positive H norm'),
        ('negative_upper', lambda d: d.update(refined_uniform_initial_u_rr_H_upper='-1'), 'refined initial upper bound is nonnegative'),
    ]
    flags = [sys.executable, '-I', '-S', '-B']+(['-O'] if sys.flags.optimize else [])
    results = []
    with tempfile.TemporaryDirectory(prefix='hdblast-refined-c2-roles-') as tmp:
        for name, mutate, expected in variants:
            data = copy.deepcopy(original)
            mutate(data)
            corrupted = Path(tmp)/(name+'.json')
            changed = (json.dumps(data, sort_keys=True)+'\n').encode()
            corrupted.write_bytes(changed)
            sha = hashlib.sha256(changed).hexdigest()
            result = subprocess.run(flags+[str(a.checker), '--repo-root', str(a.repo_root),
                '--initial-refinement', str(corrupted), '--initial-refinement-sha256', sha],
                capture_output=True, timeout=30)
            need(result.returncode != 0 and not result.stdout and expected in result.stderr.decode(),
                 'corrupted report not rejected: '+name)
            results.append({'name': name, 'rejected': True, 'guard': expected})
    need(a.initial_refinement.read_bytes() == raw and a.checker.read_bytes() == code,
         'retained inputs unchanged on readback')
    report = {'status': 'PASS_REFINED_C2_SCOPE_CONTROLS',
        'input_sha256': {'refinement': a.initial_refinement_sha256, 'checker': a.checker_sha256},
        'corrupted_report_runs_rejected': results,
        'scope': 'Eight real corrupted input-file executions through authenticated checker; no mathematical error or accuracy claim'}
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
