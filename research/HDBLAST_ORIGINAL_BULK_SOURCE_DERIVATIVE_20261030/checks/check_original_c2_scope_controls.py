#!/usr/bin/env python3
"""Reject wrong source scope/coordinate/status at the original-target gate."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


ap = argparse.ArgumentParser()
ap.add_argument('--repo-root', type=Path, required=True)
ap.add_argument('--fx-report', type=Path, required=True)
ap.add_argument('--gate', type=Path, default=Path(__file__).with_name('check_original_target_c2.py'))
ap.add_argument('--higher-capsule', type=Path)
ap.add_argument('--causal-capsule', type=Path)
ap.add_argument('--refined-response', type=Path)
a = ap.parse_args()
raw = a.fx_report.read_bytes()
actual = json.loads(raw)
need(actual['status'].startswith('PASS') and actual['scope']['actual_original_full23_F_x_bound'] is True and
     actual['source_coordinate'] == 'x=t/ell', 'baseline is complete original scaled-time source report')
options = ['--repo-root', str(a.repo_root.resolve())]
for option, value in (('--higher-capsule', a.higher_capsule), ('--causal-capsule', a.causal_capsule),
                      ('--refined-response', a.refined_response)):
    if value is not None:
        options += [option, str(value.resolve())]
runs = []
with tempfile.TemporaryDirectory(prefix='hdblast-c2-scope-') as temp:
    for injection in ('filtered_source_scope', 'physical_t_instead_of_scaled_x', 'negative_cap', 'nonpassing_status'):
        changed = copy.deepcopy(actual)
        expected = 'complete original scaled-time source derivative scope'
        if injection == 'filtered_source_scope':
            changed['scope']['actual_original_full23_F_x_bound'] = False
        elif injection == 'physical_t_instead_of_scaled_x':
            changed['source_coordinate'] = 't'
        elif injection == 'negative_cap':
            changed['genuine_scaled_F_x_H_cap'] = '-1'
            expected = 'nonnegative actual raw F_x H cap'
        else:
            changed['status'] = 'NOT_PASS'
        path = Path(temp)/(injection+'.json')
        path.write_text(json.dumps(changed, sort_keys=True))
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        for mode in ('normal', 'optimized'):
            command = [sys.executable, '-I', '-S', '-B']+(['-O'] if mode == 'optimized' else [])
            command += [str(a.gate.resolve())]+options+['--fx-report', str(path), '--fx-sha256', digest]
            result = subprocess.run(command, capture_output=True, timeout=60)
            need(result.returncode != 0 and expected in result.stderr.decode(errors='replace'),
                 'wrong source input not rejected: '+injection+'/'+mode)
            runs.append({'mode': mode, 'injection': injection, 'rejected': True})
report = {'status': 'PASS_ORIGINAL_C2_INPUT_SCOPE_CONTROLS',
          'actual_source_report_sha256': hashlib.sha256(raw).hexdigest(),
          'gate_checker_sha256': hashlib.sha256(a.gate.read_bytes()).hexdigest(),
          'rejected_run_count': len(runs), 'rejected_runs': runs,
          'scope': 'Input-role and unit guards only; mathematical bounds require the genuine complete-source derivative proof.'}
print(json.dumps(report, indent=2, sort_keys=True))
