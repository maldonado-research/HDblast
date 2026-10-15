#!/usr/bin/env python3
"""Fresh-process, mutation-specific failures of the direct 5D gauge check."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--checker', type=Path,
                    default=Path(__file__).resolve().parents[1]/'gauge/check_initial_harmonic_completion.py')
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    mutations = ['omit_lapse_quadratic', 'halve_radial_U_Uy', 'omit_radial_fy',
                 'halve_angular_quadratic', 'omit_mixed_metric_cross',
                 'omit_second_metric_lapse_value']
    cases = []
    with tempfile.TemporaryDirectory() as tmp:
        for mode in ['normal', 'optimized']:
            for name in mutations:
                target = Path(tmp)/(mode+'-'+name+'.json')
                cmd = [sys.executable, '-I', '-S', '-B']
                if mode == 'optimized':
                    cmd += ['-O']
                cmd += [str(args.checker), '--mutation', name, '--output', str(target)]
                r = subprocess.run(cmd, text=True, capture_output=True, timeout=60)
                reason = 'FAIL: nonzero initial relative-harmonic residual under '+name
                if r.returncode == 0 or reason not in r.stderr or target.exists():
                    raise ValueError('MUTATION_SPECIFIC_CLI_REJECTION_FAILED: '+mode+' '+name)
                cases.append({'mode': mode, 'mutation': name, 'returncode': r.returncode,
                              'expected_reason': reason, 'no_output': True})
    report = {'status': 'PASS_TWELVE_FRESH_GAUGE_FAILURE_CONTROLS', 'cases': cases,
              'checker_sha256': hashlib.sha256(args.checker.read_bytes()).hexdigest()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, sort_keys=True, indent=2)+'\n')


if __name__ == '__main__':
    main()
