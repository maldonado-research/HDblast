#!/usr/bin/env python3
"""Check normal/optimized equivalence and reason-specific proof mutations."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from check_evolved_scalar_signs import MUTATIONS


def need(condition, reason):
    if not condition:
        raise ValueError(reason)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    need(not args.output.exists(), 'refuse_existing_output')
    here = Path(__file__).resolve().parent
    checker = here / 'check_evolved_scalar_signs.py'
    with tempfile.TemporaryDirectory(prefix='evolved-sign-controls-') as temp:
        temp = Path(temp)
        base_reports, controls = {}, []
        for mode, flags in [('normal', []), ('optimized', ['-O'])]:
            output = temp / (mode + '.json')
            base = subprocess.run([sys.executable, *flags, str(checker), '--output', str(output)],
                                  capture_output=True, text=True)
            need(base.returncode == 0 and output.exists(), 'base_' + mode)
            base_reports[mode] = output.read_bytes()
            for mutation, expected_reason in MUTATIONS.items():
                output = temp / (mode + '-' + mutation + '.json')
                result = subprocess.run([sys.executable, *flags, str(checker), '--output', str(output),
                                         '--inject', mutation], capture_output=True, text=True)
                need(result.returncode != 0 and not output.exists(), 'reject_' + mode + '_' + mutation)
                tail = result.stderr.strip().splitlines()[-1]
                need(tail == 'ValueError: ' + expected_reason, 'reason_' + mode + '_' + mutation)
                controls.append({'mode': mode, 'mutation': mutation, 'reason': expected_reason,
                                 'rejected': True, 'output_absent': True})
        need(base_reports['normal'] == base_reports['optimized'], 'normal_optimized_byte_equivalence')
    report = {'status': 'PASS_EVOLVED_SIGN_REASON_CONTROLS',
              'normal_optimized_byte_identical': True,
              'base_report_sha256': hashlib.sha256(base_reports['normal']).hexdigest(),
              'checker_sha256': hashlib.sha256(checker.read_bytes()).hexdigest(),
              'rejections': controls}
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(report['status'])


if __name__ == '__main__':
    main()
