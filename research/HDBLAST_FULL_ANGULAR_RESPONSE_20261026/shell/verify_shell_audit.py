#!/usr/bin/env python3
"""Run only the new bounded shell checks and reason-specific negative controls."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def need(ok, label):
    if not ok:
        raise ValueError(label)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    need(not args.output_dir.exists(), 'FRESH_OUTPUT_DIRECTORY_REQUIRED')
    args.output_dir.mkdir(parents=True)
    producer = Path(__file__).resolve().with_name('check_full23_shell.py')
    before = sha(producer)
    runs, controls = [], []
    reasons = {'boundary-sign': 'FULL_BOUNDARY_SOURCE:w:2',
               'missing-TT': 'FULL_BOUNDARY_SOURCE:T:-',
               'normal-trace': 'FULL_BOUNDARY_SOURCE:K:0',
               'wrong-Dirichlet': 'EXACT_DIRICHLET_SLOTS'}
    reports = []
    for optimized in (False, True):
        mode = 'optimized' if optimized else 'normal'
        command = [sys.executable, '-I', '-S', '-B']+(['-O'] if optimized else [])
        command += [str(producer), '--repo-root', str(args.repo_root.resolve())]
        report = args.output_dir/('FULL23_SHELL_'+mode.upper()+'.json')
        process = subprocess.run(command+['--output', str(report)], capture_output=True, text=True)
        (args.output_dir/(mode+'.log')).write_text(process.stdout+process.stderr)
        need(process.returncode == 0 and report.is_file(), 'POSITIVE_RUN:'+mode)
        reports.append(report)
        runs.append({'optimized': optimized, 'report': report.name, 'sha256': sha(report)})
        for defect, reason in reasons.items():
            failed_output = args.output_dir/('REJECT_'+mode+'_'+defect+'.json')
            failed = subprocess.run(command+['--output', str(failed_output), '--inject-error', defect],
                                    capture_output=True, text=True)
            need(failed.returncode != 0 and ('ValueError: '+reason) in failed.stderr
                 and not failed_output.exists(), 'REASON_SPECIFIC_CONTROL:'+mode+':'+defect)
            (args.output_dir/('reject-'+mode+'-'+defect+'.log')).write_text(failed.stdout+failed.stderr)
            controls.append({'optimized': optimized, 'defect': defect, 'reason': reason,
                             'rejected': True, 'no_report_written': True})
    need(reports[0].read_bytes() == reports[1].read_bytes(), 'NORMAL_OPTIMIZED_BYTES')
    need(sha(producer) == before, 'NEW_PRODUCER_CHANGED_DURING_RUN')
    receipt = {'schema': 'HDBLAST_FULL23_SHELL_NEW_CHECK_EXECUTION_V1', 'status': 'PASS',
               'producer_sha256': before, 'normal_optimized_byte_identical': True,
               'positive_runs': runs, 'negative_controls': controls,
               'scope': 'Only new bounded universal polynomial checks were run; no old campaign replay.'}
    (args.output_dir/'EXECUTION_RECEIPT.json').write_text(json.dumps(receipt, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status': 'PASS', 'positive_runs': len(runs), 'negative_controls': len(controls),
                      'report_sha256': sha(reports[0])}, sort_keys=True))


if __name__ == '__main__':
    main()
