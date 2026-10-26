#!/usr/bin/env python3
"""Run the independent new full23 kernel review and four precise mutants."""
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
    parser.add_argument('--work-root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    need(not args.output_dir.exists(), 'FRESH_OUTPUT_DIRECTORY_REQUIRED')
    args.output_dir.mkdir(parents=True)
    script = Path(__file__).resolve().with_name('check_angular_ridge.py')
    code_pin = sha(script)
    reasons = {'wrong-reflection': 'POLYNOMIAL_DIAMOND:value',
               'missing-one-D-row': 'MISSING_ONE_DIRICHLET_ROW_GUARD:v:4',
               'commute-matrices': 'COUPLED_PHYSICAL_POLYGON:value',
               'omit-PN-projection': 'COUPLED_PHYSICAL_POLYGON:value'}
    runs, controls, reports = [], [], []
    for optimized in (False, True):
        mode = 'optimized' if optimized else 'normal'
        command = [sys.executable, '-I', '-S', '-B']+(['-O'] if optimized else [])
        command += [str(script), '--repo-root', str(args.repo_root.resolve()),
                    '--work-root', str(args.work_root.resolve())]
        report = args.output_dir/('ANGULAR_RIDGE_'+mode.upper()+'.json')
        process = subprocess.run(command+['--output', str(report)], capture_output=True, text=True)
        (args.output_dir/(mode+'.log')).write_text(process.stdout+process.stderr)
        need(process.returncode == 0 and report.is_file(), 'POSITIVE_RUN:'+mode)
        reports.append(report)
        runs.append({'optimized': optimized, 'report': report.name, 'sha256': sha(report)})
        for defect, reason in reasons.items():
            output = args.output_dir/('REJECT_'+mode+'_'+defect+'.json')
            failed = subprocess.run(command+['--output', str(output), '--inject-error', defect],
                                    capture_output=True, text=True)
            (args.output_dir/('reject-'+mode+'-'+defect+'.log')).write_text(failed.stdout+failed.stderr)
            need(failed.returncode != 0 and ('ValueError: '+reason) in failed.stderr
                 and not output.exists(), 'REASON_SPECIFIC_CONTROL:'+mode+':'+defect)
            controls.append({'optimized': optimized, 'defect': defect, 'reason': reason,
                             'rejected': True, 'no_report_written': True})
    need(reports[0].read_bytes() == reports[1].read_bytes(), 'NORMAL_OPTIMIZED_BYTES')
    need(sha(script) == code_pin, 'REVIEWER_CHANGED_DURING_EXECUTION')
    receipt = {'schema': 'HDBLAST_FULL23_KERNEL_INDEPENDENT_EXECUTION_V1', 'status': 'PASS',
               'reviewer_sha256': code_pin, 'normal_optimized_byte_identical': True,
               'positive_runs': runs, 'negative_controls': controls,
               'scope': 'New full23 independent exact controls only; no old campaign replay.'}
    (args.output_dir/'EXECUTION_RECEIPT.json').write_text(json.dumps(receipt, sort_keys=True, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', 'positive_runs': len(runs), 'mutant_controls': len(controls),
                      'report_sha256': sha(reports[0])}, sort_keys=True))


if __name__ == '__main__':
    main()
