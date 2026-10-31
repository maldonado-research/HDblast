#!/usr/bin/env python3
"""Strict subprocess controls for the proper-clock background certificate."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(argv):
    return subprocess.run(argv, text=True, capture_output=True)


def parsed_failure(proc):
    lines = proc.stderr.strip().splitlines()
    line = lines[-1] if lines else ''
    prefix = 'ValueError: '
    return line[len(prefix):] if line.startswith(prefix) else None


def successful_stdout(proc):
    for line in proc.stdout.splitlines():
        try:
            value = json.loads(line)
        except (ValueError, TypeError):
            continue
        if isinstance(value, dict):
            status = value.get('status', '')
            if isinstance(status, str) and status.startswith(('PASS', 'CERTIFIED')):
                return True
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    need(not args.output.exists() and not args.output.is_symlink(), 'FRESH_CONTROL_OUTPUT')
    need(all(not p.is_symlink() for p in args.output.absolute().parents), 'OUTPUT_ANCESTOR_SYMLINK')
    own_pin = sha(__file__)
    pins_path = HERE/'SOURCE_PINS.json'
    pins_pin = sha(pins_path)
    pins = json.loads(pins_path.read_text())['files']
    baseline = HERE/'BACKGROUND_NORMAL.json'
    baseline_pin = sha(baseline)
    need(baseline.read_bytes() == (HERE/'BACKGROUND_OPTIMIZED.json').read_bytes(),
         'NORMAL_OPTIMIZED_BACKGROUND_RECEIPTS_IDENTICAL')
    rejections, positives = [], []
    with tempfile.TemporaryDirectory(prefix='proper-clock-background-controls-') as scratch:
        scratch = Path(scratch)
        copyroot = scratch/'source-copy'
        for rel in pins:
            relative = Path(rel)
            need(not relative.is_absolute() and '..' not in relative.parts, 'SAFE_SOURCE_COPY_PATH')
            path = copyroot/relative
            path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(args.repo_root/relative, path)
        source_relative = Path('research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/interval_core.py')
        corrupted_source = copyroot/source_relative
        corrupted_source.write_text("raise RuntimeError('IMPORT_SENTINEL_MUST_NOT_RUN')\n"+
                                    corrupted_source.read_text())
        changed_receipt = json.loads(baseline.read_text())
        changed_receipt['compact_tube']['x_plus_r_lt_8_margin'] = '-1'
        changed_receipt_path = scratch/'changed-margin-receipt.json'
        changed_receipt_path.write_text(json.dumps(changed_receipt, sort_keys=True, indent=2)+'\n')

        def reject(case, category, optimized, argv, expected, requested_output=None):
            proc = run(argv)
            observed = parsed_failure(proc)
            no_success_stdout = not successful_stdout(proc)
            output_absent = requested_output is None or not requested_output.exists()
            need(proc.returncode != 0 and observed == expected and no_success_stdout and output_absent,
                 'CONTROL_NOT_REJECTED:'+case)
            if category == 'source_authentication':
                need('IMPORT_SENTINEL_MUST_NOT_RUN' not in proc.stderr+proc.stdout,
                     'CORRUPTED_SOURCE_EXECUTED')
            rejections.append({'case': case, 'category': category, 'optimized': optimized,
                               'subprocess_exit': proc.returncode,
                               'expected_witness': expected, 'observed_witness': observed,
                               'success_status_absent_from_stdout': no_success_stdout,
                               'success_artifact_absent': output_absent,
                               'success_artifact_kind': 'requested_output_report' if requested_output else 'stdout_success_receipt',
                               'output_report_requested': requested_output is not None})

        for optimized in (False, True):
            py = [sys.executable, '-I', '-B']+(['-O'] if optimized else [])
            mode = 'optimized' if optimized else 'normal'
            for case, category, expected in (
                    ('wrong_acceleration_sign', 'formula', 'inward launch must have negative rapidity; wrong-sign control'),
                    ('wrong_clock_unit', 'units', 'physical proper clock unit dτ/dt=ρ sechχ'),
                    ('corrupt_domain_margin', 'domain', 'declared chart tube must contain all connecting trajectories')):
                target = scratch/(case+'-'+mode+'.json')
                reject(case, category, optimized,
                       py+[str(HERE/'background_certificate.py'), '--repo-root', str(args.repo_root),
                           '--inject-error', case, '--output', str(target)], expected, target)
            target = scratch/('changed-source-'+mode+'.json')
            reject('on_disk_imported_source_mutation', 'source_authentication', optimized,
                   py+[str(HERE/'background_certificate.py'), '--repo-root', str(copyroot), '--output', str(target)],
                   'repository authentication before import: '+str(source_relative), target)
            reject('on_disk_margin_receipt_mutation', 'receipt_integrity', optimized,
                   py+[str(HERE/'check_background_certificate.py'), '--repo-root', str(args.repo_root),
                       '--receipt', str(changed_receipt_path)],
                   'background receipt differs from exact recomputation')

        # Positive relocation checks are reported separately from rejections.
        shutil.copyfile(args.repo_root/source_relative, copyroot/source_relative)
        for optimized in (False, True):
            py = [sys.executable, '-I', '-B']+(['-O'] if optimized else [])
            proc = run(py+[str(HERE/'check_background_certificate.py'), '--repo-root', str(copyroot)])
            lines = [line for line in proc.stdout.splitlines() if line.strip()]
            value = json.loads(lines[-1]) if lines else {}
            need(proc.returncode == 0 and value.get('status') == 'PASS_BACKGROUND_RECEIPT_REPLAY'
                 and value.get('sha256') == baseline_pin, 'RELOCATED_SOURCE_CLOSURE_REPLAY')
            positives.append({'case': 'relocated_authenticated_source_closure', 'optimized': optimized,
                              'subprocess_exit': proc.returncode, 'status': value['status'],
                              'canonical_background_sha256': value['sha256']})
    need(sha(__file__) == own_pin and sha(pins_path) == pins_pin and sha(baseline) == baseline_pin,
         'CONTROL_INPUTS_CHANGED')
    report = {'schema': 'proper-clock-background-strict-controls-v1',
              'status': 'PASS_BACKGROUND_PROPER_CLOCK_CONTROLS',
              'checker_sha256': own_pin, 'source_pins_sha256': pins_pin,
              'canonical_background_sha256': baseline_pin,
              'normal_optimized_receipts_identical': True,
              'distinct_rejection_controls': 5, 'subprocess_rejections': len(rejections),
              'positive_subprocess_checks': len(positives),
              'rejections': rejections, 'positive_checks': positives,
              'scope': 'Five distinct rejection controls run in both modes; two successful relocation checks are counted separately.'}
    with args.output.open('x') as handle:
        handle.write(json.dumps(report, sort_keys=True, indent=2)+'\n')
    print(json.dumps({'status': report['status'], 'subprocess_rejections': len(rejections),
                      'positive_subprocess_checks': len(positives)}, sort_keys=True))


if __name__ == '__main__':
    main()
