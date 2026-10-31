#!/usr/bin/env python3
"""Real normal/optimized runs and corrupt-code/source/background subprocesses."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', required=True, type=Path)
    ap.add_argument('--background-report', required=True, type=Path)
    ap.add_argument('--background-sha256', required=True)
    ap.add_argument('--output-dir', required=True, type=Path)
    args = ap.parse_args()
    if args.output_dir.exists():
        raise ValueError('FRESH_OUTPUT_DIR')
    args.output_dir.mkdir(parents=True)
    checker = HERE/'check_combined_tt.py'
    own = sha(checker)
    baseline_args = ['--repo-root', str(args.repo_root), '--background-report', str(args.background_report),
                     '--background-sha256', args.background_sha256]
    reports = {}
    for label, opt in [('NORMAL', []), ('OPTIMIZED', ['-O'])]:
        path = args.output_dir/('TT_'+label+'.json')
        p = subprocess.run([sys.executable, '-I', '-S', '-B', *opt, str(checker), *baseline_args,
                            '--output', str(path)], capture_output=True, text=True)
        if p.returncode:
            raise ValueError('BASELINE_FAILED:'+p.stderr[-1800:])
        reports[label] = json.loads(path.read_bytes())
    if reports['NORMAL'] != reports['OPTIMIZED']:
        raise ValueError('NORMAL_OPTIMIZED_CANONICAL_DIFFER')
    controls = []
    mutants = [
        ('wrong_Neumann_combo_sign', 'NEUMANN_COMBO_SIGN = -1', 'NEUMANN_COMBO_SIGN = 1', 'COMBINED_KERNEL_BOOST_IDENTITY'),
        ('wrong_reflected_time_sign', 'REFLECTED_COMBO_SIGN = 1', 'REFLECTED_COMBO_SIGN = -1', 'COMBINED_KERNEL_BOOST_IDENTITY'),
        ('wrong_positive_exp_weight', 'POSITIVE_WEIGHT_IS_EXP_MINUS = True', 'POSITIVE_WEIGHT_IS_EXP_MINUS = False', 'COMBINED_KERNEL_BOOST_IDENTITY'),
        ('wrong_radial_fold_sign', 'RADIAL_REFLECTED_SIGN = -1', 'RADIAL_REFLECTED_SIGN = 1', 'RADIAL_KERNEL_ENDPOINT_SIGN'),
        ('feedback_disabled', 'FEEDBACK_GAIN = 1', 'FEEDBACK_GAIN = 0', 'COMPLETE_FEEDBACK_GAIN'),
        ('extra_half_second_derivative', 'SECOND_DERIVATIVE_NORMALIZATION = 1', 'SECOND_DERIVATIVE_NORMALIZATION = Q(1, 2)', 'NO_EXTRA_HALF_SECOND_DERIVATIVE'),
        ('full_feedback_q_disabled', 'q = T*sigma', 'q = 0*sigma', 'COMPLETE_FEEDBACK_Q'),
    ]
    source = checker.read_text()
    with tempfile.TemporaryDirectory(prefix='tt-controls-', dir=args.output_dir) as tmp:
        tmp = Path(tmp)

        def run(name, script, argv, witness):
            out = tmp/(name+'.json')
            p = subprocess.run([sys.executable, '-I', '-S', '-B', str(script), *argv,
                                '--output', str(out)], capture_output=True, text=True)
            error_line = p.stderr.strip().splitlines()[-1] if p.stderr.strip() else ''
            passed = p.returncode != 0 and not out.exists() and error_line == 'ValueError: '+witness
            controls.append({'name': name, 'exit_code': p.returncode,
                             'rejected': passed, 'failed_witness': witness,
                             'success_artifact_absent': not out.exists()})
            if not passed:
                raise ValueError('CONTROL_NOT_REJECTED:'+name+':'+p.stderr[-1500:])

        for name, original, corrupt, witness in mutants:
            if source.count(original) != 1:
                raise ValueError('UNIQUE_MUTATION_ANCHOR:'+name)
            script = tmp/(name+'.py')
            script.write_text(source.replace(original, corrupt))
            run(name, script, baseline_args, witness)

        # Real byte corruption in an isolated copy; original checkout unchanged.
        cap = 'HDBLAST_POSTCONTACT_SCREEN_CURVATURE_20261031'
        scratch_repo = tmp/'corrupt_repo'
        scratch_cap = scratch_repo/'research'/cap
        shutil.copytree(args.repo_root/'research'/cap, scratch_cap)
        report = scratch_cap/'main/POSTCONTACT_NORMAL.json'
        report.write_bytes(report.read_bytes()+b' ')
        argv = baseline_args.copy()
        argv[1] = str(scratch_repo)
        run('authenticated_source_byte_corruption', checker, argv,
            cap+':MEMBER_BYTES:main/POSTCONTACT_NORMAL.json')
        # A report byte mutation must fail before any background value is used.
        bg = tmp/'corrupt_background.json'
        bg.write_bytes(args.background_report.read_bytes()+b' ')
        argv = baseline_args.copy()
        argv[3] = str(bg)
        run('background_report_byte_corruption', checker, argv, 'BACKGROUND_REPORT_PIN')
    if sha(checker) != own:
        raise ValueError('CHECKER_CHANGED')
    receipt = {'schema': 'HDBLAST_PROPER_CLOCK_TT_SUBPROCESS_CONTROLS_V1',
               'status': 'PASS', 'checker_sha256': own,
               'background_report_sha256': args.background_sha256,
               'normal_optimized_canonical_equal': True,
               'baseline_reports': {k: sha(args.output_dir/('TT_'+k+'.json')) for k in reports},
               'controls': controls, 'rejected_subprocess_count': len(controls),
               'limits': 'These controls verify formula wiring and authentication; they do not reprove inherited physical source boxes.'}
    (args.output_dir/'CONTROLS.json').write_text(json.dumps(receipt, sort_keys=True, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', 'rejected_subprocess_count': len(controls),
                      'canonical_equal': True}))


if __name__ == '__main__':
    main()
