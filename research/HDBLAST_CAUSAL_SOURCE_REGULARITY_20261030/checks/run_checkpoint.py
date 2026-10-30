#!/usr/bin/env python3
"""Reproduce this checkpoint in normal/optimized isolated standard-library Python."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def need(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    root = args.repo_root.resolve()
    cap = Path(__file__).resolve().parent.parent
    out = args.output_dir.resolve()
    need(not out.exists(), 'Use a fresh output directory')
    out.mkdir(parents=True)
    hi = root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030/reference'
    geom = root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016'
    cform = root/'research/HDBLAST_HARMONIC_CONSTRAINT_VALUES_20261030/reference/FORMULAS.json'
    gx = root/'research/HDBLAST_SHELL_DERIVATIVE_C2_20261026/checks/check_full23_shell_derivative.py'
    records = {}
    controls = []
    for mode in ('normal', 'optimized'):
        dest = out/mode
        dest.mkdir()
        flags = [sys.executable, '-I', '-S', '-B'] + (['-O'] if mode == 'optimized' else [])
        p = lambda name: dest/(name+'.json')

        def run(name, script, options=(), file_output=False):
            command = flags+[str(script)]+list(map(str, options))
            result = subprocess.run(command, cwd=root, capture_output=True, timeout=180)
            need(result.returncode == 0, name+' failed: '+result.stderr.decode(errors='replace')[-4000:])
            if not file_output:
                p(name).write_bytes(result.stdout)
            data = p(name).read_bytes()
            report = json.loads(data)
            need('PASS' in report.get('status', ''), name+' status')
            need(data == (cap/'reference'/(name+'.json')).read_bytes(), name+' reference mismatch')
            digest = hashlib.sha256(data).hexdigest()
            if name in records:
                need(records[name] == digest, name+' normal/optimized mismatch')
            else:
                records[name] = digest

        run('GX', gx, ('--repo-root', root))
        run('INITIAL_ABSOLUTE', cap/'checks/check_absolute_initial.py', ('--repo-root', root))
        run('FILTER', cap/'checks/check_causal_filter.py', ('--repo-root', root))
        run('DIFFERENTIATION', cap/'checks/check_time_differentiation.py')
        run('MASTER4', cap/'checks/check_master_fourth.py',
            ('--root', root, '--output', p('MASTER4')), True)
        run('MASTER4_AUDIT', cap/'checks/check_master4_review.py',
            ('--repo-root', root, '--producer-report', p('MASTER4')))
        run('CANDIDATE', cap/'checks/check_candidate_second.py',
            ('--repo-root', root, '--initial-absolute-report', p('INITIAL_ABSOLUTE'), '--gx-report', p('GX')))
        run('CANDIDATE_AUDIT', cap/'checks/check_candidate_bound_audit.py',
            ('--repo-root', root, '--candidate-report', p('CANDIDATE'),
             '--initial-report', p('INITIAL_ABSOLUTE'), '--gx-report', p('GX')))
        run('NATURAL', cap/'checks/check_natural_residual.py',
            ('--source-map', geom/'projection/GEOMETRIC_COEFFICIENT_MAPS.json',
             '--source-report', geom/'numeric/ACTUAL_SOURCES_NORMAL.json',
             '--time-report', hi/'TIME.json', '--constraint-formulas', cform,
             '--linear-operator', geom/'projection/inputs/linear_operator.py',
             '--jet-caps', p('CANDIDATE'), '--jet-caps-sha256', records['CANDIDATE'],
             '--higher-transfer', hi/'TRANSFER.json'))
        run('NATURAL_AUDIT', cap/'checks/check_natural_bound_audit.py',
            ('--natural-report', p('NATURAL'), '--candidate-report', p('CANDIDATE'),
             '--time-report', hi/'TIME.json', '--transfer-report', hi/'TRANSFER.json',
             '--source-report', geom/'numeric/ACTUAL_SOURCES_NORMAL.json'))
        for injection, reason in (
            ('drop_reflection_trace', 'exact reflected Krr identity'),
            ('drop_Hsecond_endpoint', 'incoming Isecond includes fx and fr endpoints'),
            ('drop_initial_force_derivatives', 'actual Uxx initial forcing derivatives retained'),
        ):
            target = dest/('corrupted-'+injection+'.json')
            result = subprocess.run(flags+[str(cap/'checks/check_master_fourth.py'),
                '--root', str(root), '--output', str(target), '--inject-error', injection],
                cwd=root, capture_output=True, timeout=180)
            need(result.returncode != 0 and reason in result.stderr.decode(errors='replace')
                 and not target.exists(), 'Wrong master formula was not rejected: '+injection)
            controls.append({'mode': mode, 'injection': injection, 'rejected': True})
    summary = {'status': 'PASS_CHECKPOINT_REPLAY', 'normal_optimized_reference_reports': len(records),
        'report_sha256': records, 'corrupted_master_runs_rejected': controls,
        'scope': 'Conditional mathematics; actual master C4 and auxiliary candidate bounds, not actual full F_x, C2 target accuracy, nonlinear realization or Big Bang evidence'}
    (out/'VERIFICATION.json').write_text(json.dumps(summary, indent=2, sort_keys=True)+'\n')
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
