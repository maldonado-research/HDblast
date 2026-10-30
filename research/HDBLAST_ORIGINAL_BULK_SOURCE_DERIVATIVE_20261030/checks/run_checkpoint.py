#!/usr/bin/env python3
"""Reproduce original-source derivative and original-target C2 estimates."""
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
    geom = root/'research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016'
    master4 = root/'research/HDBLAST_CAUSAL_SOURCE_REGULARITY_20261030/reference/MASTER4.json'
    master3 = geom/'numeric/inputs/HIGHER_JETS_NORMAL.json'
    source_map = geom/'projection/GEOMETRIC_COEFFICIENT_MAPS.json'
    identity = cap/'checks/check_bulk_derivative_identities.py'
    identity_sha = hashlib.sha256(identity.read_bytes()).hexdigest()
    records = {}
    controls = []
    for mode in ('normal', 'optimized'):
        dest = out/mode
        dest.mkdir()
        flags = [sys.executable, '-I', '-S', '-B'] + (['-O'] if mode == 'optimized' else [])
        p = lambda name: dest/(name+'.json')

        def run(name, script, options=(), file_output=False):
            result = subprocess.run(flags+[str(script)]+list(map(str, options)),
                cwd=root, capture_output=True, timeout=480)
            need(result.returncode == 0, name+' failed: '+result.stderr.decode(errors='replace')[-4000:])
            if not file_output:
                p(name).write_bytes(result.stdout)
            data = p(name).read_bytes()
            report = json.loads(data)
            need(report.get('status', '').startswith('PASS'), name+' status')
            need(data == (cap/'reference'/(name+'.json')).read_bytes(), name+' reference mismatch')
            digest = hashlib.sha256(data).hexdigest()
            if name in records:
                need(records[name] == digest, name+' normal/optimized mismatch')
            else:
                records[name] = digest

        run('GAUGE4', cap/'checks/check_actual_gauge_fourth.py', ('--repo-root', root))
        run('GAUGE4_AUDIT', cap/'checks/check_gauge4_review.py',
            ('--repo-root', root, '--producer-report', p('GAUGE4')))
        run('FX_IDENTITIES', identity, ('--source-map', source_map))
        run('FX', cap/'checks/check_original_bulk_fx.py',
            ('--repo-root', root, '--master4-report', master4,
             '--gauge4-report', p('GAUGE4'), '--gauge4-sha256', records['GAUGE4'], '--output', p('FX')), True)
        run('FX_AUDIT', cap/'checks/check_bulk_derivative_numeric.py',
            ('--fx-report', p('FX'), '--fx-sha256', records['FX'], '--source-map', source_map,
             '--identity-checker', identity, '--identity-checker-sha256', identity_sha,
             '--gauge4-report', p('GAUGE4'), '--master4-report', master4, '--master3-report', master3))
        run('C2', cap/'checks/check_original_target_c2.py',
            ('--repo-root', root, '--fx-report', p('FX'), '--fx-sha256', records['FX']))
        run('C2_AUDIT', cap/'checks/check_original_c2_audit.py',
            ('--repo-root', root, '--gate-report', p('C2'), '--fx-report', p('FX'), '--fx-sha256', records['FX']))
        run('C2_SCOPE', cap/'checks/check_original_c2_scope_controls.py',
            ('--repo-root', root, '--fx-report', p('FX')))
        for injection, reason in (
            ('normalized_taylor_order3', 'ordinary third cosh-minus-three-halves derivative polynomial'),
            ('drop_initial_P_squared', 'stabilized initial fourth jet retains noncommuting P0 squared source'),
            ('flip_radial_recovery_sign', 'PDE fourth radial recovery sign'),
        ):
            result = subprocess.run(flags+[str(cap/'checks/check_actual_gauge_fourth.py'),
                '--repo-root', str(root), '--inject-error', injection],
                cwd=root, capture_output=True, timeout=60)
            need(result.returncode != 0 and reason in result.stderr.decode(errors='replace')
                 and not result.stdout, 'Wrong gauge formula was not rejected: '+injection)
            controls.append({'mode': mode, 'injection': injection, 'rejected': True})
    summary = {'status': 'PASS_ORIGINAL_SOURCE_AND_C2_CHECKPOINT_REPLAY',
        'normal_optimized_reference_reports': len(records), 'report_sha256': records,
        'corrupted_gauge_runs_rejected': controls,
        'wrong_scope_runs_per_scope_report': 8,
        'scope': 'Conditional actual original source derivative and unchanged-candidate absolute C2 error bounds; no source averaging, relative physical accuracy, nonlinear family or Big Bang evidence'}
    (out/'VERIFICATION.json').write_text(json.dumps(summary, indent=2, sort_keys=True)+'\n')
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
