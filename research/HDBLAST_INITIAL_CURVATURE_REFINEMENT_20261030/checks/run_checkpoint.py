#!/usr/bin/env python3
"""Replay exact initial diagnostics, refinement, audits and C2 transfer."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def need(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def authenticate(cap):
    raw = (cap/'MANIFEST.json').read_bytes()
    seen = {'MANIFEST.json'}
    for row in json.loads(raw)['files']:
        rel = Path(row['path'])
        need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen, 'manifest member path')
        seen.add(str(rel))
        p = cap/rel
        need(p.is_file() and not p.is_symlink(), 'manifest member type')
        data = p.read_bytes()
        need(len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'], 'capsule bytes '+str(rel))
    need(all(not p.is_symlink() and (not p.is_file() or str(p.relative_to(cap)) in seen)
             for p in cap.rglob('*')), 'closed capsule file tree')
    return hashlib.sha256(raw).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    a = ap.parse_args()
    root = a.repo_root.resolve()
    cap = Path(__file__).resolve().parent.parent
    manifest = authenticate(cap)
    out = a.output_dir.resolve()
    need(not out.exists() and cap not in out.parents, 'use a fresh output directory outside capsule')
    out.mkdir(parents=True)
    records, counts, corruptions = {}, {}, []
    checks = cap/'checks'
    for mode in ('normal', 'optimized'):
        dest = out/mode
        dest.mkdir()
        flags = [sys.executable, '-I', '-S', '-B']+(['-O'] if mode == 'optimized' else [])
        p = lambda name: dest/(name+'.json')

        def run(name, script, options=()):
            result = subprocess.run(flags+[str(checks/script)]+list(map(str, options)),
                                    cwd=root, capture_output=True, timeout=480)
            need(result.returncode == 0, name+' failed: '+result.stderr.decode(errors='replace')[-4000:])
            data = result.stdout
            report = json.loads(data)
            need(report.get('status', '').startswith('PASS'), name+' successful status')
            need(data == (cap/'reference'/(name+'.json')).read_bytes(), name+' reference mismatch')
            p(name).write_bytes(data)
            digest = hashlib.sha256(data).hexdigest()
            need(name not in records or records[name] == digest, name+' mode mismatch')
            records[name] = digest
            if name == 'JOINS':
                counts[name] = report['checks']['group_count']
            elif name == 'JOINS_AUDIT':
                counts[name] = report['group_count']
            else:
                counts[name] = len(report.get('named_checks', []))

        run('JOINS', 'check_inner_join_error.py', ('--repo-root', root))
        run('JOINS_AUDIT', 'check_inner_join_independent.py',
            ('--repo-root', root, '--producer-report', p('JOINS')))
        run('INITIAL', 'check_refined_initial_curvature.py', ('--repo-root', root))
        run('INITIAL_AUDIT', 'check_initial_refinement_audit.py',
            ('--repo-root', root, '--refined-report', p('INITIAL'), '--refined-sha256', records['INITIAL']))
        run('C2', 'check_refined_initial_c2.py',
            ('--repo-root', root, '--initial-refinement', p('INITIAL'), '--initial-refinement-sha256', records['INITIAL']))
        run('C2_AUDIT', 'check_refined_transfer_audit.py',
            ('--repo-root', root, '--refinement', p('INITIAL'), '--refinement-sha256', records['INITIAL'],
             '--transfer', p('C2'), '--transfer-sha256', records['C2']))
        transfer = checks/'check_refined_initial_c2.py'
        run('C2_SCOPE', 'check_refined_c2_scope_controls.py',
            ('--repo-root', root, '--initial-refinement', p('INITIAL'), '--initial-refinement-sha256', records['INITIAL'],
             '--checker', transfer, '--checker-sha256', hashlib.sha256(transfer.read_bytes()).hexdigest()))
        for injection, reason in (
            ('half-mass', 'kernel_half_mass'), ('half-moment', 'kernel_half_first_moment'),
            ('curvature-factor', 'ordinary_scaled_curvature_factor'), ('join-branches', 'exact_one_sided_branches'),
            ('cutoff-derivative', 'packet_matching_derivative_0'), ('gram-weight', 'scalar_mean_positive_Gram_weight'),
        ):
            result = subprocess.run(flags+[str(checks/'check_inner_join_error.py'),
                '--repo-root', str(root), '--inject-error', injection], cwd=root, capture_output=True, timeout=60)
            need(result.returncode != 0 and not result.stdout and reason in result.stderr.decode(errors='replace'),
                 'corrupted formula not rejected: '+injection)
            corruptions.append({'mode': mode, 'injection': injection, 'rejected': True})
    need(authenticate(cap) == manifest, 'unchanged capsule readback')
    result = {'status': 'PASS_INITIAL_CURVATURE_REFINEMENT_CHECKPOINT',
        'normal_optimized_reference_reports': len(records), 'report_sha256': records,
        'named_check_counts': counts, 'corrupted_join_formula_runs_rejected': corruptions,
        'corrupted_C2_input_file_runs_per_mode': 8,
        'candidate': 'unchanged W1_quadratic; actual data and sources unchanged',
        'scope': 'Conditional coordinate initial-error diagnostics and tighter absolute C2 upper bounds; no true-error reduction, physical relative accuracy or nonlinear/cosmological realization'}
    (out/'VERIFICATION.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
