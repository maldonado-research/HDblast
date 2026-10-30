#!/usr/bin/env python3
"""Fresh exact replays of refined error bounds and explicitly conditional transfers."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction as Q
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
        need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen, 'manifest path')
        seen.add(str(rel))
        p = cap/rel
        need(p.is_file() and not p.is_symlink(), 'manifest member type')
        data = p.read_bytes()
        need(len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'], 'capsule bytes '+str(rel))
    need(all(not p.is_symlink() and (not p.is_file() or str(p.relative_to(cap)) in seen)
             for p in cap.rglob('*')), 'closed capsule tree')
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
    checks = cap/'checks'
    def run_mode(mode):
        records, corruptions = {}, []
        dest = out/mode
        dest.mkdir()
        flags = [sys.executable, '-I', '-S', '-B']+(['-O'] if mode == 'optimized' else [])
        p = lambda name: dest/(name+'.json')

        def run(name, script, options=()):
            result = subprocess.run(flags+[str(script)]+list(map(str, options)),
                                    cwd=root, capture_output=True, timeout=1200)
            need(result.returncode == 0, name+' failed: '+result.stderr.decode(errors='replace')[-4000:])
            data = result.stdout
            report = json.loads(data)
            need(report.get('status', '').startswith('PASS'), name+' successful status')
            need(data == (cap/'reference'/(name+'.json')).read_bytes(), name+' reference mismatch')
            p(name).write_bytes(data)
            digest = hashlib.sha256(data).hexdigest()
            need(name not in records or records[name] == digest, name+' normal/optimized mismatch')
            records[name] = digest
            print(json.dumps({'mode': mode, 'completed': name}), flush=True)

        run('M3', checks/'check_nominal_m3_refinement.py', ('--repo-root', root))
        run('PAIRED', checks/'check_paired_cutoff_refinement.py', ('--repo-root', root))
        inputs = ('--repo-root', root, '--nominal-report', p('M3'), '--nominal-sha256', records['M3'],
                  '--paired-report', p('PAIRED'), '--paired-sha256', records['PAIRED'])
        run('COMBINED', checks/'check_combined_initial.py', inputs)
        run('INITIAL_AUDIT', checks/'check_combined_initial_audit.py', inputs)
        combined = json.loads(p('COMBINED').read_bytes())
        audit = json.loads(p('INITIAL_AUDIT').read_bytes())
        need(Q(combined['refined_uniform_initial_u_rr_H_upper']) == Q(audit['refined_uniform_initial_u_rr_H_upper']),
             'independent complete initial cap equals combined producer')
        run('ENDPOINT', root/'research/HDBLAST_HIGHER_INITIAL_JETS_20261030/checks/check_endpoint_curvature.py',
            ('--repo-root', root))
        run('C2', checks/'check_combined_c2.py',
            ('--repo-root', root, '--initial-refinement', p('COMBINED'), '--initial-refinement-sha256', records['COMBINED']))
        norm_inputs = ('--repo-root', root, '--combined-report', p('COMBINED'), '--combined-sha256', records['COMBINED'],
                       '--c2-report', p('C2'), '--c2-sha256', records['C2'])
        run('NORMALIZATION', checks/'check_coordinate_normalization.py', norm_inputs)
        run('NORMALIZATION_AUDIT', checks/'check_transfer_normalization_audit.py',
            inputs+('--initial-report', p('COMBINED'), '--initial-sha256', records['COMBINED'],
                    '--transfer-report', p('C2'), '--transfer-sha256', records['C2'],
                    '--normalization-report', p('NORMALIZATION'), '--normalization-sha256', records['NORMALIZATION']))
        natural_inputs = ('--repo-root', root, '--c2-report', p('C2'), '--c2-sha256', records['C2'])
        run('NATURAL', checks/'check_natural_error_transfer.py',
            natural_inputs+('--assume-first-harmonic-neighbourhood', '--assume-formal-subsidiary-completion'))
        run('NATURAL_AUDIT', checks/'check_natural_error_audit.py',
            natural_inputs+('--natural-report', p('NATURAL'), '--natural-sha256', records['NATURAL']))
        transfer = checks/'check_combined_c2.py'
        run('C2_SCOPE', checks/'check_combined_c2_scope_controls.py',
            ('--repo-root', root, '--initial-refinement', p('COMBINED'), '--initial-refinement-sha256', records['COMBINED'],
             '--checker', transfer, '--checker-sha256', hashlib.sha256(transfer.read_bytes()).hexdigest()))
        normalizer = checks/'check_coordinate_normalization.py'
        run('NORMALIZATION_SCOPE', checks/'check_coordinate_scope_controls.py',
            norm_inputs+('--checker', normalizer, '--checker-sha256', hashlib.sha256(normalizer.read_bytes()).hexdigest()))
        natural = checks/'check_natural_error_transfer.py'
        run('NATURAL_SCOPE', checks/'check_natural_scope_controls.py',
            natural_inputs+('--checker', natural, '--checker-sha256', hashlib.sha256(natural.read_bytes()).hexdigest()))
        for script, cases in (
            ('check_nominal_m3_refinement.py', (
                ('missing-cover-cell', 'complete nominal parent cover'), ('wrong-branch', 'forced plain nominal branch'),
                ('missing-factorial', 'ordinary third derivative factorial'), ('missing-ell-cube', 'scaled radial third derivative ell cube'),
                ('omit-inner-jump', 'retained inner jump allowance'), ('parent-maximum', 'parent cap intersection uses minimum'))),
            ('check_paired_cutoff_refinement.py', (
                ('omit-curvature-factor', 'ordinary scaled radial curvature factor'),
                ('omit-radial-scale', 'ordinary scaled radial curvature factor'),
                ('changed-width', 'same fixed convolution widths'), ('wrong-parent-cover', 'authenticated exact ten-leaf parent cover'))),
        ):
            for injection, reason in cases:
                result = subprocess.run(flags+[str(checks/script), '--repo-root', str(root), '--inject-error', injection],
                                        cwd=root, capture_output=True, timeout=120)
                need(result.returncode != 0 and not result.stdout and reason in result.stderr.decode(errors='replace'),
                     'corrupted refinement formula not rejected: '+injection)
                corruptions.append({'mode': mode, 'checker': script, 'injection': injection, 'rejected': True})
        return records, corruptions

    # Modes have disjoint output directories and fresh child processes. All
    # scientific inputs are read-only; the two independent executions can
    # therefore use the runner's two CPU workers without shared script state.
    with ThreadPoolExecutor(max_workers=2) as pool:
        normal_future = pool.submit(run_mode, 'normal')
        optimized_future = pool.submit(run_mode, 'optimized')
        normal = normal_future.result()
        optimized = optimized_future.result()
    records = normal[0]
    need(records == optimized[0], 'complete normal/optimized report set equality')
    corruptions = normal[1]+optimized[1]
    need(authenticate(cap) == manifest, 'unchanged capsule readback')
    result = {'status': 'PASS_REFINED_ERROR_AND_CONDITIONAL_TRANSFER_CHECKPOINT',
        'normal_optimized_reference_reports': len(records), 'report_sha256': records,
        'corrupted_refinement_runs_rejected': corruptions,
        'corrupted_scope_runs_per_mode': {'C2': 8, 'coordinate_normalization': 7, 'natural_transfer': 18},
        'candidate': 'unchanged W1_quadratic; target and original sources unchanged',
        'scope': 'Global coordinate seminorm accuracy coexists with a local relative error above99.9%; natural residual transfer has additional explicit formal constraint assumptions. No physical smallness or nonlinear/cosmological realization.'}
    (out/'VERIFICATION.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
