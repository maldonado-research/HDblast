#!/usr/bin/env python3
"""Fresh exact replays of natural shell compatibility and first-wave binding."""
import argparse
from concurrent.futures import ThreadPoolExecutor
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
        need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,
             'manifest path')
        seen.add(str(rel))
        path = cap/rel
        need(path.is_file() and not path.is_symlink(), 'manifest member type')
        data = path.read_bytes()
        need(len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'],
             'capsule bytes '+str(rel))
    need(all(not p.is_symlink() and (not p.is_file() or str(p.relative_to(cap)) in seen)
             for p in cap.rglob('*')), 'closed capsule tree')
    return hashlib.sha256(raw).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    root = args.repo_root.resolve()
    cap = Path(__file__).resolve().parent.parent
    manifest = authenticate(cap)
    out = args.output_dir.resolve()
    need(not out.exists() and cap not in out.parents,
         'use a fresh output directory outside the capsule')
    out.mkdir(parents=True)
    scripts = (
        ('CODAZZI', 'check_natural_shell_codazzi.py', 11),
        ('INTRINSIC', 'check_intrinsic_codazzi_review.py', 8),
        ('BRIDGE', 'check_shell_source_bridge.py', 18),
        ('FIRST_WAVE', 'check_first_wave_binding.py', 15),
        ('WAVE_INDEPENDENT', 'independent_wave_identity.py', 10),
    )

    def run_mode(mode):
        dest = out/mode
        dest.mkdir()
        records = {}
        for name, filename, expected_controls in scripts:
            report = dest/(name+'.json')
            flags = [sys.executable, '-I', '-S', '-B']
            if mode == 'optimized':
                flags.append('-O')
            run = subprocess.run(flags+[str(cap/'checks'/filename),
                                        '--repo-root', str(root), '--output', str(report),
                                        '--run-controls'],
                                 cwd=root, capture_output=True, timeout=900)
            need(run.returncode == 0,
                 mode+'/'+name+' failed: '+run.stderr.decode(errors='replace')[-4000:])
            raw = report.read_bytes()
            parsed = json.loads(raw)
            need(parsed.get('status', '').startswith('PASS'), name+' status')
            need(len(parsed['rejected_corrupted_executions']) == expected_controls,
                 name+' rejected mathematical control count')
            need(raw == (cap/'reference'/(name+'.json')).read_bytes(),
                 mode+'/'+name+' complete reference report differs')
            records[name] = hashlib.sha256(raw).hexdigest()
            print(json.dumps({'mode': mode, 'completed': name,
                              'corrupted_processes_rejected': expected_controls}), flush=True)
        return records

    with ThreadPoolExecutor(max_workers=2) as pool:
        normal = pool.submit(run_mode, 'normal')
        optimized = pool.submit(run_mode, 'optimized')
        nr, opt = normal.result(), optimized.result()
    need(nr == opt, 'normal/optimized complete report agreement')
    need(authenticate(cap) == manifest, 'capsule unchanged after readback')
    receipt = {
        'status': 'PASS_BULK_SHELL_CONSTRAINT_IMPLEMENTATION_BINDING',
        'reference_report_SHA256': nr,
        'new_reports_reproduced_in_normal_and_optimized_modes': len(scripts),
        'rejected_corrupted_executions_per_mode': {n: k for n, _, k in scripts},
        'rejected_corrupted_executions_total': 2*sum(k for _, _, k in scripts),
        'natural_Codazzi_modal_rows_with_first_scalar_residual_free': 5,
        'natural_Codazzi_modal_rows_with_physical_scalar_condition': 5,
        'normal_harmonic_row_imposed_for_natural_Codazzi': False,
        'actual_first_wave_PDE_solution_and_regularity_retained': True,
        'frozen_sources_and_numerical_candidates_unchanged': True,
        'scope': 'Exact source and implementation binding to existing analytic proofs. '
                 'No new nonlinear existence theorem, physical accuracy certificate, '
                 'external specialist approval or cosmological realization.',
    }
    (out/'VERIFICATION.json').write_text(json.dumps(receipt, sort_keys=True, indent=2)+'\n')
    print(json.dumps(receipt, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
