#!/usr/bin/env python3
"""Fresh exact implementation joins, canonical reports and real corruptions."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile


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
             'capsule member bytes '+str(rel))
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
    # Independent scripts preserve their original command-line interfaces.
    scripts = (
        ('NATURAL_MASTER', 'check_natural_master.py', '--repo', '--run-controls',
         'wrong_formula_controls', 22),
        ('FIRST_SHELL', 'check_first_shell_master_binding.py', '--repo-root', '--controls',
         'corrupted_subprocesses', 12),
        ('NATURAL_INDEPENDENT', 'independent_natural_review.py', '--repo-root', '--run-controls',
         'corrupted_subprocess_controls', 22),
    )

    def flags(mode):
        return [sys.executable, '-I', '-S', '-B'] + (['-O'] if mode == 'optimized' else [])

    def reference(name, report):
        raw = report.read_bytes()
        need(raw == (cap/'reference'/(name+'.json')).read_bytes(),
             name+' complete canonical report differs')
        need(json.loads(raw).get('status', '').startswith('PASS'), name+' report status')
        return hashlib.sha256(raw).hexdigest()

    def run_mode(mode):
        dest = out/mode
        dest.mkdir()
        records = {}
        for name, filename, root_flag, control_flag, key, count in scripts:
            report = dest/(name+'.json')
            cmd = flags(mode)+[str(cap/'checks'/filename), root_flag, str(root),
                               '--output', str(report), control_flag]
            run = subprocess.run(cmd, cwd=root, capture_output=True, timeout=900)
            need(run.returncode == 0, mode+'/'+name+': '+run.stderr.decode(errors='replace')[-4000:])
            need(len(json.loads(report.read_bytes())[key]) == count, name+' control count')
            records[name] = reference(name, report)
            print(json.dumps({'mode': mode, 'completed': name,
                              'corrupted_processes_rejected': count}), flush=True)
        report = dest/'B0_INDEPENDENT.json'
        script = cap/'checks/review/check_full_frame_B0.py'
        run = subprocess.run(flags(mode)+[str(script), '--output', str(report)],
                             cwd=root, capture_output=True, timeout=180)
        need(run.returncode == 0, 'independent full-frame contraction: '+run.stderr.decode(errors='replace'))
        records['B0_INDEPENDENT'] = reference('B0_INDEPENDENT', report)
        with tempfile.TemporaryDirectory(prefix='hdblast-B0-control-') as tmp:
            corrupt = Path(tmp)/'must-not-exist.json'
            run = subprocess.run(flags(mode)+[str(script), '--output', str(corrupt),
                                 '--mutation', 'time_connection_sign'],
                                 cwd=root, capture_output=True, timeout=180)
            need(run.returncode != 0 and not corrupt.exists() and
                 b'ValueError: FULL_FRAME_B0_NONZERO:' in run.stderr,
                 'independent B0 control was not mathematically rejected')
        print(json.dumps({'mode': mode, 'completed': 'B0_INDEPENDENT',
                          'corrupted_processes_rejected': 1}), flush=True)
        return records

    def run_bianchi():
        dest = out/'bianchi'
        cmd = flags('normal')+[str(cap/'checks/run_bianchi_verification.py'),
                              '--repo-root', str(root), '--output-dir', str(dest), '--workers', '2']
        run = subprocess.run(cmd, cwd=root, capture_output=True, timeout=1200)
        need(run.returncode == 0, 'Bianchi replay: '+run.stderr.decode(errors='replace')[-4000:])
        report = json.loads((dest/'VERIFICATION.json').read_bytes())
        controls = report['fresh_corrupted_subprocesses']
        need(report['normal_optimized_reports_byte_identical'] and report['closure_pre_post_equal'],
             'Bianchi source/decision agreement')
        need(len(controls) == 26 and
             all(r['returncode'] != 0 and r['successful_report_absent'] and
                 r['exact_nonzero_witness'] for r in controls), 'Bianchi failed-control roster')
        need({r['mode'] for r in controls} == {'NORMAL', 'OPTIMIZED'} and
             all(sum(r['mode'] == mode for r in controls) == 13 for mode in ('NORMAL', 'OPTIMIZED')),
             'Bianchi mode roster')
        normal = reference('BULK_BIANCHI', dest/'NORMAL.json')
        optimized = reference('BULK_BIANCHI', dest/'OPTIMIZED.json')
        need(normal == optimized, 'Bianchi canonical mode agreement')
        print(json.dumps({'completed': 'BULK_BIANCHI', 'corrupted_processes_rejected': 26}), flush=True)
        return normal

    with ThreadPoolExecutor(max_workers=3) as pool:
        normal = pool.submit(run_mode, 'normal')
        optimized = pool.submit(run_mode, 'optimized')
        bianchi = pool.submit(run_bianchi)
        nr, opt, br = normal.result(), optimized.result(), bianchi.result()
    need(nr == opt, 'normal/optimized complete report agreement')
    nr['BULK_BIANCHI'] = br
    need(authenticate(cap) == manifest, 'capsule unchanged after readback')
    receipt = {
        'status': 'PASS_NATURAL_FIRST_FIELD_SHELL_AND_BULK_SOURCE_BINDING',
        'canonical_reference_report_SHA256': nr,
        'scientific_reports_reproduced_in_normal_and_optimized_modes': 5,
        'rejected_corrupted_executions_per_mode': {
            'NATURAL_MASTER': 22, 'FIRST_SHELL': 12, 'NATURAL_INDEPENDENT': 22,
            'B0_INDEPENDENT': 1, 'BULK_BIANCHI': 13},
        'rejected_corrupted_executions_total': 140,
        'natural_Bianchi_modal_rows': 8,
        'Bianchi_input_scope': 'Registered seven-profile first quadrupole ansatz; coherent ordinary metric jets through 3 and scalar jets through 2.',
        'first_residual_terms_retained_before_specialization': True,
        'actual_background_master_gauge_PDE_and_regularity_are_retained_analytic_inputs': True,
        'frozen_sources_and_numerical_candidates_unchanged': True,
        'scope': 'Exact implementation binding to existing analytic proofs. No new PDE existence theorem, quantitative error improvement, exact finite-amplitude nonlinear family, external specialist assessment or cosmological realization.',
    }
    (out/'VERIFICATION.json').write_text(json.dumps(receipt, sort_keys=True, indent=2)+'\n')
    print(json.dumps(receipt, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
