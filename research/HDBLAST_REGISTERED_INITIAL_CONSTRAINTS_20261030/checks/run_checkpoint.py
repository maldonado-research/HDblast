#!/usr/bin/env python3
"""Reproduce the exact initial-data binding and its historical proof checks."""
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
        p = cap/rel
        need(p.is_file() and not p.is_symlink(), 'manifest member type')
        data = p.read_bytes()
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

    historical = (
        ('HISTORICAL_PHYSICAL', 'HDBLAST_EXPLICIT_SECOND_INITIAL_JET_20261015',
         '1b27636663347db5b7dc3112d5f157f69f8d615cf2a95df2c1309a06e0bf8e86'),
        ('HISTORICAL_HARMONIC', 'HDBLAST_SECOND_INITIAL_NORMS_GAUGE_20261015',
         '4c0dae7de0b63acca15e5204681a859f49e5b5028340457dea5c29246bf0626c'),
    )
    records = {}
    for name, folder, pin in historical:
        script = root/'research'/folder/'reproduce.py'
        need(hashlib.sha256(script.read_bytes()).hexdigest() ==
             'd1827d61adf34075073964681e81881be83060660885121a975134be879e5002',
             'historical runner pin')
        dest = out/name.lower()
        run = subprocess.run([sys.executable, '-I', '-S', '-B', str(script),
                              '--manifest-sha256', pin, '--output-dir', str(dest)],
                             cwd=root, capture_output=True, timeout=900)
        need(run.returncode == 0,
             name+' failed: '+run.stderr.decode(errors='replace')[-4000:])
        raw = (dest/'REPRODUCTION_RECORD.json').read_bytes()
        need(raw == (cap/'reference'/(name+'.json')).read_bytes(),
             name+' complete historical reproduction record differs')
        records[name] = hashlib.sha256(raw).hexdigest()
        print(json.dumps({'completed': name}), flush=True)

    scripts = (
        ('PHYSICAL', 'check_physical_initial_binding.py'),
        ('ADM', 'check_registered_adm_binding.py'),
        ('INDEPENDENT', 'check_independent_profile_reductions.py'),
    )

    def run_mode(mode):
        dest = out/mode
        dest.mkdir()
        result = {}
        for name, filename in scripts:
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
            need(raw == (cap/'reference'/(name+'.json')).read_bytes(),
                 mode+'/'+name+' complete reference report differs')
            result[name] = hashlib.sha256(raw).hexdigest()
            print(json.dumps({'mode': mode, 'completed': name}), flush=True)
        return result

    with ThreadPoolExecutor(max_workers=2) as pool:
        normal = pool.submit(run_mode, 'normal')
        optimized = pool.submit(run_mode, 'optimized')
        nr, opt = normal.result(), optimized.result()
    need(nr == opt, 'normal/optimized complete report agreement')
    records.update(nr)
    control_fields = {'PHYSICAL': ('wrong_formula_controls', 8),
                      'ADM': ('rejected_corrupted_executions', 14),
                      'INDEPENDENT': ('mutations_rejected', 8)}
    controls_per_mode = {}
    for name, (field, expected) in control_fields.items():
        report = json.loads((out/'normal'/(name+'.json')).read_bytes())
        need(len(report[field]) == expected, name+' mathematical control count')
        controls_per_mode[name] = expected
    need(authenticate(cap) == manifest, 'capsule unchanged after readback')
    receipt = {
        'status': 'PASS_REGISTERED_INITIAL_CONSTRAINT_BINDING_CHECKPOINT',
        'reference_report_SHA256': records,
        'new_reports_reproduced_in_normal_and_optimized_modes': len(scripts),
        'historical_complete_replay_rosters': len(historical),
        'rejected_new_corrupted_executions_per_mode': controls_per_mode,
        'rejected_new_corrupted_executions_total': 2*sum(controls_per_mode.values()),
        'scope': 'Exact implementation binding to retained analytic initial physical/ADM proofs. '
                 'The evolution, physical boundary and subsidiary proof chain remain explicit '
                 'analytic inputs. No finite-amplitude nonlinear or cosmological realization.',
        'frozen_sources_unchanged': True,
    }
    (out/'VERIFICATION.json').write_text(json.dumps(receipt, sort_keys=True, indent=2)+'\n')
    print(json.dumps(receipt, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
