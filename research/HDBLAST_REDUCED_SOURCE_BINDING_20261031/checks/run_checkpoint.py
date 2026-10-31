#!/usr/bin/env python3
"""Portable exact reduced-source and actual third-jet checkpoint replay."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

ADAPTER_PIN = '4847c79f034a09f423cc020ab3bef86e03d729970084a9bb088e94caa3b9b18a'


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


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
             'capsule member '+str(rel))
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
    closure = authenticate(cap)
    out = a.output_dir.resolve()
    need(not out.exists() and cap not in out.parents, 'fresh output outside capsule required')
    out.mkdir(parents=True)
    checks = cap/'checks'
    wave_data = out/'wave-data'
    wave_data.mkdir()
    # Read-only producer output is data for the independent coefficient comparison.
    for name in ('WAVE_NORMAL.json', 'WAVE_OPTIMIZED.json'):
        shutil.copyfile(cap/'reference/WAVE.json', wave_data/name)
    shutil.copyfile(checks/'covector_wave_rows.py', wave_data/'covector_wave_rows.py')

    def flags(optimized=False):
        return [sys.executable, '-I', '-S', '-B']+(['-O'] if optimized else [])

    def execute(command):
        done = subprocess.run(command, cwd=root, capture_output=True, timeout=1500)
        need(done.returncode == 0, 'fresh replay failed: '+done.stderr.decode(errors='replace')[-5000:])

    def reference(name, path):
        raw = path.read_bytes()
        need(raw == (cap/'reference'/(name+'.json')).read_bytes(), name+' complete canonical report differs')
        need(json.loads(raw).get('status', '').startswith('PASS'), name+' status')
        return hashlib.sha256(raw).hexdigest()

    def reduced():
        dest = out/'reduced'
        execute(flags()+[str(checks/'run_reduced_verification.py'), '--repo-root', str(root),
                         '--output-dir', str(dest), '--workers', '2'])
        receipt = json.loads((dest/'VERIFICATION.json').read_bytes())
        rows = receipt['fresh_corrupted_subprocesses']
        need(receipt['normal_optimized_reports_byte_identical'] and receipt['closure_pre_post_equal'],
             'reduced mode/closure agreement')
        need(len(rows) == 40 and all(r['returncode'] != 0 and r['successful_report_absent'] and
             r['exact_nonzero_witness'] for r in rows), 'reduced corruption roster')
        nr = reference('REDUCED', dest/'NORMAL.json')
        need(nr == reference('REDUCED', dest/'OPTIMIZED.json'), 'reduced mode report')
        print(json.dumps({'completed': 'REDUCED', 'mathematical_corruptions_rejected': 40}), flush=True)
        return {'REDUCED': nr}

    def independent_roles():
        dest = out/'independent-roles'
        execute(flags()+[str(checks/'run_independent_review.py'), '--root', str(root),
                         '--wave-dir', str(wave_data), '--output-dir', str(dest)])
        receipt = json.loads((dest/'RUNTIME_RECEIPT.json').read_bytes())
        need(receipt['canonical_normal_optimized_reports_identical'] and
             receipt['rejected_corrupted_subprocess_count'] == 40 and
             len(receipt['independent_eight_wave_row_data_comparisons']) == 8,
             'independent role/wave comparison roster')
        rows = [r for r in receipt['executions'] if r['mutation'] is not None]
        need(len(rows) == 40 and all(r['returncode'] != 0 and not r['successful_output_exists'] and
             r['exact_nonzero_collected_witness'] for r in rows), 'independent role corruptions')
        nr = reference('INDEPENDENT_ROLES', dest/'SOURCE_ROLES_NORMAL.json')
        need(nr == reference('INDEPENDENT_ROLES', dest/'SOURCE_ROLES_OPTIMIZED.json'), 'independent role mode')
        print(json.dumps({'completed': 'INDEPENDENT_ROLES', 'mathematical_corruptions_rejected': 40}), flush=True)
        return {'INDEPENDENT_ROLES': nr}

    def wave_and_jets():
        digests = {}
        for optimized in (False, True):
            mode = 'optimized' if optimized else 'normal'
            dest = out/mode
            dest.mkdir()
            for name, filename, extra, count in (
                ('WAVE', 'check_covector_wave.py', ['--run-controls'], 13),
                ('ACTUAL_THIRD', 'check_actual_third_jet_adapter.py',
                 ['--repo-root', str(root), '--run-controls'], 10),
                ('INDEPENDENT_THIRD', 'check_independent_third_jets.py',
                 ['--repo-root', str(root), '--adapter', str(checks/'actual_first_third_jets.py'),
                  '--adapter-sha256', ADAPTER_PIN], 0),
            ):
                path = dest/(name+'.json')
                execute(flags(optimized)+[str(checks/filename), '--output', str(path)]+extra)
                if count:
                    need(len(json.loads(path.read_bytes())['rejected_corrupted_executions']) == count,
                         name+' corruption count')
                digest = reference(name, path)
                need(name not in digests or digests[name] == digest, name+' mode agreement')
                digests[name] = digest
                print(json.dumps({'mode': mode, 'completed': name,
                                  'mathematical_corruptions_rejected': count}), flush=True)
        witness = out/'witness'
        execute(flags()+[str(checks/'run_witness_verification.py'), '--repo-root', str(root),
                         '--output-dir', str(witness)])
        receipt = json.loads((witness/'VERIFICATION.json').read_bytes())
        rows = receipt['fresh_negative_subprocesses']
        need(receipt['baseline_reports_byte_identical'] and len(rows) == 4 and
             all(r['returncode'] != 0 and r['successful_report_absent'] and r['exact_nonzero_witness']
                 for r in rows), 'coherent witness corruption roster')
        digests['WITNESS'] = reference('WITNESS', witness/'NORMAL.json')
        need(digests['WITNESS'] == reference('WITNESS', witness/'OPTIMIZED.json'), 'witness mode report')
        print(json.dumps({'completed': 'WITNESS', 'mathematical_corruptions_rejected': 4}), flush=True)
        corrupt = out/'third-disk-corruptions'
        corrupt.mkdir()
        path = corrupt/'REAL_CORRUPTION_RECEIPTS.json'
        execute(flags()+[str(checks/'run_real_corruptions.py'), '--repo-root', str(root),
                         '--adapter', str(checks/'actual_first_third_jets.py'),
                         '--checker', str(checks/'check_independent_third_jets.py'), '--output', str(path)])
        receipt = json.loads(path.read_bytes())
        need(receipt['semantic_cases'] == 10 and receipt['source_authentication_cases'] == 6 and
             receipt['fresh_subprocesses'] == 32 and len(receipt['executions']) == 32,
             'real disk corruption roster')
        need(all(r['returncode'] != 0 and r.get('output_absent', r.get('side_effect_absent', False))
                 and r['witness'] for r in receipt['executions']), 'real disk corruption failures')
        print(json.dumps({'completed': 'THIRD_DISK_CONTROLS', 'mathematical_corruptions_rejected': 20,
                          'source_authentication_or_readback_cases_rejected': 12}), flush=True)
        return digests

    with ThreadPoolExecutor(max_workers=3) as pool:
        jobs = [pool.submit(reduced), pool.submit(independent_roles), pool.submit(wave_and_jets)]
        reports = {}
        for job in jobs:
            reports.update(job.result())
    need(len(reports) == 6 and authenticate(cap) == closure, 'final report/closed-source agreement')
    receipt = {
        'status': 'PASS_REDUCED_SOURCE_WAVE_AND_ACTUAL_THIRD_JET_BINDING',
        'canonical_reference_report_SHA256': reports,
        'scientific_reports_reproduced_in_normal_and_optimized_modes': 6,
        'mathematical_corruptions_rejected_per_mode': {
            'REDUCED': 20, 'INDEPENDENT_ROLES': 20, 'WAVE': 13, 'ACTUAL_THIRD': 10,
            'WITNESS': 2, 'INDEPENDENT_THIRD_DISK': 10},
        'mathematical_corruptions_rejected_total': 150,
        'source_authentication_import_or_readback_corruptions_rejected_total': 12,
        'corrupted_subprocesses_rejected_total': 162,
        'general_reduced_Bianchi_modal_rows': 8,
        'actual_profile_jet_slots_through_order_three': {'metric': 70, 'scalar': 10},
        'new_pure_radial_third_profile_entries': 8,
        'earlier_lower_two_and_24_time_containing_third_intervals_preserved': True,
        'third_bounds_are_coarse_uniform_magnitudes_under_inherited_actual_certificates': True,
        'old_sources_maps_numerical_candidates_and_error_certificates_unchanged': True,
        'scope': 'Exact implementation/proof joins and additive actual derivative enclosures. Analytic PDE/regularity and subsidiary premises remain; no new existence theorem, approximation error improvement, finite-amplitude nonlinear family or cosmological realization.',
    }
    (out/'VERIFICATION.json').write_text(json.dumps(receipt, sort_keys=True, indent=2)+'\n')
    print(json.dumps(receipt, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
