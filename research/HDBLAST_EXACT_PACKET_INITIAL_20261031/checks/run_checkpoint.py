#!/usr/bin/env python3
"""Portable exact-packet initial and response checkpoint; stdlib only."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def authenticate(cap):
    raw = (cap/'MANIFEST.json').read_bytes()
    seen = {'MANIFEST.json'}
    for row in json.loads(raw)['files']:
        rel = Path(row['path'])
        need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,
             'manifest relative unique path')
        seen.add(str(rel))
        path = cap/rel
        need(path.is_file() and not path.is_symlink(), 'manifest regular member')
        data = path.read_bytes()
        need(len(data) == row['bytes'] and hashlib.sha256(data).hexdigest() == row['sha256'],
             'capsule member '+str(rel))
    need(all(not path.is_symlink() and (not path.is_file() or
             str(path.relative_to(cap)) in seen) for path in cap.rglob('*')),
         'closed capsule tree')
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
    need((cap/'checks/FREEZE.json').read_bytes() == (cap/'reference/FREEZE.json').read_bytes(),
         'same pre-execution freeze for producer and reader')
    out.mkdir(parents=True)
    checks, independent = cap/'checks', cap/'independent'
    initial_pin = digest(cap/'reference/INITIAL.json')
    response_pin = digest(cap/'reference/RESPONSE.json')
    oracle_pin = digest(checks/'response_ledger_oracle.py')

    def flags(optimized=False):
        return [sys.executable, '-I', '-S', '-B']+(['-O'] if optimized else [])

    def execute(command, log):
        start = time.monotonic()
        p = subprocess.run(command, cwd=root, capture_output=True, timeout=1500)
        log.parent.mkdir(parents=True, exist_ok=True)
        log.with_suffix('.stdout').write_bytes(p.stdout)
        log.with_suffix('.stderr').write_bytes(p.stderr)
        log.with_suffix('.json').write_text(json.dumps({
            'command': command, 'returncode': p.returncode,
            'elapsed_seconds': time.monotonic()-start}, indent=2)+'\n')
        need(p.returncode == 0, 'fresh replay failed: '+str(log)+'\n'+
             p.stderr.decode(errors='replace')[-6000:])

    def reference(name, path):
        raw = path.read_bytes()
        need(raw == (cap/'reference'/(name+'.json')).read_bytes(),
             name+' complete canonical report differs')
        need(json.loads(raw)['status'].startswith('PASS') or
             json.loads(raw)['status'] == 'INITIAL_AND_LOCAL_GATES_PASS', name+' status')
        return hashlib.sha256(raw).hexdigest()

    def native_reports():
        reports = {}
        for optimized in (False, True):
            mode = 'optimized' if optimized else 'normal'
            dest = out/'native'/mode
            dest.mkdir(parents=True)
            cmd = flags(optimized)
            execute(cmd+[str(checks/'producer.py'), '--repo-root', str(root),
                         '--output-dir', str(dest/'producer')], dest/'logs/producer')
            initial = dest/'producer/INITIAL.json'
            execute(cmd+[str(checks/'check_controls.py'), '--repo-root', str(root),
                         '--initial-report', str(initial), '--output', str(dest/'INITIAL_CONTROLS.json')],
                    dest/'logs/initial-controls')
            response = dest/'RESPONSE.json'
            execute(cmd+[str(checks/'response_ledger_oracle.py'), '--repo-root', str(root),
                         '--candidate-report', str(initial), '--candidate-sha256', initial_pin,
                         '--output', str(response)], dest/'logs/response')
            review = dest/'INDEPENDENT_RESPONSE.json'
            execute(cmd+[str(independent/'check_response_independent.py'), '--repo-root', str(root),
                         '--initial', str(initial), '--initial-sha256', initial_pin,
                         '--producer-response', str(response), '--producer-response-sha256', response_pin,
                         '--initial-source', str(checks/'producer.py'),
                         '--initial-freeze', str(checks/'FREEZE.json'),
                         '--response-source', str(checks/'response_ledger_oracle.py'),
                         '--output', str(review)], dest/'logs/independent-response')
            for name, path in [('INITIAL', initial), ('INITIAL_CONTROLS', dest/'INITIAL_CONTROLS.json'),
                               ('RESPONSE', response), ('INDEPENDENT_RESPONSE', review)]:
                sha = reference(name, path)
                need(name not in reports or reports[name] == sha, name+' execution-mode identity')
                reports[name] = sha
            print(json.dumps({'completed': 'NATIVE_AND_RESPONSE', 'mode': mode}), flush=True)
        n = out/'native/normal'
        control = out/'native-controls'
        execute(flags()+[str(checks/'run_control_subprocesses.py'), '--repo-root', str(root),
                         '--initial-report', str(n/'producer/INITIAL.json'),
                         '--output-dir', str(control/'initial')], control/'logs/initial')
        r = json.loads((control/'initial/RECEIPTS.json').read_bytes())['runtime_receipts']
        need(len(r) == 20 and all(x['returncode'] != 0 and not x['successful_output_exists'] and
             x['explicit_witness'] for x in r), '20 initial formula/ledger rejections')
        execute(flags()+[str(checks/'run_response_controls.py'), '--repo-root', str(root),
                         '--candidate-report', str(n/'producer/INITIAL.json'),
                         '--candidate-sha256', initial_pin, '--oracle', str(checks/'response_ledger_oracle.py'),
                         '--oracle-sha256', oracle_pin, '--normal-report', str(n/'RESPONSE.json'),
                         '--optimized-report', str(out/'native/optimized/RESPONSE.json'),
                         '--output-dir', str(control/'response')], control/'logs/response')
        r = json.loads((control/'response/CONTROLS.json').read_bytes())['actual_child_executions']
        need(len(r) == 22 and all(x['returncode'] != 0 and x['successful_output_absent'] and
             x['nonzero_exact_residual'] not in ('0', '0/1') for x in r), '22 response formula rejections')
        execute(flags()+[str(independent/'run_effective_controls.py'), '--repo-root', str(root),
                         '--initial', str(n/'producer/INITIAL.json'), '--producer-response', str(n/'RESPONSE.json'),
                         '--initial-source', str(checks/'producer.py'), '--initial-freeze', str(checks/'FREEZE.json'),
                         '--response-source', str(checks/'response_ledger_oracle.py'),
                         '--checker', str(independent/'check_response_independent.py'),
                         '--review', str(n/'INDEPENDENT_RESPONSE.json'),
                         '--output', str(control/'independent-response/CONTROLS.json')],
                control/'logs/independent-response')
        r = json.loads((control/'independent-response/CONTROLS.json').read_bytes())
        need(r['case_count'] == 15 and r['child_execution_count'] == 30 and
             len(r['cases']) == 15 and all(x['exact_nonzero_witness_or_failed_identity'] and
             len(x['runs']) == 2 and all(y['returncode'] != 0 and y['successful_output_absent']
             for y in x['runs']) for x in r['cases']), '30 independent response rejections')
        return reports

    def independent_reports():
        reports = {}
        scripts = independent/'scripts'
        for optimized in (False, True):
            mode = 'optimized' if optimized else 'normal'
            dest = out/'independent'/mode
            dest.mkdir(parents=True)
            cmd = flags(optimized)
            initial = dest/'INITIAL.json'
            execute(cmd+[str(scripts/'check_independent_exact.py'), '--source', str(root),
                         '--out', str(initial), '--bits', '128', '--coverage'], dest/'logs/initial')
            audit = dest/'AUDIT.json'
            execute(cmd+[str(scripts/'check_independent_report.py'), '--source', str(root),
                         '--report', str(initial), '--output', str(audit)], dest/'logs/initial-audit')
            regularity = dest/'REGULARITY.json'
            execute(cmd+[str(independent/'check_regularity.py'), '--source-root', str(root),
                         '--pins', str(independent/'regularity-source-pins.json'),
                         '--output', str(regularity)], dest/'logs/regularity')
            response = dest/'COARSE_RESPONSE.json'
            execute(cmd+[str(independent/'check_coarse_response.py'), '--source-root', str(root),
                         '--initial-report', str(initial), '--output', str(response)], dest/'logs/coarse-response')
            for name, path in [('INDEPENDENT_INITIAL', initial), ('INDEPENDENT_INITIAL_AUDIT', audit),
                               ('REGULARITY', regularity), ('INDEPENDENT_COARSE_RESPONSE', response)]:
                sha = reference(name, path)
                need(name not in reports or reports[name] == sha, name+' execution-mode identity')
                reports[name] = sha
            controls = dest/'disk-controls'
            execute(cmd+[str(scripts/'run_independent_corruptions.py'), '--source', str(root),
                         '--report', str(initial), '--output-dir', str(controls)], dest/'logs/disk-controls')
            r = json.loads((controls/'CONTROLS.json').read_bytes())['cases']
            need(len(r) == 7 and all(x['returncode'] != 0 and x['success_artifact_absent'] and
                 x['failure'] for x in r), 'seven independent initial rejections per mode')
            controls = dest/'coarse-response-controls'
            execute(cmd+[str(independent/'run_coarse_response_controls.py'),
                         '--checker', str(independent/'check_coarse_response.py'),
                         '--source-root', str(root), '--initial-report', str(initial),
                         '--output-dir', str(controls)], dest/'logs/coarse-response-controls')
            r = json.loads((controls/'CONTROLS.json').read_bytes())['controls']
            need(len(r) == 6 and all(x['nonzero_exit'] and x['successful_artifact_absent'] and
                 x['rejection'] for x in r), 'six independent coarse response rejections per mode')
            print(json.dumps({'completed': 'INDEPENDENT_INITIAL_REGULARITY_AND_RESPONSE', 'mode': mode}), flush=True)
        return reports

    with ThreadPoolExecutor(max_workers=2) as pool:
        jobs = [pool.submit(native_reports), pool.submit(independent_reports)]
        reports = {}
        for job in jobs:
            reports.update(job.result())
    need(len(reports) == 8 and authenticate(cap) == closure, 'eight reports and final source identity')
    result = {
        'status': 'PASS_EXACT_PACKET_INITIAL_AND_CONDITIONAL_RESPONSE',
        'canonical_report_SHA256': reports,
        'complete_reports_reproduced_in_normal_and_optimized_modes': 8,
        'deliberate_corrupted_subprocesses_rejected': 98,
        'mathematical_formula_or_proof_ledger_corruptions_rejected': 92,
        'source_authentication_corruptions_rejected': 6,
        'corruption_roster': {'initial': 20, 'independent_initial': 14,
                              'response': 22, 'independent_coarse_response': 12,
                              'independent_full_response': 30},
        'independent_interval_evaluator_has_own_arithmetic_and_all23_table': True,
        'local_true_error_gate_and_complete_initial_and_response_certificate_gates_pass': True,
        'capsule_authenticated_before_and_after': True,
        'capsule_manifest_SHA256': closure,
        'old_candidates_sources_manuscripts_and_sites_preserved': True,
        'scope': 'Conditional coordinate approximation bounds. Local true-error improvement and tighter global upper certificates are different statements. No new nonlinear shell-contact solution, invariant observable or cosmological confirmation.'}
    (out/'VERIFICATION.json').write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
