#!/usr/bin/env python3
"""Fresh normal/optimized exact review plus meaningful corrupted executions."""
import argparse
import ast
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
CHECKER_PIN = '6c1384d380127c89149e5752da839dd280ad507f3cf9550de06d686f0cbe7cef'
WAVE_REPORT_PIN = '226117368a98a7a01255547147fcbf56d9019e1e2d29880136bd5592f8a3f549'
WAVE_HELPER_PIN = 'cb025648deda6eec4cfd77b3bb96c6cf6e0713f2b55987b3fdfff22dbb066bbf'


def need(ok, why):
    if not ok:
        raise ValueError(why)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def compare_wave_rows(independent, wave):
    own = independent['independent_arbitrary_covector_operator']['independent_eight_modal_wave_rows']
    other = wave['modal_rows']
    checks = []
    for field in ('y', '0', 'a'):
        need(len(own[field]) == len(other[field]), 'WAVE_ROW_COUNT:' + field)
        for index, (mine, theirs) in enumerate(zip(own[field], other[field])):
            left = {}
            for row in mine:
                powers = {}
                for name, exponent in row['powers'].items():
                    if name.startswith('z:'):
                        sector, degree, suffix = name[2:].split(':')
                        name = {'radial': 'f', 'temporal': 'g', 'angular': 'h'}[sector] + degree + ('_' + suffix if suffix else '')
                    powers[name] = exponent
                left[tuple(sorted(powers.items()))] = Fraction(row['coefficient'])
            right = {tuple((n, e) for n, e in pairs): Fraction(coefficient) for pairs, coefficient in theirs}
            need(left == right, 'INDEPENDENT_WAVE_ROWS_DIFFER:' + field + ':' + str(index))
            checks.append({'field': field, 'index': index, 'collected_terms': len(left), 'residual_terms': 0})
    return checks


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', type=Path, required=True)
    ap.add_argument('--wave-dir', type=Path, required=True)
    ap.add_argument('--output-dir', type=Path, required=True)
    args = ap.parse_args()
    checker = HERE / 'check_independent_source_roles.py'
    need(sha(checker) == CHECKER_PIN, 'CHECKER_PIN')
    runner = Path(__file__).resolve()
    before = {'checker': sha(checker), 'runner': sha(runner)}
    tree = ast.parse(checker.read_text())
    controls = None
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'CONTROLS' for t in node.targets):
            controls = ast.literal_eval(node.value)
    need(isinstance(controls, tuple) and len(controls) == 20, 'CONTROL_ROSTER')
    for path in (checker, runner):
        need(not any(isinstance(n, ast.Assert) for n in ast.walk(ast.parse(path.read_text()))), 'ASSERT:' + path.name)
    output = args.output_dir.resolve()
    need(not output.exists(), 'OUTPUT_DIRECTORY_MUST_BE_ABSENT')
    output.mkdir(parents=True)
    root = args.root.resolve()
    jobs = []
    for optimized in (False, True):
        mode = 'OPTIMIZED' if optimized else 'NORMAL'
        for mutation in (None,) + controls:
            destination = output / ('SOURCE_ROLES_' + mode + '.json' if mutation is None else 'REJECT_' + mode + '_' + mutation + '.json')
            command = [sys.executable, '-I', '-S', '-B'] + (['-O'] if optimized else [])
            command += [str(checker), '--root', str(root), '--output', str(destination)]
            if mutation is not None:
                command += ['--mutation', mutation]
            jobs.append((mode, mutation, destination, command))

    def execute(job):
        mode, mutation, destination, command = job
        need(not destination.exists(), 'SUBPROCESS_OUTPUT_MUST_BE_ABSENT')
        process = subprocess.run(command, capture_output=True, text=True, timeout=180)
        record = {'mode': mode, 'mutation': mutation, 'command': command,
                  'returncode': process.returncode, 'stdout': process.stdout.strip(),
                  'successful_output_exists': destination.exists()}
        if mutation is None:
            need(process.returncode == 0 and destination.is_file(), 'POSITIVE_FAILED:' + mode + ':' + process.stderr)
            record['output_sha256'] = sha(destination)
        else:
            last = process.stderr.strip().splitlines()[-1] if process.stderr.strip() else ''
            need(process.returncode != 0 and not destination.exists() and ':NONZERO:' in last,
                 'CONTROL_NOT_MATHEMATICALLY_REJECTED:' + str(mutation) + ':' + process.stderr)
            record['exact_nonzero_collected_witness'] = last.removeprefix('ValueError: ')
        return record

    # Independent read-only subprocesses have distinct absent output paths.
    with ThreadPoolExecutor(max_workers=2) as pool:
        records = list(pool.map(execute, jobs))
    normal = output / 'SOURCE_ROLES_NORMAL.json'
    optimized = output / 'SOURCE_ROLES_OPTIMIZED.json'
    need(normal.read_bytes() == optimized.read_bytes(), 'NORMAL_OPTIMIZED_REPORTS_DIFFER')
    wave_dir = args.wave_dir.resolve()
    wave_normal = wave_dir / 'WAVE_NORMAL.json'
    wave_optimized = wave_dir / 'WAVE_OPTIMIZED.json'
    wave_helper = wave_dir / 'covector_wave_rows.py'
    need(sha(wave_normal) == WAVE_REPORT_PIN and sha(wave_optimized) == WAVE_REPORT_PIN, 'WAVE_REPORT_PIN')
    need(sha(wave_helper) == WAVE_HELPER_PIN, 'WAVE_HELPER_PIN')
    row_checks = compare_wave_rows(json.loads(normal.read_text()), json.loads(wave_normal.read_text()))
    need(before == {'checker': sha(checker), 'runner': sha(runner)}, 'REVIEW_SOURCE_CHANGED')
    need(sha(wave_normal) == WAVE_REPORT_PIN and sha(wave_optimized) == WAVE_REPORT_PIN and sha(wave_helper) == WAVE_HELPER_PIN, 'WAVE_DATA_CHANGED')
    receipt = {'status': 'PASS_INDEPENDENT_SOURCE_ROLES_AND_WAVE_REVIEW',
               'review_source_authentication': before,
               'canonical_normal_optimized_reports_identical': True,
               'positive_subprocess_count': 2, 'rejected_corrupted_subprocess_count': 40,
               'executions': records,
               'independent_eight_wave_row_data_comparisons': row_checks,
               'wave_data_authentication': {'report_sha256': WAVE_REPORT_PIN, 'helper_sha256': WAVE_HELPER_PIN},
               'wave_comparison_independence': 'Producer wave JSON is read as data after own frame/quotient arithmetic computes every coefficient. Producer scientific code is never imported.'}
    path = output / 'RUNTIME_RECEIPT.json'
    with path.open('x') as handle:
        json.dump(receipt, handle, sort_keys=True, indent=2)
        handle.write('\n')
    print(json.dumps({'status': receipt['status'], 'positive': 2, 'rejected': 40,
                      'canonical_report_sha256': sha(normal), 'receipt_sha256': sha(path)}, sort_keys=True))


if __name__ == '__main__':
    main()
