#!/usr/bin/env python3
"""Replay pinned finite algebra and background identities without network or Git."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

PACKAGE = Path(__file__).resolve().parent
JOBS = {'matrix': 'finite_commuting_geometry.py',
        'rotation': 'check_rotational_constraints.py',
        'quadrupole': 'check_quadrupole_admissibility.py',
        'background': 'check_explicit_background.py'}
INJECTED_JOBS = ('rotation', 'quadrupole', 'background')
MODES = [('normal', []), ('optimized', ['-O'])]
EXPECTED_PATHS = {
    'code/finite_commuting_geometry.py', 'code/check_rotational_constraints.py',
    'code/check_quadrupole_admissibility.py', 'code/check_explicit_background.py',
    'code/ROTATIONAL_SECOND_ORDER_SHELL_OBSTRUCTION.md',
    'expected/external-controls.json',
} | {'expected/' + job + '.' + mode + '.json' for job in JOBS for mode, _ in MODES}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(path):
    data = path.read_bytes()
    return {'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


def snapshot():
    result = {}
    for path in sorted(PACKAGE.rglob('*')):
        require(not path.is_symlink(), 'Package contains a symbolic link: ' + str(path))
        if path.is_file():
            result[path.relative_to(PACKAGE).as_posix()] = digest(path)
    return result


def read_json(path):
    with path.open(encoding='utf-8') as stream:
        return json.load(stream)


def write_json(path, value):
    with path.open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write('\n')


def verify_pins(pins):
    require(pins['schema'] == 'focused-finite-algebra-input-pins-v1', 'Unsupported input pins')
    require(set(pins['files']) == EXPECTED_PATHS, 'Pinned input inventory differs')
    for relative, expected in pins['files'].items():
        path = PACKAGE / relative
        require(path.is_file() and not path.is_symlink(), 'Missing regular input: ' + relative)
        require(digest(path) == {key: expected[key] for key in ('bytes', 'sha256')},
                'Pinned bytes/hash mismatch: ' + relative)
    for subdirectory in ('code', 'expected'):
        actual = {p.relative_to(PACKAGE).as_posix() for p in (PACKAGE / subdirectory).rglob('*') if p.is_file()}
        require(actual == {p for p in EXPECTED_PATHS if p.startswith(subdirectory + '/')},
                'Unexpected files in ' + subdirectory)
    require(pins['mode_comparison_ignored_fields'] ==
            {'matrix': ['python_optimization'], 'rotation': [], 'quadrupole': [], 'background': []},
            'Undeclared mode comparison exclusions')
    require(platform.python_version() == pins['environment']['python'],
            'Python version mismatch: require ' + pins['environment']['python'])
    for package in ('sympy', 'mpmath'):
        require(importlib.metadata.version(package) == pins['environment'][package],
                package + ' version mismatch: require ' + pins['environment'][package])


def classify(name, report):
    if name == 'matrix':
        items = report['checks']
        require(report['all_pass'] is True and all(item['pass'] is True for item in items), 'Matrix check did not pass')
        require(report['check_count'] == len(items), 'Matrix count differs from its actual list')
        wrong = sum(item.get('details', {}).get('kind') == 'wrong-formula control' for item in items)
        return {'positive_checks': len(items) - wrong, 'embedded_wrong_formula_controls': wrong,
                'total_generated_checks_including_controls': len(items)}
    if name == 'rotation':
        items, controls = report['checks'], report['wrong_controls']
        require(report['scientific_status'] == 'exact_algebra_controls_passed', 'Rotational check status differs')
        require(all(item['passed'] is True for item in items), 'Rotational positive check did not pass')
        require(all(item['detected'] is True for item in controls), 'Rotational control did not detect its wrong formula')
        require(report['counts'] == {'checks': len(items), 'wrong_controls': len(controls)}, 'Rotational count differs from actual list')
    elif name == 'quadrupole':
        items, controls = report['checks'], report['controls']
        require(report['status'] == 'PASS_EXACT_QUADRUPOLE_SYMMETRY_CHECKS', 'Quadrupole check status differs')
        require(all(item['residual'] == '0' for item in items), 'Quadrupole positive residual differs')
        require(all(item['incorrect_formula_rejected'] is True for item in controls), 'Quadrupole control did not reject wrong formula')
        require(report['positive_checks'] == len(items) and report['wrong_formula_controls'] == len(controls),
                'Quadrupole count differs from actual lists')
    else:
        require(name == 'background', 'Unknown checker name')
        items, controls = report['checks'], report['controls']
        require(report['status'] == 'PASS_EXACT_BACKGROUND_CHECKS', 'Background check status differs')
        require(all(item['pass'] is True for item in items), 'Background positive check did not pass')
        require(all(item['incorrect_formula_rejected'] is True for item in controls), 'Background control did not reject wrong formula')
        require(report['positive_checks'] == len(items) and report['wrong_formula_controls'] == len(controls),
                'Background count differs from actual lists')
    return {'positive_checks': len(items), 'embedded_wrong_formula_controls': len(controls)}


def invoke(output, key, name, flags, inject=False):
    result_path = output / (key + '.json')
    require(not result_path.exists(), 'Result output already exists: ' + key)
    command = [sys.executable, '-B', *flags, 'code/' + JOBS[name], '--output', str(result_path)]
    if inject:
        command.append('--inject-failure')
    child_environment = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    proc = subprocess.run(command, cwd=PACKAGE, capture_output=True, text=True,
                          env=child_environment, timeout=600)
    for suffix, data in [('stdout', proc.stdout), ('stderr', proc.stderr)]:
        with (output / (key + '.' + suffix + '.txt')).open('x', encoding='utf-8') as stream:
            stream.write(data)
    report = read_json(result_path) if result_path.is_file() else None
    return proc, report, {'job': key, 'kind': 'external_injection' if inject else 'positive',
                          'command': command, 'exit_code': proc.returncode,
                          'report': result_path.name if result_path.exists() else None,
                          'stdout': key + '.stdout.txt', 'stderr': key + '.stderr.txt'}


def replay(output, pins, records):
    reports, counts = {}, {}
    controls = read_json(PACKAGE / 'expected/external-controls.json')
    require(set(controls) == {name + '.injected.' + mode for name in INJECTED_JOBS for mode, _ in MODES},
            'External failure contract inventory differs')
    for mode, flags in MODES:
        counts[mode] = {}
        for name in JOBS:
            key = name + '.' + mode
            proc, actual, record = invoke(output, key, name, flags)
            records.append(record)
            require(proc.returncode == 0, key + ' exited ' + str(proc.returncode))
            require(actual is not None, key + ' produced no report')
            expected = read_json(PACKAGE / ('expected/' + key + '.json'))
            require(actual == expected, key + ' differs from its complete packaged expected report')
            require(proc.stderr == '', key + ' produced unexpected stderr')
            if name == 'matrix':
                require(actual['python_optimization'] == (1 if flags else 0), 'Incorrect matrix optimization metadata')
                require(json.loads(proc.stdout) == {k: v for k, v in actual.items() if k != 'checks'}, 'Matrix stdout differs from report')
            elif name == 'rotation':
                require(proc.stdout == '', 'Unexpected rotational stdout')
            else:
                require(proc.stdout == 'PASS ' + str(len(actual['checks'])) + ' exact identities; ' + str(len(actual['controls'])) + ' wrong-formula controls\n',
                        'Unexpected ' + name + ' stdout')
            counts[mode][name] = classify(name, actual)
            reports[key] = actual
            record['full_expected_report_agreement'] = True

        for name in INJECTED_JOBS:
            key = name + '.injected.' + mode
            proc, report, record = invoke(output, key, name, flags, inject=True)
            records.append(record)
            if name in ('rotation', 'background'):
                require(proc.stderr == '', key + ' produced unexpected stderr')
                reason = report.get('error') if isinstance(report, dict) else None
            else:
                require(proc.stderr.startswith('Traceback (most recent call last):\n'), key + ' has no Python traceback')
                reason = proc.stderr.splitlines()[-1] if proc.stderr else ''
            actual = {'exit_code': proc.returncode, 'report_exists': (output / (key + '.json')).exists(),
                      'report': report, 'failure_reason': reason, 'stdout': proc.stdout}
            require(actual == controls[key], key + ' differs from its packaged expected failure contract')
            record['failure_reason'] = reason
            record['expected_failure_contract_agreement'] = True

        counts[mode]['total_positive_checks'] = sum(counts[mode][name]['positive_checks'] for name in JOBS)
        counts[mode]['total_embedded_wrong_formula_controls'] = sum(counts[mode][name]['embedded_wrong_formula_controls'] for name in JOBS)
        require(counts[mode] == pins['observed_counts_per_mode'], 'Actual counts differ from captured expected inventory: ' + mode)

    for name in JOBS:
        normal, optimized = dict(reports[name + '.normal']), dict(reports[name + '.optimized'])
        for field in pins['mode_comparison_ignored_fields'][name]:
            del normal[field]
            del optimized[field]
        require(normal == optimized, name + ' full cross-mode report comparison failed')
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, required=True, help='A new directory outside this supplement; must not exist')
    args = parser.parse_args()
    output = args.output_dir.resolve()
    require(output != PACKAGE and PACKAGE not in output.parents, 'Output directory must be outside the supplement')
    require(not output.exists(), 'Output directory must be fresh; refusing overwrite: ' + str(output))
    output.mkdir(parents=True)
    records = []
    receipt = {'schema': 'focused-finite-algebra-replay-v1', 'status': 'failed',
               'scope': 'Exact finite algebra and explicit-background identities; no general PDE existence or novelty assessment.',
               'records': records, 'python': platform.python_version(), 'runner_python_optimization': sys.flags.optimize}
    before = None
    try:
        before = snapshot()
        receipt['package_files_before'] = before
        pins = read_json(PACKAGE / 'INPUT_PINS.json')
        verify_pins(pins)
        receipt['counts_per_mode'] = replay(output, pins, records)
        positive_processes = len(JOBS) * len(MODES)
        external_injection_processes = len(INJECTED_JOBS) * len(MODES)
        require(len([r for r in records if r['kind'] == 'positive']) == positive_processes, 'Incomplete positive process inventory')
        require(len([r for r in records if r['kind'] == 'external_injection']) == external_injection_processes, 'Incomplete injection process inventory')
        after = snapshot()
        receipt['package_files_after'] = after
        require(before == after, 'Supplement changed during replay')
        verify_pins(pins)
        receipt.update(status='PASS_EXACT_AGREEMENT', positive_processes=positive_processes,
                       external_injection_processes=external_injection_processes,
                       all_full_expected_reports_agree=True, all_external_failure_contracts_agree=True,
                       all_cross_mode_reports_agree=True, unchanged_package=True,
                       mode_comparison_ignored_fields=pins['mode_comparison_ignored_fields'])
        exit_code = 0
    except Exception as error:
        receipt['error'] = type(error).__name__ + ': ' + str(error)
        if before is not None:
            try:
                receipt['package_files_after'] = snapshot()
                receipt['unchanged_package'] = before == receipt['package_files_after']
            except Exception as snapshot_error:
                receipt['snapshot_error'] = str(snapshot_error)
        exit_code = 1
    write_json(output / 'REPLAY_RECEIPT.json', receipt)
    print(json.dumps({'status': receipt['status'], 'receipt': str(output / 'REPLAY_RECEIPT.json'),
                      'error': receipt.get('error')}, sort_keys=True))
    return exit_code


if __name__ == '__main__':
    raise SystemExit(main())
