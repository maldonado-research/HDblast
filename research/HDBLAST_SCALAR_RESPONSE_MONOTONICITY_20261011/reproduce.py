#!/usr/bin/env python3
"""Replay frozen proof arithmetic without modifying this capsule.

Outputs go to a new directory outside the capsule. Every frozen input is hashed
before execution. The producer, independent checker and symbolic verifier run
normally and with Python optimization enabled, with no imports of the producer
by its independent checker.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
FROZEN = {
    'theory/JACOBI_ENDPOINT_MONOTONICITY.md': 'f1e2571a551ec122dd7755a13145bae73bd92bd7f7426bdc3002e2f52f8dbfc8',
    'bounds/evaluate_rescaled_bounds.py': 'f427cbbc3b6bac37337fc723a7d27535ad18db426d3750d570efcbab5b231aed',
    'bounds/RESCALED_BOUND_CERTIFICATE.json': '22c6cd4e505ee2c7e549420355b94a5b3c079ebb6945314f5b285c13bb67835d',
    'bounds/public_baseline/EXACT_UNIFORM_CONSTANTS.json': '95923d9824f3f367385f67392ac10391b489337138851b642aa4553377d11ba1',
    'baseline/UNIFORM_BRANCH_THEOREM.md': '544c6dbced51fbc68dccab3bbbe84210545849c6e7447e1b6bdae0fc98e76f51',
    'baseline/REVIEW_ACCEPTANCE.json': '9540c183f321312cabb25e72ef5aab423bcb27dcd6e9dedc9725bc8eab0d0846',
    'baseline/independent-math/UNIFORM_PROOF_REVIEW.md': 'ea0fa00a0f93a2923d4e29081798def400d6c4412f720880601e2b4d0afc53e5',
    'baseline/independent-math/UNIFORM_PROOF_REVIEW.json': 'cd8296e247110ee7165a10f7dc5ec586143a47378a9d39ccca6c259d1f24d902',
    'baseline/independent-math/INDEPENDENT_EXACT_UNIFORM_CONSTANTS.json': 'db0fed7832ad259da38e37bc2627baa8c50805713a1a67bc9d18fbff81032ae0',
    'independent/check_rescaled_certificate.py': '975365e061f2dec196eb8e5661e1a46e82b4f98355ecbeb655e1abb3c180474a',
    'independent/verify_identities.py': '45ecc031a1263c06312ecc7b1967b2aeac52c04fc18b38d007984aa987e171a9',
    'independent/verify_identities.json': '271bdf1310580deffa929b0a4f5c4233afd6ffa14345e2228450a8b38cd7f323',
    'review/FINAL_REVIEW_RECEIPT.json': '955c28656536daae436cbe0b0393a76d7cdbe53082c1e6387e65d2889175c692',
    'review/INDEPENDENT_PORTABLE_EXACT_CHECK.json': 'ba28687b1667282835a6a099a4173078fb2c3b3d534019999b00f868ead128f8',
    'review/INDEPENDENT_PORTABLE_EXACT_CHECK_OPTIMIZED.json': 'ba28687b1667282835a6a099a4173078fb2c3b3d534019999b00f868ead128f8',
}
COMPLETE_MANIFEST_STATUS = 'COMPLETE_SELF_CONTAINED_CAPSULE'


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fail(message):
    raise RuntimeError(message)


def safe_relative(text):
    path = Path(text)
    if path.is_absolute() or '..' in path.parts or not path.parts:
        fail('Unsafe capsule manifest path: ' + text)
    return path


def files_under(directory):
    found = {}
    for path in sorted(directory.rglob('*')):
        if path.is_symlink():
            fail('Capsule replay does not accept symbolic links: ' + str(path.relative_to(directory)))
        if path.is_file():
            found[path.relative_to(directory).as_posix()] = sha256(path)
    return found


def verify_source(allow_pending):
    manifest_path = ROOT / 'MANIFEST.json'
    manifest = json.loads(manifest_path.read_text())
    status = manifest.get('status')
    if status != COMPLETE_MANIFEST_STATUS and not allow_pending:
        fail('Capsule manifest is not final; use --allow-pending-manifest only for an explicitly provisional replay.')
    reported_pins = manifest.get('frozen_pins', {})
    if reported_pins != FROZEN:
        fail('Manifest frozen input pins differ from the replay contract')
    for name, expected in FROZEN.items():
        path = ROOT / safe_relative(name)
        if not path.is_file() or sha256(path) != expected:
            fail('Frozen input hash mismatch: ' + name)
    inventoried = {}
    for entry in manifest.get('files', []):
        name = entry['path']
        path = ROOT / safe_relative(name)
        if name in inventoried or name == 'MANIFEST.json':
            fail('Invalid manifest inventory entry: ' + name)
        if not path.is_file() or path.stat().st_size != entry['bytes'] or sha256(path) != entry['sha256']:
            fail('Manifest file mismatch: ' + name)
        inventoried[name] = entry['sha256']
    snapshot = files_under(ROOT)
    if status == COMPLETE_MANIFEST_STATUS and set(inventoried) != set(snapshot) - {'MANIFEST.json'}:
        fail('Final manifest inventory does not cover exactly the capsule files')
    if not set(FROZEN) <= set(snapshot):
        fail('A required frozen input is missing')
    return manifest, snapshot


def run_command(command, working_directory, log, expected_failure=None, environment=None):
    result = subprocess.run(command, cwd=working_directory, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
                            env=environment)
    log.write_text(result.stdout + result.stderr)
    if expected_failure is None:
        if result.returncode != 0:
            fail(f'Command exited {result.returncode}; see {log.name}')
    elif result.returncode == 0 or expected_failure not in result.stderr:
        fail('Expected negative control was not rejected: ' + log.name)
    return result.returncode


def python_command(optimized, script, *arguments):
    return [sys.executable, '-B'] + (['-O'] if optimized else []) + [str(script), *map(str, arguments)]


def check_receipt(path, expected):
    if path.read_bytes() != expected.read_bytes():
        fail('Independent exact-check receipt differs from the archived receipt')
    report = json.loads(path.read_text())
    if report['status'] != 'PASS_INDEPENDENT_EXACT_ARITHMETIC':
        fail('Independent arithmetic checker did not pass')
    if len(report['reports']) != 8 or len(report['corrupt_certificate_rejections']) != 6:
        fail('Independent checker executed an unexpected test inventory')
    if report['reports']['primary_uniform_certificate']['compared_scalar_constants'] != 56:
        fail('Independent checker did not compare all primary constants')
    return report


def check_symbolic(path, expected):
    report = json.loads(path.read_text())
    archived = json.loads(expected.read_text())
    if report['all_passed'] is not True or report['check_count'] != 36:
        fail('The 36 symbolic identity and inequality checks did not all pass')
    if any(item['passed'] is not True for item in report['checks'].values()):
        fail('A symbolic check failed')
    for key in ('checks', 'assumptions', 'derived_constants', 'source_sha256', 'sympy_version'):
        if report[key] != archived[key]:
            fail('Symbolic evidence differs from the frozen report: ' + key)
    return report


def replay_mode(optimized, work, output):
    name = 'optimized' if optimized else 'normal'
    stage = output / name
    stage.mkdir()
    proof = work / 'theory/JACOBI_ENDPOINT_MONOTONICITY.md'
    producer = work / 'bounds/evaluate_rescaled_bounds.py'
    checker = work / 'independent/check_rescaled_certificate.py'
    archived_certificate = work / 'bounds/RESCALED_BOUND_CERTIFICATE.json'
    certificate = stage / 'RESCALED_BOUND_CERTIFICATE.json'
    receipt = stage / 'INDEPENDENT_EXACT_CHECK.json'
    run_command(python_command(optimized, producer, '--proof', proof,
                    '--proof-sha256', FROZEN['theory/JACOBI_ENDPOINT_MONOTONICITY.md'],
                    '--out', certificate), work, stage / 'producer.log')
    if certificate.read_bytes() != archived_certificate.read_bytes():
        fail('Replayed producer certificate differs from the frozen certificate in ' + name)
    run_command(python_command(optimized, checker, '--certificate', certificate,
                    '--producer', producer, '--proof', proof, '--require-proof-pin',
                    '--output', receipt), work, stage / 'independent_checker.log')
    check_receipt(receipt, work / 'review/INDEPENDENT_PORTABLE_EXACT_CHECK.json')

    symbolic_stage = stage / 'symbolic'
    symbolic_stage.mkdir()
    symbolic_script = symbolic_stage / 'verify_identities.py'
    shutil.copyfile(work / 'independent/verify_identities.py', symbolic_script)
    run_command(python_command(optimized, symbolic_script), symbolic_stage,
                stage / 'symbolic_verifier.log')
    symbolic_report = symbolic_script.with_suffix('.json')
    check_symbolic(symbolic_report, work / 'independent/verify_identities.json')

    wrong_destination = stage / 'WRONG_PROOF_MUST_NOT_CREATE_CERTIFICATE.json'
    run_command(python_command(optimized, producer, '--proof', proof,
                    '--proof-sha256', '0' * 64, '--out', wrong_destination), work,
                stage / 'wrong_proof_producer.log',
                expected_failure='Runtime analytic proof SHA256 pin mismatch')
    if wrong_destination.exists():
        fail('Producer wrote a certificate after receiving a wrong proof hash')
    wrong_proof = stage / 'ALTERED_PROOF_FOR_NEGATIVE_CONTROL.md'
    wrong_proof.write_bytes(proof.read_bytes() + b'\nNegative-control alteration.\n')
    wrong_receipt = stage / 'WRONG_PROOF_MUST_NOT_CREATE_RECEIPT.json'
    run_command(python_command(optimized, checker, '--certificate', certificate,
                    '--producer', producer, '--proof', wrong_proof, '--require-proof-pin',
                    '--output', wrong_receipt), work, stage / 'wrong_proof_checker.log',
                expected_failure='Certificate analytic proof digest mismatch')
    if wrong_receipt.exists():
        fail('Independent checker wrote a receipt for an altered proof')
    return dict(mode=name, certificate_exactly_matches_archive=True,
                independent_receipt_exactly_matches_archive=True,
                independent_certificate_variants=8, rational_constants_per_variant=56,
                symbolic_polynomials_per_variant=7, corrupt_certificate_rejections=6,
                symbolic_identity_checks=36, wrong_proof_producer_rejected_without_output=True,
                wrong_proof_checker_rejected_without_output=True,
                certificate_sha256=sha256(certificate), independent_receipt_sha256=sha256(receipt),
                symbolic_report_sha256=sha256(symbolic_report))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', '--output', dest='output', type=Path, required=True,
                        help='New external directory for the copied capsule and generated evidence.')
    parser.add_argument('--allow-pending-manifest', action='store_true',
                        help='Permit a provisional replay before the entire publication capsule is frozen.')
    parser.add_argument('--figures', action='store_true',
                        help='Also regenerate the optional plot; requires Matplotlib in addition to the proof dependencies.')
    args = parser.parse_args()
    destination = args.output.resolve()
    if destination == ROOT or ROOT in destination.parents:
        fail('Replay output must be outside the source capsule')
    if destination.exists():
        fail('Replay output directory already exists; choose a new directory')
    try:
        import sympy
        import mpmath
    except ImportError as error:
        fail('Install the pinned requirements.txt dependencies before replay: ' + str(error))
    if sympy.__version__ != '1.14.0' or mpmath.__version__ != '1.3.0':
        fail('Replay requires sympy==1.14.0 and mpmath==1.3.0')
    matplotlib_version = None
    if args.figures:
        try:
            matplotlib_version = importlib.metadata.version('matplotlib')
        except importlib.metadata.PackageNotFoundError:
            fail('Optional --figures replay requires Matplotlib; omit --figures to run every proof check without plotting dependencies.')
    manifest, before = verify_source(args.allow_pending_manifest)
    destination.mkdir(parents=True)
    work = destination / 'capsule'
    shutil.copytree(ROOT, work, symlinks=False)
    if files_under(work) != before:
        fail('External capsule copy does not match the verified source snapshot')
    modes = [replay_mode(False, work, destination), replay_mode(True, work, destination)]
    normal, optimized = destination / 'normal', destination / 'optimized'
    for name in ('RESCALED_BOUND_CERTIFICATE.json', 'INDEPENDENT_EXACT_CHECK.json',
                 'symbolic/verify_identities.json'):
        if (normal / name).read_bytes() != (optimized / name).read_bytes():
            fail('Normal and optimized Python produced different evidence: ' + name)
    figure_result = dict(requested=args.figures, affects_mathematical_acceptance=False)
    if args.figures:
        figure_output = destination / 'figures'
        environment = dict(os.environ)
        environment['MPLCONFIGDIR'] = str(destination / 'figure-config')
        run_command(python_command(False, work / 'figures/make_figure.py',
                                   '--output', figure_output), work,
                    destination / 'figure_regeneration.log', environment=environment)
        for name in ('controlled-response.csv', 'figure-data.json'):
            if (figure_output / name).read_bytes() != (work / 'figures' / name).read_bytes():
                fail('Optional figure mathematical data differ from the archived data: ' + name)
        for name, signature in [('controlled-response.png', b'\x89PNG\r\n\x1a\n'),
                                ('controlled-response.pdf', b'%PDF-')]:
            data = (figure_output / name).read_bytes()
            if len(data) < 1000 or not data.startswith(signature):
                fail('Optional figure output is missing or malformed: ' + name)
        figure_result.update(matplotlib_version=matplotlib_version,
                             mathematical_data_exactly_match_archive=True,
                             generated_files={name: sha256(figure_output / name)
                                              for name in ('controlled-response.csv', 'figure-data.json',
                                                           'controlled-response.png', 'controlled-response.pdf')})
    if files_under(ROOT) != before:
        fail('Source capsule changed during replay')
    if files_under(work) != before:
        fail('The externally copied frozen capsule changed during replay')
    result = dict(status='PASS_SELF_CONTAINED_EXTERNAL_REPLAY',
                  capsule_manifest_status=manifest['status'],
                  provisional=manifest['status'] != COMPLETE_MANIFEST_STATUS,
                  source_capsule_unchanged=True, copied_frozen_capsule_unchanged=True,
                  normal_and_optimized_evidence_identical=True,
                  optional_figure_regeneration=figure_result,
                  runtime_versions=dict(python=platform.python_version(), sympy=sympy.__version__,
                                        mpmath=mpmath.__version__),
                  frozen_input_pins=FROZEN, modes=modes,
                  scope='Exact mathematical branch bounds and independent symbolic checks. No ODE trajectory, stability, cosmological-origin or external peer-review claim.')
    receipt = destination / 'REPRODUCTION_RECEIPT.json'
    receipt.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(status=result['status'], provisional=result['provisional'],
                         modes=2, symbolic_checks_per_mode=36, certificate_variants_per_mode=8,
                         corrupt_certificate_rejections_per_mode=6, receipt=str(receipt)), indent=2))


if __name__ == '__main__':
    main()
