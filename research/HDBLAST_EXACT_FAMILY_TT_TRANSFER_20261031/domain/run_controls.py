#!/usr/bin/env python3
"""Run normal/-O checks and real copied-source mutation guards privately."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--pins-sha256', required=True)
    parser.add_argument('--pins', type=Path, default=Path(__file__).with_name('SOURCE_PINS.json'))
    args = parser.parse_args()
    source_root = Path(__file__).resolve().parent
    root = args.output_dir.resolve()
    if root.exists() and any(root.iterdir()):
        raise RuntimeError('Use a new or empty output directory')
    root.mkdir(parents=True, exist_ok=True)
    repo_root = args.repo_root.resolve(strict=True)
    if root.is_relative_to(repo_root):
        raise RuntimeError('Output directory must be outside the source checkout')
    pins_path = args.pins.resolve(strict=True)
    pins = json.loads(pins_path.read_text())
    digest = hashlib.sha256(pins_path.read_bytes()).hexdigest()
    if digest != args.pins_sha256:
        raise RuntimeError('Externally supplied source manifest digest differs')
    control_root = root / 'controls'
    control_root.mkdir(exist_ok=True)
    mutations = {
        'packet_velocity': ('packet', b'Ce^{-z}q_\\ell', b'Cq_\\ell'),
        'actual_first_domain': ('first_actual_domain', b'x+r<=8', b'x+r<=7'),
        'actual_background_rho': ('background_bounds', b'"131"', b'"128"'),
        'strict_theorem_margin': ('full_theorem', b'Z>L+c_*T', b'Z>=L+c_*T'),
        'zero_fixing_boundary_flow': ('boundary_flow', b'It fixes the zero perturbation.', b'It moves the zero perturbation.'),
    }
    prepared = {}
    for name, (source, before, after) in mutations.items():
        relative = pins['sources'][source]['path']
        original_path = (repo_root / relative).resolve(strict=True)
        if not original_path.is_relative_to(repo_root):
            raise RuntimeError('Source path escaped checkout')
        original = original_path.read_bytes()
        if before not in original:
            raise RuntimeError('Mutation anchor missing: ' + name)
        changed = original.replace(before, after, 1)
        dest = control_root / (name + Path(relative).suffix)
        dest.write_bytes(changed)
        overrides = control_root / (name + '.overrides.json')
        overrides.write_text(json.dumps({source: str(dest)}, indent=2) + '\n')
        prepared[name] = overrides
    expected_errors = {
        'causal_equality': 'STRICT_FULL_CORRECTION_CAUSAL_MARGIN',
        'packet_touches_annulus': 'STRICT_PACKET_IN_ANNULUS',
        'collar_touches_annulus': 'FIXED_COLLAR_IN_ZERO_INITIAL_REGION',
        'target_outside_slab': 'TARGET_INTERIOR',
        'first_domain_escape': 'SLAB_IN_ACTUAL_FIRST_TRIANGLE',
    }
    records = []
    for mode, flags in [('NORMAL', []), ('OPTIMIZED', ['-O'])]:
        base = [sys.executable, '-I', '-B', *flags, str(source_root / 'check_domain_seed.py'),
                '--pins', str(pins_path), '--pins-sha256', digest, '--repo-root', str(repo_root)]
        cp = subprocess.run(base, capture_output=True, text=True)
        if cp.returncode != 0:
            raise RuntimeError('Nominal failure: ' + cp.stderr)
        (root / (mode + '.json')).write_text(cp.stdout)
        records.append({'mode': mode, 'case': 'nominal', 'exit': cp.returncode,
                        'stdout_sha256': hashlib.sha256(cp.stdout.encode()).hexdigest()})
        for name, overrides in prepared.items():
            cp = subprocess.run(base + ['--overrides', str(overrides)], capture_output=True, text=True)
            if cp.returncode != 1 or 'SOURCE_AUTHENTICATION_FAILED:' not in cp.stderr or cp.stdout:
                raise RuntimeError('Mutation guard failed: ' + name + ': ' + cp.stderr)
            records.append({'mode': mode, 'case': name, 'kind': 'actual_source_mutation',
                            'exit': cp.returncode, 'rejected_before_algebra': True,
                            'no_success_stdout': not cp.stdout,
                            'success_artifact_absent': not cp.stdout,
                            'reason': cp.stderr.strip()})
        for name, expected in expected_errors.items():
            cp = subprocess.run(base + ['--case', name], capture_output=True, text=True)
            if cp.returncode != 1 or cp.stderr.strip() != expected or cp.stdout:
                raise RuntimeError('Inequality guard failed: ' + name + ': ' + cp.stderr)
            records.append({'mode': mode, 'case': name, 'kind': 'strict_inequality_control',
                            'exit': cp.returncode, 'no_success_stdout': not cp.stdout,
                            'success_artifact_absent': not cp.stdout, 'reason': cp.stderr.strip()})
    if (root / 'NORMAL.json').read_bytes() != (root / 'OPTIMIZED.json').read_bytes():
        raise RuntimeError('Normal and optimized receipts differ')
    report = {'status': 'PASS_DOMAIN_SEED_CONTROLS', 'source_pins_sha256': digest,
              'normal_optimized_scientific_bytes_equal': True,
              'actual_source_mutations_rejected': 10,
              'strict_inequality_cases_rejected': 10,
              'symbolic_wrong_seed_formulas_rejected': 6,
              'success_artifact_kind': 'Checker success is stdout JSON only; rejected subprocesses emit none.',
              'runs': records,
              'scope': 'Arithmetic, seed algebra, and source authentication controls; not PDE proof certification.'}
    (root / 'CONTROL_REPORT.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: value for key, value in report.items() if key != 'runs'}, indent=2))


if __name__ == '__main__':
    main()
