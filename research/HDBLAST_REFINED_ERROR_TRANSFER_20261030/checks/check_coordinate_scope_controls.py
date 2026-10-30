#!/usr/bin/env python3
"""Execute the normalization checker with seven invalid formula/scope modes."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys


def need(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo-root', type=Path, required=True)
    ap.add_argument('--combined-report', type=Path, required=True)
    ap.add_argument('--combined-sha256', required=True)
    ap.add_argument('--c2-report', type=Path, required=True)
    ap.add_argument('--c2-sha256', required=True)
    ap.add_argument('--checker', type=Path, required=True)
    ap.add_argument('--checker-sha256', required=True)
    a = ap.parse_args()
    for p, digest in ((a.combined_report, a.combined_sha256), (a.c2_report, a.c2_sha256), (a.checker, a.checker_sha256)):
        need(hashlib.sha256(p.read_bytes()).hexdigest() == digest, 'authenticate guarded scope input')
    variants = (
        ('affine-error-denominator', 'denominator uses actual target derivative rather than affine approximation error'),
        ('wrong-Gram', 'actual scalar mean has unit Gram weight'),
        ('drop-trace', 'retained classical target initial rr trace is an explicit required premise'),
        ('pointwise-relative', 'only declared global coordinate rr and max-second seminorm ratios justified'),
        ('separate-xx-denominator', 'only declared global coordinate rr and max-second seminorm ratios justified'),
        ('full-C2-norm', 'only declared global coordinate rr and max-second seminorm ratios justified'),
        ('physical-accuracy', 'only declared global coordinate rr and max-second seminorm ratios justified'),
    )
    flags = [sys.executable, '-I', '-S', '-B']+(['-O'] if sys.flags.optimize else [])
    result = []
    for name, expected in variants:
        run = subprocess.run(flags+[str(a.checker), '--repo-root', str(a.repo_root),
            '--combined-report', str(a.combined_report), '--combined-sha256', a.combined_sha256,
            '--c2-report', str(a.c2_report), '--c2-sha256', a.c2_sha256,
            '--inject-error', name], capture_output=True, timeout=30)
        need(run.returncode != 0 and not run.stdout and expected in run.stderr.decode(errors='replace'),
             'invalid coordinate scope/formula was accepted: '+name)
        result.append({'mode': name, 'rejected': True, 'guard': expected})
    print(json.dumps({'status': 'PASS_COORDINATE_NORMALIZATION_SCOPE_CONTROLS',
        'input_sha256': {'combined': a.combined_sha256, 'C2': a.c2_sha256, 'checker': a.checker_sha256},
        'invalid_formula_or_scope_runs_rejected': result,
        'scope': 'Seven real invalid-formula/scope executions; no assertion of physical or pointwise relative accuracy.'}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
