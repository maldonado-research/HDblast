#!/usr/bin/env python3
"""Portable fresh-output CLI for the exact source-free curvature fixture."""
import argparse
import json
from pathlib import Path
import subprocess
import sys


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo-root', required=True, type=Path,
                   help='accepted for uniform CLI; the fixture reads no production sources')
    p.add_argument('--output', required=True, type=Path)
    a = p.parse_args()
    here = Path(__file__).resolve().parent
    if a.output.exists() or a.output.is_symlink() or here in a.output.resolve().parents:
        raise ValueError('FRESH_OUTPUT_OUTSIDE_CODE')
    argv = [sys.executable, '-I', '-S', '-B']+(['-O'] if sys.flags.optimize else [])
    argv += [str(here/'radial-review/check_direct_radial.py'), '--combined']
    result = subprocess.run(argv, capture_output=True)
    if result.returncode:
        raise ValueError('EXACT_DIRECT_FIXTURE_FAILURE:'+result.stderr.decode(errors='replace')[-2000:])
    report = json.loads(result.stdout)
    if report['status'] != 'PASS_RADIAL_DIRECT_COMBINED_CURVATURE_WITNESS':
        raise ValueError('EXACT_DIRECT_FIXTURE_STATUS')
    with a.output.open('xb') as handle:
        handle.write(result.stdout)
    print(json.dumps({'status': report['status']}, sort_keys=True))


if __name__ == '__main__':
    main()
