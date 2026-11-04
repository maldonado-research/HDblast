#!/usr/bin/env python3
"""Real subprocess controls reaggregate the fixed raw witness without oracles."""
import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        print(json.dumps({'status': 'REJECTED', 'witness': 'OUTPUT_ALREADY_EXISTS'}))
        return 2
    spec = importlib.util.spec_from_file_location('portable_preparation_driver', HERE/'reproduce.py')
    driver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    try:
        authentication = driver.authenticate(args.repo_root.resolve(strict=True))
        cases = []
        with tempfile.TemporaryDirectory(prefix='preparation-controls-') as tmp:
            for name, field in driver.INJECTIONS.items():
                out = Path(tmp)/name
                command = [sys.executable, '-I', '-S', '-B']
                if sys.flags.optimize:
                    command.append('-O')
                command += [str(HERE/'reproduce.py'), '--repo-root', str(args.repo_root),
                            '--output-dir', str(out), '--inject', name,
                            '--use-frozen-raw-for-control']
                run = subprocess.run(command, capture_output=True, text=True, check=False)
                expected = {'status': 'REJECTED', 'witness': 'FROZEN_SCIENCE_MISMATCH:'+field}
                try:
                    parsed = json.loads(run.stdout)
                except ValueError:
                    parsed = None
                driver.require(run.returncode == 2 and run.stderr == '' and parsed == expected
                               and len(run.stdout.splitlines()) == 1 and not out.exists(),
                               'NEGATIVE_CONTROL_PROTOCOL:'+name)
                cases.append({'case': name, 'expected': expected, 'exit_code': 2,
                              'stdout_single_JSON': True, 'stderr_empty': True,
                              'no_artifact': True,
                              'child_required_flags': ['-I', '-S', '-B'],
                              'optimized_mode_matches_parent': ('-O' in command) == bool(sys.flags.optimize),
                              'execution': 'Real CLI subprocess, authenticated raw reaggregation only; no new bump/background/first-field calls.'})
        result = {'status': 'PASS_PORTABLE_REAL_NEGATIVE_CONTROLS', 'cases': cases,
                  'case_count': len(cases),
                  'authentication_pins': {k: authentication[k] for k in
                                          ('local_pins_sha256', 'source_pins_sha256')},
                  'new_oracle_calls': {'background': 0, 'evolved_actual_first_field': 0},
                  'positive_scientific_campaign_rerun': False}
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
        print(json.dumps({'status': 'PASS', 'witness': 'REAL_NEGATIVE_CONTROLS',
                          'case_count': len(cases)}, sort_keys=True))
        return 0
    except (driver.Reject, OSError, ValueError, KeyError) as error:
        witness = error.witness if isinstance(error, driver.Reject) else 'CONTROL_INPUT_FAILURE'
        print(json.dumps({'status': 'REJECTED', 'witness': witness}, sort_keys=True))
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
