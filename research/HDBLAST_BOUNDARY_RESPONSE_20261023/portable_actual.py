#!/usr/bin/env python3
"""Relocate only authenticated input roots of the unmodified actual producer."""
import argparse
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
PRODUCER_SHA256 = '6fd86cd06582c7182ffd82ad26eac19f22087a8a234ae8d19e2eae5fc7e8b5bc'


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ValueError('MODULE_LOADER:' + name)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    if not (sys.flags.isolated and sys.flags.no_site and sys.flags.dont_write_bytecode):
        raise ValueError('USE_PYTHON_I_S_B')
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--manifest-sha256', required=True)
    ap.add_argument('--source-root', type=Path, required=True)
    ap.add_argument('--previous-source-root', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    auth = load('boundary_capsule_auth', ROOT / 'source_auth.py')
    before = auth.authenticate(args.manifest_sha256)
    auth.need(args.output.is_absolute() and not args.output.exists() and
              not args.output.resolve().is_relative_to(ROOT), 'FRESH_EXTERNAL_OUTPUT_REQUIRED')
    dep = load('boundary_dependency_chain', ROOT / 'dependency_chain.py')
    old = args.source_root.resolve(); previous = args.previous_source_root.resolve()
    old_before = dep.inventory(old, dep.GEOMETRIC)
    previous_before = dep.inventory(previous, dep.PREVIOUS)
    path = ROOT / 'compute/certify_boundary_response_v2.py'
    auth.need(auth.sha(path.read_bytes()) == PRODUCER_SHA256, 'PRODUCER_EXTERNAL_PIN')
    producer = load('actual_boundary_producer', path)
    # These are the only producer globals changed. Its numerical code and
    # __file__ identity are untouched, and its own input authentication runs.
    producer.OLD = old
    producer.NEW = previous
    producer.FOLLOWON = ROOT / 'inputs/followon'
    report = producer.run(16)
    auth.need(auth.sha(path.read_bytes()) == PRODUCER_SHA256, 'PRODUCER_CHANGED')
    auth.need(dep.inventory(old, dep.GEOMETRIC) == old_before and
              dep.inventory(previous, dep.PREVIOUS) == previous_before, 'DEPENDENCY_CHANGED')
    auth.need(auth.authenticate(args.manifest_sha256) == before, 'SOURCE_CHANGED')
    with args.output.open('x') as handle:
        handle.write(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': report['status'], 'sha256': auth.sha(args.output.read_bytes())}, sort_keys=True))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError) as error:
        sys.exit('PORTABLE_ACTUAL_FAILED:' + str(error))
