"""Bounded scientific child. Authenticate all bytes before importing methods."""
import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import resource
import sys
import time

from registration_guard import authenticate, require, SourceAuthorization


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', required=True)
    parser.add_argument('--receipt', required=True)
    parser.add_argument('--receipt-sha256', required=True)
    parser.add_argument('--registration-sha256', required=True)
    parser.add_argument('--route', choices=('primary', 'independent'), required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--fabricated-only', action='store_true')
    args = parser.parse_args()
    require(os.getuid() != 0, 'Scientific child must be nonroot')
    root = Path(args.root).absolute()
    # No numerical source or third-party method import has occurred yet.
    verified = authenticate(root, args.receipt, args.receipt_sha256,
                            args.registration_sha256)
    require(sys.version_info[:3] == (3, 12, 14), 'Registered Python3.12.14 required')
    if not args.fabricated_only:
        require(not verified['receipt'].get('fabricated_fixture', False) and
                not verified['registration'].get('fabricated_fixture', False),
                'Fabricated readback cannot authorize a real source')
    require(Path(__file__).resolve() == root/'execution/run_registered.py',
            'Execute the frozen child, not an external copy')
    require(sys.modules['registration_guard'].__file__ ==
            str(root/'execution/registration_guard.py'), 'Unfrozen guard import')
    output = Path(args.output).absolute()
    require(not output.exists() and root not in output.parents,
            'Exclusive external output filename required')
    limits = verified['contract']['resources']
    require(resource.getrlimit(resource.RLIMIT_AS)[0] ==
            limits['each_route_peak_rss_kib']*1024 and
            resource.getrlimit(resource.RLIMIT_CPU)[0] ==
            limits['each_route_wall_seconds'], 'Outer resource limits absent/different')
    resource.setrlimit(resource.RLIMIT_NPROC, (0, 0))
    sys.dont_write_bytecode = True
    sys.set_int_max_str_digits(10000)
    # Freeze-verified paths only; strip inherited PYTHONPATH/import injection.
    sys.path.insert(0, str(root/'protocol'))
    sys.path.insert(0, str(root/args.route))
    auth = SourceAuthorization(verified, args.route, output.parent/'SOURCE_ATTEMPTS.jsonl')
    started = time.monotonic()
    if args.route == 'primary':
        import flint
        require(flint.__version__ == '0.9.0', 'python-flint pin differs')
        require(flint.__FLINT_VERSION__ == '3.6.0', 'FLINT pin differs')
        flint.ctx.threads = 1
        require(flint.ctx.threads == 1, 'FLINT thread pin differs')
    method = importlib.import_module('route')
    require(Path(method.__file__).resolve() == root/args.route/'route.py',
            'Wrong method import path')
    payload, method_evidence = (method.run_fabricated() if args.fabricated_only
                                else method.run(auth))
    contract = importlib.import_module('operator_probe_contract')
    semantic = contract.validate_payload(payload)
    if args.fabricated_only:
        require(not auth.events, 'Real callback invoked during fabricated entry')
        chronology = {'route': args.route, 'source_bundle_constructions': 0,
                      'archive_arrays_decoded': 0, 'events': [],
                      'fabricated_only': True}
    else:
        chronology = auth.complete()
    # Rehash all frozen sources/proofs after the run; no source repair in place.
    authenticate(root, args.receipt, args.receipt_sha256, args.registration_sha256)
    envelope = {'schema_version': 1, 'route': args.route, 'payload': payload,
        'method_evidence': method_evidence, 'semantic_validation': semantic,
        'execution': {'uid': os.getuid(), 'python': sys.version,
            'optimization': sys.flags.optimize,
            'wall_seconds_internal': time.monotonic()-started,
            'peak_rss_kib_internal': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'registration_sha256': args.registration_sha256,
            'remote_receipt_sha256': args.receipt_sha256,
            'freeze_commit': verified['receipt']['freeze_commit'],
            'source_chronology': chronology}}
    with output.open('x') as handle:
        json.dump(envelope, handle, sort_keys=True, separators=(',', ':'))
        handle.write('\n')
    require(output.stat().st_size > 0, 'Empty output')
    print(json.dumps({'status': ('COMPLETED_FABRICATED_ROUTE' if args.fabricated_only
                                 else 'COMPLETED_REGISTERED_SOURCE_ROUTE'),
        'route': args.route, 'output_bytes': output.stat().st_size,
        'output_sha256': hashlib.sha256(output.read_bytes()).hexdigest(),
        'semantic_validation': semantic}), flush=True)


if __name__ == '__main__':
    main()
