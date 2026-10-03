#!/usr/bin/env python3
"""Registered independent source-free saved-ledger diagnostic; no mode evolution."""
from __future__ import annotations
import argparse
from decimal import Decimal
from fractions import Fraction
import gc
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import re
import resource
import sys
import time
import zipfile
import mpmath
import numpy as np
from exact_binary import require
from ledger_core import (MEMBERS, EPS_RATIO, PI_RATIO, CUTOFFS, SETTINGS,
                         NAMES, evaluate_arrays)

SOURCES = ('positive_B', 'signed_uB')
PRECISIONS = (80, 100)
GATES = {'profile': '0.0000002', 'D_cont': '0.0000002', 'E_flow': '0.0000002',
         'decomposition': '0.000000000001', 'precision_gap': '0.000000000001',
         'serialization': '0.000000000001', 'attribution': '0.000002'}
CONFIGURATION = {
    'schema_version': 1, 'route': 'independent_mpf_bare_bilinear_fourier',
    'settings': SETTINGS, 'gates': GATES, 'sources': list(SOURCES), 'cutoffs': list(CUTOFFS),
    'decimal_precisions': [80, 100], 'phase_refresh_samples': 64,
    'profile_interval': ['-2.5', '-1.5'], 'profile_mask': 'inclusive',
    'epsilon_exact_binary_ratio': list(EPS_RATIO), 'pi_exact_binary_ratio': list(PI_RATIO),
    'analytic_background': 'L=-1/exact represented binary eta at each MP precision',
    'canonical_ledger': 'per-node initial-mode antiderivative with direct MPexp final phase',
    'recurrence_ledger_gate_each_precision': '0.000000000001',
    'simpson': 'original full3-column LD globalprefix; subtract in LD then exact binary conversion',
    'double_step_simpson': 'fine whole global F[::2],2dt,global endpoint prefixes and LD subtraction',
    'direct_reset_simpson': 'separate diagnostic on inclusive interval slice in LD',
    'endpoint_increment': 'MP subtraction of separately exact-binary converted retained R endpoints',
    'wall_seconds': 900, 'peak_rss_kib': 262144,
    'resource_scope': 'all12 cases and both precision passes in one standalone independent process',
    'physical_consistency_failure_exit': 0, 'internal_failure_exit': 2,
    'phase_projection': 'none; retain unrestricted nonzero Re(c)',
    'dtype_abi': 'native little-endian16-byte x87LD nmant63 and32-byte complexLD',
    'source_free_guard': 'saved source6jets,forcing4jets andbaseline contacts allzero over profile mask',
    'gap_universe': 'every scalar/profile scientific decimal and every streamed endpoint-flow decimal',
    'publication_or_old_calibration_status_change': False,
}
SCIENTIFIC_PATTERN = re.compile(r'-?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?\Z')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            digest.update(block)
    return digest.hexdigest()


def json_read(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, 'Duplicate JSON key')
            result[key] = value
        return result
    def invalid(value):
        raise RuntimeError('Nonfinite JSON literal: ' + value)
    return json.loads(Path(path).read_text(), object_pairs_hook=unique, parse_constant=invalid)


def exact_decimal(text):
    require(isinstance(text, str) and SCIENTIFIC_PATTERN.fullmatch(text) is not None,
            'Scientific values must be finite decimal strings')
    return Fraction(Decimal(text))


def runtime_check():
    info = np.finfo(np.longdouble)
    require(sys.version_info[:2] == (3, 12), 'Pinned Python 3.12 major/minor required')
    require(np.__version__ == '2.2.6' and mpmath.__version__ == '1.3.0',
            'Pinned NumPy2.2.6/mpmath1.3.0 required')
    require(sys.byteorder == 'little' and info.nmant == 63 and info.maxexp == 16384
            and info.minexp == -16382 and np.dtype(np.longdouble).itemsize == 16
            and np.dtype(np.clongdouble).itemsize == 32,
            'Saved native little-endian x87 binary80 ABI required')
    require(np.longdouble('0.0001').as_integer_ratio() == EPS_RATIO,
            'Represented producer epsilon differs')
    require(np.arccos(np.longdouble(-1)).as_integer_ratio() == PI_RATIO,
            'Represented producer pi differs')
    require(sys.flags.optimize == 0, 'The physical saved-data diagnostic uses normal Python')
    return {'python': sys.version.split()[0], 'numpy': np.__version__,
            'mpmath': mpmath.__version__, 'longdouble_nmant': int(info.nmant),
            'longdouble_itemsize': np.dtype(np.longdouble).itemsize,
            'complex_longdouble_itemsize': np.dtype(np.clongdouble).itemsize,
            'byteorder': sys.byteorder}


def frozen_guard(root, registration_sha256, freeze_commit):
    """Check explicit registration SHA before importing the shared integrity code."""
    require(re.fullmatch('[0-9a-f]{64}', registration_sha256) is not None,
            'Explicit SHA256 must be64 lower-case hex digits')
    require(re.fullmatch('[0-9a-f]{40}', freeze_commit) is not None,
            'Explicit public freeze commit must be40 lower-case hex digits')
    registration_path = root / 'FULL_REGISTRATION.json'
    require(not registration_path.is_symlink() and registration_path.is_file()
            and sha(registration_path) == registration_sha256, 'Registration bytes differ')
    registration = json_read(registration_path)
    require(registration['schema_version'] == 1 and isinstance(registration['files'], dict),
            'Unsupported registration schema')
    helper_path = root / 'code' / 'ledger_integrity.py'
    require(not helper_path.is_symlink() and helper_path.is_file()
            and registration['files'].get('code/ledger_integrity.py') == sha(helper_path),
            'Shared integrity helper must be registered before import')
    # Every route dependency must be inside this checkpoint and listed in both pins.
    route = Path(__file__).resolve().parent
    require(route == (root / 'independent').resolve(), 'Route source must come from frozen checkpoint')
    manifest_path = route / 'MANIFEST.json'
    require(registration['files'].get('independent/MANIFEST.json') == sha(manifest_path),
            'Independent source manifest is not registered')
    manifest = json_read(manifest_path)
    require(manifest['schema_version'] == 1 and isinstance(manifest['files'], dict),
            'Unknown independent manifest schema')
    for name in ('diagnostic_independent.py', 'ledger_core.py', 'profile_kernel.py', 'exact_binary.py',
                 'CONFIGURATION.json'):
        require(not (route / name).is_symlink() and manifest['files'].get(name) == sha(route / name)
                and registration['files'].get('independent/' + name) == sha(route / name),
                'Imported route dependency is not exactly source-pinned: ' + name)
    for module_name in ('ledger_core', 'profile_kernel', 'exact_binary'):
        require(Path(sys.modules[module_name].__file__).resolve() == route / (module_name + '.py'),
                'Loaded route module came from a different source path')
    configuration = json_read(route / 'CONFIGURATION.json')
    require(configuration == CONFIGURATION,
            'Unknown independent numerical configuration')
    spec = importlib.util.spec_from_file_location('registered_ledger_integrity', helper_path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    evidence = helper.verify_frozen(root, registration_sha256, freeze_commit)
    return evidence


def input_metadata(root):
    """Only metadata, raw byte hashes, CRC and NPY headers; no array payload decoding."""
    path = root / 'INPUT_MANIFEST.json'
    require(not path.is_symlink() and path.is_file(), 'Input manifest must be a regular file')
    manifest = json_read(path)
    require(manifest['schema_version'] == 1 and manifest['scientific_calibration_status'] == 'FAIL_UNCHANGED',
            'Input lineage/failed experiment status differs')
    require(tuple(manifest['source_constants']['epsilon']['exact_represented_binary80_ratio']) == EPS_RATIO
            and tuple(manifest['source_constants']['pi']['exact_represented_binary80_ratio']) == PI_RATIO,
            'Input manifest producer constants differ')
    require(manifest['source_constants']['model'] == {
        'H': 1, 'incoming_state': 'BD', 'mass_squared': 2, 'phi': 'fixed', 'xi': 0},
        'Inherited fixed physical model differs')
    require(manifest['source_constants']['source_defined_cutoffs'] == list(CUTOFFS)
            and manifest['source_constants']['source_defined_quantity_names'] == list(NAMES),
            'Inherited cutoff or quantity schema differs')
    expected_names = {f'metric_modes_{source}_{setting}.npz'
                      for source in SOURCES for setting in SETTINGS}
    require(len(manifest['capsules']) == 4
            and {c['capsule_path'] for c in manifest['capsules']} == expected_names,
            'Frozen four input capsules differ')
    capsules = {}
    for capsule in manifest['capsules']:
        name = capsule['capsule_path']
        setting = 'coarse' if name.endswith('_coarse.npz') else 'fine'
        spec = SETTINGS[setting]
        n, nt = spec['nodes'], spec['history_samples']
        expected_headers = {
            'k': ((n,), '<f16'), 'momentum_weights': ((n,), '<f16'),
            **{key: ((n,), '<c32') for key in ('u_4', 'w_4', 'u_5', 'w_5')},
            'observation_eta': ((6,), '<f16'), 'history_eta': ((nt,), '<f16'),
            'history_values': ((nt, 3, 9), '<f16'),
            'history_ledger_integrand': ((nt, 3), '<f16'),
            'history_baseline_contact': ((nt, 3), '<f16'),
            'history_source_jet': ((nt, 6), '<f16'),
            'history_forcing_jet': ((nt, 4), '<f16'),
            'history_geometry': ((nt, 10), '<f16'),
            'history_cutoffs': ((3,), '<i8'),
            'history_quantity_names': ((9,), '<U8'),
            'history_geometry_names': ((10,), '<U8'),
        }
        require(Path(name).name == name, 'Capsule path escapes input directory')
        capsule_path = root / name
        require(not capsule_path.is_symlink() and capsule_path.is_file()
                and capsule_path.stat().st_size == capsule['capsule_bytes']
                and sha(capsule_path) == capsule['capsule_sha256'], 'Capsule byte identity failed')
        expected_members = {member + '.npy' for member in MEMBERS}
        require(len(capsule['members']) == len(MEMBERS)
                and {m['name'] for m in capsule['members']} == expected_members,
                'Input manifest member inventory differs')
        with zipfile.ZipFile(capsule_path) as archive:
            require(len(archive.namelist()) == len(MEMBERS)
                    and set(archive.namelist()) == expected_members and archive.testzip() is None,
                    'Capsule member/CRC identity failed')
            for member in capsule['members']:
                raw = archive.read(member['name'])
                require(len(raw) == member['bytes']
                        and hashlib.sha256(raw).hexdigest() == member['sha256'], 'NPY member hash differs')
                stream = io.BytesIO(raw)
                version = np.lib.format.read_magic(stream)
                shape, order, dtype = np.lib.format._read_array_header(stream, version)
                require(list(shape) == member['header']['shape']
                        and dtype.str == member['header']['dtype_descriptor']
                        and bool(order) == member['header']['fortran_order'], 'NPY header schema differs')
                expected_shape, expected_dtype = expected_headers[member['name'][:-4]]
                require(shape == expected_shape and dtype.str == expected_dtype and not order,
                        'Frozen NPY shape/dtype/order rejected before payload decode')
        capsules[name] = capsule_path
    return manifest, capsules


def compare_values(first, second):
    """Exact decimal rational comparison; numerical JSON values are rejected."""
    if isinstance(first, list):
        require(isinstance(second, list) and len(first) == len(second), 'Precision profile lengths differ')
        return max((compare_values(a, b) for a, b in zip(first, second)), default=Fraction(0))
    return abs(exact_decimal(first) - exact_decimal(second))


def compare_rows(first, second):
    require(len(first) == len(second) == 3, 'Three canonical cutoff rows required per level')
    maximum = Fraction(0)
    for a, b, K in zip(first, second, CUTOFFS):
        require(a.keys() == b.keys() and a['K'] == b['K'] == K, 'Precision row inventory differs')
        for key in a:
            if key != 'K':
                maximum = max(maximum, compare_values(a[key], b[key]))
    return maximum


def compare_flow_files(path80, path100):
    maximum = Fraction(0)
    count = 0
    with path80.open() as a, path100.open() as b:
        while True:
            one, two = a.readline(), b.readline()
            require(bool(one) == bool(two), 'Precision streamed flow lengths differ')
            if not one:
                break
            first, second = json.loads(one), json.loads(two)
            require(first.keys() == second.keys() and first['node_index'] == second['node_index'] == count,
                    'Precision streamed flow node inventory differs')
            for key in first:
                if key != 'node_index':
                    maximum = max(maximum, compare_values(first[key], second[key]))
            count += 1
    return maximum, count


def run(root, output, started):
    def resource_check():
        require(time.monotonic() - started <= 900, 'Shared all-case900s wall budget exceeded')
        require(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss <= 262144,
                'Shared all-case262144KiB memory budget exceeded')
    manifest, capsules = input_metadata(root / 'inputs')
    records, science_failures, internal_failures = [], [], []
    maximum_gap = Fraction(0)
    streamed_files = []
    progress_path = output / 'progress.jsonl'
    with progress_path.open('x') as progress:
        def note(event, **fields):
            resource_check()
            progress.write(json.dumps({'event': event, 'elapsed_seconds': time.monotonic() - started,
                                       'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                                       **fields}, sort_keys=True) + '\n')
            progress.flush()
        note('immutable_input_bytes_verified_before_array_decode')
        for source in SOURCES:
            for setting in SETTINGS:
                name = f'metric_modes_{source}_{setting}.npz'
                with np.load(capsules[name], allow_pickle=False) as archive:
                    arrays = {member: archive[member] for member in MEMBERS}
                levels, flow_paths = {}, {}
                for dps in PRECISIONS:
                    flow_path = output / f'flow_{source}_{setting}_{dps}.jsonl'
                    with flow_path.open('x') as flow:
                        def sink(record):
                            flow.write(json.dumps(record, separators=(',', ':')) + '\n')
                        def block_progress():
                            resource_check()
                        levels[str(dps)] = evaluate_arrays(arrays, setting, dps, sink, block_progress)
                    flow_paths[dps] = flow_path
                    streamed_files.append({'path': flow_path.name, 'bytes': flow_path.stat().st_size,
                                           'sha256': sha(flow_path)})
                    note('precision_pass_complete', source=source, setting=setting, dps=dps)
                gap = compare_rows(levels['80'], levels['100'])
                flow_gap, count = compare_flow_files(flow_paths[80], flow_paths[100])
                require(count == SETTINGS[setting]['nodes'], 'Not all raw endpoint flow defects retained')
                gap = max(gap, flow_gap)
                maximum_gap = max(maximum_gap, gap)
                if gap > Fraction(GATES['precision_gap']):
                    internal_failures.append({'source': source, 'setting': setting,
                                              'gate': 'precision_gap', 'exact_gap': str(gap)})
                for row in levels['100']:
                    for key, gate in (('R_profile_max_error', 'profile'), ('P_profile_max_error', 'profile'),
                                      ('D_cont', 'D_cont'), ('E_flow', 'E_flow')):
                        if abs(exact_decimal(row[key])) > Fraction(GATES[gate]):
                            science_failures.append({'source': source, 'setting': setting, 'K': row['K'],
                                                     'gate': key, 'value': row[key], 'threshold': GATES[gate]})
                for dps in PRECISIONS:
                    for row in levels[str(dps)]:
                        if abs(exact_decimal(row['decomposition_error'])) > Fraction(GATES['decomposition']):
                            internal_failures.append({'source': source, 'setting': setting, 'K': row['K'],
                                                      'dps': dps, 'gate': 'decomposition_error',
                                                      'value': row['decomposition_error']})
                        if abs(exact_decimal(row['I_recurrence_minus_direct'])) > Fraction('1e-12'):
                            internal_failures.append({'source': source, 'setting': setting, 'K': row['K'],
                                                      'dps': dps, 'gate': 'direct_recurrence_ledger',
                                                      'value': row['I_recurrence_minus_direct']})
                records.append({'source': source, 'setting': setting, 'levels': levels,
                                'precision_gap_exact_rational': str(gap),
                                'flow_precision_gap_exact_rational': str(flow_gap)})
                note('archive_complete', source=source, setting=setting)
                del arrays
                gc.collect()
        resource_check()
    attribution = []
    if not science_failures and not internal_failures:
        for record in records:
            for row in record['levels']['100']:
                defect = abs(exact_decimal(row['D_S']))
                if defect > Fraction(GATES['attribution']) \
                   and abs(exact_decimal(row['D_cont'])) <= defect / 10 \
                   and abs(exact_decimal(row['E_Q'])) >= 9 * defect / 10:
                    attribution.append({'source': record['source'], 'setting': record['setting'], 'K': row['K']})
    classification = 'CONSISTENCY_FAILURE' if science_failures or internal_failures else \
        'LEDGER_ERROR_DEMONSTRATED' if attribution else 'NO_GATE_SCALE_ATTRIBUTION'
    return {'schema_version': 1, 'route': 'independent_mpf_bare_bilinear_fourier',
            'status': 'FATAL_DIAGNOSTIC_FAILURE' if internal_failures else 'COMPLETED_SAVED_DATA_DIAGNOSTIC',
            'classification': classification, 'old_metric_status': 'FAIL',
            'fixed_cases': 12, 'attribution_cases': attribution,
            'failures': science_failures, 'fatal_failures': internal_failures,
            'records': records, 'flow_files': streamed_files,
            'precision_gap_max_exact_rational': str(maximum_gap),
            'input_manifest_sha256': sha(root / 'inputs' / 'INPUT_MANIFEST.json'),
            'input_lineage_public_commit': manifest['published_lineage_commit'],
            'resources': {'seconds': time.monotonic() - started,
                          'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                          'scope': 'all twelve cases and both decimal precision passes in this independent route',
                          'limit_seconds': 900, 'limit_peak_rss_kib': 262144},
            'arithmetic': {'levels_decimal_dps': [80, 100], 'rounding': 'MPF nearest',
                           'phase_refresh_samples': 64, 'shared_lower_precision_coefficients': False,
                           'gap_universe': 'every non-K scalar/profile in records.levels and every non-node_index decimal in all streamed endpoint flow files',
                           'serialization': 'exact Fraction(Decimal(string)) compared with exact MPF dyadic; absolute bound1e-12',
                           'empirical_arithmetic_check_not_interval_certificate': True},
            'gates': GATES}


def main():
    started = time.monotonic()
    parser = argparse.ArgumentParser()
    parser.add_argument('--checkpoint-root', required=True, type=Path)
    parser.add_argument('--registration-sha256', required=True)
    parser.add_argument('--freeze-commit', required=True)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    require(not args.checkpoint_root.is_symlink(), 'Checkpoint root symlink is forbidden')
    root = args.checkpoint_root.resolve()
    require(not args.output_dir.exists(), 'Fresh output directory required; old evidence is immutable')
    require(not args.output_dir.resolve().is_relative_to(root),
            'Diagnostic output must be outside immutable checkpoint')
    runtime = runtime_check()
    provenance = frozen_guard(root, args.registration_sha256, args.freeze_commit)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    try:
        result = run(root, args.output_dir, started)
        provenance_after = frozen_guard(root, args.registration_sha256, args.freeze_commit)
        require(provenance == provenance_after, 'Frozen provenance changed during diagnostic')
        result['provenance'] = {'registration_sha256': args.registration_sha256,
                                'public_freeze_commit': args.freeze_commit,
                                'shared_integrity_receipt_before': provenance,
                                'shared_integrity_receipt_after': provenance_after}
        result['runtime'] = runtime
        target = args.output_dir / 'diagnostic.json'
        target.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
        print(json.dumps({'status': result['status'], 'classification': result['classification'],
                          'output': str(target), 'sha256': sha(target)}))
        return 2 if result['fatal_failures'] else 0
    except Exception as error:
        failure = {'schema_version': 1, 'status': 'FAIL_INTERNAL_DIAGNOSTIC',
                   'old_metric_status': 'FAIL', 'error_type': type(error).__name__,
                   'error': str(error), 'registration_sha256': args.registration_sha256,
                   'public_freeze_commit': args.freeze_commit,
                   'resources': {'seconds': time.monotonic() - started,
                                 'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}}
        (args.output_dir / 'FAILED_DIAGNOSTIC.json').write_text(json.dumps(failure, indent=2) + '\n')
        raise


if __name__ == '__main__':
    sys.exit(main())
