"""Synthetic-only provenance, input schema and exact-decimal comparison guards."""
from __future__ import annotations
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import zipfile
import numpy as np
from diagnostic_independent import (input_metadata, exact_decimal, compare_rows,
                                    frozen_guard, sha, CUTOFFS, SOURCES)
from exact_binary import require
from ledger_core import validate_arrays, MEMBERS, NAMES, GEOMETRY_NAMES, EPS_RATIO, PI_RATIO

LD, CD = np.longdouble, np.clongdouble


def rejects(function, message):
    try:
        function()
    except (RuntimeError, ValueError, TypeError, FileNotFoundError):
        return
    raise RuntimeError('Mutation accepted: ' + message)


def fabricated_arrays(setting='coarse'):
    n, nt, denominator = (8192, 577, 128) if setting == 'coarse' else (16384, 1153, 256)
    node_denominator = 32 if setting == 'coarse' else 64
    return {
        'k': (np.arange(n, dtype=LD) + LD('0.5')) / LD(node_denominator),
        'momentum_weights': np.full(n, LD(1) / LD(node_denominator), dtype=LD),
        **{name: np.zeros(n, dtype=CD) for name in ('u_4', 'w_4', 'u_5', 'w_5')},
        'observation_eta': np.array(['-5.5', '-4.5', '-4', '-3.5', '-2.5', '-1.5'], dtype=LD),
        'history_eta': LD(-6) + np.arange(nt, dtype=LD) / LD(denominator),
        'history_values': np.zeros((nt, 3, 9), dtype=LD),
        'history_ledger_integrand': np.zeros((nt, 3), dtype=LD),
        'history_baseline_contact': np.zeros((nt, 3), dtype=LD),
        'history_source_jet': np.zeros((nt, 6), dtype=LD),
        'history_forcing_jet': np.zeros((nt, 4), dtype=LD),
        'history_geometry': np.zeros((nt, 10), dtype=LD),
        'history_geometry_names': np.array(GEOMETRY_NAMES, dtype='<U8'),
        'history_cutoffs': np.array(CUTOFFS, dtype=np.int64),
        'history_quantity_names': np.array(NAMES, dtype='<U8'),
    }


def fabricated_capsules(root):
    # Fabricated arrays test only byte/header guards, never retained physical data.
    capsules = []
    for source in SOURCES:
        for setting in ('coarse', 'fine'):
            members, raw_members = [], {}
            arrays = fabricated_arrays(setting)
            for member in MEMBERS:
                stream = io.BytesIO()
                np.lib.format.write_array(stream, arrays[member], allow_pickle=False)
                raw = stream.getvalue()
                raw_members[member + '.npy'] = raw
                header = io.BytesIO(raw)
                version = np.lib.format.read_magic(header)
                shape, order, dtype = np.lib.format._read_array_header(header, version)
                members.append({'name': member + '.npy', 'bytes': len(raw),
                                'sha256': hashlib.sha256(raw).hexdigest(),
                                'header': {'shape': list(shape), 'dtype_descriptor': dtype.str,
                                           'fortran_order': bool(order)}})
            name = f'metric_modes_{source}_{setting}.npz'
            path = root / name
            with zipfile.ZipFile(path, 'w', compression=zipfile.ZIP_STORED) as archive:
                for member, raw in raw_members.items():
                    archive.writestr(member, raw)
            capsules.append({'capsule_path': name, 'capsule_bytes': path.stat().st_size,
                             'capsule_sha256': sha(path), 'members': members})
    manifest = {'schema_version': 1, 'scientific_calibration_status': 'FAIL_UNCHANGED',
                'source_constants': {'epsilon': {'exact_represented_binary80_ratio': list(EPS_RATIO)},
                                     'pi': {'exact_represented_binary80_ratio': list(PI_RATIO)},
                                     'model': {'H':1,'incoming_state':'BD','mass_squared':2,'phi':'fixed','xi':0},
                                     'source_defined_cutoffs':list(CUTOFFS),
                                     'source_defined_quantity_names':list(NAMES)},
                'capsules': capsules}
    (root / 'INPUT_MANIFEST.json').write_text(json.dumps(manifest))
    return manifest


def run():
    passed = []
    arrays = fabricated_arrays()
    indices, boundaries, dt = validate_arrays(arrays, 'coarse')
    require(len(indices) == 129 and indices[0] == 448 and indices[-1] == 576,
            'Inclusive mask changed')
    require(boundaries == (2048, 4096, 8192), 'Cutoff prefixes changed')
    passed.append('native_int64_cutoffs_and_inclusive_mask')
    for key in ('k', 'u_4'):
        original = arrays[key]
        arrays[key] = original.astype(np.float64 if key == 'k' else np.complex128)
        rejects(lambda: validate_arrays(arrays, 'coarse'), 'binary64 array ' + key)
        arrays[key] = original
    passed.append('binary64_array_mutations_rejected')
    for key in ('history_source_jet', 'history_forcing_jet', 'history_baseline_contact'):
        arrays[key][indices[0], 0] = LD(1)
        rejects(lambda: validate_arrays(arrays, 'coarse'), 'nonzero source/contact ' + key)
        arrays[key][indices[0], 0] = LD(0)
    passed.append('nonzero_source_forcing_baseline_contact_mutations_rejected')
    original = arrays['history_cutoffs']
    arrays['history_cutoffs'] = original.astype(LD)
    rejects(lambda: validate_arrays(arrays, 'coarse'), 'LD rather than int64 cutoff labels')
    arrays['history_cutoffs'] = original
    passed.append('wrong_cutoff_dtype_mutation_rejected')
    original = arrays['history_eta'].copy()
    arrays['history_eta'][indices[5]] += LD(1) / LD(1024)
    rejects(lambda: validate_arrays(arrays, 'coarse'), 'nonuniform retained history')
    arrays['history_eta'] = original
    passed.append('nonuniform_history_mutation_rejected')
    for value in (1.0, 1, 'NaN', 'Infinity', '1/3', ' 1', '', True):
        rejects(lambda value=value: exact_decimal(value), 'non-decimal-science value')
    passed.append('numeric_JSON_nonfinite_malformed_decimal_mutations_rejected')
    rows = [{'K': K, 'I_ab': '0.125', 'R': ['1.000', '-0.500']} for K in CUTOFFS]
    other = json.loads(json.dumps(rows))
    other[0]['R'][0] = '1.000000000002'
    require(compare_rows(rows, other) > exact_decimal('1e-12'), 'Precision gap not exact')
    other[0]['R'][0] = 1.0
    rejects(lambda: compare_rows(rows, other), 'numeric scientific field')
    passed.append('exact_decimal_precision_gap_and_schema_guard')
    with tempfile.TemporaryDirectory(prefix='ledger-independent-synthetic-') as directory:
        root = Path(directory)
        manifest = fabricated_capsules(root)
        _, capsules = input_metadata(root)
        require(len(capsules) == 4, 'Synthetic immutable capsule check failed')
        passed.append('fabricated_capsule_hash_CRC_NPY_header_preflight')
        path = root / 'INPUT_MANIFEST.json'
        manifest['capsules'][0]['capsule_sha256'] = '0' * 64
        path.write_text(json.dumps(manifest))
        rejects(lambda: input_metadata(root), 'capsule SHA mutation')
        manifest['capsules'][0]['capsule_sha256'] = sha(root / manifest['capsules'][0]['capsule_path'])
        passed.append('capsule_SHA_mutation_rejected')
        manifest['capsules'][0]['members'][0]['sha256'] = '0' * 64
        path.write_text(json.dumps(manifest))
        rejects(lambda: input_metadata(root), 'raw NPY member SHA mutation')
        passed.append('NPY_member_SHA_mutation_rejected')
        (root / 'FULL_REGISTRATION.json').write_text('{"schema_version":1,"files":{}}')
        rejects(lambda: frozen_guard(root, '0' * 64, '1' * 40), 'registration SHA mismatch')
        passed.append('registration_SHA_guard_before_any_array_decoder')
    return {'status': 'PASS_SYNTHETIC_ONLY', 'physical_values_loaded': False,
            'checks_passed': len(passed), 'checks': passed, 'optimized': not __debug__}


if __name__ == '__main__':
    result = run()
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    if path:
        require(not path.exists(), 'Guard receipt must be fresh')
        path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
