#!/usr/bin/env python3
"""Copy authenticated NPY member bytes into deterministic uncompressed capsules.

This is pre-registration file preparation, not a scientific calculation.
It does not deserialize, interpret, compare or evaluate ndarray payload values.
Only hashes, CRCs, lengths and NPY header metadata are inspected. The recorded
EPS/PI exact ratios are source scalar-constant metadata, not archive values.
"""
import ast
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import resource
import shutil
import struct
import subprocess
import tempfile
import traceback
import zipfile
import zlib

REPO = Path('/workspace/HDblast')
STAGE = Path(__file__).resolve().parent
PUBLISHED_COMMIT = '71d00cc423e9049c8166b7ee7afbefbac28d8a18'
CHECKPOINT = 'research/HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE_MEMORY_FOLLOWUP'
ARCHIVE_DIRECTORY = CHECKPOINT + '/outputs/fresh/independent'
MEMBERS = ('k', 'momentum_weights', 'u_4', 'w_4', 'u_5', 'w_5', 'observation_eta',
           'history_eta', 'history_values', 'history_ledger_integrand',
           'history_baseline_contact', 'history_source_jet', 'history_forcing_jet',
           'history_geometry', 'history_cutoffs', 'history_quantity_names', 'history_geometry_names')
ARCHIVES = tuple(f'metric_modes_{s}_{g}.npz' for s in ('positive_B', 'signed_uB') for g in ('coarse', 'fine'))
CHUNK_BYTES = 65536
CODE_AND_EVIDENCE = (
    ('independent/forced_metric.py', 'reference_code/forced_metric.py.txt'),
    ('independent/stream_npz.py', 'reference_code/stream_npz.py.txt'),
    ('review/saved_data_diagnostics/raw_metric_audit_v2.py', 'reference_code/raw_metric_audit_v2.py.txt'),
    ('review/saved_data_diagnostics/verify_c_metadata_v2.py', 'reference_code/verify_c_metadata_v2.py.txt'),
    ('review/saved_data_diagnostics/C_METADATA_EXACT.json', 'evidence/C_METADATA_EXACT.json'),
    ('review/saved_data_diagnostics/C_METADATA_EXACT_OPTIMIZED.json', 'evidence/C_METADATA_EXACT_OPTIMIZED.json'),
    ('review/saved_data_diagnostics/NORMAL_OPTIMIZED_INVARIANCE.json', 'evidence/NORMAL_OPTIMIZED_INVARIANCE.json'),
)
PIN_ONLY_EVIDENCE = (
    'review/saved_data_diagnostics/memory-record-normal-001/FAILED_SAVED_DATA_DIAGNOSTIC.json',
    'review/saved_data_diagnostics/memory-record-optimized-001/FAILED_SAVED_DATA_DIAGNOSTIC.json',
)

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def stream_digest(stream):
    h = hashlib.sha256(); length = 0
    for chunk in iter(lambda: stream.read(CHUNK_BYTES), b''):
        h.update(chunk); length += len(chunk)
    return h.hexdigest(), length

def sha(path):
    with path.open('rb') as stream:
        return stream_digest(stream)[0]

def git_pin(path):
    process = subprocess.Popen(['git', '-C', str(REPO), 'cat-file', 'blob', PUBLISHED_COMMIT + ':' + path], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    digest, length = stream_digest(process.stdout)
    error = process.stderr.read().decode()
    require(process.wait() == 0, 'Pinned Git object unavailable: ' + path + ': ' + error)
    return digest, length

def authenticate(path):
    local = REPO / path
    published_sha, published_bytes = git_pin(path)
    actual_sha = sha(local)
    require(actual_sha == published_sha and local.stat().st_size == published_bytes, 'Tracked publication blob mismatch: ' + path)
    return {'repository_path': path, 'published_commit': PUBLISHED_COMMIT,
            'sha256': actual_sha, 'bytes': published_bytes,
            'local_file_vs_pinned_Git_blob': 'IDENTICAL'}

def header(archive, member):
    with archive.open(member, 'r') as stream:
        magic = stream.read(8)
        require(magic[:6] == b'\x93NUMPY', 'Bad NPY magic: ' + member)
        version = tuple(magic[6:8])
        require(version in ((1, 0), (2, 0), (3, 0)), 'Unsupported NPY version')
        length_bytes = 2 if version == (1, 0) else 4
        raw_length = stream.read(length_bytes)
        require(len(raw_length) == length_bytes, 'Truncated NPY header length')
        length = int.from_bytes(raw_length, 'little')
        require(length <= 65536, 'Unbounded NPY header')
        raw_header = stream.read(length)
        require(len(raw_header) == length, 'Truncated NPY header')
        metadata = ast.literal_eval(raw_header.decode('utf8' if version == (3, 0) else 'latin1').strip())
        require(set(metadata) == {'descr', 'fortran_order', 'shape'}, 'Unexpected NPY header fields')
        dtype, shape, fortran = metadata['descr'], metadata['shape'], metadata['fortran_order']
        require(isinstance(dtype, str) and isinstance(shape, tuple) and all(type(d) is int and d >= 0 for d in shape) and type(fortran) is bool, 'Unexpected dtype/shape/order header')
        match = re.fullmatch(r'[<>|=]([fciuUSb])(\d+)', dtype)
        require(match is not None, 'Object/structured/unsupported dtype forbidden')
        kind, units = match.group(1), int(match.group(2))
        item_bytes = units * 4 if kind == 'U' else units
        offset = 8 + length_bytes + length
        expected_bytes = offset + math.prod(shape) * item_bytes
        require(expected_bytes == archive.getinfo(member).file_size, 'NPY header/payload length mismatch')
        return {'npy_version': list(version), 'dtype_descriptor': dtype, 'shape': list(shape),
                'fortran_order': fortran, 'header_bytes_including_magic': offset,
                'header_sha256': hashlib.sha256(magic + raw_length + raw_header).hexdigest(),
                'payload_bytes': expected_bytes - offset, 'array_payload_interpreted': False}

def deterministic_info(name):
    info = zipfile.ZipInfo(name, (1980, 1, 1, 0, 0, 0))
    info.compress_type = zipfile.ZIP_STORED
    info.create_system = 3
    info.external_attr = 0o100644 << 16
    info.internal_attr = 0
    info.comment = b''; info.extra = b''
    return info

def build_capsule(source, output):
    selected = []
    with zipfile.ZipFile(source, 'r') as before, zipfile.ZipFile(output, 'x', compression=zipfile.ZIP_STORED, allowZip64=True) as after:
        require(len(before.namelist()) == 631 and len(set(before.namelist())) == 631, 'Original archive member inventory differs')
        for base in MEMBERS:
            member = base + '.npy'
            source_info = before.getinfo(member)
            metadata = header(before, member)
            digest = hashlib.sha256(); crc = 0; length = 0
            with before.open(member, 'r') as src, after.open(deterministic_info(member), 'w', force_zip64=True) as dst:
                for chunk in iter(lambda: src.read(CHUNK_BYTES), b''):
                    dst.write(chunk); digest.update(chunk)
                    crc = zlib.crc32(chunk, crc); length += len(chunk)
            require(length == source_info.file_size and crc == source_info.CRC, 'Selected original member CRC/size mismatch')
            selected.append({'name': member, 'sha256': digest.hexdigest(), 'bytes': length,
                             'crc32_hex': f'{crc:08x}', 'original_compression': source_info.compress_type,
                             'capsule_compression': zipfile.ZIP_STORED, 'header': metadata})
    return selected

def verify_capsule(path, selected):
    with zipfile.ZipFile(path, 'r') as archive:
        require(archive.namelist() == [base + '.npy' for base in MEMBERS], 'Capsule member order/inventory differs')
        for record in selected:
            member = record['name']; info = archive.getinfo(member)
            require(info.compress_type == zipfile.ZIP_STORED and info.date_time == (1980, 1, 1, 0, 0, 0), 'Capsule ZIP contract differs')
            require(header(archive, member) == record['header'], 'Header bytes/metadata changed')
            with archive.open(member, 'r') as stream:
                digest, length = stream_digest(stream)
            require(digest == record['sha256'] and length == record['bytes'] and info.CRC == int(record['crc32_hex'], 16), 'Copied NPY member bytes/CRC differ')
    return {'raw_member_sha256_identity': 'PASS', 'CRC_verified_members': len(selected),
            'header_only_verified_members': len(selected), 'array_values_loaded_or_interpreted': False}

def main():
    require(not (STAGE / 'INPUT_MANIFEST.json').exists(), 'Refusing overwrite')
    original_results_pin = authenticate(ARCHIVE_DIRECTORY + '/results.json')
    shutil.copyfile(REPO / (ARCHIVE_DIRECTORY + '/results.json'), STAGE / 'original_results.json')
    require(sha(STAGE / 'original_results.json') == original_results_pin['sha256'], 'Copied results bytes changed')
    references = []
    for relative, destination in CODE_AND_EVIDENCE:
        pin = authenticate(CHECKPOINT + '/' + relative)
        target = STAGE / destination
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(REPO / (CHECKPOINT + '/' + relative), target)
        require(sha(target) == pin['sha256'], 'Reference copy changed')
        references.append(dict(pin, capsule_path=destination, executed=False))
    pin_only = [authenticate(CHECKPOINT + '/' + relative) for relative in PIN_ONLY_EVIDENCE]
    archives = []
    for name in ARCHIVES:
        original = REPO / ARCHIVE_DIRECTORY / name
        original_pin = authenticate(ARCHIVE_DIRECTORY + '/' + name)
        target = STAGE / name
        selected = build_capsule(original, target)
        verification = verify_capsule(target, selected)
        target_sha = sha(target)
        with tempfile.TemporaryDirectory(prefix='raw-rebuild-', dir=STAGE) as temp:
            rebuilt = Path(temp) / name
            rebuilt_selected = build_capsule(original, rebuilt)
            require(rebuilt_selected == selected and sha(rebuilt) == target_sha and rebuilt.stat().st_size == target.stat().st_size, 'Deterministic ZIP rebuild differs')
        archives.append({'capsule_path': name, 'capsule_sha256': target_sha, 'capsule_bytes': target.stat().st_size,
                         'original': original_pin, 'original_members': 631, 'selected_members': 17,
                         'selected_member_bytes': sum(r['bytes'] for r in selected), 'members': selected,
                         'verification': verification, 'deterministic_rebuild_same_sha256': True,
                         'physical_producer_run_status': 'FAILED_IN_PRIOR_REGISTERED_RUN_RETAINED'})
    files = {p.relative_to(STAGE).as_posix(): {'sha256': sha(p), 'bytes': p.stat().st_size}
             for p in sorted(STAGE.rglob('*')) if p.is_file() and p.name not in ('INPUT_MANIFEST.json', 'EXTRACTION_FAILURE.json')}
    manifest = {'schema_version': 1, 'status': 'PASS_RAW_BYTE_INPUT_EXTRACTION_ONLY',
                'recorded_utc': datetime.now(timezone.utc).isoformat(),
                'purpose': 'Self-contained raw input capsules for a separately prospectively registered twelve-case saved-ledger diagnostic.',
                'extraction_is_not_calculation': True, 'new_physical_evaluations': 0,
                'array_payload_values_loaded_interpreted_compared_or_evaluated': False,
                'source_or_reader_functions_executed': False,
                'allowed_pre_registration_operations': ['file copying', 'streamed SHA-256', 'CRC validation', 'NPY header dtype/shape/size inspection', 'deterministic ZIP rebuilding', 'source scalar-constant metadata recording'],
                'scientific_calibration_status': 'FAIL_UNCHANGED',
                'published_lineage_commit': PUBLISHED_COMMIT,
                'physical_run_public_freeze_commit': '19fde76912af6e2f30d6f55e066b27d88cc7e34b',
                'physical_run_registration_sha256': '1c5bc9b21b34b1d036b7a7d12ccdb6766ac2ba7835ea04182dba868843866607',
                'authentication_scope': 'Local input bytes match pinned published-commit Git blob objects; no remote action or network query performed by extractor.',
                'subset_zip_contract': {'compression': 'ZIP_STORED', 'compression_code': 0,
                                        'member_order': [m + '.npy' for m in MEMBERS],
                                        'timestamp': [1980, 1, 1, 0, 0, 0], 'create_system': 3,
                                        'external_attr': 0o100644 << 16, 'force_zip64_member_headers': True,
                                        'archive_comment': '', 'member_comments': '',
                                        'copy_chunk_bytes': CHUNK_BYTES, 'npy_member_bytes_reencoded': False,
                                        'longdouble_padding_bytes_preserved': True},
                'source_constants': {'source': 'reference_code/forced_metric.py.txt',
                                     'model': {'H': 1, 'mass_squared': 2, 'xi': 0, 'phi': 'fixed', 'incoming_state': 'BD'},
                                     'epsilon': {'source_expression': 'np.longdouble("0.0001")',
                                                 'intended_decimal_ratio': [1, 10000],
                                                 'exact_represented_binary80_ratio': [3777893186295716171, 37778931862957161709568]},
                                     'pi': {'source_expression': 'np.arccos(np.longdouble(-1))',
                                            'exact_represented_binary80_ratio': [14488038916154245685, 4611686018427387904],
                                            'interpretation': 'Exact ratio of the computed finite-precision constant, not an exact rational value of mathematical pi.'},
                                     'constant_capture_runtime': {'python': '3.12.14', 'numpy': '2.2.6', 'longdouble_nmant': 63},
                                     'momentum_measure_source_expression': 'k*k/(2*pi*pi)',
                                     'source_defined_cutoffs': [64, 128, 256],
                                     'source_defined_observation_eta': [-5.5, -4.5, -4.0, -3.5, -2.5, -1.5],
                                     'requested_source_observation_indices': [4, 5],
                                     'source_defined_quantity_names': ['q', 'q_prime', 'q_second', 'rho', 'p', 'Q0', 'rho0', 'p0', 'current']},
                'original_results': dict(original_results_pin, capsule_path='original_results.json', copied_without_numerical_reanalysis=True),
                'corrected_V2_reader': next(r for r in references if r['capsule_path'] == 'reference_code/raw_metric_audit_v2.py.txt'),
                'reference_code_and_small_receipts': references, 'corrected_V2_full_evidence_pin_only': pin_only,
                'capsules': archives, 'total_selected_NPY_members': 68,
                'total_header_only_verified_members': 68, 'total_CRC_verified_members': 68,
                'all_four_deterministic_rebuilds_identical': True,
                'extractor_peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                'files': files}
    (STAGE / 'INPUT_MANIFEST.json').write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': manifest['status'], 'capsules': [(r['capsule_path'], r['capsule_bytes'], r['capsule_sha256']) for r in archives],
                      'selected_members': 68, 'manifest_sha256': sha(STAGE / 'INPUT_MANIFEST.json'),
                      'extractor_peak_rss_kib': manifest['extractor_peak_rss_kib']}))

if __name__ == '__main__':
    try:
        main()
    except Exception:
        (STAGE / 'EXTRACTION_FAILURE.json').write_text(json.dumps({'status': 'FAIL_RAW_BYTE_EXTRACTION_ONLY', 'new_physical_evaluations': 0, 'array_values_evaluated': False, 'traceback': traceback.format_exc()}, indent=2) + '\n')
        raise
