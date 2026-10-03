#!/usr/bin/env python3
"""Post-run storage audit; one saved ndarray pair at a time, no physics."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import resource
import sys
import numpy as np

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), 'Refusing overwrite')
    comparison = json.loads(args.comparison.read_text())
    require(comparison['status'] == 'PASS_EXACT_SAVED_ARRAY_INVARIANCE', 'Exact-value comparison must complete first')
    # The ABI contract is deliberately specific. x87 binary80 has a 64-bit
    # significand and 15-bit exponent/sign, stored in the first 10 bytes of
    # little-endian 16-byte long-double lanes; the other six bytes are padding.
    require(sys.byteorder == 'little' and platform.machine() == 'x86_64', 'Unsupported byte-lane ABI')
    require(np.dtype(np.longdouble).itemsize == 16 and np.finfo(np.longdouble).nmant == 63 and np.finfo(np.longdouble).nexp == 15, 'Unsupported long-double representation')
    require(np.dtype(np.clongdouble).itemsize == 32, 'Unsupported complex-long-double representation')
    fields = []
    value_byte_mismatches = padding_byte_differences = 0
    for record in comparison['archives']:
        before, after = Path(record['prior_file']), Path(record['current_file'])
        require(sha(before) == record['prior_sha256'] and sha(after) == record['current_sha256'], 'Archives changed after exact comparison')
        with np.load(before, allow_pickle=False) as old, np.load(after, allow_pickle=False) as new:
            for field in record['fields']:
                if field['array_storage_bytes_equal']:
                    continue
                name = field['name']
                a, b = old[name], new[name]
                require(a.dtype == b.dtype and a.shape == b.shape and np.array_equal(a, b), 'Saved values changed')
                require(a.dtype.str in ('<f16', '<c32'), 'Storage difference outside registered long-double lanes')
                av = np.frombuffer(a.tobytes(order='C'), dtype=np.uint8).reshape(-1, 16)
                bv = np.frombuffer(b.tobytes(order='C'), dtype=np.uint8).reshape(-1, 16)
                value_count = int(np.count_nonzero(av[:, :10] != bv[:, :10]))
                padding_count = int(np.count_nonzero(av[:, 10:] != bv[:, 10:]))
                fields.append({'source': record['source'], 'setting': record['setting'], 'name': name,
                               'dtype': a.dtype.str, 'longdouble_lanes': int(av.shape[0]),
                               'binary80_value_bytes_different': value_count,
                               'unused_padding_bytes_different': padding_count,
                               'storage_difference_is_padding_only': value_count == 0 and padding_count > 0})
                value_byte_mismatches += value_count
                padding_byte_differences += padding_count
                del av, bv, a, b
    status = 'PASS_BINARY80_VALUE_BYTES_IDENTICAL_PADDING_ONLY' if value_byte_mismatches == 0 and all(f['storage_difference_is_padding_only'] for f in fields) else 'FAIL_STORAGE_DIFFERENCE_NOT_PADDING_ONLY'
    report = {'status': status, 'new_physical_evaluations': 0,
              'scientific_status': 'FAIL', 'metric_calibration_passed': False,
              'scope': 'Storage-only diagnosis of all nonidentical ndarray byte strings from the exact saved-array comparison.',
              'abi': {'machine': platform.machine(), 'byteorder': sys.byteorder, 'lane_bytes': 16,
                      'binary80_value_offsets': list(range(10)), 'unused_padding_offsets': list(range(10, 16)),
                      'longdouble_nmant': np.finfo(np.longdouble).nmant, 'longdouble_nexp': np.finfo(np.longdouble).nexp},
              'audited_storage_unequal_arrays': len(fields), 'binary80_value_bytes_different': value_byte_mismatches,
              'unused_padding_bytes_different': padding_byte_differences,
              'comparison_report_sha256': sha(args.comparison), 'audit_source_sha256': sha(Path(__file__)),
              'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'recorded_utc': datetime.now(timezone.utc).isoformat(), 'fields': fields}
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: report[k] for k in ('status', 'audited_storage_unequal_arrays', 'binary80_value_bytes_different', 'unused_padding_bytes_different', 'peak_rss_kib')}))
    require(status.startswith('PASS_'), 'All storage differences preserved; padding-only audit failed')

if __name__ == '__main__':
    main()
