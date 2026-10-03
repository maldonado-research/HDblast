"""Synthetic-only timing including exact conversion, direct phases and flow I/O."""
from __future__ import annotations
import json
from pathlib import Path
import resource
import sys
import tempfile
import time
import numpy as np
from ledger_core import calculate_profiles, decimal
from diagnostic_independent import compare_flow_files, exact_decimal, compare_values
from exact_binary import require


def run():
    LD, CD = np.longdouble, np.clongdouble
    n, samples = 1024, 257
    # Fabricated spectrum and coefficients, deliberately distinct from retained data.
    k = (np.arange(n, dtype=LD) + LD(1)) / LD(47)
    weights = np.full(n, LD(1) / LD(991), dtype=LD)
    real = (np.arange(n, dtype=LD) % LD(19) + LD(1)) / LD(1000003)
    imag = -(np.arange(n, dtype=LD) % LD(17) + LD(1)) / LD(700003)
    ua = real.astype(CD) + CD(1j) * imag
    wa = (imag / LD(13)).astype(CD) + CD(1j) * (real / LD(17))
    ub = ua + CD(LD(1) / LD(123457))
    wb = wa - CD(1j) * (LD(1) / LD(765431))
    times = LD('-2.5') + np.arange(samples, dtype=LD) / LD(256)
    rows = []
    encoded_levels = {}
    started = time.monotonic()
    with tempfile.TemporaryDirectory(prefix='ledger-synthetic-benchmark-') as directory:
        paths = {}
        for dps in (80, 100):
            path = Path(directory) / f'flow_{dps}.jsonl'
            pass_started = time.monotonic()
            with path.open('x') as flow:
                def sink(record):
                    flow.write(json.dumps(record, separators=(',', ':')) + '\n')
                ctx, _, results = calculate_profiles(k, weights, ua, wa, ub, wb, times,
                                                     (256, 512, 1024), dps, sink)
                encoded = [{key: [decimal(ctx, value) for value in values]
                            if isinstance(values, list) else decimal(ctx, values)
                            for key, values in result.items()} for result in results]
                require(len(json.dumps(encoded)) > 1000, 'Synthetic serialization missing')
                encoded_levels[dps] = encoded
            elapsed = time.monotonic() - pass_started
            paths[dps] = path
            rows.append({'dps': dps, 'synthetic_nodes': n, 'samples': samples,
                         'pairs': n * samples, 'seconds': elapsed,
                         'flow_jsonl_bytes': path.stat().st_size,
                         'projected_10534912_pair_seconds': elapsed * 10534912 / (n * samples)})
        comparison_started = time.monotonic()
        gap, count = compare_flow_files(paths[80], paths[100])
        for first, second in zip(encoded_levels[80], encoded_levels[100]):
            for key in first:
                gap = max(gap, compare_values(first[key], second[key]))
        require(count == n and gap <= exact_decimal('1e-12'), 'Synthetic precision gap failed')
        comparison_seconds = time.monotonic() - comparison_started
    return {'status': 'PASS_SYNTHETIC_ONLY_TIMING_AND_IO', 'physical_values_loaded': False,
            'rows': rows, 'comparison_seconds': comparison_seconds,
            'synthetic_flow_precision_gap_exact_rational': str(gap),
            'total_seconds': time.monotonic() - started,
            'projected_both_precision_seconds_including_node_comparison':
                sum(row['projected_10534912_pair_seconds'] for row in rows)
                + comparison_seconds * (49152 / n),
            'peak_rss_kib': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'projection_is_not_full_retained_data_resource_proof': True}


if __name__ == '__main__':
    result = run()
    path = Path(sys.argv[1])
    require(not path.exists(), 'Timing receipt must be fresh')
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
