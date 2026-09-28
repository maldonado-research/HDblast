#!/usr/bin/env python3
"""Audit script 5: compare re-run JSON outputs (verify/rerun/) with the audited JSON files, number by number.
Writes verify/RERUN_COMPARISON.json with the maximum relative difference per file (runtime fields excluded)."""
import json, math, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
def walk(a, b, path, acc):
    if isinstance(a, dict):
        for k in a:
            if k in ('runtime_s',) or k not in b: 
                if k not in ('runtime_s',) and k not in b: acc['missing'].append(path + '/' + k)
                continue
            walk(a[k], b[k], path + '/' + k, acc)
    elif isinstance(a, list):
        if len(a) != len(b): acc['length_mismatch'].append(path); return
        for i, (x, y) in enumerate(zip(a, b)): walk(x, y, path + '[%d]' % i, acc)
    elif isinstance(a, (int, float)) and not isinstance(a, bool) and isinstance(b, (int, float)):
        if a == b: acc['n_equal'] += 1; return
        d = abs(a - b)/max(abs(a), abs(b), 1e-300)
        acc['n_diff'] += 1
        if d > acc['max_rel'][0]: acc['max_rel'] = (d, path, a, b)
    elif a != b:
        acc['other_diff'].append(path)
out = {}
for f in ('PREHEATING_RESULTS.json', 'CONTROLS.json', 'CORRECTION_FIT.json', 'SUMMARY.json'):
    p = HERE/'rerun'/f
    if not p.exists(): out[f] = 'not re-run'; continue
    acc = {'n_equal': 0, 'n_diff': 0, 'max_rel': (0.0, None, None, None), 'missing': [], 'length_mismatch': [], 'other_diff': []}
    walk(json.loads((HERE.parent/f).read_text()), json.loads(p.read_text()), '', acc)
    acc['missing'] = acc['missing'][:20]; acc['other_diff'] = acc['other_diff'][:20]
    out[f] = acc
(HERE/'RERUN_COMPARISON.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
print(json.dumps(out, indent=1, default=float))
