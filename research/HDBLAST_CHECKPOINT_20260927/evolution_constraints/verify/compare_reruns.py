#!/usr/bin/env python3
"""Compare verify/rerun/*.json|npz with the workstream's saved runs (records, snapshots, final state)."""
import json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
PAIRS = {'base_h4e-4': 'final/base_h4e-4', 'proj_h4e-4': 'final/proj_h4e-4', 'k10_h4e-4': 'final/k10_h4e-4',
         'k10p_h4e-4': 'final/k10p_h4e-4', 'o6p_h4e-4': 'final/o6p_h4e-4', 'km5_h4e-4': 'ctrl/o4km5_h4e-4',
         'proj_h2e-4': 'final/proj_h2e-4', 'short_base_h5e-5': 'short/base_h5e-5', 'short_proj_h5e-5': 'short/proj_h5e-5'}
out = {}
for mine, theirs in PAIRS.items():
    fm, ft = HERE / f'rerun/{mine}.json', ROOT / f'runs/{theirs}.json'
    if not fm.exists(): out[mine] = 'not run'; continue
    a, b = json.loads(fm.read_text()), json.loads(ft.read_text())
    keys = ['H_max', 'M_max', 'H_L2', 'Cplus_w_max', 'Cminus_w_max', 'shell_phi_dev', 'max_abs_f']
    rel = 0.
    for ra, rb in zip(a['records'], b['records']):
        assert abs(ra['time'] - rb['time']) < 1e-12
        for k in keys:
            rel = max(rel, abs(ra[k] - rb[k]) / max(abs(rb[k]), 1e-300))
    A, B = np.load(str(fm)[:-5] + '.npz'), np.load(ROOT / f'runs/{theirs}.npz')
    fin = float(np.abs(A['final'] - B['final']).max() / np.abs(B['final']).max())
    snaps = {t: a['snapshots'][t]['H_max'] for t in a['snapshots']}
    out[mine] = dict(records=len(a['records']), max_relative_record_difference=rel, final_state_relative_max_difference=fin,
                     bit_identical_final=bool(np.array_equal(A['final'], B['final'])), H_max_snapshots_rerun=snaps,
                     H_max_snapshots_saved={t: b['snapshots'][t]['H_max'] for t in b['snapshots']},
                     projection_iterations=(len(a['projection']['iterations']) if a.get('projection') else None),
                     runtime_rerun=a['runtime_seconds'], runtime_saved=b['runtime_seconds'])
(HERE / 'compare_reruns.json').write_text(json.dumps(out, indent=1) + '\n')
for k, v in out.items():
    print(k, v if isinstance(v, str) else {x: v[x] for x in ('max_relative_record_difference', 'final_state_relative_max_difference', 'bit_identical_final', 'projection_iterations')})
