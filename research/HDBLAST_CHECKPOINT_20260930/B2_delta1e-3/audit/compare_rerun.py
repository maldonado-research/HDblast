#!/usr/bin/env python3
"""AUDIT: compare the auditor's rerun of main_Y1_dc1e-2_S1 (same command, no checkpoint) with the producer's saved output.
Output audit/AUDIT_RERUN.json."""
import json
from pathlib import Path
import numpy as np
H = Path(__file__).resolve().parent.parent
a = np.load(H/'audit/rerun/rerun_main_Y1_dc1e-2_S1_timeseries.npz'); b = np.load(H/'runs/main/main_Y1_dc1e-2_S1_timeseries.npz')
ca, cb = list(a['cols']), list(b['cols']); ra, rb = a['rec'], b['rec']
sa = json.load(open(H/'audit/rerun/rerun_main_Y1_dc1e-2_S1_summary.json')); sb = json.load(open(H/'runs/main/main_Y1_dc1e-2_S1_summary.json'))
n = min(len(ra), len(rb)); out = dict(n_records=[len(ra), len(rb)], rerun=dict(T_end=sa['T_end'], H0tau_end=sa['H0tau_end'], stop=sa['stop_reason']),
    producer=dict(T_end=sb['T_end'], H0tau_end=sb['H0tau_end'], stop=sb['stop_reason']), max_abs_diff={})
for c in ('T', 'H0tau', 'phi_b', 'H_over_H0', 'lapse', 'R', 'Wy', 'rad'):
    x, y = ra[:n, ca.index(c)], rb[:n, cb.index(c)]; m = np.isfinite(x) & np.isfinite(y)
    out['max_abs_diff'][c] = float(np.max(np.abs(x[m] - y[m])))
print(json.dumps(out, indent=1)); (H/'audit/AUDIT_RERUN.json').write_text(json.dumps(out, indent=1) + '\n')
