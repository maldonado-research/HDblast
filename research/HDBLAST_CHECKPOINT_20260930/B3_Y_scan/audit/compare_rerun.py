#!/usr/bin/env python3
"""Audit: compare the audit re-run of main_Y2_dc1e-2_dzf5e-4 (same code hash, same args except no --save_T) with the
producer's saved run, and apply the audit's own plateau finder.  Output audit/AUDIT_RERUN.json.  Numerical."""
import json, sys
from pathlib import Path
import numpy as np
H = Path(__file__).resolve().parent; sys.path.insert(0, str(H))
import importlib.util
spec = importlib.util.spec_from_file_location('ab', H/'audit_b3.py')
src = (H/'audit_b3.py').read_text().split("out = dict(status='numerical (audit)'")[0]
ns = {'__file__': str(H/'audit_b3.py')}; exec(src, ns)
s0, t0 = ns['series'](H.parent/'runs/main/main_Y2_dc1e-2_dzf5e-4_summary.json')
s1, t1 = ns['series'](H/'rerun/rerun_Y2_dc1e-2_dzf5e-4_summary.json')
n = min(len(t0['T']), len(t1['T']))
d = {k: float(np.nanmax(np.abs(t0[k][:n] - t1[k][:n]))) for k in ['H0tau', 'H_over_H0', 'Wy', 'R', 'phi_b']}
out = dict(n_common_records=n, max_abs_diff=d, code_sha_equal=s0['code_sha256'] == s1['code_sha256'],
           plateau_saved=ns['plateau'](t0), plateau_rerun=ns['plateau'](t1), stop_saved=s0['stop_reason'], stop_rerun=s1['stop_reason'])
for k in ['plateau_saved', 'plateau_rerun']: out[k].pop('idx', None)
(H/'AUDIT_RERUN.json').write_text(json.dumps(out, indent=1, default=float) + '\n'); print(json.dumps(out, indent=1, default=float))
