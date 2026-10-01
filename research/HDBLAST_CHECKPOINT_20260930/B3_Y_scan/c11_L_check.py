#!/usr/bin/env python3
"""C11: L = 10 (B3) versus L = 16 (A1) at equal spacing, Y = 3, dc = 1e-2.  Max relative difference of W a^4 and H up to
H0 tau = 20 (registered tolerance 1e-3).  Writes C11_L_CHECK.json."""
import json, sys, hashlib
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import a1_analyze as A
A1 = Path('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum/runs/main')
out = dict(status='numerical', tolerance=1e-3, pairs={})
for dz in ['1e-3', '5e-4']:
    fb = HERE/'runs/main'/('main_Y3_dc1e-2_dzf%s_summary.json' % dz); fa = A1/('main_dstar_Y3_dc1e-2_dzf%s_summary.json' % dz)
    if not fb.exists(): continue
    sb, tb = A.load(fb); sa, ta = A.load(fa); db = A.derived(sb, tb); da = A.derived(sa, ta)
    taus = np.linspace(0.5, 20, 400)
    def at(d, k):
        t = d['H0tau']; keep = np.concatenate(([True], np.diff(t) > 1e-12)); return np.interp(taus, t[keep], d[k][keep])
    rel = lambda k: float(np.max(np.abs(at(db, k) - at(da, k))/np.maximum(np.abs(at(da, k)), 1e-300)))
    # W a^4 changes sign early (W ~ 0 crossing); use the relative difference where |W a^4| > 1% of its max on the window
    wa, wb = at(da, 'Wa4'), at(db, 'Wa4'); m = np.abs(wa) > 0.01*np.abs(wa).max()
    out['pairs']['dzf' + dz] = dict(max_rel_dH=rel('H_over_H0'), max_rel_dWa4=float(np.max(np.abs(wb[m] - wa[m])/np.abs(wa[m]))),
                                     max_abs_dphi=float(np.max(np.abs(at(db, 'phi_b') - at(da, 'phi_b')))), xc=[sb['params']['xc'], sa['params']['xc']], L=[sb['params']['L'], sa['params']['L']])
    out['pairs']['dzf' + dz]['PASS'] = bool(out['pairs']['dzf' + dz]['max_rel_dH'] <= 1e-3 and out['pairs']['dzf' + dz]['max_rel_dWa4'] <= 1e-3)
out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(HERE/'C11_L_CHECK.json').write_text(json.dumps(out, indent=1) + '\n'); print(json.dumps(out, indent=1))
