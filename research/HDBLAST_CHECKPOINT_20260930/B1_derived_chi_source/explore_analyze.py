#!/usr/bin/env python3
"""Exploratory cell (third dated note in REGISTRATION.md): applies the registered per-run rule (b1_analyze.classify) to
runs/explore and reports run-end and last-reliable values side by side with the y = 1, phi* = 0.5, G = 100 grid cells
(b_cut and b_need).  Informative only; not used for any registered verdict.  Output: B1_EXPLORE.json."""
import json, glob, hashlib
from pathlib import Path
import numpy as np
import b1_analyze as BA
HERE = Path(__file__).resolve().parent
out = dict(status='numerical, exploratory (not registered)', runs={})
files = sorted(glob.glob(str(HERE/'runs/explore/*_summary.json'))) + [str(HERE/('runs/main/main_ps0.5_G100_y1_%s_%s_summary.json' % (b, z))) for b in ('bcut', 'bneed') for z in ('dzf1e-3', 'dzf5e-4')]
for f in files:
    s, ts = BA.load(f); d = BA.derived(ts); c = BA.classify(s, d)
    i = len(d['H0tau']) - 1
    while i > 0 and not np.isfinite(d['r'][i]): i -= 1
    end = dict(H0tau=float(d['H0tau'][i]), r=float(d['r'][i]), Omega_r=float(d['Omega_r'][i]), Omega_chi=float(d['Omega_chi'][i]),
               Omega_vac=float(d['Omega_vac'][i]), ln_a=float(d['ln_a'][i]))
    out['runs'][s['tag']] = dict(b=s['matter']['b'], lambda_c_full=BA.cutoff(s, d)['lambda_c_full'], classification=c['classification'],
                                 reliable_end_H0tau=c['H0tau_end'], reliable_end_reason=c['reliable_end_reason'],
                                 max_rel_H_near_shell=c.get('max_rel_H_near_shell'), max_rel_M_near_shell=c.get('max_rel_M_near_shell'),
                                 max_Rj=c['max_Rj'], last_reliable=c['last'], run_end=end, stop=s['stop_reason'])
    o = out['runs'][s['tag']]
    print('%-40s b=%-8.3g lam_c=%-6.3g rel_end=%5.2f maxM=%.3g Rj=%.3g | reliable-end r=%-9.4g Om_r=%-8.3g | run-end tau=%.2f r=%-9.4g Om_r=%.3g Om_chi=%.3g' % (
        s['tag'], o['b'], o['lambda_c_full'], o['reliable_end_H0tau'], o['max_rel_M_near_shell'] or -1, o['max_Rj'], o['last_reliable']['r'], o['last_reliable']['Omega_r'],
        end['H0tau'], end['r'], end['Omega_r'], end['Omega_chi']))
out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(HERE/'B1_EXPLORE.json').write_text(json.dumps(BA.clean(out), indent=1) + '\n')
