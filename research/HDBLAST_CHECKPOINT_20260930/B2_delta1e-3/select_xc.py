#!/usr/bin/env python3
"""x_c selection per REGISTRATION.md section 4 + dated note 1 (amendment 1).
Chain per (Y, dc): registered pre-run (rho grid dz_f 4e-4, kappa 0) -> (a) S1 kappa 0 -> (a') S1 kappa 10 -> (b) x_c = inf.
Writes runs/pre/xc_Y<Y>_dc<dc>.json with key xc_suggest (read by evolve_a1.py --xc_from) and the chain used.  Output also XC_SELECTION.json."""
import json, sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
P = HERE/'runs/pre'
out = {}
for Y, dc in [('1', '1e-2'), ('1', '1e-4'), ('0', '1e-2'), ('0', '1e-4')]:
    chain = ['pre_Y%s_dc%s' % (Y, dc), 'preS1_Y%s_dc%s' % (Y, dc), 'preS1k10_Y%s_dc%s' % (Y, dc)]
    used, steps, xc = None, [], None
    for tag in chain:
        f = P/(tag + '_summary.json')
        if not f.exists():
            steps.append(dict(tag=tag, status='not run')); continue
        s = json.load(open(f))
        steps.append(dict(tag=tag, stop_reason=s['stop_reason'], xc_suggest=s.get('xc_suggest'), T_end=s['T_end'], H0tau_end=s['H0tau_end']))
        if s.get('xc_suggest') is not None:
            used, xc = tag, s['xc_suggest']; break
    rec = dict(Y=float(Y), dc=float(dc), xc_suggest=xc, chosen_from=used if used else '(b) x_c = inf (no pre-run met the criterion)', chain=steps)
    (P/('xc_Y%s_dc%s.json' % (Y, dc))).write_text(json.dumps(rec, indent=1) + '\n')
    out['Y=%s dc=%s' % (Y, dc)] = rec
(HERE/'XC_SELECTION.json').write_text(json.dumps(dict(status='numerical', rule='REGISTRATION.md section 4 + dated note 1, amendment 1', selection=out), indent=1) + '\n')
print(json.dumps(out, indent=1))
