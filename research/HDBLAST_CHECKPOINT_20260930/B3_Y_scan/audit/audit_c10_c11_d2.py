#!/usr/bin/env python3
"""Audit: (1) recompute C11 (L=10 B3 Y=3 runs vs L=16 A1 runs) over the WHOLE common reliable range, not only H0tau<=20;
(2) recompute the C10 drifts with an independent settled reference (median W a^4 over H0tau 13-14 of the full 5e-4 run);
(3) re-run the D2 eigenvalue model problem and also report the RK4 check on the imaginary axis (limit 2*sqrt(2)) for
the interior wave modes.  Output audit/AUDIT_C10_C11_D2.json.  Numerical."""
import json, sys, math
from pathlib import Path
import numpy as np
B3 = Path(__file__).resolve().parent.parent; sys.path.insert(0, str(B3))
A1 = Path('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum/runs/main')
def raw(f):
    z = np.load(str(f).replace('_summary.json', '_timeseries.npz')); c = [str(x) for x in z['cols']]
    t = {k: z['rec'][:, i] for i, k in enumerate(c)}; t['Wa4'] = t['Wy']*np.exp(4*t['ln_a']); return t
def at(t, k, x):
    keep = np.concatenate(([True], np.diff(t['H0tau']) > 1e-12)); return np.interp(x, t['H0tau'][keep], t[k][keep])
out = {}
c11 = {}
for dz in ['1e-3', '5e-4']:
    b = raw(B3/'runs/main'/('main_Y3_dc1e-2_dzf%s_summary.json' % dz)); a = raw(A1/('main_dstar_Y3_dc1e-2_dzf%s_summary.json' % dz))
    for lim in [20.0, 30.0]:
        x = np.linspace(0.5, lim, 400)
        c11['dzf%s_to_%g' % (dz, lim)] = dict(max_rel_dWa4=float(np.max(np.abs(at(b, 'Wa4', x)/at(a, 'Wa4', x) - 1))),
                                              max_abs_dH=float(np.max(np.abs(at(b, 'H_over_H0', x) - at(a, 'H_over_H0', x)))))
out['C11_recomputed'] = c11
full = raw(B3/'runs/main/main_Y3_dc1e-2_dzf5e-4_summary.json'); c10 = raw(B3/'runs/ctl/ctl_Y3_dc1e-2_restart1e-3to5e-4_summary.json')
co = raw(B3/'runs/main/main_Y3_dc1e-2_dzf1e-3_summary.json')
# ctl is restart-only: splice parent (1e-3) records before restart
m = co['H0tau'] < c10['H0tau'][0] - 1e-9
spl = {k: np.concatenate([co[k][m], c10[k]]) for k in c10}
ref = float(np.median(full['Wa4'][(full['H0tau'] >= 13) & (full['H0tau'] <= 14)]))
rows = []
for x in [18, 20, 22, 24]:
    rows.append(dict(H0tau=x, restart=float(at(spl, 'Wa4', x)/ref - 1), full5e4=float(at(full, 'Wa4', x)/ref - 1), coarse=float(at(co, 'Wa4', x)/ref - 1)))
out['C10_recomputed'] = dict(ref=ref, rows=rows)
import evolve_a1 as E
d2 = []
for Y in [2, 3, 5]:
    n = 400; dz = 1e-3; D1, D2, g1, g2 = E.build_operators(np.full(n, dz), np.zeros(n), 4)
    M = np.zeros((2*n, 2*n)); M[:n, n:] = np.eye(n); M[n:, :n] = D2.toarray(); M[n:, 2*n - 1] += g2*(-Y/2); M[0, :] = 0; M[n, :] = 0
    ev = np.linalg.eigvals(M)*dz
    # RK4 amplification |1+z+z^2/2+z^3/6+z^4/24| for each eigenvalue at cfl 0.5 and 0.25
    amp = lambda z: np.abs(1 + z + z**2/2 + z**3/6 + z**4/24)
    d2.append(dict(Y=Y, max_abs=float(np.abs(ev).max()), most_negative_real=float(ev.real.min()), max_imag=float(np.abs(ev.imag).max()),
                   max_amp_cfl0p5=float(amp(0.5*ev).max()), max_amp_cfl0p25=float(amp(0.25*ev).max())))
out['D2_rerun'] = d2
(Path(__file__).parent/'AUDIT_C10_C11_D2.json').write_text(json.dumps(out, indent=1) + '\n'); print(json.dumps(out, indent=1))
