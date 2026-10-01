#!/usr/bin/env python3
"""Third-note trigger for Y = 3: apparent order of the late W a^4 drift between the 5e-4 run and the 2.5e-4 hybrid (and, if
present, between the 2.5e-4 and 1.25e-4 hybrids), relative to the settled value of the finest run; writes TERTIARY_ORDER.json"""
import json, sys, math, hashlib
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import a1_analyze as A, b3_analyze as B
Y = sys.argv[1] if len(sys.argv) > 1 else '3'
sfx = '_cfl0.25' if Y == '5' else ''
runs = [('5e-4', HERE/'runs/main'/('main_Y%s_dc1e-2_dzf5e-4%s%s_summary.json' % (Y, sfx, '_ext' if Y == '5' else ''))),
        ('2.5e-4', HERE/'runs/fine'/('fine_Y%s_dc1e-2_dzf2.5e-4%s_summary.json' % (Y, sfx))),
        ('1.25e-4', HERE/'runs/fine'/('fine_Y%s_dc1e-2_dzf1.25e-4%s_summary.json' % (Y, sfx)))]
D = {}
for k, f in runs:
    if f.exists(): s, t = B.load_hybrid(f); D[k] = A.derived(s, t)
ks = list(D)
fin = D[ks[-1]]; tau = fin['H0tau']; lna = fin['ln_a']
sl = np.gradient(np.log(np.abs(fin['Wa4'])), lna)
i_set = next(i for i in range(len(tau)) if tau[i] > (11 if Y == '3' else 13) and abs(sl[i]) < 0.05)
W0 = float(np.median(fin['Wa4'][(tau >= tau[i_set]) & (tau <= tau[i_set] + 1)]))
out = dict(status='numerical', Wa4_settled=W0, pairs={})
def at(d, k, t):
    tt = d['H0tau']; keep = np.concatenate(([True], np.diff(tt) > 1e-12)); return float(np.interp(t, tt[keep], d[k][keep]))
for a, b in zip(ks[:-1], ks[1:]):
    te = min(D[a]['H0tau'][A.reliable_end(D[a])], D[b]['H0tau'][A.reliable_end(D[b])])
    rows = []
    for t in [20, 22, 24, 26, te]:
        if t > te: continue
        ea, eb = at(D[a], 'Wa4', t)/W0 - 1, at(D[b], 'Wa4', t)/W0 - 1
        rows.append(dict(H0tau=float(t), drift_coarse=ea, drift_fine=eb, order=(math.log2(ea/eb) if eb != 0 and ea/eb > 0 else None)))
    o = rows[-1]['order']
    out['pairs'][a + '/' + b] = dict(samples=rows, order_at_last_common_reliable=o, consistent_3_5=bool(o is not None and 3 <= o <= 5))
out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(HERE/('TERTIARY_ORDER.json' if Y == '3' else 'TERTIARY_ORDER_Y%s.json' % Y)).write_text(json.dumps(out, indent=1) + '\n'); print(json.dumps(out, indent=1))
