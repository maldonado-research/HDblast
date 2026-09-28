#!/usr/bin/env python3
"""Audit script 2b: repeat the trajectory-like out-of-sample toy (mix3) with a shorter span.
With span 1.4 the cubic Delta(s) has a second zero at s=-1.46, just outside the window, so the vacuum was set in a
non-adiabatic region (a defect of the audit toy, not of the audited code).  Span 1.0 avoids it.
Three q values allow a check of the O(1/q^2) remainder.  Writes verify/INDEP_MIX3.json."""
import json
from pathlib import Path
import indep_modes as I
HERE = Path(__file__).resolve().parent
p = dict(v=0.5786840415164352, A=0.2863880601008191, B=-1.0512888681550976, h0=0.8878608255431975, h1=-0.19628251354630122)
out = {'params': p, 'closed_form': I.closed_form(**p), 'runs': {}}
for span in (1.0, 0.8):
    rows = []
    for qt in (2400.0, 4800.0, 9600.0):
        r = I.number(I.Toy(span=span, **p), qt/p['v'], eps=1e-4)
        rows.append({'q': r['q'], 'D_q': r['rel']*r['q'], 'edges': r['edges'], 'wronskian_dev': r['wronskian_dev']})
    (q1, d1), (q2, d2), (q3, d3) = [(x['q'], x['D_q']) for x in rows]
    R12 = (q2*d2 - q1*d1)/(q2 - q1); R23 = (q3*d3 - q2*d2)/(q3 - q2)
    out['runs'][str(span)] = {'rows': rows, 'richardson_12': R12, 'richardson_23': R23,
                              'dev_closed_23': R23 - out['closed_form']}
    print(span, out['runs'][str(span)], flush=True)
(HERE/'INDEP_MIX3.json').write_text(json.dumps(out, indent=1) + '\n')
