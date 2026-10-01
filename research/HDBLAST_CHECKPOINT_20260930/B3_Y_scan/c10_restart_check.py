#!/usr/bin/env python3
"""C10 (restart validation): hybrid run 'restart from the dz 1e-3 run onto a 5e-4 grid' versus the full 5e-4 run and the
1e-3 run, Y = 3, dc = 1e-2.  Late W a^4 drift (relative to the settled value of the full 5e-4 run) at H0 tau = 18..24.
PASS if |drift_restart| <= 2 |drift_full5e-4| (or both < 1e-3) and |drift_1e-3| >= 5 |drift_restart| at every sample.
Writes C10_RESTART_CHECK.json."""
import json, sys, hashlib
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import a1_analyze as A, b3_analyze as B
sR, tR = B.load_hybrid(HERE/'runs/ctl/ctl_Y3_dc1e-2_restart1e-3to5e-4_summary.json')
s5, t5 = A.load(HERE/'runs/main/main_Y3_dc1e-2_dzf5e-4_summary.json'); s1, t1 = A.load(HERE/'runs/main/main_Y3_dc1e-2_dzf1e-3_summary.json')
dR, d5, d1 = A.derived(sR, tR), A.derived(s5, t5), A.derived(s1, t1)
trig = json.load(open(HERE/'SECONDARY_TRIGGER.json'))['cases']['main_Y3_dc1e-2']; W0 = trig['Wa4_settled']
def at(d, k, t):
    tt = d['H0tau']; keep = np.concatenate(([True], np.diff(tt) > 1e-12)); return float(np.interp(t, tt[keep], d[k][keep]))
rows = []; ok = True
for t in [18, 20, 22, 24]:
    eR, e5, e1 = at(dR, 'Wa4', t)/W0 - 1, at(d5, 'Wa4', t)/W0 - 1, at(d1, 'Wa4', t)/W0 - 1
    good = (abs(eR) <= 2*abs(e5) or max(abs(eR), abs(e5)) < 1e-3) and abs(e1) >= 5*abs(eR)
    ok &= good; rows.append(dict(H0tau=t, drift_restart_1e3_to_5e4=eR, drift_full_5e4=e5, drift_1e3=e1, ok=bool(good)))
out = dict(status='numerical', restart_info=sR.get('restart_info'), Wa4_settled_ref=W0, samples=rows, PASS=bool(ok),
           script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(HERE/'C10_RESTART_CHECK.json').write_text(json.dumps(out, indent=1) + '\n'); print(json.dumps(out, indent=1))
