#!/usr/bin/env python3
"""REGISTRATION section 4 trigger for the secondary (restart) pair: for each (Y, dc) with an unreliable primary pair, the
apparent convergence order of the late W a^4 drift between dz_fine 1e-3 and 5e-4, measured against the settled value of the
5e-4 run (first time after friction-off with |dln(W a^4)/dln a| < 0.05, median over the next 1 H0^-1), at the last time where
both runs are reliable (min of the two reliable ends; the coarse run has no reliable data beyond its own).
Writes SECONDARY_TRIGGER.json."""
import json, sys, math, glob, hashlib
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import a1_analyze as A
out = dict(status='numerical', cases={})
for fc in sorted([f for f in glob.glob(str(HERE/'runs/main/main_Y*_dzf1e-3*_summary.json')) if '_ext_' not in f]):
    ff = fc.replace('dzf1e-3', 'dzf5e-4')
    if Path(ff.replace('_summary.json', '_ext_summary.json')).exists(): ff = ff.replace('_summary.json', '_ext_summary.json')   # fourth dated note
    if not Path(ff).exists(): continue
    import b3_analyze as B0
    sc, tc = A.load(fc); sf, tf = B0.load_hybrid(ff); dc_ = A.derived(sc, tc); df = A.derived(sf, tf)
    rc, rf = A.classify(sc, dc_), A.classify(sf, df)
    rc['tag'] = 'c'; rf['tag'] = 'f'
    import b3_analyze as B
    prim_rel = B.pair_reliability(rf, rc)['reliable']
    key = Path(fc).name.replace('_summary.json', '').replace('_dzf1e-3', '')
    tau = df['H0tau']; lna = df['ln_a']; Y = sf['params']['Y']; H = df['H_over_H0']
    with np.errstate(divide='ignore', invalid='ignore'):
        src = np.nan_to_num(np.abs(Y*df['v_over_H0']**2/(4*H*df['R'])), nan=1.0)
    iH = int(np.argmax(H <= 0)) if np.any(H <= 0) else len(H)
    i_off = next((i for i in range(iH) if (src[i:iH] < 1e-3).all()), None)
    if i_off is None: out['cases'][key] = dict(note='friction never switched off before H = 0'); continue
    sl = np.gradient(np.log(np.abs(df['Wa4'])), lna)
    i_set = next((i for i in range(i_off, iH) if abs(sl[i]) < 0.05), None)
    if i_set is None: out['cases'][key] = dict(note='W a^4 never settled before H = 0'); continue
    m = (tau >= tau[i_set]) & (tau <= tau[i_set] + 1); W0 = float(np.median(df['Wa4'][m]))
    t_eval = min(dc_['H0tau'][A.reliable_end(dc_)], tau[A.reliable_end(df)])
    if rf.get('plateau'): t_eval = min(t_eval, rf['plateau']['H0tau'])
    def at(d, k):
        t = d['H0tau']; keep = np.concatenate(([True], np.diff(t) > 1e-12)); return float(np.interp(t_eval, t[keep], d[k][keep]))
    ec, ef = at(dc_, 'Wa4')/W0 - 1, at(df, 'Wa4')/W0 - 1
    order = math.log2(ec/ef) if ef != 0 and ec/ef > 0 else None
    out['cases'][key] = dict(primary_classes=[rc['classification'], rf['classification']], t_eval_H0tau=float(t_eval), Wa4_settled=W0,
                             drift_1e3=ec, drift_5e4=ef, order=order, primary_reliable=prim_rel, order_in_3_5=bool(order is not None and 3 <= order <= 5),
                             trigger=bool((not prim_rel) and order is not None and 3 <= order <= 5))
out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(HERE/'SECONDARY_TRIGGER.json').write_text(json.dumps(out, indent=1) + '\n'); print(json.dumps(out, indent=1))
