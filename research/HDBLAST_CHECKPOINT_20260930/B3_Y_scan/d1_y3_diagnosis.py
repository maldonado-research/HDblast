#!/usr/bin/env python3
"""B3 step 1: why did the A1 Y = 3 runs not reach the registered plateau?  Reads the saved A1 time series (read-only) and
writes D1_Y3_DIAGNOSIS.json.  Numerical diagnostics only (no classification change)."""
import json, sys, math, hashlib
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import a1_analyze as A
BASE = Path('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum/runs/main')
out = dict(status='numerical', source='A1 runs/main/main_dstar_Y3_* (read-only)', runs={}, resolution_pairs={})
def at(d, key, tau0):
    tau = d['H0tau']; keep = np.concatenate(([True], np.diff(tau) > 1e-12)); return float(np.interp(tau0, tau[keep], d[key][keep]))
D = {}
for dc in ['1e-2', '1e-4']:
    for dz in ['1e-3', '5e-4']:
        tag = 'main_dstar_Y3_dc%s_dzf%s' % (dc, dz)
        s, ts = A.load(BASE/(tag + '_summary.json')); d = A.derived(s, ts); D[tag] = (s, d)
        Y = s['params']['Y']; rb = s['rho_b']
        tau = d['H0tau']; H = d['H_over_H0']; lna = d['ln_a']; v = d['v_over_H0']; R = d['R']
        with np.errstate(divide='ignore', invalid='ignore'):
            src_ratio = (Y*v*v)/(4*H*R)          # (source of R a^4)/(dilution), H0 units cancel (Y v^2 / (4 H R), same units)
        # friction off: first time after which src_ratio < 1e-3 stays
        ok = np.nan_to_num(np.abs(src_ratio), nan=1.0) < 1e-3
        i_off = None
        for i in range(len(ok)):
            if ok[i] and ok[i:np.argmax(H <= 0) if np.any(H <= 0) else len(ok)].all(): i_off = i; break
        iH0 = int(np.argmax(H <= 0)) if np.any(H <= 0) else None
        ie = A.reliable_end(d)
        # settled W a^4: median over tau in [i_off, i_off + 3 H0^-1]
        t_off = float(tau[i_off])
        # settled: first record after friction-off where the local logarithmic slope |d ln(W a^4)/d ln a| < 0.05
        sl_ = np.gradient(np.log(np.abs(d['Wa4'])), lna)
        i_set = next(i for i in range(i_off, len(tau)) if abs(sl_[i]) < 0.05)
        t_set = float(tau[i_set]); m = (tau >= t_set) & (tau <= t_set + 1.0)
        W_set = float(np.median(d['Wa4'][m]))
        drift = np.abs(d['Wa4']/W_set - 1)
        i5 = next((i for i in range(i_off, len(tau)) if tau[i] > t_set and drift[i] > 0.05), None)
        # Weyl identity residual restricted to the late window (friction off -> reliable end)
        sl = slice(i_off, ie + 1)
        dd = {k: (v_[sl] if isinstance(v_, np.ndarray) and v_.ndim == 1 and len(v_) == len(tau) else v_) for k, v_ in d.items()}
        wi = A.weyl_identity(s, dd)
        out['runs'][tag] = dict(xc=s['params']['xc'], stop_reason=s['stop_reason'], H0tau_end_run=float(tau[-1]),
            reliable_end_H0tau=float(tau[ie]), reliable_end_reason='near-shell constraint residual > 0.05 (checked every 10th record)' if ie < len(tau) - 1 else 'end',
            H_zero_H0tau=float(tau[iH0]) if iH0 else None, ln_a_max=float(lna.max()),
            friction_off=dict(criterion='Y v^2/(4 H R) < 1e-3 for all later times with H > 0', H0tau=t_off, ln_a=float(lna[i_off]), lapse=float(d['lapse'][i_off])),
            e_folds_between_friction_off_and_turnaround=float(lna.max() - lna[i_off]),
            Wa4_settled=W_set, r_settled=float(np.median(d['r'][m])), Omega_r_settled=float(np.median(d['Omega_r'][m])), Omega_vac_settled=float(np.median(d['Omega_vac'][m])),
            first_5pct_drift_of_Wa4=dict(H0tau=float(tau[i5]) if i5 else None, ln_a=float(lna[i5]) if i5 else None, lapse=float(d['lapse'][i5]) if i5 else None,
                                         e_folds_after_settling=float(lna[i5] - lna[i_set]) if i5 else None),
            settled=dict(criterion='first time after friction-off with |dln(W a^4)/dln a| < 0.05; W a^4 = median over the next 1 H0^-1', H0tau=t_set, ln_a=float(lna[i_set]), lapse=float(d['lapse'][i_set]), e_folds_settled_to_turnaround=float(lna.max() - lna[i_set])),
            Ra4_change_after_friction_off=float(np.max(np.abs(d['Ra4'][i_off:ie]/d['Ra4'][i_off] - 1))),
            weyl_identity_late_window=wi, vac_end=float(d['vac'][ie]))
    a, b = 'main_dstar_Y3_dc%s_dzf1e-3' % dc, 'main_dstar_Y3_dc%s_dzf5e-4' % dc
    t0 = out['runs'][b]['settled']['H0tau'] - 2; W0 = out['runs'][b]['Wa4_settled']
    rows = []
    for tt in np.arange(t0 + 2, min(D[a][1]['H0tau'][A.reliable_end(D[a][1])], D[b][1]['H0tau'][A.reliable_end(D[b][1])]), 2.0):
        ea = at(D[a][1], 'Wa4', tt)/W0 - 1; eb = at(D[b][1], 'Wa4', tt)/W0 - 1
        rows.append(dict(H0tau=float(tt), lapse=at(D[b][1], 'lapse', tt), drift_dzf1e3=ea, drift_dzf5e4=eb,
                         ratio=(ea/eb if abs(eb) > 1e-3 else None), order=(math.log2(abs(ea/eb)) if abs(eb) > 1e-3 and ea/eb > 0 else None)))
    out['resolution_pairs']['dc' + dc] = rows
out['conclusions'] = [
 'Friction is OFF well before the runs end (Y v^2/(4HR) < 1e-3 from H0tau ~ 8.5 (dc 1e-2) / 13.5 (dc 1e-4)); R a^4 is then constant to < 1e-3.',
 'W a^4 settles after the source switches off (r ~ 0.010) and then drifts upward; the drift is resolution dependent (ratio between dzf 1e-3 and 5e-4 ~ 2^3.5-2^4, i.e. the 4th-order discretisation error) and grows with the shell lapse of the bounded chart (lapse 10 -> 50).  The Weyl transport identity, satisfied to ~1e-5 in median, fails in the late window (p95 0.1-0.3): the drift is numerical, not physical (the identity source terms are all proportional to v ~ 1e-7 there).',
 'Physically available plateau window: between friction switch-off and turnaround (H = 0, caused by the tuned negative residual vacuum ~ -3.5e-4 H0^2) there are only ~0.75 e-folds, so the registered 0.5 e-fold window must end close to the turnaround, where the lapse is ~20-30 and the late drift exceeds 5% at dzf 1e-3 (and is borderline at 5e-4).',
 "The runs' 'non-finite' stop happens after the turnaround (H < 0) and is irrelevant; the A1 'H0tau 31-45' is the reliable end (near-shell constraint residual > 0.05)."]
out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(HERE/'D1_Y3_DIAGNOSIS.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
print(json.dumps(out, indent=1, default=float))
