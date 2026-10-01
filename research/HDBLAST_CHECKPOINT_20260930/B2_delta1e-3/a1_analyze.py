#!/usr/bin/env python3
"""A1 analysis: applies the registered rule (REGISTRATION.md section 3) and the calibration/control criteria (section 5)
to the saved run time series.  Output: A1_RESULTS.json (and a printed table).  Numerical; every number traces to
runs/**/<tag>_timeseries.npz + <tag>_summary.json produced by evolve_a1.py."""
import json, glob, math, hashlib, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import static_w as S

PLUS_BRANCH_D01 = dict(phi_b=0.9914346595842017, H_over_H0=math.sqrt(0.40616846475555957))   # M2 (27 Sept), delta = 0.1, registered c
PILOT = Path('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/mechanisms/pilot_5d/runs')

def load(path_summary):
    s = json.load(open(path_summary))
    d = np.load(str(path_summary).replace('_summary.json', '_timeseries.npz'))
    cols = list(d['cols']); rec = d['rec']
    ts = {c: rec[:, i] for i, c in enumerate(cols)}
    ts['modes'] = d['modes']
    return s, ts

def derived(s, ts):
    H = ts['H_over_H0']; lna = ts['ln_a']
    out = dict(ts)
    out['Wa4'] = ts['Wy']*np.exp(4*lna)
    out['Ra4'] = ts['R']*np.exp(4*lna)
    with np.errstate(divide='ignore', invalid='ignore'):
        out['Omega_r'] = ts['rad']/H**2
        out['r'] = np.where(ts['rad'] > 0, ts['Wy']/ts['rad'], np.nan)
        out['Omega_vac'] = ts['vac']/H**2
        out['Omega_W'] = ts['Wy']/H**2
    out['closure'] = ts['vac'] + ts['rad'] + ts['kin'] + ts['fric'] + ts['Wy'] - H**2
    return out

def plateau_index(d, dlna=0.5, tol=0.05, need_R=True):
    """first record p such that over ln a in [ln a_p - dlna, ln a_p]: H > 0, and max-min of W a^4 and R a^4 <= tol*|value at p|"""
    lna = d['ln_a']; H = d['H_over_H0']
    Wa4 = d['Wa4']; Ra4 = d['Ra4']
    n = len(lna)
    j0 = 0
    for p in range(n):
        if lna[p] - lna[0] < dlna: continue
        while j0 < p and lna[j0] < lna[p] - dlna: j0 += 1
        if lna[p] < lna[j0] + dlna - 0.02: continue     # ln a not monotone enough to span the window
        win = slice(j0, p + 1)
        if np.any(H[win] <= 0): continue
        w = Wa4[win]; okW = (w.max() - w.min()) <= tol*abs(Wa4[p]) if Wa4[p] != 0 else False
        if need_R:
            rr = Ra4[win]; okR = Ra4[p] > 0 and (rr.max() - rr.min()) <= tol*Ra4[p]
        else:
            okR = True
        if okW and okR: return p, j0
    return None, None

def reliable_end(d, thr=0.05):
    """index of the last record before the near-shell relative constraint residual (checked every 10th record) exceeds thr"""
    hn = d.get('Hmax_near'); mn = d.get('Mmax_near')
    n = len(d['T'])
    if hn is None: return n - 1
    bad = np.flatnonzero((np.nan_to_num(hn, nan=0.0) > thr) | (np.nan_to_num(mn, nan=0.0) > thr))
    return (int(bad[0]) - 1) if len(bad) else n - 1

def truncate(d, iend):
    return {k: (v[:iend + 1] if isinstance(v, np.ndarray) and v.ndim == 1 and len(v) == len(d['T']) else v) for k, v in d.items()}

def classify(s, d_full):
    iend = reliable_end(d_full)
    d = truncate(d_full, iend)
    Y = s['params']['Y']
    H = d['H_over_H0']; tau = d['H0tau']
    res = dict(tag=s['tag'], Y=Y, dc=s['params']['dc'], dzf=s['params']['dzf'], d=s['params']['d'], xc=s['params']['xc'],
               stop_reason=s['stop_reason'], H0tau_end_run=float(d_full['H0tau'][-1]), H0tau_end=float(tau[-1]),
               reliable_end_reason=('near-shell constraint > 0.05' if iend < len(d_full['T']) - 1 else 'end of run'),
               T_end=s['T_end'], runtime_s=s['runtime_s'], chart=s['params'].get('chart_kind'),
               vac_end=float(d['vac'][-1]), phi_b_end=float(d['phi_b'][-1]))
    # recollapse
    neg = np.flatnonzero(H < -0.05)
    good = (d['Omega_r'] >= 0.9) & (np.abs(d['r']) <= 0.1)
    first_good = np.flatnonzero(good)
    res['recollapse_H0tau'] = float(tau[neg[0]]) if len(neg) else None
    i0 = np.flatnonzero(H <= 0)
    res['H_zero_H0tau'] = float(np.interp(0.0, -H[max(i0[0]-1, 0):i0[0]+1], tau[max(i0[0]-1, 0):i0[0]+1])) if len(i0) and i0[0] > 0 else None
    res['max_Omega_r'] = float(np.nanmax(d['Omega_r'])) if Y > 0 else 0.0
    im = int(np.nanargmax(d['Omega_r'])) if Y > 0 else 0
    res['at_max_Omega_r'] = dict(H0tau=float(tau[im]), r=float(d['r'][im]) if Y > 0 else None, H_over_H0=float(H[im]),
                                 Omega_vac=float(d['Omega_vac'][im]), Omega_W=float(d['Omega_W'][im]), phi_b=float(d['phi_b'][im]))
    p, j0 = plateau_index(d, need_R=(Y > 0))
    cls = None
    if p is not None:
        res['plateau'] = dict(H0tau=float(tau[p]), window_start_H0tau=float(tau[j0]), r=float(d['r'][p]) if Y > 0 else None,
                              Omega_r=float(d['Omega_r'][p]) if Y > 0 else 0.0, Wa4=float(d['Wa4'][p]), Ra4=float(d['Ra4'][p]),
                              R_over_sigma=None, H_over_H0=float(H[p]), phi_b=float(d['phi_b'][p]), Omega_vac=float(d['Omega_vac'][p]))
    else:
        res['plateau'] = None
    # rule (registration section 3); recollapse checked first in time order
    if len(neg) and (len(first_good) == 0 or first_good[0] > neg[0]) and (p is None or p > neg[0]):
        cls = 'FAIL-no-radiation-era (recollapse)'
        res['classification_H0tau'] = float(tau[neg[0]])
    elif p is not None:
        om = res['plateau']['Omega_r']; r = res['plateau']['r']
        res['classification_H0tau'] = float(tau[p])
        if Y == 0: cls = 'FAIL-no-radiation-era (Y=0 baseline: no radiation)'
        elif abs(r) > 0.1: cls = 'FAIL-Weyl'
        elif om < 0.5: cls = 'FAIL-no-radiation-era (Omega_r < 0.5 at plateau)'
        elif om >= 0.9: cls = 'PASS-combined' if abs(r) <= 0.03 else 'PASS-conservative'
        else: cls = 'INCONCLUSIVE (0.5 <= Omega_r < 0.9 with |r| <= 0.1)'
    else:
        cls = 'INCONCLUSIVE (no plateau within the evolution)'
        res['classification_H0tau'] = float(tau[-1])
    res['classification'] = cls
    # constraint residuals near the shell up to the classification time
    tc = res['classification_H0tau']
    m = (tau <= tc + 1e-9) & np.isfinite(d.get('Hmax_near', np.full_like(tau, np.nan)))
    if 'Hmax_near' in d and m.any():
        res['max_rel_H_near_shell'] = float(np.nanmax(d['Hmax_near'][m])); res['max_rel_M_near_shell'] = float(np.nanmax(d['Mmax_near'][m]))
        res['max_rel_H_domain'] = float(np.nanmax(d['Hmax'][m])); res['max_rel_M_domain'] = float(np.nanmax(d['Mmax'][m]))
    res['closure_max_abs'] = float(np.nanmax(np.abs(d['closure'])))
    return res

def weyl_identity(s, d, variant='ok'):
    """V2-style Weyl transport identity with matter in H0 units; returns median/p95 relative residual (FD on records)."""
    p = s['params']; ten = S.Tension(p['delta'], p['c'], p['d']); Y = p['Y']; rb = s['rho_b']
    tau = d['H0tau']; keep = np.concatenate(([True], np.diff(tau) > 1e-9))
    tau = tau[keep]; phi = d['phi_b'][keep]; h = d['H_over_H0'][keep]; vh = d['v_over_H0'][keep]; R = d['R'][keep]; Wy = d['Wy'][keep]
    k = rb*(ten.s(phi) + (R if variant != 'noR' else 0))/6
    wv = -(rb*ten.s1(phi) + Y*vh)/2
    vdot = np.gradient(vh, tau, edge_order=2)
    wdot = -(rb*ten.s2(phi)*vh + Y*vdot)/2
    dW = np.gradient(Wy, tau, edge_order=2)
    fac = 3 if variant == '3H' else 4
    rhs = (4*k*wv*vh - 4*h*vh*vh - vh*vdot + wv*wdot - rb**2*S.U1(phi)*vh)/6
    lhs = dW + fac*h*Wy
    scale = (np.abs(4*k*wv*vh) + np.abs(4*h*vh*vh) + np.abs(vh*vdot) + np.abs(wv*wdot) + np.abs(rb**2*S.U1(phi)*vh))/6 + np.abs(dW) + np.abs(fac*h*Wy) + 1e-300
    r = (np.abs(lhs - rhs)/scale)[3:-3]
    return dict(median=float(np.median(r)), p95=float(np.percentile(r, 95)))

def ledger(s, d):
    p = s['params']; Y = p['Y']; rb = s['rho_b']
    if Y == 0: return None
    tau = d['H0tau']; keep = np.concatenate(([True], np.diff(tau) > 1e-9))
    t = tau[keep]*rb; R = d['R'][keep]; H = d['H_over_H0'][keep]/rb; v = d['v_over_H0'][keep]/rb
    dR = np.gradient(R, t, edge_order=2)
    L = dR + 4*H*R - Y*v*v; sc = np.abs(dR) + np.abs(4*H*R) + np.abs(Y*v*v) + 1e-300
    r = (np.abs(L)/sc)[3:-3]
    return dict(median=float(np.median(r)), p95=float(np.percentile(r, 95)))

def growth_rate(s, d, target):
    """fit ln|phi_b - phi_b(static bg)| against H0 tau in the window 30 x initial deviation < dev < 2e-3, tau > 1"""
    dev = np.abs(d['phi_b'] - s['phi_b_static']); tau = d['H0tau']
    m = (dev > 30*dev[0]) & (dev < 2e-3) & (tau > 1.0)
    if m.sum() < 20: return None
    sl = np.polyfit(tau[m], np.log(dev[m]), 1)[0]
    # sliding windows of width 1.5
    sw = []
    t0 = tau[m][0]
    while t0 + 1.5 <= tau[m][-1] + 1e-9:
        mm = m & (tau >= t0) & (tau <= t0 + 1.5)
        if mm.sum() > 20: sw.append(float(np.polyfit(tau[mm], np.log(dev[mm]), 1)[0]))
        t0 += 0.5
    return dict(rate=float(sl), window=[float(tau[m][0]), float(tau[m][-1])], sliding=sw, target=target, diff=float(sl - target))

def interp_at(d, key, tau_pts):
    tau = d['H0tau']; keep = np.concatenate(([True], np.diff(tau) > 1e-12))
    return np.interp(tau_pts, tau[keep], d[key][keep])

def compare(dA, dB, keys, tau_pts):
    return {k: float(np.max(np.abs(interp_at(dA, k, tau_pts) - interp_at(dB, k, tau_pts)))) for k in keys}

def main():
    out = dict(status='numerical', registration='REGISTRATION.md', runs={}, calibration={}, controls={}, convergence={}, aggregate={})
    summaries = sorted(glob.glob(str(HERE/'runs/**/*_summary.json'), recursive=True))
    data = {}
    for f in summaries:
        s, ts = load(f); d = derived(s, ts); tag = Path(f).name.replace('_summary.json', '')
        data[tag] = (s, d)
    # ---------------- calibration
    cal = {}
    if 'cal_reg_Y0_dc1e-2_xcinf' in data:
        s, d = data['cal_reg_Y0_dc1e-2_xcinf']
        P = np.load(PILOT/'pilot_Y0_dc1e-2_timeseries.npz')['rec']
        pts = np.array([1.0, 2.0, 3.0, 4.0])
        tp = P[:, 8]; keep = np.concatenate(([True], np.diff(tp) > 1e-12))
        ph_p = np.interp(pts, tp[keep], P[keep, 1]); H_p = np.interp(pts, tp[keep], P[keep, 2])
        ph = interp_at(d, 'phi_b', pts); H = interp_at(d, 'H_over_H0', pts)
        cal['C1_old_chart_vs_pilot'] = dict(H0tau=pts.tolist(), phi_b_new=ph.tolist(), phi_b_pilot=ph_p.tolist(), H_new=H.tolist(), H_pilot=H_p.tolist(),
                                            max_abs_dphi=float(np.max(np.abs(ph - ph_p))), max_abs_dH=float(np.max(np.abs(H - H_p))),
                                            PASS=bool(np.max(np.abs(ph - ph_p)) <= 1e-3 and np.max(np.abs(H - H_p)) <= 1e-3))
    for tag, target in [('cal_growth_d0.1', 1.62702), ('cal_growth_d0.001', 1.65719)]:
        if tag in data:
            g = growth_rate(*data[tag], target)
            tol = 1e-3 if 'd0.1' in tag else 2e-4
            cal['C2_' + tag] = dict(g or {}, tol=tol, PASS=bool(g is not None and abs(g['diff']) <= tol))
    # chart independence
    ci = [t for t in data if t.startswith('cal_reg_Y0_dc1e-2_xc') and not t.endswith('xcinf')]
    if 'cal_reg_Y0_dc1e-2_xcinf' in data and ci:
        dinf = data['cal_reg_Y0_dc1e-2_xcinf'][1]
        pre = np.linspace(0.5, 4.3, 39)
        c3 = dict(before_freeze={})
        for t in ci:
            c3['before_freeze'][t] = compare(dinf, data[t][1], ['phi_b', 'H_over_H0', 'Wy'], pre)
        c3['after_freeze'] = {}
        for i_, a in enumerate(sorted(ci)):
            for b in sorted(ci)[i_ + 1:]:
                da = truncate(data[a][1], reliable_end(data[a][1])); db = truncate(data[b][1], reliable_end(data[b][1]))
                e = min(da['H0tau'][-1], db['H0tau'][-1])
                if e > 4.6:
                    post = np.linspace(4.3, e - 0.05, 60)
                    c3['after_freeze'][a + ' vs ' + b] = dict(H0tau_range=[4.3, float(e - 0.05)], diffs=compare(da, db, ['phi_b', 'H_over_H0', 'Wy'], post))
        bf = max(max(v.values()) for v in c3['before_freeze'].values())
        c3['PASS_before'] = bool(bf <= 1e-3)
        if c3['after_freeze']: c3['PASS_after'] = bool(max(max(v['diffs'].values()) for v in c3['after_freeze'].values()) <= 1e-2)
        cal['C3_chart_independence'] = c3
    # known end state
    for t in ci:
        s, d = data[t]; d = truncate(d, reliable_end(d))
        cal.setdefault('C4_end_state', {})[t] = dict(H0tau_end=float(d['H0tau'][-1]), H_over_H0_end=float(d['H_over_H0'][-1]),
                                                     phi_b_end=float(d['phi_b'][-1]), Wy_end=float(d['Wy'][-1]),
                                                     plus_branch=PLUS_BRANCH_D01,
                                                     dH=float(d['H_over_H0'][-1] - PLUS_BRANCH_D01['H_over_H0']),
                                                     dphi=float(d['phi_b'][-1] - PLUS_BRANCH_D01['phi_b']))
    out['calibration'] = cal
    # ---------------- tuned runs
    for tag, (s, d) in data.items():
        if not tag.startswith('main_') and not tag.startswith('ctl_'): continue
        r = classify(s, d)
        r['weyl_identity'] = weyl_identity(s, d); r['weyl_identity_control_3H'] = weyl_identity(s, d, '3H')
        r['weyl_identity_control_noR'] = weyl_identity(s, d, 'noR') if s['params']['Y'] > 0 else None
        r['radiation_ledger'] = ledger(s, d)
        (out['runs'] if tag.startswith('main_') else out['controls'])[tag] = r
    # ---------------- convergence pairs (dzf 1e-3 vs 5e-4) and aggregate
    main = out['runs']
    groups = {}
    for tag, r in main.items():
        key = (r['Y'], r['dc'], round(r['d'], 6))
        groups.setdefault(key, {})[r['dzf']] = tag
    agg = {}
    for key, g in sorted(groups.items()):
        entry = dict(Y=key[0], dc=key[1], d=key[2], runs=g)
        if len(g) >= 2:
            fine = main[g[min(g)]]; coarse = main[g[max(g)]]
            same = fine['classification'].split(' ')[0] == coarse['classification'].split(' ')[0]
            agree = None
            if fine.get('plateau') and coarse.get('plateau') and fine['plateau']['r'] is not None and coarse['plateau']['r'] is not None:
                rf, rc = fine['plateau']['r'], coarse['plateau']['r']
                agree = abs(rf - rc) <= max(0.2*abs(rf), 0.02)
            elif 'recollapse' in fine['classification'] and 'recollapse' in coarse['classification']:
                agree = abs(fine['recollapse_H0tau'] - coarse['recollapse_H0tau']) <= 0.05
            cons_ok = fine.get('max_rel_H_near_shell', 1) < 0.05 and fine.get('max_rel_M_near_shell', 1) < 0.05
            entry.update(fine=fine['classification'], coarse=coarse['classification'], same_class=same, value_agreement=agree,
                         constraints_ok=cons_ok, reliable=bool(same and agree and cons_ok),
                         class_reliable=(fine['classification'] if (same and agree and cons_ok) else 'UNRELIABLE'))
        else:
            entry.update(reliable=False, class_reliable='UNRELIABLE (single resolution)')
        agg['Y=%g dc=%g d=%.4f' % key] = entry
    out['convergence'] = agg
    # aggregate verdict at delta = 0.1 (registered rule)
    verdict = 'INCONCLUSIVE'
    Ys = sorted(set(k[0] for k in groups if k[0] > 0))
    per_Y = {}
    for Y in Ys:
        cls = [agg['Y=%g dc=%g d=%.4f' % k]['class_reliable'] for k in groups if k[0] == Y]
        per_Y[Y] = cls
    if any(len(c) >= 2 and all(x.startswith('PASS') for x in c) for c in per_Y.values()): verdict = 'PASS'
    elif per_Y and all(len(c) >= 2 and all(x.startswith('FAIL') for x in c) for c in per_Y.values()): verdict = 'FAIL'
    out['aggregate'] = dict(per_Y=per_Y, verdict_delta_0p1=verdict, Y_values=Ys)
    out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    def clean(v):
        if isinstance(v, dict): return {str(k): clean(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)): return [clean(x) for x in v]
        if isinstance(v, (np.floating, float)): return None if not math.isfinite(float(v)) else float(v)
        if isinstance(v, (np.integer,)): return int(v)
        if isinstance(v, np.bool_): return bool(v)
        return v
    (HERE/'A1_RESULTS.json').write_text(json.dumps(clean(out), indent=1) + '\n')
    for tag, r in sorted(main.items()):
        pl = r.get('plateau') or {}
        print('%-40s %-45s r=%s Om_r=%s recoll=%s Hzero=%s maxOm=%.3f end=%.2f %s' % (tag, r['classification'], pl.get('r'), pl.get('Omega_r'),
              r['recollapse_H0tau'], r['H_zero_H0tau'], r['max_Omega_r'], r['H0tau_end'], r['stop_reason']))
    print(json.dumps(out['aggregate'], indent=1))
    print(json.dumps(clean(out['calibration']), indent=1)[:3000])

if __name__ == '__main__':
    main()
