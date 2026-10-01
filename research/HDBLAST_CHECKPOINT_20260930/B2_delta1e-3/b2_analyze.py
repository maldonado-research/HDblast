#!/usr/bin/env python3
"""B2 analysis (copy of the A1 a1_analyze.py; classification, plateau rule, Weyl identity and ledger functions unchanged).
Applies the rule of REGISTRATION.md (= A1 section 3) and the B2 calibration/control criteria to runs/**.  Output: B2_RESULTS.json.
Numerical; every number traces to runs/**/<tag>_timeseries.npz + <tag>_summary.json produced by evolve_a1.py (B2 copy)."""
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


def growth_rate_chat14(s, d, target=1.65719, key='H0tau'):
    """B2 registered C2 estimator: log-slope of |phi_b - phi_b(static)| against H0 tau on [a, a+1], a in {2, 2.5, 3, 3.5, 4},
    records with |dev| < 1e-3 only; rate = mean, spread = max - min."""
    dev = np.abs(d['phi_b'] - s['phi_b_static']); t = d[key]
    rates = {}
    for a in [2.0, 2.5, 3.0, 3.5, 4.0]:
        m = (t >= a) & (t <= a + 1) & (dev < 1e-3) & (dev > 0)
        rates['%.1f-%.1f' % (a, a + 1)] = float(np.polyfit(t[m], np.log(dev[m]), 1)[0]) if m.sum() > 20 else None
    v = [x for x in rates.values() if x is not None]
    if len(v) < 5: return dict(windows=rates, rate=None, spread=None, PASS=False, note='fewer than 5 windows in the linear regime')
    rate = float(np.mean(v)); spread = float(max(v) - min(v))
    return dict(windows=rates, rate=rate, spread=spread, diff=rate - target, target=target,
                PASS=bool(abs(rate - target) <= 2e-4 and spread <= 2e-4))

def growth_rate_late(s, d, target=1.65719, key='H0tau', starts=(4.5, 5.0, 5.5, 6.0)):
    """SUPPLEMENTARY (dated note 1, post-registration): same fit as growth_rate_chat14 on later windows [a, a+1]."""
    dev = np.abs(d['phi_b'] - s['phi_b_static']); t = d[key]
    rates = {}
    for a in starts:
        m = (t >= a) & (t <= a + 1) & (dev < 1e-3) & (dev > 0)
        rates['%.1f-%.1f' % (a, a + 1)] = float(np.polyfit(t[m], np.log(dev[m]), 1)[0]) if m.sum() > 20 else None
    v = [x for x in rates.values() if x is not None]
    if len(v) < len(starts): return dict(windows=rates, rate=None, spread=None, PASS=False, note='fewer windows than required in the linear regime')
    rate = float(np.mean(v)); spread = float(max(v) - min(v))
    return dict(windows=rates, rate=rate, spread=spread, diff=rate - target, target=target,
                PASS=bool(abs(rate - target) <= 2e-4 and spread <= 2e-4))

def transient_check(s, d, target=1.65719):
    """window rates on [a, a+1], a = 1.0, 1.5, ..., 6.0 (|dev| < 1e-3): monotone approach and geometric shrinking of the gaps"""
    dev = np.abs(d['phi_b'] - s['phi_b_static']); t = d['H0tau']
    ws = []
    for a in np.arange(1.0, 6.01, 0.5):
        m = (t >= a) & (t <= a + 1) & (dev < 1e-3) & (dev > 0)
        if m.sum() > 20: ws.append((float(a), float(np.polyfit(t[m], np.log(dev[m]), 1)[0])))
    r = np.array([w[1] for w in ws])
    gaps = np.abs(np.diff(r))
    ratios = (gaps[1:]/gaps[:-1]).tolist() if len(gaps) > 1 else []
    mono = bool(np.all(np.diff(r) > 0) or np.all(np.diff(r) < 0))
    return dict(windows={'%.1f-%.1f' % (a, a + 1): v for a, v in ws}, monotone=mono, gap_ratios=ratios,
                last_window_rate=float(r[-1]) if len(r) else None,
                aitken_extrapolated=(float(r[-1] - gaps[-1]**2/(gaps[-1] - gaps[-2])*np.sign(r[-1]-r[-2])) if len(gaps) > 1 and gaps[-1] != gaps[-2] else None))

def main():
    out = dict(status='numerical', registration='REGISTRATION.md', calibration={}, pre_runs={}, runs={}, controls={}, exploratory={}, convergence={}, aggregate={})
    summaries = sorted(glob.glob(str(HERE/'runs/**/*_summary.json'), recursive=True))
    data = {}
    for f in summaries:
        s, ts = load(f); d = derived(s, ts); tag = Path(f).name.replace('_summary.json', '')
        data[tag] = (s, d, Path(f).parent.name)
    # ---------------- C2
    cal = {}
    for tag, (s, d, grp) in data.items():
        if grp not in ('cal', 'diag'): continue
        g = growth_rate_chat14(s, d); gt = growth_rate_chat14(s, d, key='T')
        a1 = growth_rate(s, d, 1.65719)
        dev = np.abs(d['phi_b'] - s['phi_b_static'])
        early = {}
        for a in [1.0, 1.5]:
            m = (d['H0tau'] >= a) & (d['H0tau'] <= a + 1) & (dev < 1e-3)
            early['%.1f-%.1f' % (a, a + 1)] = float(np.polyfit(d['H0tau'][m], np.log(dev[m]), 1)[0]) if m.sum() > 20 else None
        cal[tag] = dict(registered_estimator=g, chart_time_version=gt, A1_estimator=a1, early_windows=early,
                        late_estimator_supplementary=growth_rate_late(s, d), transient_check=transient_check(s, d), T_final=s['T_end'],
                        dc=s['params']['dc'], dzf=s['params']['dzf'], grid=s['params'].get('grid_kind'), table_dx=s['params'].get('table_dx'),
                        points_per_wall=s.get('points_per_wall'), n_points=s['n_points'], runtime_s=s['runtime_s'])
    c2 = {k: v for k, v in cal.items() if k in ('C2a_rho1.7e-4', 'C2b_rho9.5e-5')}
    cal['C2_verdict'] = dict(runs=sorted(c2), PASS=bool(len(c2) == 2 and all(v['registered_estimator']['PASS'] for v in c2.values())),
                             note='registered estimator on C2a (T_final 6.5) and C2b, exactly as registered')
    c2s = {k: v for k, v in cal.items() if k in ('C2a8_rho1.7e-4', 'C2b_rho9.5e-5')}
    cal['C2_supplementary_late_verdict'] = dict(runs=sorted(c2s), label='post-registration (dated note 1, amendment 2)',
        PASS=bool(len(c2s) == 2 and all(v['late_estimator_supplementary']['PASS'] for v in c2s.values())))
    if 'C2a_rho1.7e-4' in cal and 'C2c_rho1.7e-4_dc3e-8' in cal:
        ra, rc = cal['C2a_rho1.7e-4']['registered_estimator']['rate'], cal['C2c_rho1.7e-4_dc3e-8']['registered_estimator']['rate']
        cal['C2c_seed_independence'] = dict(diff=(None if ra is None or rc is None else rc - ra), expected_within=1e-4)
    out['calibration'] = cal
    # ---------------- pre-runs
    for tag, (s, d, grp) in data.items():
        if grp == 'pre':
            out['pre_runs'][tag] = dict(xc_suggest=s.get('xc_suggest'), stop_reason=s['stop_reason'], T_end=s['T_end'], H0tau_end=s['H0tau_end'],
                                        grid=s['params'].get('grid_kind'), dzf=s['params']['dzf'], dzcap=s['params'].get('dzcap'), runtime_s=s['runtime_s'])
    # ---------------- tuned runs
    for tag, (s, d, grp) in data.items():
        if grp not in ('main', 'ctl', 'explore', 'res'): continue
        r = classify(s, d)
        r['weyl_identity'] = weyl_identity(s, d); r['weyl_identity_control_3H'] = weyl_identity(s, d, '3H')
        r['weyl_identity_control_noR'] = weyl_identity(s, d, 'noR') if s['params']['Y'] > 0 else None
        r['radiation_ledger'] = ledger(s, d)
        r['points_per_wall'] = s.get('points_per_wall'); r['n_points'] = s['n_points']; r['dzcap'] = s['params'].get('dzcap')
        r['kappa'] = s['params']['kappa']
        if r.get('plateau') and s['params']['Y'] > 0:
            p = r['plateau']
            # R/sigma at the plateau (flag for the low-energy identification)
            i = int(np.argmin(np.abs(d['H0tau'] - p['H0tau'])))
            ten = S.Tension(s['params']['delta'], s['params']['c'], s['params']['d'])
            p['R_over_sigma'] = float(d['R'][i]/ten.s(d['phi_b'][i]))
        {'main': out['runs'], 'ctl': out['controls'], 'explore': out['exploratory'], 'res': out['exploratory']}[grp][tag] = r
        r['chart_s_max'] = s['params'].get('s_max', 20.0)
    main_ = out['runs']
    groups = {}
    for tag, r in main_.items():
        key = (r['Y'], r['dc'])
        groups.setdefault(key, {})[r['dzf']] = tag
    agg = {}
    for key, g in sorted(groups.items()):
        entry = dict(Y=key[0], dc=key[1], runs=g)
        if len(g) >= 2:
            fine = main_[g[min(g)]]; coarse = main_[g[max(g)]]
            same = fine['classification'].split(' ')[0] == coarse['classification'].split(' ')[0]
            agree = None
            if fine.get('plateau') and coarse.get('plateau') and fine['plateau']['r'] is not None and coarse['plateau']['r'] is not None:
                rf, rc = fine['plateau']['r'], coarse['plateau']['r']
                agree = abs(rf - rc) <= max(0.2*abs(rf), 0.02)
            elif 'recollapse' in fine['classification'] and 'recollapse' in coarse['classification']:
                agree = abs(fine['recollapse_H0tau'] - coarse['recollapse_H0tau']) <= 0.05
            cons_ok = fine.get('max_rel_H_near_shell', 1) < 0.05 and fine.get('max_rel_M_near_shell', 1) < 0.05
            decisive = not fine['classification'].startswith('INCONCLUSIVE')
            entry.update(fine=fine['classification'], coarse=coarse['classification'], same_class=same, value_agreement=agree,
                         constraints_ok=cons_ok, reliable=bool(same and agree and cons_ok and decisive),
                         class_reliable=(fine['classification'] if (same and agree and cons_ok and decisive) else
                                         ('INCONCLUSIVE' if not decisive else 'UNRELIABLE')))
        else:
            entry.update(reliable=False, class_reliable='UNRELIABLE (single resolution)')
        agg['Y=%g dc=%g' % key] = entry
    out['convergence'] = agg
    y1 = [agg[k]['class_reliable'] for k in agg if agg[k]['Y'] == 1]
    verdict = 'INCONCLUSIVE'
    if len(y1) == 2 and all(x.startswith('PASS') for x in y1): verdict = 'PASS'
    elif len(y1) == 2 and all(x.startswith('FAIL') for x in y1): verdict = 'FAIL'
    out['aggregate'] = {'Y1_classes': y1, 'verdict_delta_1e-3': verdict, 'C2_pass': cal['C2_verdict']['PASS']}
    out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    def clean(v):
        if isinstance(v, dict): return {str(k): clean(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)): return [clean(x) for x in v]
        if isinstance(v, (np.floating, float)): return None if not math.isfinite(float(v)) else float(v)
        if isinstance(v, (np.integer,)): return int(v)
        if isinstance(v, np.bool_): return bool(v)
        return v
    (HERE/'B2_RESULTS.json').write_text(json.dumps(clean(out), indent=1) + '\n')
    for tag, v in sorted(cal.items()):
        if 'registered_estimator' in v:
            print('%-34s rate=%s spread=%s PASS=%s  A1-est=%s early=%s' % (tag, v['registered_estimator'].get('rate'), v['registered_estimator'].get('spread'),
                  v['registered_estimator'].get('PASS'), (v['A1_estimator'] or {}).get('rate'), v['early_windows']))
    print(json.dumps(clean(out['pre_runs']), indent=0))
    for tag, r in sorted(list(main_.items()) + list(out['controls'].items()) + list(out['exploratory'].items())):
        pl = r.get('plateau') or {}
        print('%-32s %-45s r=%s Om_r=%s recoll=%s maxOm=%.3f end=%.2f rel_end=%s %s' % (tag, r['classification'], pl.get('r'), pl.get('Omega_r'),
              r['recollapse_H0tau'], r['max_Omega_r'], r['H0tau_end'], r['reliable_end_reason'], r['stop_reason']))
    print(json.dumps(clean(out['aggregate']), indent=1))

if __name__ == '__main__':
    main()
