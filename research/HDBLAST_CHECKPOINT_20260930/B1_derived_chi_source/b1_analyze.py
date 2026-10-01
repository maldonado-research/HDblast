#!/usr/bin/env python3
"""B1 analysis: applies the registered rule (REGISTRATION.md section 2-4, the A1 rule with R = decay radiation) and the
calibration/ledger criteria to the saved runs.  Output: B1_RESULTS.json and a printed table.  Numerical; every number traces
to runs/**/<tag>_timeseries.npz + <tag>_summary.json produced by evolve_b1.py, and to B1_REDUCED.json (conditional)."""
import json, glob, math, hashlib, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
A1 = Path('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum')

def load(path_summary):
    s = json.load(open(path_summary))
    d = np.load(str(path_summary).replace('_summary.json', '_timeseries.npz'))
    cols = list(d['cols']); rec = d['rec']
    ts = {str(c): rec[:, i] for i, c in enumerate(cols)}
    return s, ts

def derived(ts):
    H = ts['H_over_H0']; lna = ts['ln_a']
    out = dict(ts)
    chi = ts.get('chi', np.zeros_like(H))
    out['chi'] = chi
    out['Wa4'] = ts['Wy']*np.exp(4*lna)
    out['Ra4'] = ts['R']*np.exp(4*lna)
    with np.errstate(divide='ignore', invalid='ignore'):
        out['Omega_r'] = ts['rad']/H**2
        out['r'] = np.where(ts['rad'] > 0, ts['Wy']/ts['rad'], np.nan)
        out['Omega_vac'] = ts['vac']/H**2
        out['Omega_W'] = ts['Wy']/H**2
        out['Omega_chi'] = chi/H**2
        out['Rj'] = np.abs(ts.get('J', np.zeros_like(H)))/np.abs(ts['sigma1']) if 'sigma1' in ts else np.zeros_like(H)
    out['closure'] = ts['vac'] + ts['rad'] + chi + ts['kin'] + ts['fric'] + ts['Wy'] - H**2
    return out

def plateau_index(d, dlna=0.5, tol=0.05, need_R=True):
    """identical to A1 a1_analyze.plateau_index"""
    lna = d['ln_a']; H = d['H_over_H0']; Wa4 = d['Wa4']; Ra4 = d['Ra4']
    n = len(lna); j0 = 0
    for p in range(n):
        if lna[p] - lna[0] < dlna: continue
        while j0 < p and lna[j0] < lna[p] - dlna: j0 += 1
        if lna[p] < lna[j0] + dlna - 0.02: continue
        win = slice(j0, p + 1)
        if np.any(H[win] <= 0): continue
        w = Wa4[win]; okW = (w.max() - w.min()) <= tol*abs(Wa4[p]) if Wa4[p] != 0 else False
        if need_R:
            rr = Ra4[win]; okR = Ra4[p] > 0 and (rr.max() - rr.min()) <= tol*Ra4[p]
        else: okR = True
        if okW and okR: return p, j0
    return None, None

def reliable_end(d, thr=0.05):
    hn = d.get('Hmax_near'); mn = d.get('Mmax_near')
    n = len(d['T'])
    if hn is None: return n - 1
    bad = np.flatnonzero((np.nan_to_num(hn, nan=0.0) > thr) | (np.nan_to_num(mn, nan=0.0) > thr))
    return (int(bad[0]) - 1) if len(bad) else n - 1

def truncate(d, iend):
    return {k: (v[:iend + 1] if isinstance(v, np.ndarray) and v.ndim == 1 and len(v) == len(d['T']) else v) for k, v in d.items()}

def has_radiation(s):
    m = s.get('matter', {}); p = s['params']
    if m.get('kind') == 'chi': return m.get('G', 0) > 0 and m.get('b', 0) > 0
    return p.get('Y', 0) > 0

def classify(s, d_full):
    """A1 rule (verbatim logic of a1_analyze.classify) with 'radiation present' in place of Y > 0"""
    iend = reliable_end(d_full); d = truncate(d_full, iend)
    rad_on = has_radiation(s)
    H = d['H_over_H0']; tau = d['H0tau']
    res = dict(tag=s['tag'], stop_reason=s['stop_reason'], H0tau_end_run=float(d_full['H0tau'][-1]), H0tau_end=float(tau[-1]),
               reliable_end_reason=('near-shell constraint > 0.05' if iend < len(d_full['T']) - 1 else 'end of run'),
               T_end=s['T_end'], runtime_s=s['runtime_s'], xc=s['params']['xc'], dzf=s['params']['dzf'],
               vac_end=float(d['vac'][-1]), phi_b_end=float(d['phi_b'][-1]))
    neg = np.flatnonzero(H < -0.05)
    good = (d['Omega_r'] >= 0.9) & (np.abs(d['r']) <= 0.1)
    first_good = np.flatnonzero(good)
    res['recollapse_H0tau'] = float(tau[neg[0]]) if len(neg) else None
    i0 = np.flatnonzero(H <= 0)
    res['H_zero_H0tau'] = float(np.interp(0.0, -H[max(i0[0]-1, 0):i0[0]+1], tau[max(i0[0]-1, 0):i0[0]+1])) if len(i0) and i0[0] > 0 else None
    if rad_on and np.any(np.isfinite(d['Omega_r'])):
        im = int(np.nanargmax(np.where(np.isfinite(d['Omega_r']), d['Omega_r'], -1)))
        res['max_Omega_r'] = float(d['Omega_r'][im])
        res['at_max_Omega_r'] = dict(H0tau=float(tau[im]), r=float(d['r'][im]), H_over_H0=float(H[im]), Omega_vac=float(d['Omega_vac'][im]),
                                     Omega_W=float(d['Omega_W'][im]), Omega_chi=float(d['Omega_chi'][im]), phi_b=float(d['phi_b'][im]))
    p, j0 = plateau_index(d, need_R=rad_on)
    if p is not None:
        res['plateau'] = dict(H0tau=float(tau[p]), window_start_H0tau=float(tau[j0]), r=float(d['r'][p]) if rad_on else None,
                              Omega_r=float(d['Omega_r'][p]) if rad_on else 0.0, Omega_vac=float(d['Omega_vac'][p]), Omega_chi=float(d['Omega_chi'][p]),
                              Wa4=float(d['Wa4'][p]), Ra4=float(d['Ra4'][p]), R_over_sigma=float(d['R'][p]/d['sigma'][p]) if 'sigma' in d else None,
                              H_over_H0=float(H[p]), phi_b=float(d['phi_b'][p]), ln_a=float(d['ln_a'][p]))
    else:
        res['plateau'] = None
    if len(neg) and (len(first_good) == 0 or first_good[0] > neg[0]) and (p is None or p > neg[0]):
        cls = 'FAIL-no-radiation-era (recollapse)'; res['classification_H0tau'] = float(tau[neg[0]])
    elif p is not None:
        om = res['plateau']['Omega_r']; r = res['plateau']['r']; res['classification_H0tau'] = float(tau[p])
        if not rad_on: cls = 'FAIL-no-radiation-era (baseline: no radiation)'
        elif abs(r) > 0.1: cls = 'FAIL-Weyl'
        elif om < 0.5: cls = 'FAIL-no-radiation-era (Omega_r < 0.5 at plateau)'
        elif om >= 0.9: cls = 'PASS-combined' if abs(r) <= 0.03 else 'PASS-conservative'
        else: cls = 'INCONCLUSIVE (0.5 <= Omega_r < 0.9 with |r| <= 0.1)'
    else:
        cls = 'INCONCLUSIVE (no plateau within the evolution)'; res['classification_H0tau'] = float(tau[-1])
    res['classification'] = cls
    tc = res['classification_H0tau']
    m = (tau <= tc + 1e-9) & np.isfinite(d.get('Hmax_near', np.full_like(tau, np.nan)))
    if 'Hmax_near' in d and m.any():
        res['max_rel_H_near_shell'] = float(np.nanmax(d['Hmax_near'][m])); res['max_rel_M_near_shell'] = float(np.nanmax(d['Mmax_near'][m]))
    res['closure_max_abs'] = float(np.nanmax(np.abs(d['closure'])))
    mm = tau <= tc + 1e-9
    res['max_Rj'] = float(np.nanmax(d['Rj'][mm])) if 'sigma1' in d else None
    res['Rj_at_class_time'] = float(d['Rj'][mm][-1]) if 'sigma1' in d else None
    # latest-time diagnostics (informative, also for runs without a plateau)
    res['last'] = dict(H0tau=float(tau[-1]), r=float(d['r'][-1]) if rad_on else None, Omega_r=float(d['Omega_r'][-1]) if rad_on else None,
                       Omega_chi=float(d['Omega_chi'][-1]), Omega_W=float(d['Omega_W'][-1]), Omega_vac=float(d['Omega_vac'][-1]),
                       H_over_H0=float(H[-1]), Wa4=float(d['Wa4'][-1]), Ra4=float(d['Ra4'][-1]), ln_a=float(d['ln_a'][-1]))
    return res

def _fd(x, t): return np.gradient(x, t, edge_order=2)

def weyl_identity(s, d, variant='ok'):
    """W_dot + 4HW = (1/6)[4 k_s w v - 4H v^2 - v v_dot + w w_dot - U' v] in H0 units, k_s = (sigma + rho)/6, w = -(sigma' + J)/2
    (J = kappa5^2 (j + j_prod) + Y v, all from the records; w_dot by finite differences)"""
    rb = s['rho_b']
    tau = d['H0tau']; keep = np.concatenate(([True], np.diff(tau) > 1e-9))
    g = lambda k: d[k][keep]
    tau = tau[keep]; h = g('H_over_H0'); vh = g('v_over_H0'); Wy = g('Wy')
    rho = g('R') + (g('X') if 'X' in d else 0)
    J = g('J') if 'J' in d else 0*h
    k = rb*(g('sigma') + (rho if variant != 'noX' else g('R')))/6
    wv = -rb*(g('sigma1') + (J if variant != 'noJ' else 0))/2
    vdot = _fd(vh, tau); wdot = _fd(wv, tau); dW = _fd(Wy, tau)
    fac = 3 if variant == '3H' else 4
    U1 = g('U1')
    rhs = (4*k*wv*vh - 4*h*vh*vh - vh*vdot + wv*wdot - rb**2*U1*vh)/6
    lhs = dW + fac*h*Wy
    scale = (np.abs(4*k*wv*vh) + np.abs(4*h*vh*vh) + np.abs(vh*vdot) + np.abs(wv*wdot) + np.abs(rb**2*U1*vh))/6 + np.abs(dW) + np.abs(fac*h*Wy) + 1e-300
    r = (np.abs(lhs - rhs)/scale)[3:-3]
    return dict(median=float(np.median(r)), p95=float(np.percentile(r, 95)))

def ledgers(s, d, sign=1.0):
    """(a) total shell matter: d(X+R)/dtau + 3H(X+R+P) - J0 v = 0;  (b) radiation: dR/dtau + 4HR - b Qhat/rb^2 = 0;
    (c) chi: dX/dtau + 3H(X + P - R/3) - J0 v + b Qhat/rb^2 = 0.  Model units; FD on records.  sign=-1: wrong-sign source."""
    rb = s['rho_b']; b = s.get('matter', {}).get('b', 0.0)
    tau = d['H0tau']; keep = np.concatenate(([True], np.diff(tau) > 1e-9))
    g = lambda k: d[k][keep]
    t = tau[keep]*rb; H = g('H_over_H0')/rb; v = g('v_over_H0')/rb
    X, R, P, J0 = g('X'), g('R'), g('P'), g('J0')
    Q = b*g('Q_hat')/rb**2 if 'Q_hat' in d else 0*R
    out = {}
    for name, dens, press, src in [('total', X + R, P, sign*J0*v), ('radiation', R, R/3, Q + 0*R), ('chi', X, P - R/3, sign*J0*v - Q)]:
        dd = _fd(dens, t)
        extra = (H*R*0 if name != 'radiation' else 0*R)
        L = dd + (4*H*R if name == 'radiation' else 3*H*(dens + press)) - src
        sc = np.abs(dd) + np.abs(3*H*(dens + press)) + np.abs(src) + 1e-300
        m = np.isfinite(L) & (sc > 1e-14*np.nanmax(sc))
        rr = (np.abs(L)/sc)[m][3:-3]
        out[name] = dict(median=float(np.median(rr)) if len(rr) else None, p95=float(np.percentile(rr, 95)) if len(rr) else None, n=int(len(rr)))
    return out

def cutoff(s, d):
    m = s.get('matter', {})
    if m.get('kind') != 'chi' or not m.get('G'): return None
    G, ps, b = m['G'], m['phistar'], m['b']
    ms = G*np.abs(d['phi_b'] - ps)
    qs = [e['q'] for e in m.get('events', [])]
    lam_full = max(float(ms.max()), max([math.sqrt(q) for q in qs], default=0.0))
    ev = m.get('events', [])
    post = ms[d['H0tau'] >= ev[0]['s_star']] if ev else ms
    lam_post = max(float(post.max()), max([math.sqrt(q) for q in qs], default=0.0))
    return dict(Lambda_hat_full=lam_full, Lambda_hat_post=lam_post, lambda_c_full=lam_full*b**(1/3), lambda_c_post=lam_post*b**(1/3),
                n_crossings=len(ev), events=[{k: e.get(k) for k in ('s_star', 'q', 'D_over_q', 'number_factor', 'screen_pass', 'G_min_screen',
                                                                    'D_outside_validity', 'kappa5sq_rho_jump_over_H0', 'massless_energy_proxy_over_H0', 'ramp_H0tau')} for e in ev])

def interp_T(dA, dB, keys):
    """max |A - B| on the common coordinate-time records (identical T sampling expected)"""
    TA, TB = dA['T'], dB['T']; n = min(len(TA), len(TB))
    out = dict(n_common=n, max_dT=float(np.max(np.abs(TA[:n] - TB[:n]))))
    for k in keys: out[k] = float(np.nanmax(np.abs(dA[k][:n] - dB[k][:n])))
    return out

def clean(v):
    if isinstance(v, dict): return {str(k): clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)): return [clean(x) for x in v]
    if isinstance(v, (np.floating, float)): return None if not math.isfinite(float(v)) else float(v)
    if isinstance(v, (np.integer,)): return int(v)
    if isinstance(v, np.bool_): return bool(v)
    return v

def main():
    out = dict(status='numerical (5D runs) / conditional (reduced estimate); MODEL CHANGE B1 (derived chi source) on the A1 tuned model d = d*(0.1)',
               registration='REGISTRATION.md', calibration={}, runs={}, controls={}, convergence={}, aggregate={})
    data = {}
    for f in sorted(glob.glob(str(HERE/'runs/**/*_summary.json'), recursive=True)):
        s, ts = load(f); tag = Path(f).name.replace('_summary.json', '')
        data[tag] = (s, derived(ts))
    # ---------------- calibration K1 / K2 against the A1 archive
    cal = {}
    pairs = [('K1_chioff_Y0_dc1e-2_dzf1e-3', 'main_dstar_Y0_dc1e-2_dzf1e-3', 1e-8), ('K2_toy_Y1_dc1e-2_dzf1e-3', 'main_dstar_Y1_dc1e-2_dzf1e-3', 1e-6)]
    for tag, ref, tol in pairs:
        if tag not in data: continue
        s, d = data[tag]; sa, ta = load(A1/'runs/main'/(ref + '_summary.json')); da = derived(ta)
        cmp_ = interp_T(d, da, ['phi_b', 'H_over_H0', 'Wy', 'R'])
        ca = classify(sa, da); cb = classify(s, d)
        rA = (ca.get('plateau') or {}).get('r'); rB = (cb.get('plateau') or {}).get('r')
        entry = dict(ref=ref, diffs=cmp_, stop_B1=s['stop_reason'], stop_A1=sa['stop_reason'], H0tau_end_B1=s['H0tau_end'], H0tau_end_A1=sa['H0tau_end'],
                     class_B1=cb['classification'], class_A1=ca['classification'], r_B1=rB, r_A1=rA, tol_fields=tol)
        same_class = cb['classification'].split(' ')[0] == ca['classification'].split(' ')[0]
        fields_ok = max(cmp_['phi_b'], cmp_['H_over_H0'], cmp_['Wy']) <= tol
        if tag.startswith('K1'):
            entry['PASS'] = bool(fields_ok and s['stop_reason'] == sa['stop_reason'])
        else:
            r_ok = (rA is not None and rB is not None and abs(rA - rB) <= max(0.2*abs(rA), 0.02))
            entry['PASS_formal_A1_tolerance'] = bool(same_class and r_ok)
            entry['fields_within_1e-6'] = bool(fields_ok)
        cal[tag] = entry
    out['calibration'] = cal
    # ---------------- per-run classification, ledgers, Weyl identity
    for tag, (s, d) in data.items():
        if not (tag.startswith('main_') or tag.startswith('ctl_')): continue
        r = classify(s, d)
        r['matter'] = {k: s.get('matter', {}).get(k) for k in ('kind', 'G', 'phistar', 'y', 'b', 'mode', 'n_cohorts')}
        r['matter']['ramp'] = s['params'].get('matter', {}).get('ramp'); r['matter']['apply_D'] = s['params'].get('matter', {}).get('apply_D')
        r['cutoff'] = cutoff(s, d)
        r['weyl_identity'] = weyl_identity(s, d); r['weyl_identity_control_3H'] = weyl_identity(s, d, '3H')
        r['weyl_identity_control_noJ'] = weyl_identity(s, d, 'noJ')
        r['ledgers'] = ledgers(s, d); r['ledgers_control_wrong_sign'] = ledgers(s, d, sign=-1.0)
        (out['runs'] if tag.startswith('main_') else out['controls'])[tag] = r
    # ---------------- convergence pairs and reliability
    main = out['runs']; groups = {}
    for tag, r in main.items():
        cell = tag.rsplit('_dzf', 1)[0]
        groups.setdefault(cell, {})[r['dzf']] = tag
    agg = {}
    for cell, g in sorted(groups.items()):
        entry = dict(runs=g)
        if len(g) >= 2:
            fine = main[g[min(g)]]; coarse = main[g[max(g)]]
            same = fine['classification'].split(' ')[0] == coarse['classification'].split(' ')[0]
            agree = None
            if fine.get('plateau') and coarse.get('plateau') and fine['plateau']['r'] is not None and coarse['plateau']['r'] is not None:
                rf, rc = fine['plateau']['r'], coarse['plateau']['r']; agree = abs(rf - rc) <= max(0.2*abs(rf), 0.02)
            elif 'recollapse' in fine['classification'] and 'recollapse' in coarse['classification']:
                agree = abs(fine['recollapse_H0tau'] - coarse['recollapse_H0tau']) <= 0.05
            cons_ok = fine.get('max_rel_H_near_shell', 1) < 0.05 and fine.get('max_rel_M_near_shell', 1) < 0.05
            rel = bool(same and agree and cons_ok)
            entry.update(fine=fine['classification'], coarse=coarse['classification'], same_class=same, value_agreement=agree,
                         constraints_ok=cons_ok, reliable=rel, class_reliable=(fine['classification'] if rel else 'UNRELIABLE'),
                         r_fine=(fine.get('plateau') or {}).get('r'), r_coarse=(coarse.get('plateau') or {}).get('r'),
                         Omega_r_fine=(fine.get('plateau') or {}).get('Omega_r'), max_Rj_fine=fine.get('max_Rj'),
                         lambda_c_full=(fine.get('cutoff') or {}).get('lambda_c_full'))
        else:
            entry.update(reliable=False, class_reliable='UNRELIABLE (single resolution)',
                         single=main[list(g.values())[0]]['classification'])
        agg[cell] = entry
    out['convergence'] = agg
    # ---------------- aggregate verdicts (registration section 4)
    cells_cut = {c: e for c, e in agg.items() if c.endswith('_bcut')}
    cells_need = {c: e for c, e in agg.items() if c.endswith('_bneed')}
    red = json.load(open(HERE/'B1_REDUCED.json')) if (HERE/'B1_REDUCED.json').exists() else None
    red_fail_all = None
    if red:
        prim = {k: v for k, v in red['cases'].items() if k.startswith('main_dstar_Y0_dc1e-2_dzf5e-4|')}
        red_fail_all = bool(len(prim) == 8 and all((v['cls_at_b_cut_full_lowE']['cls'].startswith('FAIL')) for v in prim.values()))
    expected_cut = ['ps%s_G%s_y%s_bcut' % (a, b, c) for a in ('0.5', '0.9') for b in ('100', '1000') for c in ('1', '0.1')]
    done_cut = [c for c in expected_cut if ('main_' + c) in cells_cut]
    classes_cut = {c: cells_cut['main_' + c]['class_reliable'] for c in done_cut}
    if any(v.startswith('PASS') for v in classes_cut.values()): v_cut = 'B1-PASS at lambda_c = 1'
    elif done_cut and all(v.startswith('FAIL') for v in classes_cut.values()) and red_fail_all: v_cut = 'B1-FAIL at lambda_c = 1'
    else: v_cut = 'INCONCLUSIVE at lambda_c = 1'
    classes_need = {c: e['class_reliable'] for c, e in cells_need.items()}
    out['aggregate'] = dict(verdict_lambda_c_1=v_cut, b_cut_cells_completed=done_cut, b_cut_cells_not_run=[c for c in expected_cut if c not in done_cut],
                            b_cut_classes=classes_cut, reduced_all_eight_fail_at_b_cut=red_fail_all, b_need_classes=classes_need,
                            any_reliable_PASS_at_b_need=any(v.startswith('PASS') for v in classes_need.values()))
    if red:
        out['reduced_summary'] = {k.split('|', 1)[1]: dict(q=v['q'], D_over_q=v['D_over_q'], Lambda_hat_full=v['Lambda_hat_full'], Lambda_hat_post=v['Lambda_hat_post'],
                                                            b_cut_full=v['b_cut_full'], class_at_b_cut=v['cls_at_b_cut_full_lowE']['cls'], r_at_b_cut=v['cls_at_b_cut_full_lowE']['r'],
                                                            class_at_b_cut_post=v['cls_at_b_cut_post_lowE']['cls'], r_at_b_cut_post=v['cls_at_b_cut_post_lowE']['r'],
                                                            b_need_r0p1=v.get('b_need_r0p1'), lambda_c_full_r0p1=v.get('lambda_c_full_r0p1'), lambda_c_post_r0p1=v.get('lambda_c_post_r0p1'),
                                                            b_need_r0p03=v.get('b_need_r0p03'), lambda_c_full_r0p03=v.get('lambda_c_full_r0p03'),
                                                            b_interp_r0p1=v.get('b_interp_r0p1'), b_interp_r0p03=v.get('b_interp_r0p03'),
                                                            Rj_max_at_b_need=(v.get('backreaction_at_b_need_r0p1') or {}).get('max_Rj'),
                                                            Rj_end_at_b_need=(v.get('backreaction_at_b_need_r0p1') or {}).get('Rj_end'),
                                                            Rj_max_at_b_cut=(v.get('backreaction_at_b_cut_full') or {}).get('max_Rj'))
                                  for k, v in red['cases'].items() if k.startswith('main_dstar_Y0_dc1e-2_dzf5e-4|')}
        out['reduced_grid_spread'] = {}
        for k, v in red['cases'].items():
            if not k.startswith('main_dstar_Y0_dc1e-2_dzf5e-4|'): continue
            key = k.split('|', 1)[1]; alt = [red['cases'].get(t + '|' + key) for t in ('main_dstar_Y0_dc1e-2_dzf1e-3', 'main_dstar_Y0_dc1e-4_dzf5e-4')]
            r0 = v['cls_at_b_cut_full_lowE']['r']
            out['reduced_grid_spread'][key] = dict(r_at_b_cut=[r0] + [a['cls_at_b_cut_full_lowE']['r'] for a in alt if a],
                                                   b_interp_r0p1=[v.get('b_interp_r0p1')] + [a.get('b_interp_r0p1') for a in alt if a])
    # production-model controls vs the coarse main run
    base = main.get('main_ps0.5_G100_y1_bneed_dzf1e-3')
    if base:
        pc = {}
        for tag, r in out['controls'].items():
            rb_ = (base.get('plateau') or {}).get('r'); rc = (r.get('plateau') or {}).get('r')
            pc[tag] = dict(r=rc, r_base=rb_, class_=r['classification'], rel_change=(abs(rc/rb_ - 1) if (rb_ and rc) else None),
                           last_r=r['last']['r'], last_r_base=base['last']['r'])
        out['production_model_controls'] = pc
    out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (HERE/'B1_RESULTS.json').write_text(json.dumps(clean(out), indent=1) + '\n')
    for tag, r in sorted(main.items()):
        pl = r.get('plateau') or {}
        print('%-44s %-40s r=%-9s Om_r=%-7s Rj=%-8s lam_c=%-7s end=%6.2f last_r=%s Wid=%.1e/%.1e led=%s' % (
            tag, r['classification'][:40], ('%.4g' % pl['r']) if pl.get('r') is not None else '-', ('%.3f' % pl['Omega_r']) if pl.get('Omega_r') is not None else '-',
            ('%.3g' % r['max_Rj']) if r.get('max_Rj') is not None else '-', ('%.3g' % r['cutoff']['lambda_c_full']) if r.get('cutoff') else '-',
            r['H0tau_end'], ('%.3g' % r['last']['r']) if r['last']['r'] is not None else '-', r['weyl_identity']['median'], r['weyl_identity_control_3H']['median'],
            ('%.1e' % r['ledgers']['total']['median']) if r['ledgers']['total']['median'] is not None else '-'))
    print(json.dumps(clean(out['aggregate']), indent=1))
    print(json.dumps(clean(out['calibration']), indent=1))

if __name__ == '__main__':
    main()
