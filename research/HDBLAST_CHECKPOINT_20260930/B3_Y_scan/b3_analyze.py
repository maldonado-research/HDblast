#!/usr/bin/env python3
"""B3 analysis: r(Y) scan of the tuned model (d = d*(0.1)) with the friction closure.  Applies the A1 per-run rule
UNCHANGED (a1_analyze.classify, copied verbatim from A1) and the A1 reliability criterion to resolution pairs, as registered
in REGISTRATION.md.  Adds: (i) Richardson-extrapolated plateau values (supplementary, not decisive), (ii) the targeted C5
Weyl-identity controls, (iii) a fit of r(Y).  Output B3_RESULTS.json.  Numerical."""
import json, glob, math, hashlib, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent; sys.path.insert(0, str(HERE))
import a1_analyze as A
import static_w as S
A1RUNS = Path('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum/runs/main')

def load_hybrid(f):
    """restart run: parent records with T < T_restart, then the restart run's own records (REGISTRATION section 4)"""
    s, ts = A.load(f)
    ri = s.get('restart_info')
    if not ri: return s, ts
    st = Path(ri['file']); parent_tag = st.name.split('_state_T')[0]
    pf = st.parent.parent/(parent_tag + '_summary.json')
    sp, tp = A.load(pf)
    m = tp['T'] < ri['T'] - 1e-9
    out = {}
    for k in ts:
        if k == 'modes': out[k] = np.concatenate([tp[k][m], ts[k]]); continue
        if k in tp: out[k] = np.concatenate([tp[k][m], ts[k]])
    s = dict(s); s['hybrid_parent'] = parent_tag
    return s, out

def load_all():
    files = sorted(glob.glob(str(HERE/'runs/main/*_summary.json'))) + sorted(glob.glob(str(HERE/'runs/fine/*_summary.json')))
    # A1 runs reused (registered in REGISTRATION.md section 3): Y = 0.3, 1, 3 at dzf 1e-3 / 5e-4
    files += sorted(glob.glob(str(A1RUNS/'main_dstar_Y*_summary.json')))
    data = {}
    # a wall-limited run continued on the same grid ('<tag>_ext', fourth dated note) replaces its parent
    files = [f for f in files if not Path(f.replace('_summary.json', '_ext_summary.json')).exists()]
    for f in files:
        s, ts = load_hybrid(f); d = A.derived(s, ts)
        tag = Path(f).name.replace('_summary.json', '')
        if str(A1RUNS) in f: tag = 'A1_' + tag
        if s['params']['Y'] == 0: continue
        data[tag] = (s, d, 'A1' if str(A1RUNS) in f else 'B3')
    return data

def weyl_terms(s, d, iend):
    p = s['params']; ten = S.Tension(p['delta'], p['c'], p['d']); Y = p['Y']; rb = s['rho_b']
    tau = d['H0tau'][:iend + 1]; keep = np.concatenate(([True], np.diff(tau) > 1e-9))
    g = lambda k: d[k][:iend + 1][keep]
    tau = tau[keep]; phi = g('phi_b'); h = g('H_over_H0'); vh = g('v_over_H0'); R = g('R'); Wy = g('Wy')
    vdot = np.gradient(vh, tau, edge_order=2); dW = np.gradient(Wy, tau, edge_order=2)
    return dict(tau=tau, phi=phi, h=h, vh=vh, R=R, Wy=Wy, vdot=vdot, dW=dW, ten=ten, Y=Y, rb=rb)

def weyl_resid(T, variant='ok'):
    """V2-form identity dW/dtau + 4 H W = (4 k w v - 4 H v^2 - v vdot + w wdot - U' v)/6 (H0 units); variants are wrong-factor controls"""
    ten, Y, rb = T['ten'], T['Y'], T['rb']; phi, h, vh, R, Wy, vdot, dW = T['phi'], T['h'], T['vh'], T['R'], T['Wy'], T['vdot'], T['dW']
    Rk = {'noR': 0*R, 'R3': 3*R, 'Rneg': -R}.get(variant, R)
    k = rb*(ten.s(phi) + Rk)/6
    Yw = 0.0 if variant == 'noYw' else Y
    wv = -(rb*ten.s1(phi) + Yw*vh)/2
    wdot = -(rb*ten.s2(phi)*vh + Yw*vdot)/2
    fac = 3 if variant == '3H' else 4
    terms = [4*k*wv*vh, -4*h*vh*vh, -vh*vdot, wv*wdot, -rb**2*S.U1(phi)*vh]
    rhs = sum(terms)/6
    lhs = dW + fac*h*Wy
    scale = sum(np.abs(t) for t in terms)/6 + np.abs(dW) + np.abs(fac*h*Wy) + 1e-300
    wR = np.abs(4*rb*R/6*wv*vh)/6/scale            # weight of the R part of the k-term in the identity
    return np.abs(lhs - rhs)/scale, wR

def c5_controls(s, d, iend, thetaR=1e-3):
    T = weyl_terms(s, d, iend)
    ok, wR = weyl_resid(T, 'ok')
    sl = slice(3, -3)
    out = dict(window='records up to the classification time (or reliable end)', all={}, R_window={})
    for v in ['ok', '3H', 'noR', 'R3', 'Rneg', 'noYw']:
        r, _ = weyl_resid(T, v)
        out['all'][v] = float(np.median(r[sl]))
    m = np.zeros_like(wR, dtype=bool); m[sl] = wR[sl] >= thetaR
    out['R_window'] = dict(thetaR=thetaR, n_records=int(m.sum()), H0tau_range=[float(T['tau'][m].min()), float(T['tau'][m].max())] if m.any() else None)
    if m.any():
        for v in ['ok', 'noR', 'R3', 'Rneg']:
            r, _ = weyl_resid(T, v); out['R_window'][v] = float(np.median(r[m]))
    ok_all = out['all']['ok']
    out['ratios_all'] = {v: out['all'][v]/ok_all for v in out['all'] if v != 'ok'}
    if m.any(): out['ratios_R_window'] = {v: out['R_window'][v]/out['R_window']['ok'] for v in ['noR', 'R3', 'Rneg']}
    out['median_weight_R_term'] = float(np.median(wR[sl])); out['p90_weight_R_term'] = float(np.percentile(wR[sl], 90))
    return out

def pair_reliability(fine, coarse):
    """A1 criterion (registration section 3 (i)-(ii)), identical to a1_analyze.main"""
    same = fine['classification'].split(' ')[0] == coarse['classification'].split(' ')[0]
    agree = None
    if fine.get('plateau') and coarse.get('plateau') and fine['plateau']['r'] is not None and coarse['plateau']['r'] is not None:
        rf, rc = fine['plateau']['r'], coarse['plateau']['r']; agree = abs(rf - rc) <= max(0.2*abs(rf), 0.02)
    elif 'recollapse' in fine['classification'] and 'recollapse' in coarse['classification']:
        agree = abs(fine['recollapse_H0tau'] - coarse['recollapse_H0tau']) <= 0.05
    cons_ok = fine.get('max_rel_H_near_shell', 1) < 0.05 and fine.get('max_rel_M_near_shell', 1) < 0.05
    rel = bool(same and agree and cons_ok)
    return dict(fine=fine['tag'], coarse=coarse['tag'], fine_class=fine['classification'], coarse_class=coarse['classification'],
                same_class=same, value_agreement=agree, constraints_ok=cons_ok, reliable=rel,
                class_reliable=fine['classification'] if rel else 'UNRELIABLE')

def richardson(dF, dC, p=4.0):
    """supplementary: plateau of the Richardson-extrapolated series X_RE = X_f + (X_f - X_c)/(2^p - 1) on the fine run's tau grid"""
    tauF = dF['H0tau']; tauC = dC['H0tau']; keep = np.concatenate(([True], np.diff(tauC) > 1e-12))
    ie = min(A.reliable_end(dF), np.searchsorted(tauF, tauC[A.reliable_end(dC)]) - 1)
    re = {}
    for k in ['H_over_H0', 'ln_a', 'Wy', 'R', 'rad', 'vac', 'kin', 'fric', 'phi_b', 'v_over_H0', 'lapse', 'T']:
        xc = np.interp(tauF[:ie + 1], tauC[keep], dC[k][keep])
        re[k] = dF[k][:ie + 1] + (dF[k][:ie + 1] - xc)/(2**p - 1)
    re['H0tau'] = tauF[:ie + 1]
    for k in ['Hmax_near', 'Mmax_near', 'Hmax', 'Mmax']: re[k] = dF[k][:ie + 1]
    re['Wa4'] = re['Wy']*np.exp(4*re['ln_a']); re['Ra4'] = re['R']*np.exp(4*re['ln_a'])
    H = re['H_over_H0']
    with np.errstate(divide='ignore', invalid='ignore'):
        re['Omega_r'] = re['rad']/H**2; re['r'] = np.where(re['rad'] > 0, re['Wy']/re['rad'], np.nan)
        re['Omega_vac'] = re['vac']/H**2; re['Omega_W'] = re['Wy']/H**2
    re['closure'] = 0*H
    return re

def main():
    data = load_all()
    out = dict(status='numerical', registration='REGISTRATION.md', runs={}, pairs={}, per_Y={}, richardson={}, c5={}, fit={})
    for tag, (s, d, src) in data.items():
        r = A.classify(s, d); r['source'] = src
        # extras (not part of the rule): e-folds from plateau to turnaround; a^4-weighted kinetic integral
        lna = d['ln_a']; H = d['H_over_H0']; tau = d['H0tau']
        iH = int(np.argmax(H <= 0)) if np.any(H <= 0) else None
        r['turnaround'] = dict(H0tau=float(tau[iH]) if iH else None, ln_a_max=float(lna.max()), reached=iH is not None)
        if r.get('plateau'):
            ip = int(np.searchsorted(tau, r['plateau']['H0tau']))
            r['plateau']['ln_a'] = float(lna[ip])
            r['plateau']['efolds_to_turnaround'] = float(lna.max() - lna[ip]) if iH else None
            om_W = float(d['Omega_W'][ip]); om_v = r['plateau']['Omega_vac']
            r['plateau']['Omega_W'] = om_W
            r['plateau']['efolds_to_turnaround_estimate'] = float(0.25*math.log((r['plateau']['Omega_r'] + om_W)/abs(om_v))) if om_v < 0 else None
            r['plateau']['R_over_sigma'] = float(d['R'][ip]/S.Tension(s['params']['delta'], s['params']['c'], s['params']['d']).s(d['phi_b'][ip]))
            keep = np.concatenate(([True], np.diff(tau[:ip + 1]) > 1e-12))
            t_ = tau[:ip + 1][keep]; v_ = d['v_over_H0'][:ip + 1][keep]; a4 = np.exp(4*lna[:ip + 1][keep])
            I1 = float(np.trapezoid(v_**2*a4, t_))
            r['plateau']['I_v2a4'] = I1; r['plateau']['Wa4_over_I_v2a4'] = float(r['plateau']['Wa4']/I1)
            r['plateau']['Ra4_over_Y_I_v2a4'] = float(r['plateau']['Ra4']/(s['params']['Y']*I1/s['rho_b']))
        r['weyl_identity_A1form'] = A.weyl_identity(s, d); r['radiation_ledger'] = A.ledger(s, d)
        out['runs'][tag] = r
    # group by (Y, dc): dzf -> tag
    groups = {}
    for tag, r in out['runs'].items():
        if r['source'] == 'A1': continue
        groups.setdefault((r['Y'], r['dc']), {})[r['dzf']] = tag
    for (Y, dc), g in sorted(groups.items()):
        key = 'Y=%g dc=%g' % (Y, dc); entry = dict(Y=Y, dc=dc, runs={('%g' % k): v for k, v in g.items()})
        dz = sorted(g)
        pr = {}
        for a, b in zip(dz[:-1], dz[1:]):
            pr['%g/%g' % (b, a)] = pair_reliability(out['runs'][g[a]], out['runs'][g[b]])
        entry['pairs'] = pr
        # registered order: primary pair (1e-3, 5e-4); secondary (5e-4, 2.5e-4) per REGISTRATION section 4
        prim = pr.get('0.001/0.0005'); sec = pr.get('0.0005/0.00025'); ter = pr.get('0.00025/0.000125') if Y in (3, 5) else None   # third and fifth dated notes
        if prim and prim['reliable']: entry['class'] = prim['class_reliable']; entry['class_from'] = 'primary pair 1e-3/5e-4'; best = g[5e-4]
        elif sec and sec['reliable']: entry['class'] = sec['class_reliable']; entry['class_from'] = 'secondary pair 5e-4/2.5e-4'; best = g[2.5e-4]
        elif ter and ter['reliable']: entry['class'] = ter['class_reliable']; entry['class_from'] = 'tertiary pair 2.5e-4/1.25e-4 (third/fifth dated note)'; best = g[1.25e-4]
        else: entry['class'] = 'UNRELIABLE'; entry['class_from'] = None; best = g[min(dz)]
        b = out['runs'][best]
        entry['finest_run'] = best
        pl = b.get('plateau') or {}
        entry['plateau_finest'] = dict(r=pl.get('r'), Omega_r=pl.get('Omega_r'), Omega_vac=pl.get('Omega_vac'), H0tau=pl.get('H0tau'),
                                       H_over_H0=pl.get('H_over_H0')) if pl else None
        entry['at_max_Omega_r_finest'] = b['at_max_Omega_r']
        # Richardson (supplementary) on the finest available pair
        if len(dz) >= 2:
            f_, c_ = dz[0], dz[1]
            sF, dF, _ = data[g[f_]]; sC, dC, _ = data[g[c_]]
            rich = {}
            for p in [3.5, 4.0, 4.5]:
                re = richardson(dF, dC, p)
                pi, j0 = A.plateau_index(re)
                rich['p=%g' % p] = dict(plateau=None if pi is None else dict(H0tau=float(re['H0tau'][pi]), r=float(re['r'][pi]), Omega_r=float(re['Omega_r'][pi]),
                                       Omega_vac=float(re['Omega_vac'][pi]), H_over_H0=float(re['H_over_H0'][pi]), window_start_H0tau=float(re['H0tau'][j0])))
            out['richardson'][key] = dict(pair='%g/%g' % (c_, f_), **rich)
        # C5 targeted controls on the finest run up to its classification time
        sB, dB, _ = data[best]
        ic = int(np.searchsorted(dB['H0tau'], b['classification_H0tau']))
        out['c5'][key] = c5_controls(sB, dB, min(ic, len(dB['H0tau']) - 1))
        out['pairs'][key] = entry
    # A1 context (A1 runs re-analysed with the same code; primary pair only)
    gA = {}
    for tag, r in out['runs'].items():
        if r['source'] == 'A1': gA.setdefault((r['Y'], r['dc']), {})[r['dzf']] = tag
    out['A1_context'] = {}
    for (Y, dc), g in sorted(gA.items()):
        if 1e-3 in g and 5e-4 in g:
            pr = pair_reliability(out['runs'][g[5e-4]], out['runs'][g[1e-3]]); pl = out['runs'][g[5e-4]].get('plateau') or {}
            sB, dB, _ = data[g[5e-4]]; bB = out['runs'][g[5e-4]]
            ic = int(np.searchsorted(dB['H0tau'], bB['classification_H0tau']))
            out['c5']['A1 Y=%g dc=%g' % (Y, dc)] = c5_controls(sB, dB, min(ic, len(dB['H0tau']) - 1))
            out['A1_context']['Y=%g dc=%g' % (Y, dc)] = dict(Y=Y, dc=dc, class_=pr['class_reliable'], r=pl.get('r'), Omega_r=pl.get('Omega_r'), Omega_vac=pl.get('Omega_vac'), H0tau=pl.get('H0tau'))
    # per Y summary (dc = 1e-2 primary; dc = 1e-4 repeat where run)
    for key, e in out['pairs'].items():
        out['per_Y'].setdefault('%g' % e['Y'], {})['dc=%g' % e['dc']] = dict(class_=e['class'], class_from=e['class_from'], plateau=e['plateau_finest'],
                                                                             richardson_p4=(out['richardson'].get(key) or {}).get('p=4'))
    # fit r(Y) ~ a Y^-p over reliable plateaus (dc = 1e-2)
    pts = [(e['Y'], e['plateau_finest']['r']) for e in out['pairs'].values() if e['dc'] == 0.01 and e['class'] != 'UNRELIABLE' and e['plateau_finest'] and e['plateau_finest']['r']]
    pts += [(e['Y'], e['r']) for e in out['A1_context'].values() if e['dc'] == 0.01 and e['class_'] != 'UNRELIABLE' and e['r'] and e['Y'] not in [p_[0] for p_ in pts]]
    pts = sorted(pts)
    if len(pts) >= 2:
        Ys = np.array([p[0] for p in pts]); rs = np.array([p[1] for p in pts])
        A_ = np.vstack([np.ones_like(Ys), np.log(Ys)]).T; co, *_ = np.linalg.lstsq(A_, np.log(rs), rcond=None)
        out['fit']['power_law_reliable'] = dict(note='includes A1 Y = 0.3, 1 (dc 1e-2) as context points', local_slopes=[float(-math.log(rs[i + 1]/rs[i])/math.log(Ys[i + 1]/Ys[i])) for i in range(len(Ys) - 1)], points=pts, a=float(math.exp(co[0])), p=float(-co[1]), max_rel_dev=float(np.max(np.abs(np.exp(A_ @ co)/rs - 1))))
    # trend: r = (W a^4)/(rad a^4); R a^4 = Y * int v^2 a^4 dtau exactly (ledger), so r ~ (18/sigma) (W a^4 / int v^2 a^4) / Y
    tr = []
    for key, e in out['pairs'].items():
        if e['dc'] != 0.01 or e['class'] == 'UNRELIABLE' or not e['plateau_finest']: continue
        pl = out['runs'][e['finest_run']]['plateau']; tr.append((e['Y'], pl['r'], pl['Wa4_over_I_v2a4'], 'B3'))
    for key, e in out['A1_context'].items():
        if e['dc'] != 0.01 or e['class_'] == 'UNRELIABLE' or any(abs(t[0] - e['Y']) < 1e-9 for t in tr): continue
        tag = [t for t, r in out['runs'].items() if r['source'] == 'A1' and r['Y'] == e['Y'] and r['dc'] == 0.01 and r['dzf'] == 5e-4][0]
        pl = out['runs'][tag]['plateau']; tr.append((e['Y'], pl['r'], pl['Wa4_over_I_v2a4'], 'A1'))
    tr.sort()
    if len(tr) >= 2:
        Ys = np.array([t[0] for t in tr]); q = np.array([t[2] for t in tr]); rr = np.array([t[1] for t in tr])
        cq = np.polyfit(np.log(Ys), np.log(q), 1); cr = np.polyfit(np.log(Ys), np.log(rr), 1)
        out['fit']['trend'] = dict(points=[dict(Y=t[0], r=t[1], Wa4_over_int_v2a4=t[2], source=t[3], radiation_share_1_over_1_plus_r=1/(1 + t[1])) for t in tr],
                                   Wa4_over_int_v2a4_power=float(cq[0]), r_power=float(cr[0]),
                                   Y_at_r_0p03_from_fit=float(math.exp((math.log(0.03) - cr[1])/cr[0])),
                                   explanation='R a^4 = Y int v^2 a^4 dtau exactly (ledger; checked: ratio 1.000), so r = (18/sigma)(W a^4/int v^2 a^4)/Y at low energy. The Weyl yield per unit a^4-weighted kinetic integral itself falls as Y^q (q fitted), so r ~ Y^(q-1).')
    out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    def clean(v):
        if isinstance(v, dict): return {str(k): clean(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)): return [clean(x) for x in v]
        if isinstance(v, (np.floating, float)): return None if not math.isfinite(float(v)) else float(v)
        if isinstance(v, (np.integer,)): return int(v)
        if isinstance(v, np.bool_): return bool(v)
        return v
    (HERE/'B3_RESULTS.json').write_text(json.dumps(clean(out), indent=1) + '\n')
    for tag, r in sorted(out['runs'].items(), key=lambda kv: (kv[1]['Y'], kv[1]['dc'], -kv[1]['dzf'])):
        pl = r.get('plateau') or {}
        print('%-44s %-46s r=%-8s Om_r=%-7s Om_vac=%-7s tp=%-6s end=%.1f relend=%.1f' % (tag, r['classification'][:46], ('%.4f' % pl['r']) if pl else '-', ('%.3f' % pl['Omega_r']) if pl else '-',
              ('%.3f' % pl['Omega_vac']) if pl else '-', ('%.1f' % pl['H0tau']) if pl else '-', r['H0tau_end_run'], r['H0tau_end']))
    for k, e in out['pairs'].items(): print(k, e['class'], e['class_from'], {kk: (v['same_class'], v['value_agreement'], v['constraints_ok']) for kk, v in e['pairs'].items()})
    print(json.dumps(clean(out['richardson']), indent=0)[:3000]); print(json.dumps(clean(out['fit'])))
    for k, v in out['c5'].items(): print(k, 'ratios_all', {a: round(b, 1) for a, b in v['ratios_all'].items()}, 'R_window', v['R_window'].get('n_records'), {a: round(b, 1) for a, b in (v.get('ratios_R_window') or {}).items()}, 'ok', v['all']['ok'])

if __name__ == '__main__':
    main()
