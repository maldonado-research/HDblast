#!/usr/bin/env python3
"""Independent audit checks for B3 (r(Y) scan).  Reads only the saved run outputs (npz/json); does NOT import the producer's
analysis code (own loader, own hybrid splice, own plateau finder).  Output: audit/AUDIT_CHECKS.json.  Numerical."""
import json, math, glob
from pathlib import Path
import numpy as np
B3 = Path(__file__).resolve().parent.parent
A1 = Path('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum/runs/main')

def raw(summ):
    s = json.load(open(summ)); z = np.load(str(summ).replace('_summary.json', '_timeseries.npz'))
    c = [str(x) for x in z['cols']]; rec = z['rec']
    return s, {k: rec[:, i].copy() for i, k in enumerate(c)}

def series(summ):
    """own hybrid splice: parent records strictly before the restart T, then the child's records"""
    s, t = raw(summ); ri = s.get('restart_info')
    if ri:
        st = Path(ri['file']); par = st.parent.parent/(st.name.split('_state_T')[0] + '_summary.json')
        sp, tp = series(par); m = tp['T'] < ri['T'] - 1e-9
        t = {k: np.concatenate([tp[k][m], t[k]]) for k in t if k in tp}
    return s, t

def rel_end(t, thr=0.05):
    bad = np.flatnonzero((np.nan_to_num(t['Hmax_near']) > thr) | (np.nan_to_num(t['Mmax_near']) > thr))
    return int(bad[0]) - 1 if len(bad) else len(t['T']) - 1

def plateau(t, dlna=0.5, tol=0.05):
    """own implementation: first record p with ln a_p - ln a_0 >= dlna such that on {j: ln a_j in [ln a_p - dlna, ln a_p], j<=p}
    (contiguous trailing block) H>0 and the spreads of W a^4 and R a^4 are <= tol*value at p.  Also returns the minimum
    achieved max(spreadW, spreadR)/value over all candidate p (closeness to passing)."""
    ie = rel_end(t); lna = t['ln_a'][:ie + 1]; H = t['H_over_H0'][:ie + 1]
    Wa4 = t['Wy'][:ie + 1]*np.exp(4*lna); Ra4 = t['R'][:ie + 1]*np.exp(4*lna)
    best = (np.inf, None)
    for p in range(len(lna)):
        if lna[p] - lna[0] < dlna: continue
        j = p
        while j > 0 and lna[j - 1] >= lna[p] - dlna: j -= 1
        if lna[p] - lna[j] < dlna - 0.02: continue
        if np.any(H[j:p + 1] <= 0): continue
        sW = np.ptp(Wa4[j:p + 1])/abs(Wa4[p]); sR = np.ptp(Ra4[j:p + 1])/Ra4[p] if Ra4[p] > 0 else np.inf
        m = max(sW, sR)
        if m < best[0]: best = (m, float(t['H0tau'][p]))
        if m <= tol:
            return dict(H0tau=float(t['H0tau'][p]), r=float(t['Wy'][p]/t['rad'][p]), Omega_r=float(t['rad'][p]/H[p]**2),
                        Omega_vac=float(t['vac'][p]/H[p]**2), ln_a=float(lna[p]), idx=p, min_spread=None)
    return dict(H0tau=None, min_spread=float(best[0]), at_H0tau=best[1], reliable_end_H0tau=float(t['H0tau'][ie]))

def at(t, key, tau):
    keep = np.concatenate(([True], np.diff(t['H0tau']) > 1e-12))
    return float(np.interp(tau, t['H0tau'][keep], t[key][keep]))

def r_at(t, tau): return at(t, 'Wy', tau)/at(t, 'rad', tau)

out = dict(status='numerical (audit)', runs={}, cross={}, ancestor={}, fit={}, ledger={}, turnaround={})
files = {Path(f).name.replace('_summary.json', ''): f for f in glob.glob(str(B3/'runs/main/*_summary.json')) + glob.glob(str(B3/'runs/fine/*_summary.json'))}
files.update({'A1_' + Path(f).name.replace('_summary.json', ''): f for f in glob.glob(str(A1/'main_dstar_Y*_summary.json'))})
S = {}
for tag, f in sorted(files.items()):
    s, t = series(f)
    if s['params']['Y'] == 0: continue
    S[tag] = (s, t); pl = plateau(t)
    out['runs'][tag] = dict(Y=s['params']['Y'], dc=s['params']['dc'], dzf=s['params']['dzf'], plateau=pl)
    if pl.get('H0tau') is not None:
        lna = t['ln_a']; H = t['H_over_H0']; ie = rel_end(t)
        out['runs'][tag]['efolds_to_lnamax_within_run'] = float(lna.max() - pl['ln_a'])
        out['runs'][tag]['turnaround_reached'] = bool(np.any(H <= 0))
        # ledger: R a^4 (at plateau) vs Y * int v^2 a^4 dtau (dimensionless: rho_b factors)
        p = pl['idx']; keep = np.concatenate(([True], np.diff(t['H0tau'][:p + 1]) > 1e-12))
        tau = t['H0tau'][:p + 1][keep]; v = t['v_over_H0'][:p + 1][keep]; a4 = np.exp(4*lna[:p + 1][keep])
        Ra4 = t['R'][p]*np.exp(4*lna[p]); R0a40 = t['R'][0]*np.exp(4*lna[0])
        I = np.trapezoid(v**2*a4, tau)
        out['runs'][tag]['ledger_Ra4_over_Y_I_over_rhob'] = float((Ra4 - R0a40)/(s['params']['Y']*I/s['rho_b']))

# cross-resolution r at the decisive plateau times
for key, fin, others in [('Y=2 dc=1e-2', 'fine_Y2_dc1e-2_dzf2.5e-4', ['main_Y2_dc1e-2_dzf5e-4', 'main_Y2_dc1e-2_dzf1e-3']),
                         ('Y=2 dc=1e-4', 'fine_Y2_dc1e-4_dzf2.5e-4', ['main_Y2_dc1e-4_dzf5e-4', 'main_Y2_dc1e-4_dzf1e-3']),
                         ('Y=3 dc=1e-2', 'fine_Y3_dc1e-2_dzf1.25e-4', ['fine_Y3_dc1e-2_dzf2.5e-4', 'main_Y3_dc1e-2_dzf5e-4', 'main_Y3_dc1e-2_dzf1e-3']),
                         ('Y=1.5 dc=1e-2', 'main_Y1.5_dc1e-2_dzf5e-4', ['main_Y1.5_dc1e-2_dzf1e-3']),
                         ('Y=5 dc=1e-2', 'fine_Y5_dc1e-2_dzf1.25e-4_cfl0.25', ['fine_Y5_dc1e-2_dzf2.5e-4_cfl0.25', 'main_Y5_dc1e-2_dzf5e-4_cfl0.25_ext', 'main_Y5_dc1e-2_dzf1e-3_cfl0.25'])]:
    tp = out['runs'][fin]['plateau']['H0tau']
    row = dict(plateau_H0tau=tp, r={fin: r_at(S[fin][1], tp)})
    for o in others:
        row['r'][o] = r_at(S[o][1], tp) if tp <= S[o][1]['H0tau'][rel_end(S[o][1])] else None
        row.setdefault('Wa4_rel_diff_vs_finest', {})[o] = (at(S[o][1], 'Wy', tp)*math.exp(4*at(S[o][1], 'ln_a', tp)) /
                                                          (at(S[fin][1], 'Wy', tp)*math.exp(4*at(S[fin][1], 'ln_a', tp))) - 1)
    out['cross'][key] = row

# common-ancestor error: hybrids share the 5e-4 history up to T_restart; estimate the 5e-4 error there from the primary pair
for key, fin, f5, f10 in [('Y=2 dc=1e-2', 'fine_Y2_dc1e-2_dzf2.5e-4', 'main_Y2_dc1e-2_dzf5e-4', 'main_Y2_dc1e-2_dzf1e-3'),
                          ('Y=2 dc=1e-4', 'fine_Y2_dc1e-4_dzf2.5e-4', 'main_Y2_dc1e-4_dzf5e-4', 'main_Y2_dc1e-4_dzf1e-3'),
                          ('Y=3 dc=1e-2', 'fine_Y3_dc1e-2_dzf1.25e-4', 'main_Y3_dc1e-2_dzf5e-4', 'main_Y3_dc1e-2_dzf1e-3')]:
    ri = S[fin][0]['restart_info']; tr = ri['H0tau']
    d = {}
    for k in ['Wy', 'R', 'phi_b', 'H_over_H0']:
        x5 = at(S[f5][1], k, tr); x10 = at(S[f10][1], k, tr)
        d[k] = dict(rel_diff_1e3_vs_5e4=(x10 - x5)/abs(x5), est_rel_err_5e4_p4=(x10 - x5)/abs(x5)/15)
    out['ancestor'][key] = dict(restart_T=ri['T'], restart_H0tau=tr, diffs=d)

# power law, recomputed; and local crossing of r = 0.03 by log-log interpolation between neighbours
pts = sorted([(0.3, out['runs']['A1_main_dstar_Y0.3_dc1e-2_dzf5e-4']['plateau']['r']), (0.5, out['runs']['main_Y0.5_dc1e-2_dzf5e-4']['plateau']['r']),
              (0.7, out['runs']['main_Y0.7_dc1e-2_dzf5e-4']['plateau']['r']), (1.0, out['runs']['A1_main_dstar_Y1_dc1e-2_dzf5e-4']['plateau']['r']),
              (1.5, out['runs']['main_Y1.5_dc1e-2_dzf5e-4']['plateau']['r']), (2.0, out['runs']['fine_Y2_dc1e-2_dzf2.5e-4']['plateau']['r']),
              (3.0, out['runs']['fine_Y3_dc1e-2_dzf1.25e-4']['plateau']['r'])])
Ys = np.log([p[0] for p in pts]); rs = np.log([p[1] for p in pts]); co = np.polyfit(Ys, rs, 1)
out['fit'] = dict(points=pts, p=float(-co[0]), a=float(math.exp(co[1])), max_rel_dev=float(np.max(np.abs(np.exp(np.polyval(co, Ys) - rs) - 1))),
                  Y_r003_global_fit=float(math.exp((math.log(0.03) - co[1])/co[0])),
                  Y_r003_local_interp_1p5_2=float(math.exp(np.interp(math.log(0.03), rs[[5, 4]], Ys[[5, 4]]))),
                  fit_without_Y3=float(-np.polyfit(Ys[:-1], rs[:-1], 1)[0]))
(Path(__file__).parent/'AUDIT_CHECKS.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
for k, v in out['runs'].items(): print(k, v['plateau'], v.get('efolds_to_lnamax_within_run'), v.get('ledger_Ra4_over_Y_I_over_rhob'))
print(json.dumps(out['cross'], indent=1)); print(json.dumps(out['ancestor'], indent=1)); print(json.dumps(out['fit'], indent=1))
