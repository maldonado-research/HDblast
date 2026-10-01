#!/usr/bin/env python3
"""B1 cross-check (numerical/conditional): fixed-background reduced shell ODE versus the 5D runs, over their common range.

For every 5D chi run in runs/main (evolve_b1.py), the same chi cohort (q and number factor taken from the 5D run's own
production event, instant insertion at the 5D crossing time) is evolved on the ARCHIVED A1 tuned Y = 0 background
(main_dstar_Y0_dc1e-2_dzf5e-4: phi_b(s), ln a(s), W(s)) with the shell_matter.ChiGas kinematics (free gas, per-mode
time-dilated decay Gamma = y^2 m/(8 pi) into radiation).  Compared on the 5D record times s = H0 tau up to the earlier of
(i) the 5D run's last record, (ii) the archive's end:
  * R_hat (decay radiation per unit b, H0^4 units): 5D value R * r_b / b versus reduced R_hat;
  * r = W/rad, 5D (own W, own sigma) versus reduced (archived W, archived sigma);
  * Wa4 drift in the 5D run after the decay (W a^4 at the last record over its value at the first time X/R < 0.01).
Purpose: (a) at b_cut (test-field regime, R_j << 1) it checks the fixed-background reduction against the full 5D solution;
(b) at b_need it measures how much the scalar backreaction changes R and r.  Output: B1_REDUCED_VS_5D.json."""
import json, glob, math, hashlib, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import static_w as S
import shell_matter as SM
A1 = Path('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum/runs/main')
TEN = S.Tension(0.1, S.C_REG, -3.106933495673783)

def load_npz(p):
    d = np.load(p); t = {str(c): d['rec'][:, i] for i, c in enumerate(list(d['cols']))}
    keep = np.concatenate(([True], np.diff(t['H0tau']) > 1e-9))
    return {k: v[keep] for k, v in t.items()}

def main():
    bg_tag = 'main_dstar_Y0_dc1e-2_dzf5e-4'
    sb = json.load(open(A1/(bg_tag + '_summary.json'))); tb = load_npz(A1/(bg_tag + '_timeseries.npz'))
    rb = sb['rho_b']
    sp_lna, sp_ph, sp_W = CubicSpline(tb['H0tau'], tb['ln_a']), CubicSpline(tb['H0tau'], tb['phi_b']), CubicSpline(tb['H0tau'], tb['Wy'])
    s_end_bg = float(tb['H0tau'][-1])
    out = dict(status='numerical cross-check; reduced side is conditional (fixed archived background)', background=bg_tag,
               background_H0tau_end=s_end_bg, runs={})
    for f in sorted(glob.glob(str(HERE/'runs/main/main_*_summary.json'))):
        s = json.load(open(f)); m = s.get('matter', {})
        if m.get('kind') != 'chi' or not m.get('events'): continue
        tag = s['tag']; ev = m['events'][0]
        t = load_npz(f.replace('_summary.json', '_timeseries.npz'))
        G, ps, y, b = m['G'], m['phistar'], m['y'], m['b']
        q, factor, s0 = ev['q'], ev['number_factor'], ev['s_star']
        kap, Nfin = SM.instant_nodes(q, 24, 0.0, factor)
        la0 = float(sp_lna(s0)); nk = len(Nfin)
        def rhs(sv, Y):
            N = Y[:nk]; la = float(sp_lna(sv)); p = float(sp_ph(sv)); ar = math.exp(la - la0)
            mm = G*abs(p - ps); pk = kap/ar; om = np.sqrt(pk*pk + mm*mm); Gam = y*y*mm/(8*math.pi)
            dN = -Gam*(mm/om)*N
            return np.concatenate([dN, [Gam*mm*float((N/ar**3).sum())*math.exp(4*la)]])
        s_stop = min(float(t['H0tau'][-1]), s_end_bg)
        sel = (t['H0tau'] > s0 + 0.05) & (t['H0tau'] <= s_stop)
        if sel.sum() < 5: continue
        se = t['H0tau'][sel]
        sol = solve_ivp(rhs, (s0, se[-1]), np.concatenate([Nfin, [0.0]]), t_eval=se, method='LSODA', rtol=1e-9, atol=1e-14)
        la = sp_lna(se); R_red = sol.y[nk]*np.exp(-4*la)
        sig_bg = np.array([TEN.s(float(p)) for p in sp_ph(se)])*rb
        rad_red = sig_bg*b*R_red/18 + (b*R_red)**2/36
        r_red = sp_W(se)/rad_red
        R5 = t['R'][sel]*rb/b
        r5 = t['Wy'][sel]/t['rad'][sel]
        # compare where the radiation is established (R_red above 10% of its max)
        est = R_red > 0.1*R_red.max()
        ratio_R = R5[est]/R_red[est]; ratio_r = r5[est]/r_red[est]
        # 5D W a^4 drift after the decay
        X = t['X']; Rr = t['R']; Wa4 = t['Wy']*np.exp(4*t['ln_a'])
        post = np.flatnonzero((t['H0tau'] > s0) & (X < 0.01*np.maximum(Rr, 1e-300)) & (Rr > 0))
        drift = float(Wa4[-1]/Wa4[post[0]]) if len(post) else None
        Ra4 = Rr*np.exp(4*t['ln_a'])
        driftR = float(Ra4[-1]/Ra4[post[0]]) if len(post) else None
        out['runs'][tag] = dict(G=G, phistar=ps, y=y, b=b, common_H0tau=[float(se[0]), float(se[-1])], n=int(est.sum()),
                                R5_over_Rred=dict(min=float(ratio_R.min()), max=float(ratio_R.max()), last=float(ratio_R[-1])),
                                r5_over_rred=dict(min=float(ratio_r.min()), max=float(ratio_r.max()), last=float(ratio_r[-1])),
                                r_5D_last=float(r5[-1]), r_red_last=float(r_red[-1]), H0tau_last=float(se[-1]),
                                Omega_r_5D_last=float(t['rad'][sel][-1]/t['H_over_H0'][sel][-1]**2),
                                Wa4_drift_after_decay=drift, Ra4_drift_after_decay=driftR,
                                H0tau_decay_done=(float(t['H0tau'][post[0]]) if len(post) else None),
                                dlna_after_decay=(float(t['ln_a'][-1] - t['ln_a'][post[0]]) if len(post) else None))
        o = out['runs'][tag]
        print('%-42s R5/Rred=[%.3f,%.3f] last %.3f  r5/rred last %.3f  r5=%.4g rred=%.4g at %.2f  Wa4 drift %s Ra4 drift %s dlna %s' % (
            tag, o['R5_over_Rred']['min'], o['R5_over_Rred']['max'], o['R5_over_Rred']['last'], o['r5_over_rred']['last'], o['r_5D_last'],
            o['r_red_last'], o['H0tau_last'], drift, driftR, o['dlna_after_decay']))
    out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (HERE/'B1_REDUCED_VS_5D.json').write_text(json.dumps(out, indent=1) + '\n')

if __name__ == '__main__':
    main()
