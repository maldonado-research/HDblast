#!/usr/bin/env python3
"""B1 step 3a (CONDITIONAL): reduced shell-ODE estimate of the derived chi channel on the ARCHIVED A1 tuned trajectory
(d = d*(0.1), Y = 0, fixed background: no scalar backreaction, archived Weyl term).  Output: B1_REDUCED.json.

For each (phi*, G, y): instant-preheating cohort at the archived crossing (N_k = (1 + D/q) exp(-pi k^2/q), D from centred
quartic fits of the archived phi_b(s), h(s), two windows), free chi gas with m = G |phi_b(s) - phi*|, per-mode time-dilated
decay Gamma = y^2 m/(8 pi) into radiation (R_hat' = -4 h R_hat + Gamma m n).  All matter densities are per unit b
(kappa5^2 rho/H0 = b rho_hat).
  s <= s_end (archive):  phi_b(s), h(s), ln a(s), W(s) from the archive (fixed background).
  s >  s_end:            phi frozen at phi_end, W a^4 and vac frozen at their end values, and
      variant 'fixed':   h^2 = vac + W                     (strictly fixed background; r only, Omega_r not meaningful)
      variant 'lowE':    h^2 = vac + W + sig_f b (R+X)/18 + b^2 (R+X)^2/36   (the chi/radiation energy enters H only)
The A1 plateau rule (REGISTRATION section 3 of A1, copied in a1_rules below) is applied to the reduced time series.
b_need(r) = smallest b on a log grid for which the reduced run is PASS at that |r| threshold ('lowE' variant); the M5 cutoff
b_cut = (lambda_c / Lambda_hat)^3 with Lambda_hat = max(m_chi over the trajectory, sqrt(q)) ('full') or post-crossing only
('post').  lambda_c needed = Lambda_hat b_need^(1/3).  Backreaction R_j = b |j_hat| / |r_b sigma'| along the archive.
"""
import json, math, hashlib, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import static_w as S
import shell_matter as SM

A1 = Path('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum/runs/main')
D_STAR = -3.106933495673783
TEN = S.Tension(0.1, S.C_REG, D_STAR)

def load(tag):
    s = json.load(open(A1/(tag + '_summary.json')))
    d = np.load(A1/(tag + '_timeseries.npz'))
    t = {c: d['rec'][:, i] for i, c in enumerate(list(d['cols']))}
    keep = np.concatenate(([True], np.diff(t['H0tau']) > 1e-9))
    t = {k: v[keep] for k, v in t.items()}
    return s, t

def plateau_index(lna, H, Wa4, Ra4, dlna=0.5, tol=0.05):
    """identical rule to A1 a1_analyze.plateau_index (need_R=True)"""
    n = len(lna); j0 = 0
    for p in range(n):
        if lna[p] - lna[0] < dlna: continue
        while j0 < p and lna[j0] < lna[p] - dlna: j0 += 1
        if lna[p] < lna[j0] + dlna - 0.02: continue
        win = slice(j0, p + 1)
        if np.any(H[win] <= 0): continue
        w = Wa4[win]; okW = (w.max() - w.min()) <= tol*abs(Wa4[p]) if Wa4[p] != 0 else False
        rr = Ra4[win]; okR = Ra4[p] > 0 and (rr.max() - rr.min()) <= tol*Ra4[p]
        if okW and okR: return p
    return None

def classify(tau, lna, H, W, rad):
    """A1 per-run rule (plateau, recollapse, Omega_r, |r|) on a reduced series"""
    with np.errstate(divide='ignore', invalid='ignore'):
        Om = rad/H**2; r = np.where(rad > 0, W/rad, np.nan)
    Wa4 = W*np.exp(4*lna); Ra4 = np.where(rad > 0, rad*np.exp(4*lna), 0.0)
    neg = np.flatnonzero(H < -0.05)
    good = np.flatnonzero((Om >= 0.9) & (np.abs(r) <= 0.1))
    p = plateau_index(lna, H, Wa4, Ra4)
    if len(neg) and (len(good) == 0 or good[0] > neg[0]) and (p is None or p > neg[0]):
        return dict(cls='FAIL-no-radiation-era (recollapse)', H0tau=float(tau[neg[0]]), r=None, Omega_r=None)
    if p is None:
        return dict(cls='INCONCLUSIVE (no plateau)', H0tau=float(tau[-1]), r=None, Omega_r=None)
    rr, om = float(r[p]), float(Om[p])
    if abs(rr) > 0.1: c = 'FAIL-Weyl'
    elif om < 0.5: c = 'FAIL-no-radiation-era (Omega_r < 0.5 at plateau)'
    elif om >= 0.9: c = 'PASS-combined' if abs(rr) <= 0.03 else 'PASS-conservative'
    else: c = 'INCONCLUSIVE (0.5 <= Omega_r < 0.9)'
    return dict(cls=c, H0tau=float(tau[p]), r=rr, Omega_r=om, lna=float(lna[p]))

def crossing(t, ps, half=0.1):
    tau, ph = t['H0tau'], t['phi_b']
    i = np.flatnonzero((ph[:-1] - ps)*(ph[1:] - ps) < 0)
    if not len(i): return None
    k = i[0]
    s_star = tau[k] + (ps - ph[k])/(ph[k + 1] - ph[k])*(tau[k + 1] - tau[k])
    sel = np.abs(tau - s_star) <= half
    x = tau[sel] - s_star
    cph = np.polynomial.polynomial.polyfit(x, ph[sel], 4)
    ch = np.polynomial.polynomial.polyfit(x, t['H_over_H0'][sel], 3)
    v, A, B = cph[1], 2*cph[2], 6*cph[3]
    h0, h1 = ch[0], ch[1]
    lna_star = float(np.interp(s_star, tau, t['ln_a']))
    return dict(s_star=float(s_star), v=float(v), A=float(A), B=float(B), h0=float(h0), h1=float(h1), lna_star=lna_star,
                D=float(SM.D_correction(v, A, B, h0, h1)), G_min_screen=float(SM.local_screen(v, A, h0, h1)), n_fit=int(sel.sum()))

def run_case(t, ps, G, y, b_list, nk=24, s_max=400.0):
    cr = crossing(t, ps); cr2 = crossing(t, ps, half=0.05)
    q = G*abs(cr['v']); D = cr['D']
    factor = 1 + D/q if abs(D/q) <= 0.3 else 1 + math.copysign(0.3, D)
    kap, Nfin = SM.instant_nodes(q, nk, 0.0, factor)
    tau, lna, h, ph, W, vac = t['H0tau'], t['ln_a'], t['H_over_H0'], t['phi_b'], t['Wy'], t['vac']
    s_end = tau[-1]
    sp_lna = CubicSpline(tau, lna); sp_h = CubicSpline(tau, h); sp_ph = CubicSpline(tau, ph); sp_W = CubicSpline(tau, W)
    phi_end, W_end, vac_end, lna_end = ph[-1], W[-1], vac[-1], lna[-1]
    Wa4_end = W_end*math.exp(4*lna_end)
    sig_f = TEN.s(phi_end)*t['rb']
    m_of = lambda p: G*abs(p - ps)
    # production budget quantities (per unit b)
    out = dict(phistar=ps, G=G, y=y, q=q, crossing=cr, crossing_half0p05=cr2, D_over_q=D/q, number_factor=factor,
               n_hat_star=float(Nfin.sum()), massless_energy_proxy=q**2/(4*math.pi**4)*factor)
    ms = G*np.abs(ph - ps)
    post = tau >= cr['s_star']
    Lam_full = max(float(ms.max()), math.sqrt(q)); Lam_post = max(float(ms[post].max()), math.sqrt(q))
    out.update(Lambda_hat_full=Lam_full, Lambda_hat_post=Lam_post,
               b_cut_full=(1/Lam_full)**3, b_cut_post=(1/Lam_post)**3)
    lna_star = cr['lna_star']
    def matter(N, ln_a, p):
        ar = math.exp(ln_a - lna_star)
        n = N/ar**3; pk = kap/ar; m = m_of(p); om = np.sqrt(pk*pk + m*m)
        return n, pk, om, m
    def integrate(b, variant):
        # state: N (nk), Rhat*a^4 comoving (Rc), ln_a (only used after s_end)
        def f(s, Y):
            N = Y[:nk]; Rc = Y[nk]; la = Y[nk + 1]
            if s <= s_end:
                la = float(sp_lna(s)); p = float(sp_ph(s))
            else:
                p = phi_end
            n, pk, om, m = matter(N, la, p)
            Gam = y*y*m/(8*math.pi)
            dN = -Gam*(m/om)*N
            dRc = Gam*m*float(n.sum())*math.exp(4*la)
            if s <= s_end:
                dla = float(sp_h(s))
            else:
                Wv = Wa4_end*math.exp(-4*la); X = float((n*om).sum()); R = Rc*math.exp(-4*la)
                h2 = vac_end + Wv + (sig_f*b*(R + X)/18 + (b*(R + X))**2/36 if variant == 'lowE' else 0.0)
                dla = math.sqrt(h2) if h2 > 0 else 0.0
            return np.concatenate([dN, [dRc, dla]])
        def ev_collapse(s, Y):
            if s <= s_end: return 1.0
            N = Y[:nk]; la = Y[nk + 1]; n, pk, om, m = matter(N, la, phi_end)
            X = float((n*om).sum()); R = Y[nk]*math.exp(-4*la)
            return vac_end + Wa4_end*math.exp(-4*la) + (sig_f*b*(R + X)/18 + (b*(R + X))**2/36 if variant == 'lowE' else 0.0)
        ev_collapse.terminal = True; ev_collapse.direction = -1
        Y0 = np.concatenate([Nfin, [0.0, lna_star]])
        s_eval = np.concatenate([tau[tau > cr['s_star']], np.linspace(s_end, s_max, 4000)[1:]])
        atol = np.concatenate([1e-13*Nfin.max()*np.ones(nk), [1e-13*float(Nfin.sum())*G*math.exp(4*lna_star)], [1e-12]])
        sol = solve_ivp(f, (cr['s_star'], s_max), Y0, method='LSODA', t_eval=s_eval, events=ev_collapse, rtol=1e-9, atol=atol)
        ss = sol.t; la = np.where(ss <= s_end, sp_lna(np.minimum(ss, s_end)), sol.y[nk + 1])
        # overwrite ln a after s_end by the integrated value, before by the archive
        la = np.array([float(sp_lna(s)) if s <= s_end else sol.y[nk + 1][i] for i, s in enumerate(ss)])
        Rc = sol.y[nk]; R = Rc*np.exp(-4*la)
        X = np.array([float((matter(sol.y[:nk, i], la[i], float(sp_ph(s)) if s <= s_end else phi_end)[0] *
                             matter(sol.y[:nk, i], la[i], float(sp_ph(s)) if s <= s_end else phi_end)[2]).sum()) for i, s in enumerate(ss)])
        sig = np.array([TEN.s(float(sp_ph(s)) if s <= s_end else phi_end)*t['rb'] for s in ss])
        Wv = np.where(ss <= s_end, sp_W(np.minimum(ss, s_end)), Wa4_end*np.exp(-4*la))
        rad = sig*b*R/18 + (b*R)**2/36
        chi = sig*b*X/18 + ((b*X)**2 + 2*b*R*b*X)/36
        if variant == 'lowE':
            h2 = np.where(ss <= s_end, sp_h(np.minimum(ss, s_end))**2, vac_end + Wv + rad + chi)
        else:
            h2 = np.where(ss <= s_end, sp_h(np.minimum(ss, s_end))**2, vac_end + Wv)
        # after a terminal collapse event, append H < 0 marker
        H = np.sqrt(np.maximum(h2, 0))
        tau_c = None
        if sol.status == 1 and len(sol.t_events[0]):
            tau_c = float(sol.t_events[0][0])
            ss = np.append(ss, tau_c + 1e-6); la = np.append(la, la[-1]); H = np.append(H, -1.0)
            Wv = np.append(Wv, Wv[-1]); rad = np.append(rad, rad[-1]); X = np.append(X, X[-1])
        c = classify(ss, la - la[0] + lna_star, H, Wv, rad)
        c.update(turnaround_H0tau=tau_c, chi_left_fraction_at_end=float(X[-1]/max(X.max(), 1e-300)),
                 Ra4_final=float((rad*np.exp(4*la))[-1]), Wa4_final=float((Wv*np.exp(4*la))[-1]))
        return c, ss, X
    res = {}
    for variant in ['lowE', 'fixed']:
        rows = []
        for b in b_list:
            c, ss, X = integrate(b, variant)
            c['b'] = b; rows.append(c)
        res[variant] = rows
    out['scan'] = res
    # b_need
    for variant in ['lowE']:
        for thr, key in [(0.1, 'b_need_r0p1'), (0.03, 'b_need_r0p03')]:
            ok = [r['b'] for r in res[variant] if r['cls'].startswith('PASS') and r['r'] is not None and abs(r['r']) <= thr]
            out[key] = min(ok) if ok else None
            out[key.replace('b_need', 'lambda_c_full')] = (Lam_full*min(ok)**(1/3)) if ok else None
            out[key.replace('b_need', 'lambda_c_post')] = (Lam_post*min(ok)**(1/3)) if ok else None
    # informative: log-interpolated b at which |r| = 0.1 and 0.03 among plateau rows (lowE)
    pl = [(r['b'], abs(r['r'])) for r in res['lowE'] if r['r'] is not None]
    for thr, key in [(0.1, 'b_interp_r0p1'), (0.03, 'b_interp_r0p03')]:
        out[key] = None
        for (b1, r1), (b2, r2) in zip(pl[:-1], pl[1:]):
            if r1 > thr >= r2:
                out[key] = float(math.exp(math.log(b1) + (math.log(thr) - math.log(r1))*(math.log(b2) - math.log(b1))/(math.log(r2) - math.log(r1))))
                break
    # r at the cutoff (fixed variant: r ~ 1/b exactly in the absence of backreaction)
    fx = [r for r in res['fixed'] if r['r'] is not None]
    out['cls_at_b_cut_full_lowE'] = integrate(out['b_cut_full'], 'lowE')[0]
    out['cls_at_b_cut_post_lowE'] = integrate(out['b_cut_post'], 'lowE')[0]
    return out

def backreaction(t, ps, G, y, b, nk=24):
    """R_j(s) = b |j_hat| / |r_b sigma'(phi)| and the production-work share along the archive (fixed background)"""
    cr = crossing(t, ps); q = G*abs(cr['v']); D = cr['D']
    factor = 1 + D/q if abs(D/q) <= 0.3 else 1 + math.copysign(0.3, D)
    kap, Nfin = SM.instant_nodes(q, nk, 0.0, factor)
    tau, lna, ph = t['H0tau'], t['ln_a'], t['phi_b']
    sel = tau > cr['s_star'] + 3/math.sqrt(q)
    N = Nfin.copy(); out = []
    s_prev = cr['s_star']
    for i in np.flatnonzero(sel):
        s = tau[i]; ar = math.exp(lna[i] - cr['lna_star']); m = G*abs(ph[i] - ps)
        pk = kap/ar; om = np.sqrt(pk*pk + m*m)
        Gam = y*y*m/(8*math.pi); N = N*np.exp(-Gam*(m/om)*(s - s_prev)); s_prev = s
        n = N/ar**3
        j = G*G*(ph[i] - ps)*float((n/om).sum())
        out.append((s, b*abs(j)/abs(t['rb']*TEN.s1(ph[i]))))
    a = np.array(out)
    return dict(max_Rj=float(a[:, 1].max()), H0tau_at_max=float(a[np.argmax(a[:, 1]), 0]), Rj_end=float(a[-1, 1]), H0tau_end=float(a[-1, 0]))

def main():
    b_list = [float(x) for x in np.logspace(-10, -0.5, 39)]
    out = dict(status='conditional (fixed archived background: A1 tuned Y=0 run; no scalar backreaction; archived Weyl term; '
                      'extrapolation beyond the archive with phi frozen, W a^4 and vac frozen)',
               background_runs={}, cases={})
    for tag in ['main_dstar_Y0_dc1e-2_dzf5e-4', 'main_dstar_Y0_dc1e-2_dzf1e-3', 'main_dstar_Y0_dc1e-4_dzf5e-4']:
        s, t = load(tag); t['rb'] = s['rho_b']
        out['background_runs'][tag] = dict(H0tau_end=float(t['H0tau'][-1]), phi_end=float(t['phi_b'][-1]), W_end=float(t['Wy'][-1]),
                                           vac_end=float(t['vac'][-1]), lna_end=float(t['ln_a'][-1]),
                                           Wa4_end=float(t['Wy'][-1]*math.exp(4*t['ln_a'][-1])))
        for ps in [0.5, 0.9]:
            for G in [100.0, 1000.0]:
                for y in [0.1, 1.0]:
                    key = '%s|phistar=%g|G=%g|y=%g' % (tag, ps, G, y)
                    c = run_case(t, ps, G, y, b_list)
                    c['backreaction_at_b_need_r0p1'] = backreaction(t, ps, G, y, c['b_need_r0p1']) if c.get('b_need_r0p1') else None
                    c['backreaction_at_b_cut_full'] = backreaction(t, ps, G, y, c['b_cut_full'])
                    out['cases'][key] = c
                    print('%-34s ps=%.1f G=%5g y=%.1f q=%7.1f D/q=%+.4f Lam=%.1f/%.1f b_cut=%.2e cls_cut=%s r_cut=%s | b_need(0.1)=%s lam_c=%s b_need(0.03)=%s lam_c=%s Rj=%s' % (
                        tag[-14:], ps, G, y, c['q'], c['D_over_q'], c['Lambda_hat_full'], c['Lambda_hat_post'], c['b_cut_full'],
                        c['cls_at_b_cut_full_lowE']['cls'][:18], c['cls_at_b_cut_full_lowE']['r'],
                        c.get('b_need_r0p1'), c.get('lambda_c_full_r0p1'), c.get('b_need_r0p03'), c.get('lambda_c_full_r0p03'),
                        (c['backreaction_at_b_need_r0p1'] or {}).get('max_Rj')), flush=True)
    out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    def clean(v):
        if isinstance(v, dict): return {str(k): clean(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)): return [clean(x) for x in v]
        if isinstance(v, (np.floating, float)): return None if not math.isfinite(float(v)) else float(v)
        if isinstance(v, (np.integer,)): return int(v)
        if isinstance(v, np.bool_): return bool(v)
        return v
    (HERE/'B1_REDUCED.json').write_text(json.dumps(clean(out), indent=1) + '\n')

if __name__ == '__main__':
    main()
