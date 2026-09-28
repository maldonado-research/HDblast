#!/usr/bin/env python3
"""Audit script 2: an independent chi mode-function solver (does NOT import ph_lib).

Differences from the audited code: complex-valued ODE state; vacuum and read-out basis built from exact
analytic derivatives with the second-order adiabatic frequency and a finite-difference W'; windows set by a fixed
adiabaticity epsilon (two values compared); 64-node Gauss-Legendre in x = kappa*sqrt(pi/q) on [0, 6];
Richardson extrapolation in q for the O(1/q) coefficient D.

Checks:
 (a) flat linear crossing vs exact n_k = exp(-pi(k^2+mu0^2)/q)  (calibration of THIS solver)
 (b) D on single-parameter toys in the audited fit set (reproduction)
 (c) D on OUT-OF-SAMPLE toys with all five structures on at once and different v (test of the closed form)
 (d) N on the archived trajectory (phi*=0.5, 0.75 at G=1e3) with an independent spline loader, vs the audited JSON
Writes verify/INDEP_MODES.json.
"""
import io, json, math, sys, time, zipfile
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import make_interp_spline

HERE = Path(__file__).resolve().parent
PI = math.pi


class Toy:
    """phi-phi* = v s + A s^2/2 + B s^3/6, h = h0 + h1 s, ln a = h0 s + h1 s^2/2 (a=1 at s=0)."""
    def __init__(self, v, A=0.0, B=0.0, h0=0.0, h1=0.0, span=1.4):
        self.v, self.A, self.B, self.h0, self.h1 = v, A, B, h0, h1
        self.lo, self.hi = -span, span
    def d(self, s):   # (Delta, Delta', Delta'')
        v, A, B = self.v, self.A, self.B
        return v*s + A*s*s/2 + B*s**3/6, v + A*s + B*s*s/2, A + B*s
    def hh(self, s):  # (h, h', h'', h''')
        return self.h0 + self.h1*s, self.h1, 0.0, 0.0
    def lna(self, s): return self.h0*s + self.h1*s*s/2
    def s_cross(self): return 0.0


class Archived:
    """Independent loader for the archived Chat14 trajectory: phi(s), h(s) as separate quintic splines,
    ln a from Gauss quadrature of h between nodes (no antiderivative object reused from the audited code)."""
    ZIP = Path('/home/user/unified-theory-maldonado/new-files/latest-work/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/'
               'HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/source_audit/inputs/CHAT14_COMPLETED_REFERENCE.zip')
    def __init__(self, tag, phistar, t_cut=8.5):
        with zipfile.ZipFile(self.ZIP) as z:
            pre = 'HDBLAST_CHAT14_REGISTERED_DETUNING_ROLLOFF_20260922/runs/'
            rec = np.load(io.BytesIO(z.read(pre + tag + '_timeseries.npz')))['rec']
            meta = json.loads(z.read(pre + tag + '_summary.json'))
        t_cut = min(t_cut, meta['L']*(1 - meta['taper']))
        rec = rec[rec[:, 0] < t_cut]
        self.s = rec[:, 9]; self.phistar = phistar
        self.fphi = make_interp_spline(self.s, rec[:, 1], k=5)
        self.fh = make_interp_spline(self.s, rec[:, 2], k=5)
        self.lo, self.hi = float(self.s[0]), float(self.s[-1])
        from scipy.optimize import brentq
        y = rec[:, 1] - phistar; i = int(np.flatnonzero(np.sign(y[:-1]) != np.sign(y[1:]))[0])
        self.sc = brentq(lambda x: float(self.fphi(x)) - phistar, self.s[i], self.s[i+1], xtol=1e-14)
        # ln a relative to the crossing, tabulated with 8-point Gauss per interval, then interpolated
        xg, wg = np.polynomial.legendre.leggauss(8)
        seg = np.zeros(len(self.s))
        for j in range(1, len(self.s)):
            a_, b_ = self.s[j-1], self.s[j]
            seg[j] = seg[j-1] + 0.5*(b_-a_)*np.sum(wg*self.fh(0.5*(b_-a_)*xg + 0.5*(a_+b_)))
        self.flna = make_interp_spline(self.s, seg, k=5)
        self.lna0 = float(self.flna(self.sc))
    def d(self, s):
        return float(self.fphi(s)) - self.phistar, float(self.fphi(s, 1)), float(self.fphi(s, 2))
    def hh(self, s): return tuple(float(self.fh(s, n)) for n in range(4))
    def lna(self, s): return float(self.flna(s)) - self.lna0
    def s_cross(self): return self.sc
    def v_cross(self): return float(self.fphi(self.sc, 1))


def omega2(bg, s, kap, G, mu0):
    D, D1, D2 = bg.d(s); h, h1, h2, h3 = bg.hh(s)
    k2 = kap*kap*math.exp(-2*bg.lna(s))
    O2 = k2 + mu0*mu0 + G*G*D*D - 2.25*h*h - 1.5*h1
    dO2 = -2*h*k2 + 2*G*G*D*D1 - 4.5*h*h1 - 1.5*h2
    ddO2 = (4*h*h - 2*h1)*k2 + 2*G*G*(D1*D1 + D*D2) - 4.5*(h1*h1 + h*h2) - 1.5*h3
    return O2, dO2, ddO2


def W2(bg, s, kap, G, mu0):
    O2, dO2, ddO2 = omega2(bg, s, kap, G, mu0)
    # W^2 = O^2 - (1/2)(O''/O) + (3/4)(O'/O)^2 written in terms of O2 = Omega^2:
    # O'/O = dO2/(2 O2);  O''/O = ddO2/(2 O2) - (dO2)^2/(4 O2^2)
    r1 = dO2/(2*O2); r2 = ddO2/(2*O2) - dO2**2/(4*O2**2)
    return O2 - 0.5*r2 + 0.75*r1**2


def basis(bg, s, kap, G, mu0, eps=1e-6):
    W = np.sqrt(W2(bg, s, kap, G, mu0))
    Wp = (np.sqrt(W2(bg, s + eps, kap, G, mu0)) - np.sqrt(W2(bg, s - eps, kap, G, mu0)))/(2*eps)
    f = 1/np.sqrt(2*W); fp = (-1j*W - Wp/(2*W))*f
    return f, fp


def adiab(bg, s, G, mu0):
    O2, dO2, ddO2 = omega2(bg, s, 0.0, G, mu0)
    if O2 <= 0: return np.inf
    return max(abs(dO2)/(2*O2**1.5), math.sqrt(abs(ddO2)/(2*O2*O2)))


def edges(bg, G, mu0, eps):
    sc = bg.s_cross(); out = []
    for lim in (bg.lo, bg.hi):
        xs = np.linspace(sc, lim, 3001)
        A = np.array([adiab(bg, x, G, mu0) for x in xs])
        bad = np.flatnonzero(A > eps)
        j = bad[-1] + 1 if len(bad) else 1
        out.append(float(xs[min(j, len(xs)-1)]))
    return out


def number(bg, G, mu0=0.0, eps=1e-4, nk=64, xmax=6.0, rtol=1e-11):
    v = abs(bg.d(bg.s_cross())[1]); q = G*v
    x, w = np.polynomial.legendre.leggauss(nk)
    x = 0.5*xmax*(x + 1); w = 0.5*xmax*w
    kap = x*math.sqrt(q/PI); wk = w*math.sqrt(q/PI)
    sa, sb = edges(bg, G, mu0, eps)
    f, fp = basis(bg, sa, kap, G, mu0)
    y0 = np.concatenate([f, fp]).astype(complex)
    n = len(kap)
    def rhs(s, y):
        O2 = omega2(bg, s, kap, G, mu0)[0]
        return np.concatenate([y[n:], -O2*y[:n]])
    sol = solve_ivp(rhs, (sa, sb), y0, method='DOP853', rtol=rtol, atol=1e-15)
    X, P = sol.y[:n, -1], sol.y[n:, -1]
    g, gp = basis(bg, sb, kap, G, mu0)
    # X = alpha g + beta g*;  Wronskian g gp* - g* gp = i
    beta = (g*P - gp*X)/(g*np.conj(gp) - np.conj(g)*gp)
    alpha = (np.conj(g)*P - np.conj(gp)*X)/(np.conj(g)*gp - g*np.conj(gp))
    nk_ = np.abs(beta)**2
    N = float(np.sum(wk*kap**2*nk_))/(2*PI**2)
    N0 = q**1.5/(8*PI**3)*math.exp(-PI*mu0**2/q)
    return {'q': q, 'N': N, 'N_instant': N0, 'rel': N/N0 - 1, 'edges': [sa, sb], 'nfev': int(sol.nfev),
            'wronskian_dev': float(np.max(np.abs(np.abs(alpha)**2 - np.abs(beta)**2 - 1))),
            'max_rel_nk_vs_exact_where_gt_1e-6': float(np.max(np.abs(nk_/np.exp(-PI*(kap**2 + mu0**2)/q) - 1)[np.exp(-PI*(kap**2+mu0**2)/q) > 1e-6]))}


def closed_form(v, A=0.0, B=0.0, h0=0.0, h1=0.0):
    c = {'hh': 9*PI/4 + 15/(2*PI), 'h1': 3*PI/2 - 15/(8*PI), 'AA': 45/(32*PI) - PI/8, 'B': PI/8 - 45/(96*PI), 'hA': 45/(8*PI)}
    return c['hh']*h0*h0 + c['h1']*h1 + c['AA']*(A/v)**2 + c['B']*B/v + c['hA']*h0*A/v


def D_richardson(params, q_targets=(2400.0, 9600.0), eps=1e-4):
    vals = []
    for qt in q_targets:
        r = number(Toy(**params), qt/abs(params['v']), eps=eps)
        vals.append(r)
    (q1, d1), (q2, d2) = [(r['q'], r['rel']*r['q']) for r in vals]
    return {'params': params, 'D_at_q': [d1, d2], 'q': [q1, q2], 'D_extrap': (q2*d2 - q1*d1)/(q2 - q1),
            'D_closed_form': closed_form(**params), 'wronskian_dev': max(r['wronskian_dev'] for r in vals),
            'edges': [r['edges'] for r in vals]}


def main():
    t0 = time.time(); out = {}
    # (a) calibration of this solver
    cal = {}
    for q in (10.0, 100.0, 1160.0):
        for mu0 in (0.0, 0.5*math.sqrt(q)):
            span = max(20.0, 3*math.sqrt(1/q)*40)
            r = number(Toy(v=0.5, span=span), q/0.5, mu0=mu0, eps=1e-4)
            cal['q=%g|mu0=%.3g' % (q, mu0)] = {k: r[k] for k in ('rel', 'max_rel_nk_vs_exact_where_gt_1e-6', 'wronskian_dev')}
    out['calibration_flat_linear'] = cal
    print('calibration', json.dumps(cal), flush=True)
    # (b) in-sample toys (audited fit set, v=0.6)
    ins = {'flat': dict(v=0.6), 'A=0.3': dict(v=0.6, A=0.3), 'B=0.4': dict(v=0.6, B=0.4), 'h0=0.5': dict(v=0.6, h0=0.5),
           'h1=0.3': dict(v=0.6, h1=0.3), 'h0=0.8,A=0.4': dict(v=0.6, h0=0.8, A=0.4)}
    out['in_sample'] = {}
    for k, p in ins.items():
        out['in_sample'][k] = D_richardson(p)
        print('in', k, out['in_sample'][k]['D_extrap'], out['in_sample'][k]['D_closed_form'], flush=True)
    # (c) out-of-sample toys: all structures at once, different v
    oos = {'mix1': dict(v=0.4, A=0.25, B=-0.5, h0=0.7, h1=-0.2),
           'mix2': dict(v=0.8, A=-0.6, B=0.9, h0=0.3, h1=0.25),
           'mix3_trajlike': dict(v=0.5786840415164352, A=0.2863880601008191, B=-1.0512888681550976,
                                 h0=0.8878608255431975, h1=-0.19628251354630122)}
    out['out_of_sample'] = {}
    for k, p in oos.items():
        out['out_of_sample'][k] = D_richardson(p)
        print('oos', k, out['out_of_sample'][k]['D_extrap'], out['out_of_sample'][k]['D_closed_form'], flush=True)
    # window-epsilon sensitivity on one mixed toy
    out['eps_sensitivity_mix1'] = D_richardson(oos['mix1'], eps=1e-3)['D_extrap']
    # (d) archived trajectory, independent loader
    R = json.loads((HERE.parent/'PREHEATING_RESULTS.json').read_text())
    traj = {}
    for ps in (0.5, 0.75):
        bg = Archived('v2_t1e3_plus_c4_finecoarse', ps)
        for G in (1e3, 1e4):
            r = number(bg, G, eps=1e-3)
            ref = R['primary']['%s|%g' % (ps, G)]
            traj['%s|%g' % (ps, G)] = {'N_indep': r['N'], 'N_audited': ref['analytic']['N_comoving_numeric'],
                                      'rel_diff': r['N']/ref['analytic']['N_comoving_numeric'] - 1,
                                      's_cross_indep': bg.sc, 's_cross_audited': ref['s_star'], 'q_indep': r['q'], 'q_audited': ref['q'],
                                      'edges': r['edges'], 'wronskian_dev': r['wronskian_dev']}
            print('traj', ps, G, traj['%s|%g' % (ps, G)], flush=True)
    out['archived_trajectory'] = traj
    out['runtime_s'] = round(time.time() - t0, 1)
    (HERE/'INDEP_MODES.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
    print('done', out['runtime_s'])


if __name__ == '__main__':
    main()
