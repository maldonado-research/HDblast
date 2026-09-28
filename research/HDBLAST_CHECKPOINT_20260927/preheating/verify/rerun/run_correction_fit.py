#!/usr/bin/env python3
"""O(1/q) correction to the instant-preheating number on a curved, expanding crossing.  Writes CORRECTION_FIT.json.

Local toy backgrounds around the crossing (s measured from the crossing, a(0)=1):
    phi - phi* = v s + A s^2/2 + B s^3/6 ,     h(s) = h0 + h1 s ,    ln a = h0 s + h1 s^2/2
The produced comoving number is written N = q^{3/2}/(8 pi^3) (1 + D/q + O(q^-2)), q = G|v|.
We measure D numerically for single-parameter toys and fit
    D = c_hh h0^2 + c_h1 h1 + c_AA (A/v)^2 + c_B (B/v) + c_hA h0 (A/v).
A complex-turning-point (Dykhne) expansion of the exponent alone suggests, for the flat pieces,
c_AA = 45/(32 pi), c_B = -45/(96 pi); the effective-mass shift -(9/4)h^2 - (3/2)h' alone would give pi*9/4 and pi*3/2 for c_hh, c_h1.
These analytic numbers are hypotheses to be tested here, not assumed.
Finally the fitted law is used to predict D on the archived trajectory crossings (compared with CONTROLS.json).
"""
import json, math, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ph_lib as L


class PolyToy:
    def __init__(self, v, A=0.0, B=0.0, h0=0.0, h1=0.0, span=1.4):
        self.v, self.A, self.B, self.h0, self.h1 = v, A, B, h0, h1
        self.s_min, self.s_max, self.tag = -span, span, 'poly_toy'
    def bg(self, s, nu=0):
        v, A, B, h0, h1 = self.v, self.A, self.B, self.h0, self.h1
        s = float(s)
        if nu == 0: return np.array([v*s + A*s*s/2 + B*s**3/6, h0 + h1*s])
        if nu == 1: return np.array([v + A*s + B*s*s/2, h1])
        if nu == 2: return np.array([A + B*s, 0.0])
        if nu == 3: return np.array([B, 0.0])
        return np.array([0.0, 0.0])
    def lna(self, s): return self.h0*s + self.h1*s*s/2
    def crossing(self, phistar): return 0.0


def D_of(toy, G):
    c = L.run_case(toy, G, 0.0, tol=1e-3)
    return c['analytic']['rel_diff_number']*c['q'], c


def measure(params, Gs=(4e3, 1.6e4)):
    """Richardson: D(q) = D + E/q  ->  D = (q2 D2 - q1 D1)/(q2 - q1)."""
    vals = []
    for G in Gs:
        toy = PolyToy(**params)
        d, c = D_of(toy, G)
        vals.append((c['q'], d, c['window']))
    (q1, d1, w1), (q2, d2, w2) = vals
    Dinf = (q2*d2 - q1*d1)/(q2 - q1)
    return {'params': params, 'q': [q1, q2], 'D_at_q': [d1, d2], 'D_extrapolated': Dinf,
            'windows': [w1, w2]}


def main():
    v = 0.6
    runs = {}
    runs['flat_linear'] = measure(dict(v=v))
    for A in (0.3, -0.3, 0.6):
        runs['A=%g' % A] = measure(dict(v=v, A=A))
    for B in (0.4, -0.4):
        runs['B=%g' % B] = measure(dict(v=v, B=B))
    for h0 in (0.5, 1.0):
        runs['h0=%g' % h0] = measure(dict(v=v, h0=h0))
    for h1 in (-0.3, 0.3):
        runs['h1=%g' % h1] = measure(dict(v=v, h1=h1))
    runs['h0=0.8,A=0.4'] = measure(dict(v=v, h0=0.8, A=0.4))
    runs['h0=0.8,A=-0.4'] = measure(dict(v=v, h0=0.8, A=-0.4))
    # least squares on the regressors
    keys = [k for k in runs if k != 'flat_linear']
    X = []; y = []
    for k in keys:
        p = runs[k]['params']
        a = p.get('A', 0)/v; b = p.get('B', 0)/v; h0 = p.get('h0', 0); h1 = p.get('h1', 0)
        X.append([h0*h0, h1, a*a, b, h0*a]); y.append(runs[k]['D_extrapolated'])
    X = np.array(X); y = np.array(y)
    coef, res, rank, sv = np.linalg.lstsq(X, y, rcond=None)
    names = ['c_hh', 'c_h1', 'c_AA', 'c_B', 'c_hA']
    fit = dict(zip(names, coef.tolist()))
    resid = (X @ coef - y).tolist()
    # odd-in-A check: flat A=+0.3 and A=-0.3 should agree if no linear-in-A term
    out = {'status': 'NUMERICAL', 'v': v, 'runs': runs, 'fit': fit, 'fit_residuals': resid,
           'fit_max_abs_residual': float(np.max(np.abs(resid))),
           'flat_linear_D_should_be_0': runs['flat_linear']['D_extrapolated'],
           'analytic_hypotheses': {'c_AA_exponent_only': 45/(32*math.pi), 'c_B_exponent_only': -45/(96*math.pi),
                                   'c_hh_mass_shift_only': 9*math.pi/4, 'c_h1_mass_shift_only': 3*math.pi/2}}
    # closed-form identification (conjectural for c_AA, c_B, c_hA and part of c_hh; see README)
    closed = {'c_hh': 9*math.pi/4 + 15/(2*math.pi), 'c_h1': 3*math.pi/2 - 15/(8*math.pi),
              'c_AA': 45/(32*math.pi) - math.pi/8, 'c_B': math.pi/8 - 45/(96*math.pi), 'c_hA': 45/(8*math.pi)}
    out['closed_form_candidates'] = closed
    out['closed_form_minus_fit'] = {k: closed[k] - fit[k] for k in closed}
    ccoef = np.array([closed[k] for k in names])
    out['closed_form_max_abs_residual_on_toys'] = float(np.max(np.abs(X @ ccoef - y)))
    # wrong-formula control: exponent-only (Dykhne) flat coefficients
    wrong = dict(closed); wrong['c_AA'] = 45/(32*math.pi); wrong['c_B'] = -45/(96*math.pi)
    out['wrong_exponent_only_max_abs_residual_on_toys'] = float(np.max(np.abs(X @ np.array([wrong[k] for k in names]) - y)))
    # prediction on the archived trajectory
    tr = L.Trajectory('v2_t1e3_plus_c4_finecoarse')
    ctrl = json.loads((HERE/'CONTROLS.json').read_text())['asymptotic_scaling']
    pred = {}
    for ps in (0.25, 0.5, 0.75, 0.9):
        s0 = tr.crossing(ps)
        p1, h0 = tr.bg(s0, 0)[1], tr.bg(s0, 0)[1]
        vv = float(tr.bg(s0, 1)[0]); A = float(tr.bg(s0, 2)[0]); B = float(tr.bg(s0, 3)[0])
        h0 = float(tr.bg(s0, 0)[1]); h1 = float(tr.bg(s0, 1)[1])
        reg = np.array([h0*h0, h1, (A/vv)**2, B/vv, h0*A/vv])
        Dp = float(reg @ coef)
        entry = {'v': vv, 'A': A, 'B': B, 'h0': h0, 'h1': h1, 'D_predicted_by_fit': Dp,
                 'D_predicted_by_closed_form': float(reg @ ccoef)}
        if str(ps) in ctrl:
            Gq = ctrl[str(ps)]['G'][-1]
            entry['D_measured_on_trajectory'] = ctrl[str(ps)]['G_times_rel_diff'][-1]*abs(vv)   # (rel*G)*|v| = rel*q
            entry['rel_dev_prediction'] = Dp/entry['D_measured_on_trajectory'] - 1
        pred[str(ps)] = entry
    out['trajectory_prediction'] = pred
    (HERE/'CORRECTION_FIT.json').write_text(json.dumps(out, indent=1, allow_nan=False, default=float) + '\n')
    print(json.dumps({'closed_minus_fit': out['closed_form_minus_fit'], 'closed_resid': out['closed_form_max_abs_residual_on_toys'],
                      'wrong_resid': out['wrong_exponent_only_max_abs_residual_on_toys'], 'fit': fit, 'max_resid': out['fit_max_abs_residual'], 'flat': out['flat_linear_D_should_be_0'],
                      'hyp': out['analytic_hypotheses'], 'pred': pred}, indent=1))


if __name__ == '__main__':
    main()
