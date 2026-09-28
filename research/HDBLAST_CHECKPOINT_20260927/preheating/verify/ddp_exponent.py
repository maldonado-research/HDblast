#!/usr/bin/env python3
"""Audit script 4: independent test of the 'exponent-only' (complex-turning-point) coefficients quoted in the
preheating README (section 5):  c_AA_exp = 45/(32 pi), c_B_exp = -45/(96 pi), redshift parts +15/(4 pi) h0^2 and
-15/(8 pi) h1.  Mass-shift terms are left out of Omega^2 here (they are exact and treated separately in sympy_checks).

For Omega^2(s) = kappa^2 exp(-2 ln a(s)) + G^2 Delta(s)^2 the exponent is E(kappa) = 4 Im int_0^{s_c} Omega ds with
s_c the complex turning point near i kappa/q.  N_DDP = (1/2pi^2) int kappa^2 exp(-E) dkappa and
D_exp = lim q (N_DDP/N_0 - 1), from Richardson extrapolation over two q values.  Precision: mpmath, 30 digits.
Writes verify/DDP_EXPONENT.json."""
import json, math
from pathlib import Path
import mpmath as mp
import numpy as np

mp.mp.dps = 30
HERE = Path(__file__).resolve().parent


def E_of(kap, G, v, A, B, h0, h1):
    kap = mp.mpf(kap)
    def O2(s):
        d = v*s + A*s**2/2 + B*s**3/6
        return kap**2*mp.exp(-2*(h0*s + h1*s**2/2)) + G**2*d**2
    q = G*v
    sc = mp.findroot(O2, mp.mpc(0, kap/q))
    I = sc*mp.quad(lambda t: mp.sqrt(O2(t*sc)), [0, 1])
    return 4*mp.im(I), sc


def D_exp(params, qs=(2.0e4, 8.0e4), nk=48, xmax=6.0):
    v = params['v']; vals = []
    x, w = np.polynomial.legendre.leggauss(nk)
    x = 0.5*xmax*(x + 1); w = 0.5*xmax*w
    for q in qs:
        G = q/v
        sq = math.sqrt(q/math.pi)
        tot = mp.mpf(0); tot0 = mp.mpf(0)
        for xi, wi in zip(x, w):
            kap = xi*sq
            E, _ = E_of(kap, G, v, params.get('A', 0), params.get('B', 0), params.get('h0', 0), params.get('h1', 0))
            tot += wi*kap**2*mp.exp(-E); tot0 += wi*kap**2*mp.exp(-mp.pi*kap**2/q)
        vals.append((q, float((tot/tot0 - 1)*q)))
    (q1, d1), (q2, d2) = vals
    return {'params': params, 'D_at_q': [d1, d2], 'q': [q1, q2], 'D_extrap': (q2*d2 - q1*d1)/(q2 - q1)}


def main():
    v = 0.6
    pred = {'A': lambda A: 45/(32*math.pi)*(A/v)**2, 'B': lambda B: -45/(96*math.pi)*B/v,
            'h0': lambda h: 15/(4*math.pi)*h*h, 'h1': lambda h: -15/(8*math.pi)*h}
    out = {'note': 'exponent only; mass-shift contributions excluded', 'runs': {}}
    # control: flat linear must give 0
    out['runs']['flat'] = D_exp({'v': v}); out['runs']['flat']['predicted'] = 0.0
    for key, val in (('A', 0.3), ('A', -0.3), ('B', 0.4), ('h0', 0.5), ('h1', 0.3), ('h1', -0.3)):
        r = D_exp({'v': v, key: val}); r['predicted'] = pred[key](val)
        out['runs']['%s=%g' % (key, val)] = r
        print(key, val, r['D_extrap'], r['predicted'], flush=True)
    # cross term h0*A (the README does not state whether 45/(8 pi) is exponent-only; recorded for information)
    r = D_exp({'v': v, 'h0': 0.8, 'A': 0.4})
    r['sum_of_single_predictions'] = pred['h0'](0.8) + pred['A'](0.4)
    r['implied_cross_coefficient'] = (r['D_extrap'] - r['sum_of_single_predictions'])/(0.8*0.4/v)
    r['closed_form_cross_coefficient_45_over_8pi'] = 45/(8*math.pi)
    out['runs']['h0=0.8,A=0.4'] = r
    print('cross', r, flush=True)
    for k, r in out['runs'].items():
        if 'predicted' in r: r['abs_dev'] = r['D_extrap'] - r['predicted']
    (HERE/'DDP_EXPONENT.json').write_text(json.dumps(out, indent=1, default=float) + '\n')
    print(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
