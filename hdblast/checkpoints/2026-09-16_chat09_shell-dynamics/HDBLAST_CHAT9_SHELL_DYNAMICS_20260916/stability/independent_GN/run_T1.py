"""T1: difference of neighbouring exact bulk solutions (vary phi_h) = exact mu2=0, Y=const GN-gauge perturbation
(chi = d phi/d alpha, psi = d ln rho/d alpha, E = 0).  For constant Y the momentum constraint (a 4-gradient) and the
trace-free evolution equation (coefficient of nabla nabla Y) are vacuous; HC, scalar eq., trace-evolution eq. and yy eq. must hold."""
import numpy as np, json
from gn_lib import *
import sys
res = {}
for (al, yb, n, rel) in [(8.785514640315845e-07, 8.228201550106348, 2000, 1e-3), (8.785514640315845e-07, 8.228201550106348, 2000, 1e-4), (0.3, 2.5, 1500, 1e-5), (0.3, 2.5, 1500, 1e-6)]:
    Cp = bg_coeffs(background(-1 + al*(1+rel), yb, n)); Cm = bg_coeffs(background(-1 + al*(1-rel), yb, n)); C0 = bg_coeffs(background(-1 + al, yb, n))
    D = lambda f: (f(Cp) - f(Cm))/(2*rel)       # d/d ln(alpha)
    chi = D(lambda C: C['phi']); dchi = D(lambda C: C['s']); ddchi = D(lambda C: C['spp'])
    psi = D(lambda C: np.log(C['rho'])); P = D(lambda C: C['H']); dP = D(lambda C: -C['ir2'] - C['s']**2/3)
    H, ir2, s, u1, u2 = C0['H'], C0['ir2'], C0['s'], C0['U1'], C0['U2']
    f = lambda T: np.abs(sum(T))/(sum(abs(t) for t in T) + 1e-300)
    hc = f([12*H*P, 12*psi*ir2, -s*dchi, u1*chi])
    sc = f([ddchi, 4*H*dchi, 4*P*s, -u2*chi])
    ev = f([dP, 8*H*P, 6*psi*ir2, (2/3)*u1*chi])
    yy = f([-4*dP, -8*H*P, -2*s*dchi, -(2/3)*u1*chi])
    sel = C0['y'] > 0.05
    print('  phi range', C0['phi'][0], C0['phi'][-1], ' profile y: HC SC EV YY'); [print('   %.2f %.1e %.1e %.1e %.1e' % (C0['y'][i], hc[i], sc[i], ev[i], yy[i])) for i in [np.argmin(abs(C0['y']-q)) for q in (0.1, 0.5, 1, 2, yb-0.03)]]
    res['alpha=%g,rel=%g' % (al, rel)] = dict(HC=float(hc[sel].max()), SC=float(sc[sel].max()), EVpsi=float(ev[sel].max()), YY=float(yy[sel].max()))
    print('alpha', al, 'rel step', rel, res['alpha=%g,rel=%g' % (al, rel)])
json.dump(res, open('T1_result.json', 'w'), indent=1)
