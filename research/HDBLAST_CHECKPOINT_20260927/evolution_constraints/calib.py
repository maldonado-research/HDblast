#!/usr/bin/env python3
"""Order calibration with smooth, corner-free bulk data (known-limit control).

Zero-bump reference background (epsilon=0) plus a Gaussian scalar pulse
f=AMP*exp(-((z-Z0)/W)^2) placed deep in the bulk.  The warp a=b is obtained by
integrating the *continuum* Hamiltonian constraint (A_t=1, B_t=phi_t=0, a=b)
from z_R (where a=a_z=0) towards the cone with DOP853 (rtol 1e-12), using the
dense reference profile for the background.  Thus the shell data are exactly
the static reference (no corner incompatibility) and the discrete initial H is
pure O(h^4) truncation.  The outer C2 taper of the wrapper is applied to a,b,f.
Both the shell and the taper are causally disconnected from the diagnostic
window used for the order estimate (see analyze.py) for t<=1.
"""
import argparse, json, math, time, hashlib
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import lab

p = argparse.ArgumentParser()
p.add_argument('--hmin', type=float, required=True)
p.add_argument('--amp', type=float, default=.01); p.add_argument('--z0', type=float, default=-1.4)
p.add_argument('--w', type=float, default=.05); p.add_argument('--zr', type=float, default=-1.1)
p.add_argument('--kappa', type=float, default=0.); p.add_argument('--tf', type=float, default=1.)
p.add_argument('--out', required=True)
a = p.parse_args(); L_ = 3.; stretch = 2.
t0 = time.monotonic()
s, st, meta = lab.initial_data(0., a.hmin, L_, stretch)
assert not np.any(st)
sol, _ = lab.E.B.run(0., rtol=8e-14, transport=True)
zb = sol.y[4, -1]
def bg(z):
    tg = z + zb; i = int(np.clip(np.searchsorted(sol.y[4], tg) - 1, 0, len(sol.t) - 2))
    y = brentq(lambda yy: sol.sol(yy)[4] - tg, sol.t[i], sol.t[i + 1], xtol=1e-300, rtol=4 * np.finfo(float).eps)
    r, ry, psi, sy, _, _ = sol.sol(y)
    return r, ry, psi - 1, r * sy
U = lambda ph: .5 * (ph * ph - 1) ** 2 - (2 / 3) * (1 - ph + ph ** 3 / 3) ** 2
gauss = lambda z: a.amp * np.exp(-((z - a.z0) / a.w) ** 2)
dgauss = lambda z: -2 * (z - a.z0) / a.w ** 2 * gauss(z)
def rhs(z, y):
    aa, az = y; r, hc, ph, phz = bg(z); f = gauss(z); fz = dgauss(z)
    su = r * r * (math.exp(2 * aa) * U(ph + f) - U(ph))
    azz = (-2 * su - 12 * hc * az - 6 * az * az - 2 * phz * fz - fz * fz) / 6
    return [az, azz]
zl = s.z[0]
ode = solve_ivp(rhs, (a.zr, zl), [0., 0.], method='DOP853', rtol=1e-12, atol=1e-18, dense_output=True)
if not ode.success: raise RuntimeError(ode.message)
state = np.zeros((6, s.n)); left = s.z < a.zr
state[0, left] = ode.sol(s.z[left])[0]; state[1] = state[0]; state[2] = gauss(s.z)
x = np.clip((s.z + .99 * L_) / (.14 * L_), 0, 1); taper = x ** 3 * (10 - 15 * x + 6 * x * x)
state[:3] *= taper; state[:, :2] = 0
Lb = lab.Lab(s, kappa=a.kappa)
init, H0, M0 = Lb.diagnostics(state, 0.)
print(json.dumps({'nodes': s.n, 'ode_nfev': int(ode.nfev), 'initial': init}), flush=True)
res, snaps, final = lab.evolve(Lb, state, tf=a.tf)
here = Path(__file__).resolve().parent
doc = dict(status='FLOATING_POINT_ORDER_CALIBRATION_RUN', settings=vars(a), nodes=s.n, dt=res['dt'], steps=res['steps'],
           stop_reason=res['reason'], initial=init, snapshots={str(k): v['row'] for k, v in snaps.items()}, records=res['rows'],
           runtime_seconds=time.monotonic() - t0,
           sha256={f: hashlib.sha256((here / f).read_bytes()).hexdigest() for f in ['lab.py', 'calib.py']})
out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
out.with_suffix('.json').write_text(json.dumps(doc, indent=1) + '\n')
arrs = dict(z=s.z, background=np.array([Lb.rho, Lb.phi, Lb.hc, Lb.phz]), initial=state, H_initial=H0, final=final)
for k, v in snaps.items():
    arrs['state_%.2f' % k] = v['state']; arrs['H_%.2f' % k] = v['H']; arrs['M_%.2f' % k] = v['M']
np.savez_compressed(str(out) + '.npz', **arrs)
print(json.dumps({'done': str(out), 'snap': doc['snapshots']}), flush=True)
