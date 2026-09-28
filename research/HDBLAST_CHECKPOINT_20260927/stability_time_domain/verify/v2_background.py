"""Audit V2: independent re-solution of the static +1 branch (and access to the registered original
shell) for use by the independent spectral check.  Different algorithm from td_core.plus_background:
first-order system (rho, eta, eta_y, z) with rho_y from the Hamiltonian first integral, the metric
junction located as an integration EVENT, and a 1-D bracketing root (brentq) in log10(-eta_h) for the
scalar junction.  Writes verify/v2_background.json when run as a script."""
import json, math, sys, hashlib
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

C = 2/1.0357712571566784 - 4/3
UETA = np.array([-2/27, 0., 14/9, 50/27, -1/6, -4/9, -2/27])     # U(1+eta), verified in v1


def Ud(eta, k=0):
    c = UETA.copy()
    for _ in range(k): c = np.arange(1, len(c))*c[1:]
    return np.polynomial.polynomial.polyval(eta, c)


def sig(eta, d, c=C):   # sigma(phi=1+eta)
    p = 1 + eta; return 2*(1 - p + p**3/3) + d*(1 + c*p)
def sig1(eta, d, c=C):
    p = 1 + eta; return 2*(p*p - 1) + d*c


def shoot(lx, d, y0=2e-5, rtol=1e-13, dense=False):
    eh = -10**lx; u0, u1, u2 = Ud(eh), Ud(eh, 1), Ud(eh, 2)
    A_ = -u0/36; B_ = u0*u0/4320 - u1*u1/750; Cc = u1/10; D_ = u1*(u2/280 + u0/630)
    ini = [y0 + A_*y0**3 + B_*y0**5, eh + Cc*y0**2 + D_*y0**4, 2*Cc*y0 + 4*D_*y0**3, 0.]
    def rhs(y, v):
        r, e, p, _ = v
        rp = math.sqrt(max(1 + r*r*(p*p/12 - Ud(e)/6), 0.))
        return [rp, p, Ud(e, 1) - 4*rp*p/r, 1/r]
    def ev(y, v):
        r, e, p, _ = v; rp = math.sqrt(max(1 + r*r*(p*p/12 - Ud(e)/6), 0.))
        return rp/r - sig(e, d)/6
    ev.terminal = True; ev.direction = -1
    sol = solve_ivp(rhs, (y0, 200.), ini, method='LSODA', rtol=rtol, atol=[1e-15, 1e-40, 1e-40, 1e-15],
                    events=ev, dense_output=dense, max_step=0.1)
    if sol.status != 1: raise RuntimeError('no metric-junction event for lx=%g' % lx)
    r, e, p, zz = sol.y_events[0][0]
    return sol, dict(y_b=float(sol.t_events[0][0]), rho_b=float(r), eta_b=float(e), p_b=float(p), z_end=float(zz),
                     scalar_residual=float(p + sig1(e, d)/2))


def solve_plus(d):
    f = lambda lx: shoot(lx, d)[1]['scalar_residual']
    # bracket around the archived eta_h scale; scan if needed
    grid = np.linspace(-40, -3, 75)
    vals = []
    for lx in grid:
        try: vals.append(f(lx))
        except Exception: vals.append(np.nan)
    vals = np.array(vals); roots = []
    for i in range(len(grid) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and vals[i]*vals[i + 1] < 0:
            roots.append(brentq(f, grid[i], grid[i + 1], xtol=1e-14, rtol=1e-15))
    return roots, dict(scan=grid.tolist(), residual=vals.tolist())


def sampler(d, lx):
    sol, info = shoot(lx, d, dense=True)
    yb = info['y_b']; zend = info['z_end']
    yy = np.linspace(sol.t[0], yb, 400001); vv = sol.sol(yy); zz = vv[3] - zend
    def sample(z):
        z = np.asarray(z, float)
        ys = np.interp(z, zz, yy)
        for _ in range(8):
            v = sol.sol(ys); ys = ys - (v[3] - zend - z)*v[0]
        v = sol.sol(ys); r, e, p = v[0], v[1], v[2]
        hc = np.sqrt(1 + r*r*(p*p/12 - Ud(e)/6))
        return dict(rho=r, eta=e, phi=1 + e, hc=hc, phiz=r*p, zmin=float(zz[0]))
    return sample, info


if __name__ == '__main__':
    out = {}
    arch = json.loads((Path(__file__).resolve().parents[1]/'frozen_input'/'PLUS_BRANCH_RESULTS.json').read_text())
    for d in [1e-3, 3e-3, 1e-2]:
        roots, scan = solve_plus(d)
        rows = []
        for lx in roots:
            _, info = shoot(lx, d)
            rows.append(dict(eta_h=-10**lx, **info, H2=1/info['rho_b']**2, phi_b=1 + info['eta_b']))
        a = min(arch['rows'], key=lambda q: abs(q['delta'] - d))['refined']
        best = min(rows, key=lambda r: abs(r['rho_b'] - a['rho_b'])) if rows else None
        out[str(d)] = dict(n_roots_in_scan=len(rows), roots=rows,
                           archived=dict(rho_b=a['rho_b'], phi_b=a['phi_b'], eta_h=a['eta_h'], H2=a['H2']),
                           rel_diff_rho_b=(best['rho_b']/a['rho_b'] - 1) if best else None,
                           diff_phi_b=(best['phi_b'] - a['phi_b']) if best else None)
    # H/H0 check with H0^2 = the original registered shell's H^2 (Chat 9 value 1.6093000853360052e-4)
    H0sq = 1.6093000853360052e-4
    out['H_over_H0_at_1e-3'] = math.sqrt(out['0.001']['roots'][0]['H2']/H0sq) if out['0.001']['roots'] else None
    out['checkpoint_quoted'] = dict(phi_b=0.9999159473169134, rho_b=129.9247628497, H2=5.92401479433e-5, H_over_H0=0.606721732)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps({k: (v if not isinstance(v, dict) or 'roots' not in v else dict(n=v['n_roots_in_scan'], rel_rho=v['rel_diff_rho_b'], dphi=v['diff_phi_b'], roots=[(r['eta_h'], r['rho_b'], r['phi_b']) for r in v['roots']])) for k, v in out.items()}, indent=1))
