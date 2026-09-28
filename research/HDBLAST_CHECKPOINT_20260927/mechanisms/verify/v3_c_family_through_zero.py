"""Audit V3: does the initial-shell family really 'end near c -> 0+' (workstream claim M7), or does the
continuation stop because its parameterisation phi_h = -1 + 10^x cannot represent phi_h < -1?

V1/D7 shows that at c = 0 the constant bulk phi = -1 satisfies both static junctions, so the family crosses the
trivial solution there and the linear response phi_h + 1 changes sign with c.  This script re-uses the SAME shell
ODE and junction residual as the workstream (copied module rolloff5d_v1.py, unmodified) but parameterises the cone
value as phi_h = -1 + s*10^x with s = +1 (control, reproduces M7) or s = -1 (new side, c < 0), at delta = 0.1.
It then (b) solves the '+1' static endpoint at delta = 0.1 for d = 0 and d = d*(0.1) with the workstream's own
generalised plus-branch solver (m8, imported from the audit's copy) to find phi_b of the tuned endpoint.
Output: V3_C_FAMILY_THROUGH_ZERO.json.  Time-boxed.
"""
import json, math, sys, time, warnings
import numpy as np
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE/'rerun')); sys.path.insert(0, str(HERE/'rerun/pilot_5d'))
import rolloff5d_v1 as R_
warnings.filterwarnings('ignore')
C = 2/1.0357712571566784 - 4/3

def solve_signed(ten, guess, sgn):
    def res(x):
        rho, phi, s = R_.integrate_y(-1 + sgn*10**x[0], x[1])
        rp = math.sqrt(1 + rho*rho*(s*s/12 - R_.U(phi)/6))
        return np.array([rp/rho - ten.s(phi)/6, s + ten.s1(phi)/2])
    x = np.array(guess, float)
    for _ in range(60):
        F = res(x); J = np.zeros((2, 2))
        for j in range(2):
            xp = x.copy(); xp[j] += 1e-7; xm = x.copy(); xm[j] -= 1e-7; J[:, j] = (res(xp) - res(xm))/2e-7
        dx = np.linalg.solve(J, -F); lam = 1.0
        while lam > 1e-4 and np.linalg.norm(res(x + lam*dx)) > np.linalg.norm(F): lam /= 2
        x = x + lam*dx
        if np.linalg.norm(lam*dx) < 1e-13: break
    return x, res(x)

def row(td, c, x, r, sgn):
    ph = -1 + sgn*10**x[0]; rho, phi, s = R_.integrate_y(ph, x[1])
    Hv2 = td*(1+c)/27 + td**2*((1+c)**2/36 - c*c/384)
    return dict(c=float(c), phi_h=float(ph), y_b=float(x[1]), phi_b=float(phi), rho_b=float(rho), H0sq=float(1/rho**2),
                H0sq_over_delta=float(1/rho**2/td), fH_series=float(Hv2*rho**2), max_junction_residual=float(np.max(np.abs(r))))

def main():
    t0 = time.time(); td = 0.1; out = dict(delta=td)
    # control: positive side at c = 0.0026572 (M7 last converged row: phi_h+1 = 1.565e-5, phi_b = -0.99486)
    x, r = solve_signed(R_.Tension(td, 0.002657237869143796, 0), (math.log10(1.5651473774580538e-05), 8.2282 + 0.9*math.log(1e-3/td)), +1)
    out['control_positive_side_reproduces_M7_last_row'] = row(td, 0.002657237869143796, x, r, +1)
    yb0 = x[1]
    # negative side: continue from c = -0.0027 down to c* = -0.99307 (series root at delta = 0.1)
    rows = []; c = -0.002657237869143796; x = np.array([math.log10(1.5651473774580538e-05), yb0])
    targets = list(-np.array([0.002657, 0.01, 0.03, 0.06, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99, 0.99307]))
    last_ok = None
    for ct in targets:
        if time.time() - t0 > 420: out['stopped'] = 'time box at c=%g' % ct; break
        # sub-step from last converged c to ct
        cc = last_ok['c'] if last_ok else ct; steps = [ct] if last_ok is None else list(np.linspace(cc, ct, 4)[1:])
        ok = True
        for cs in steps:
            try:
                xn, rn = solve_signed(R_.Tension(td, cs, 0), x, -1)
                if not np.all(np.isfinite(rn)) or np.max(np.abs(rn)) > 1e-9: raise RuntimeError('no root %g' % np.max(np.abs(rn)))
                x = xn; last_ok = row(td, cs, xn, rn, -1)
            except Exception as e:
                ok = False; out.setdefault('failures', []).append(dict(c=float(cs), error=str(e)[:120])); break
        if not ok: break
        rows.append(last_ok); print('c=%.5f phi_h=%.6e phi_b=%.5f H0^2/d=%.5f fH=%.4f res=%.1e [%.0fs]' % (
            last_ok['c'], last_ok['phi_h'], last_ok['phi_b'], last_ok['H0sq_over_delta'], last_ok['fH_series'], last_ok['max_junction_residual'], time.time()-t0), flush=True)
    out['negative_side_rows'] = rows
    out['reached_c_star'] = bool(rows and rows[-1]['c'] <= -0.99307 + 1e-9)
    # (b) static '+1' endpoint at delta = 0.1 with d = 0 and d = d*
    try:
        import m8_quadratic_tension_tuning as M8
        ends = {}
        g = None
        for dq in [0.0, -1.0, -2.0, -2.6, -3.0, M8.dstar(0.1)]:
            if time.time() - t0 > 900: ends['stopped'] = 'time box'; break
            try:
                r_, g = M8.plus_branch(0.1, dq, guess=g); ends['d=%.6f' % dq] = r_
                print('plus-branch delta=0.1 d=%.4f phi_b=%.5f H2=%.4e series=%.4e ok=%s' % (dq, r_['phi_b'], r_['H2'], r_['series'], r_['ok']), flush=True)
            except Exception as e:
                ends['d=%.6f' % dq] = dict(error=str(e)[:200])
        out['plus_branch_endpoint_delta_0p1'] = ends
        out['control_plus_branch_d0_H2_vs_M2'] = dict(M2_H2=0.0066143*1.0, note='M2 row delta=0.1: Hvac2/delta = 0.066143')
    except Exception as e:
        out['plus_branch_endpoint_delta_0p1'] = dict(error=str(e)[:300])
    out['runtime_s'] = time.time() - t0
    (HERE/'V3_C_FAMILY_THROUGH_ZERO.json').write_text(json.dumps(out, indent=1) + '\n')

if __name__ == '__main__':
    main()
