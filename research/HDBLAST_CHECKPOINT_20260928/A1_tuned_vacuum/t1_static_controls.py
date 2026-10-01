#!/usr/bin/env python3
"""A1 controls for the static shell in the invariant w and for the chart-mapped reference solution.
Output: T1_STATIC_CONTROLS.json.  Numerical (floating point); tolerances recorded per check.

 1. Shell root vs the archived pilot solver (rolloff5d_v1.solve_shell, read-only import): phi_h, rho_b, phi_b for
    delta = 0.1 (d = 0 and d = d*), delta = 1e-3 (d = 0 and d = d*), seed tension c + dc.
 2. Exterior profile vs rolloff5d_v1.static_on_grid on a z grid.
 3. First-order constraint along the tables (both sides of w = 0).
 4. Series/table continuity at |w| = w_start.
 5. The chart-mapped reference satisfies the conformal-gauge PDEs and constraints: residuals from its own analytic
    first/second derivatives (independent of the finite-difference solver) at random points, including points beyond
    the vertex light cone (w < 0), for xc = 3.3 and xc = infinity.  Control: perturbing l'' by 1% must give O(1e-2) residuals.
 6. H0 normalisation: static H = 1/rho_b and the H^2 shell identity with Weyl = 0 at the static shell.
"""
import json, math, sys, hashlib, time
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
PIL = '/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/mechanisms/pilot_5d'
sys.path.insert(0, PIL)
import static_w as S
import rolloff5d_v1 as R_          # read-only import of the archived pilot base code

M8 = json.load(open('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/mechanisms/M8_QUADRATIC_TENSION_TUNING.json'))
DSTAR = {0.1: M8['d_star']['0.1'], 0.001: M8['d_star']['0.001'], 0.01: M8['d_star']['0.01']}
out = dict(status='numerical', checks=[], d_star_used=DSTAR)
def rec(name, val, tol, extra=None):
    ok = bool(abs(val) <= tol)
    out['checks'].append(dict(name=name, value=float(val), tol=tol, result='PASS' if ok else 'FAIL', **(extra or {})))
    print('%s %-80s %.3e (tol %.1e)' % ('PASS' if ok else 'FAIL', name, val, tol))
    return ok

shells = {}
for delta in [0.1, 0.001]:
    for dq in [0.0, DSTAR[delta]]:
        for dc in [0.0, 1e-2]:
            if delta == 0.001 and dc == 1e-2: dc = 1e-4
            c = S.C_REG + dc
            g0 = (math.log10(0.2207*delta**1.8), 8.2282 + 0.9*math.log(1e-3/delta))
            t0 = time.time()
            ph_h, y_b, r0 = R_.solve_shell(R_.Tension(delta, c, dq), guess=g0)
            Sref = R_.static_on_grid(ph_h, y_b, np.array([-0.3, -0.1, -0.03, -0.01, 0.0]), hz=2.5e-5)
            sh = S.StaticShell(S.Tension(delta, c, dq), phh_guess=ph_h)
            key = 'delta=%g d=%.6f dc=%g' % (delta, dq, dc)
            shells[key] = sh
            rec(key + ': phi_h+1 relative to pilot solver', (sh.phh + 1)/(ph_h + 1) - 1, 1e-6)
            rec(key + ': rho_b relative to pilot solver', sh.rhob/Sref['rho'][-1] - 1, 1e-8)
            tolp = 1e-8 if delta > 0.01 else 1e-7   # the pilot solver (RK4 in y, dv = 2e-4) is itself accurate to ~1e-8 only at delta = 1e-3
            rec(key + ': phi_b vs pilot solver', sh.phb - Sref['phi'][-1], tolp)
            rec(key + ': A-junction residual', sh.junction_residuals[0], 1e-11)
            rec(key + ': scalar junction residual', sh.junction_residuals[1], 1e-9)
            rec(key + ': first-order constraint at shell (relative to 1)', sh.constraint_at_shell, 1e-10)
            rho, h, ph, phz = sh.profile_z(np.array([-0.3, -0.1, -0.03, -0.01, 0.0]))
            rec(key + ': exterior profile rho vs static_on_grid (max rel)', np.max(np.abs(rho/Sref['rho'] - 1)), 1e-8)
            rec(key + ': exterior profile phi vs static_on_grid (max abs)', np.max(np.abs(ph - Sref['phi'])), tolp)
            rec(key + ': exterior A_z = rho_z/rho vs static_on_grid Hc (max rel)', np.max(np.abs(h/Sref['Hc'] - 1)), tolp)
            for sgn in (1, -1):
                tb = sh.tab[sgn]; w = sgn*np.exp(tb['x']); l, lam, p, psi = tb['Y']
                cons = 2*lam + lam**2 - (psi**2/3 - w*np.exp(l)*S.U(p)/6)
                scl = 2*np.abs(lam) + lam**2 + psi**2/3 + np.abs(w*np.exp(l)*S.U(p))/6 + 1e-300
                big = np.abs(w) > 1e-4
                rec(key + ': first-order constraint max relative residual (|w|>1e-4) side %+d' % sgn, np.max(np.abs(cons[big])/scl[big]), 1e-10)
                rec(key + ': first-order constraint max absolute residual (|w|<=1e-4) side %+d' % sgn, np.max(np.abs(cons[~big])), 1e-12)
            ws = np.array([sh.w_start*0.999999, sh.w_start*1.000001, -sh.w_start*0.999999, -sh.w_start*1.000001])
            F = sh.fields_w(ws)
            rec(key + ': series/table continuity of lhat at |w|=w_start', max(abs(F[1][0] - F[1][1]), abs(F[1][2] - F[1][3])), 1e-9)
            rec(key + ': series/table continuity of lhat_w at |w|=w_start', max(abs(F[2][0] - F[2][1]), abs(F[2][2] - F[2][3]))/max(1e-30, abs(F[2][0])), 1e-3)
            if dc == 0.0:
                sh2 = S.StaticShell(S.Tension(delta, c, dq), phh_guess=ph_h, rtol=1e-12, atol=1e-15)
                rec(key + ': self-convergence phi_b (rtol 1e-12 vs 1e-13)', sh2.phb - sh.phb, 1e-8)
                rec(key + ': self-convergence rho_b relative (rtol 1e-12 vs 1e-13)', sh2.rhob/sh.rhob - 1, 1e-8)
            out.setdefault('shells', {})[key] = dict(phi_h=sh.phh, w_b=sh.wb, c0=sh.c0, rho_b=sh.rhob, phi_b=sh.phb, H0=1/sh.rhob,
                                                     pilot_phi_h=ph_h, pilot_y_b=y_b, seconds=round(time.time() - t0, 1))

# 5. PDE residuals of the mapped reference (analytic derivatives only)
def pde_res(ref):
    e2B = np.exp(2*ref['B']); Uv = S.U(ref['phi']); U1v = S.U1(ref['phi'])
    rA = ref['A_TT'] - (ref['A_ZZ'] - 3*ref['A_T']**2 + 3*ref['A_Z']**2 + (2/3)*e2B*Uv)
    rB = ref['B_TT'] - (ref['B_ZZ'] + 3*ref['A_T']**2 - 3*ref['A_Z']**2 - ref['phi_T']**2/2 + ref['phi_Z']**2/2 - e2B*Uv/3)
    rF = ref['phi_TT'] - (ref['phi_ZZ'] - 3*ref['A_T']*ref['phi_T'] + 3*ref['A_Z']*ref['phi_Z'] - e2B*U1v)
    M = -3*ref['A_TZ'] - 3*ref['A_T']*ref['A_Z'] + 3*ref['A_T']*ref['B_Z'] + 3*ref['A_Z']*ref['B_T'] - ref['phi_T']*ref['phi_Z']
    H = -2*Uv*e2B + 6*ref['A_T']**2 + 6*ref['A_T']*ref['B_T'] - 12*ref['A_Z']**2 + 6*ref['A_Z']*ref['B_Z'] - 6*ref['A_ZZ'] - ref['phi_T']**2 - ref['phi_Z']**2
    scale = 1 + np.abs(ref['A_ZZ']) + 3*ref['A_T']**2 + 3*ref['A_Z']**2 + np.abs(e2B*Uv) + np.abs(ref['phi_ZZ']) + ref['phi_Z']**2
    return {k: float(np.max(np.abs(v)/scale)) for k, v in dict(A=rA, B=rB, phi=rF, M=M, H=H).items()}
rng = np.random.default_rng(1)
for key in ['delta=0.1 d=-3.106933 dc=0.01', 'delta=0.001 d=-3.194242 dc=0.0001', 'delta=0.1 d=0.000000 dc=0']:
    sh = shells[key]
    for xc, kind in [(3.3, 'bounded'), (6.0, 'bounded'), (3.3, 'asinh'), (3.3, 'softplus'), (math.inf, 'asinh')]:
        ch = S.Chart(xc, kind)
        Tm = min(ch.F_inf - 0.3, 5.0) if math.isfinite(ch.F_inf) else 5.0
        T = rng.uniform(0, Tm + (3.0 if ch.kind in ('asinh', 'bounded') else 0.0), 4000); Z = -rng.uniform(0, 12.0, 4000)
        ok = (T + Z) < (ch.F_inf - 0.05 if math.isfinite(ch.F_inf) else 1e9)
        T, Z = T[ok], Z[ok]
        # keep points inside the tables
        ref = S.reference(sh, ch, T, Z)
        inside = np.abs(ref['w']) < 24*sh.wb
        ref = S.reference(sh, ch, T[inside], Z[inside])
        r = pde_res(ref)
        nneg = int(np.sum(ref['w'] < 0))
        for k, v in r.items():
            rec('%s xc=%s %s: reference PDE/constraint residual %s (n=%d, %d beyond light cone)' % (key, xc, ch.kind, k, len(ref['w']), nneg), v, 5e-9)
# control: corrupt l'' by 1 %
sh = shells['delta=0.1 d=-3.106933 dc=0.01']; ch = S.Chart(3.3)
orig = sh.fields_w
def bad(w):
    o = orig(w); o[2] = o[2]*1.01; return o
sh.fields_w = bad
T = rng.uniform(0, 3, 2000); Z = -rng.uniform(0, 1, 2000)
r = pde_res(S.reference(sh, ch, T, Z))
sh.fields_w = orig
ctrl_ok = r['A'] > 1e-5
out['checks'].append(dict(name="CONTROL: reference with l'' x 1.01 must give A-equation residual > 1e-5", value=r['A'], tol=1e-5, result='PASS' if ctrl_ok else 'FAIL'))
print('PASS' if ctrl_ok else 'FAIL', 'control corrupt lpp', r['A'])
out['control_corrupt_lpp_residuals'] = r

# 6. static shell H identity: H = 1/rho_b, Weyl = H^2 - sigma^2/36 + sigma'^2/48 - U/6 (v = 0, R = 0)
for key, sh in shells.items():
    ten = sh.ten; H = 1/sh.rhob; p = sh.phb
    Wy = H*H - ten.s(p)**2/36 + ten.s1(p)**2/48 - S.U(p)/6
    rec(key + ': static Weyl scalar from the shell identity (units H0^2)', Wy/H**2, 1e-9)

out['n_pass'] = sum(c['result'] == 'PASS' for c in out['checks']); out['n_checks'] = len(out['checks'])
out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(HERE/'T1_STATIC_CONTROLS.json').write_text(json.dumps(out, indent=1) + '\n')
print(out['n_pass'], '/', out['n_checks'])
