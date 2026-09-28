"""Leading-order fluctuation spectrum around the phi=+1 critical point with a dS4-sliced shell.

Background at leading order: pure AdS5 (k=1/9) in dS slicing, R = k rho = sinh u, 0<u<=u_b, H = k/sinh(u_b).
Mode functions f(u) Y(x) with (Box_dS - 2) Y = M^2 Y (tensor) or Box_dS Y = M^2 Y (scalar); M^2 = m^2/H^2.
   f'' + 4 coth(u) f' + (M^2/sinh^2 u - m5^2/k^2) f = 0
   tensor: m5 = 0,  Neumann f'(u_b)=0;   bulk scalar at phi=1: m5^2/k^2 = 252, Robin f'(u_b) = -18 f(u_b)
   (from the linearised scalar junction  f_y = -(sigma''/2) f,  sigma''(1)=4; mixing with the metric neglected:
    it enters through the background gradient phi'_b = O(delta) and is NOT included here).
Exact reduction: f = sinh(u)^(-3/2) w(cosh u), w = associated Legendre with nu(nu+1)=m5^2/k^2+15/4, mu^2 = 9/4 - M^2.
Parts:  A exact symbolic reduction;  B mpmath identity checks;  C tensor spectrum (exact condition);
        D decoupled bulk-scalar bound-state scan;  E calibration/controls with direct ODE shooting;
        F tensor sector on the full nonlinear branch (float, scipy) at two detunings.
Writes SPECTRUM_LEADING_ORDER.json
"""
import json, hashlib, time, math
from pathlib import Path
import sympy as sp
import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
out = dict(status='mixed: exact (A,C) / numerical (B,D,E,F); leading order in delta, scalar-metric mixing neglected')
t0 = time.time()

# ---------------- A. exact reduction to the associated Legendre equation
u, M2, m2 = sp.symbols('u M2 m2')
W0, W1, W2 = sp.symbols('W0 W1 W2')   # w, dw/dx, d2w/dx2 evaluated at x = cosh u
sh, ch = sp.sinh(u), sp.cosh(u)
# f = sh^(-3/2) w(cosh u);  d/du w = sh W1,  d2/du2 w = sh^2 W2 + ch W1
g = sh**sp.Rational(-3, 2)
gp = sp.diff(g, u); gpp = sp.diff(g, u, 2)
f0 = g * W0
f1 = gp * W0 + g * sh * W1
f2 = gpp * W0 + 2 * gp * sh * W1 + g * (sh**2 * W2 + ch * W1)
ode = f2 + 4 * ch / sh * f1 + (M2 / sh**2 - m2) * f0
leg = (1 - ch**2) * W2 - 2 * ch * W1 + (m2 + sp.Rational(15, 4) - (sp.Rational(9, 4) - M2) / (1 - ch**2)) * W0
ratio = sp.simplify((ode / g + leg).rewrite(sp.exp))
out['A_reduction_to_Legendre'] = dict(
    statement='f = sinh^{-3/2}(u) w(cosh u) turns the mode equation into (1-x^2)w\'\'-2xw\'+[nu(nu+1)-mu^2/(1-x^2)]w=0 with nu(nu+1)=m5^2/k^2+15/4, mu^2=9/4-M^2',
    residual=str(ratio), passed=bool(ratio == 0),
    nu_tensor='3/2', nu_scalar='31/2 (nu+1/2 = sqrt(252+4) = 16)')
assert out['A_reduction_to_Legendre']['passed']

# ---------------- B. mpmath identity checks (Legendre of type 3, x>1)
mp.mp.dps = 30
def P(nu_, mu_, x_):
    return mp.legenp(nu_, -mu_, x_, type=3)
def logderiv_exact(nu_, mu_, ub):
    x_ = mp.cosh(ub); s_ = mp.sinh(ub)
    # (x^2-1) dP^{-mu}_nu/dx = nu x P^{-mu}_nu - (nu-mu) P^{-mu}_{nu-1}
    return -mp.mpf(3) / 2 * x_ / s_ + (nu_ * x_ * P(nu_, mu_, x_) - (nu_ - mu_) * P(nu_ - 1, mu_, x_)) / (s_ * P(nu_, mu_, x_))
idchecks = []
for nu_ in [mp.mpf(3) / 2, mp.mpf(31) / 2]:
    for mu_ in [mp.mpf('0.3'), mp.mpf('1.1'), mp.mpf('2.7')]:
        for ub in [mp.mpf('0.4'), mp.mpf('2.0')]:
            fnum = lambda uu: mp.sinh(uu)**(-mp.mpf(3) / 2) * P(nu_, mu_, mp.cosh(uu))
            ld_num = mp.diff(fnum, ub) / fnum(ub)
            ld_ex = logderiv_exact(nu_, mu_, ub)
            # ODE residual of f
            M2v = mp.mpf(9) / 4 - mu_**2; m2v = nu_ * (nu_ + 1) - mp.mpf(15) / 4
            res = mp.diff(fnum, ub, 2) + 4 * mp.coth(ub) * mp.diff(fnum, ub) + (M2v / mp.sinh(ub)**2 - m2v) * fnum(ub)
            idchecks.append(dict(nu=float(nu_), mu=float(mu_), u=float(ub), logderiv_identity_err=float(abs(ld_num - ld_ex)),
                                 ode_rel_residual=float(abs(res) / (abs(fnum(ub)) * (1 + abs(m2v))))))
out['B_identity_checks'] = dict(max_logderiv_err=max(c['logderiv_identity_err'] for c in idchecks),
                                max_ode_rel_residual=max(c['ode_rel_residual'] for c in idchecks), rows=idchecks,
                                dps=30, passed=max(c['logderiv_identity_err'] for c in idchecks) < 1e-15)
assert out['B_identity_checks']['passed']

# ---------------- C. tensor: quantization condition (mu - 3/2) P^{-mu}_{1/2}(cosh u_b) = 0
# positivity of P^{-mu}_{1/2}(x>1) for real mu>0 (proof in README via Pfaff transformation); numerical sweep:
minP = None; rows = []
for ub in [0.01, 0.1, 0.5, 1, 2, 3.364, 5, 8]:
    for mu_ in np.linspace(0.01, 6, 60):
        v = P(mp.mpf(1) / 2, mp.mpf(mu_), mp.cosh(ub)) * mp.gamma(1 + mu_) / mp.tanh(mp.mpf(ub) / 2)**mu_
        minP = v if minP is None else min(minP, v)
    # full tensor condition F(M^2) for M^2 in (0, 9/4): count sign changes of f'(u_b)
    Fvals = [logderiv_exact(mp.mpf(3) / 2, mp.sqrt(mp.mpf(9) / 4 - m2v_), mp.mpf(ub)) for m2v_ in np.linspace(0.02, 2.24, 60)]
    rows.append(dict(u_b=ub, H_over_k=float(1 / mp.sinh(ub)), sign_changes_in_0_to_9over4=int(sum(1 for a, b in zip(Fvals, Fvals[1:]) if a * b < 0)),
                     max_F=float(max(Fvals)), min_F=float(min(Fvals))))
out['C_tensor'] = dict(
    exact_condition='f\'(u_b)=0  <=>  (mu - 3/2) P^{-mu}_{1/2}(cosh u_b) = 0 ;  P^{-mu}_{1/2}(x)>0 for x>1, mu>0',
    result='bound states: only mu=3/2 (M^2=0, f=const, the massless 4D graviton) for EVERY u_b; continuum M^2 >= 9/4, i.e. m >= (3/2) H',
    min_normalised_P_minus_mu_half=float(minP), scan=rows,
    literature='Garriga & Sasaki, https://arxiv.org/abs/hep-th/9912118 (gap m=(3/2)H); this project re-derives it in closed form')
assert minP > 0 and all(r['sign_changes_in_0_to_9over4'] == 0 for r in rows)

# ---------------- D. decoupled bulk scalar at phi=1: bound states?
nuS = mp.mpf(31) / 2; b = 18
def FS(M2v_, ub):
    mu_ = mp.sqrt(mp.mpf(9) / 4 - M2v_)
    return logderiv_exact(nuS, mu_, ub) + b
scan = []
ubs = [0.02, 0.03, 0.05, 0.07, 0.08, 0.09, 0.1, 0.12, 0.15, 0.2, 0.3, 0.5, 1, 2, 3.364, 5]
for ub in ubs:
    M2grid = np.concatenate([np.linspace(0.001, 2.249, 90)])
    vals = [FS(mp.mpf(m), mp.mpf(ub)) for m in M2grid]
    roots = []
    for (ma, va), (mb, vb) in zip(zip(M2grid, vals), zip(M2grid[1:], vals[1:])):
        if va * vb < 0:
            roots.append(float(mp.findroot(lambda m: FS(m, mp.mpf(ub)), (mp.mpf(ma), mp.mpf(mb)), solver='anderson')))
    scan.append(dict(u_b=ub, H_over_k=float(1 / mp.sinh(ub)), bound_state_M2=roots, min_F_plus_b=float(min(vals)), max_F_plus_b=float(max(vals)),
                     F_plus_b_at_M2_0=float(FS(mp.mpf('1e-12'), mp.mpf(ub)))))
# critical u_b where a bound state first appears at threshold M^2 -> 9/4 and where it reaches M^2 -> 0
def crit(M2v_):
    return float(mp.findroot(lambda ub: FS(mp.mpf(M2v_), ub), (mp.mpf('0.01'), mp.mpf('0.3')), solver='anderson'))
try:
    ub_thr = crit('2.2499'); ub_zero = None
except Exception as e:
    ub_thr = repr(e)
flat = [dict(u_b=ub, F_at_M2=[float(logderiv_exact(nuS, mp.sqrt(mp.mpf(9) / 4 - m), mp.mpf(ub))) for m in [0, 1, 2]]) for ub in [3, 5, 8]]
out['D_bulk_scalar_decoupled'] = dict(
    exact_condition='(14 cosh u_b + 18 sinh u_b) P^{-mu}_{31/2}(cosh u_b) = (31/2 - mu) P^{-mu}_{29/2}(cosh u_b)',
    no_zero_or_tachyonic_mode='exact: for M^2<=0 the regular f is positive and increasing (maximum principle), so f\'/f>0 > -18',
    scan=scan, u_b_where_bound_state_reaches_threshold=ub_thr,
    large_u_b_log_derivative=flat,
    interpretation='for u_b above the threshold value (H/k below ~ 1/sinh(u_thr)) the decoupled bulk scalar has NO discrete mode below the '
                   'universal continuum M^2>=9/4; f\'/f -> 14 = Delta-4 while the brane demands -18 = -Delta (sigma\'\'/(2k) = 3W\'\'/W = Delta).')

# ---------------- E. calibration: direct ODE shooting in pure AdS vs exact Legendre (float)
def shoot_logderiv(M2v_, ub, m5, bg=None):
    mu_ = math.sqrt(9 / 4 - M2v_); u0 = 1e-4
    # background: R=sinh u unless bg (callable returning R,R') given
    def rhs(uu, s):
        R, Rp, e, ep, fz, fp = s
        if bg is None:
            Rr, Rpp = math.sinh(uu), math.cosh(uu); Rdd = Rr; edd = 0.0
        else:
            Rr, Rpp = R, Rp
            e2 = e * e
            V = 1 - 21 * e2 - 25 * e2 * e + 2.25 * e2 * e2 + 6 * e2 * e2 * e + e2 ** 3
            Fe = 252 * e + 450 * e2 - 54 * e2 * e - 180 * e2 * e2 - 36 * e2 * e2 * e
            Rdd = R * (V - ep * ep / 4); edd = -4 * Rp / R * ep + Fe
        return [Rpp if bg is None else Rp, Rdd, ep, edd, fp, -4 * Rpp / Rr * fp - (M2v_ / Rr ** 2 - m5) * fz]
    eta_h = bg or 0.0
    s0 = [math.sinh(u0), math.cosh(u0), eta_h, 0.0, u0 ** (mu_ - 1.5), (mu_ - 1.5) * u0 ** (mu_ - 2.5)]
    sol = solve_ivp(rhs, (u0, ub), s0, method='DOP853', rtol=1e-12, atol=1e-30)
    R, Rp, e, ep, fz, fp = sol.y[:, -1]
    return fp / fz, (R, Rp, e, ep)
calib = []
for ub in [0.5, 2.0, 3.364]:
    for M2v_ in [0.3, 1.2, 2.0]:
        ex = float(logderiv_exact(mp.mpf(3) / 2, mp.sqrt(mp.mpf(9) / 4 - mp.mpf(M2v_)), mp.mpf(ub)))
        nm, _ = shoot_logderiv(M2v_, ub, 0.0)
        exs = float(logderiv_exact(nuS, mp.sqrt(mp.mpf(9) / 4 - mp.mpf(M2v_)), mp.mpf(ub)))
        nms, _ = shoot_logderiv(M2v_, ub, 252.0)
        # wrong-formula control: Legendre degree 1/2 instead of 3/2 for the tensor
        wrong = float(logderiv_exact(mp.mpf(1) / 2, mp.sqrt(mp.mpf(9) / 4 - mp.mpf(M2v_)), mp.mpf(ub)))
        calib.append(dict(u_b=ub, M2=M2v_, tensor_exact=ex, tensor_shoot=nm, scalar_exact=exs, scalar_shoot=nms,
                          tensor_abs_err=abs(ex - nm), scalar_rel_err=abs(exs - nms) / abs(exs), wrong_degree_control_err=abs(wrong - nm)))
out['E_calibration'] = dict(rows=calib, max_tensor_abs_err=max(r['tensor_abs_err'] for r in calib),
                            max_scalar_rel_err=max(r['scalar_rel_err'] for r in calib),
                            min_wrong_control_err=min(r['wrong_degree_control_err'] for r in calib),
                            note='float shooting from u0=1e-4 with leading Frobenius behaviour u^(mu-3/2); errors O(u0^2)')

# ---------------- F. tensor sector on the full nonlinear branch (float)
bvp = json.loads((HERE / 'runs' / 'BVP_SCAN_reg.json').read_text())
full = []
for row in bvp['rows']:
    if row['delta'] not in ('0.001', '0.1', '0.01'):
        continue
    eta_h = float(mp.mpf(row['eta_h'])); ub = float(mp.mpf(row['u_b']))
    Fv = []; Mg = np.linspace(0.02, 2.24, 45)
    for m in Mg:
        ld, bgend = shoot_logderiv(m, ub, 0.0, bg=eta_h if eta_h != 0 else 1e-300)
        Fv.append(ld)
    zero_mode_residual = 'h=const solves h\'\'+4(R\'/R)h\'+M^2 h/R^2=0 at M^2=0 with h\'(u_b)=0 for ANY background (exact)'
    full.append(dict(delta=row['delta'], u_b=ub, eta_h=eta_h, H2_from_mp=row['H2'], R_b_float=bgend[0], R_b_mp=row['R_b'],
                     sign_changes=int(sum(1 for a, b in zip(Fv, Fv[1:]) if a * b < 0)), min_F=float(min(Fv)), max_F=float(max(Fv)),
                     F_sign='negative' if max(Fv) < 0 else 'mixed' if min(Fv) < 0 else 'positive'))
out['F_tensor_full_branch'] = dict(rows=full, zero_mode='exact for any background: h=const (M^2=0), normalisable because the doubled bulk has finite warped volume',
                                   note='float shooting; no tensor bound state in 0<M^2<9/4 detected at the listed detunings')
out['runtime_s'] = time.time() - t0
out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(HERE / 'SPECTRUM_LEADING_ORDER.json').write_text(json.dumps(out, indent=2, default=str) + '\n')
print(json.dumps({k: (v if k not in ('B_identity_checks', 'E_calibration') else {kk: vv for kk, vv in v.items() if kk != 'rows'}) for k, v in out.items()}, indent=1, default=str)[:6000])
