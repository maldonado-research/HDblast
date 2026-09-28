"""Cheap deliberate controls for the gauge-invariant operator.

C1 known limit (decoupled AdS scalar): on the +1 branch at delta=3e-4 the backreaction is O(1e-8), so the
   mu2=0 regular solution X must be proportional to the Gegenbauer polynomial C_14^(2)(cosh(y/9))
   (massive scalar m^2 = U''(1) = 28/9 in AdS5 of radius 9).  Compare log-derivatives at the shell.
C2 detection power (scan down to mu2 = -1e6): flip the sign of the shell curvature term (sigma''/2 -> -sigma''/2) in the shell
   condition on the +1 background.  If this wrong operator produces tachyons, the scan is demonstrably
   able to see unstable modes on this background.
C3 wrong formula on the calibration: drop the gravitational term 3 lam Z/rho^2 in M on the ORIGINAL
   shell; the calibrated eigenvalue must move away from -7.7179.
C4 small-detuning limit: the original shell's roots at delta = 1e-3 and 3e-3 extrapolated linearly to
   delta -> 0 must approach the closed-form 4D value -4(3c^2-4c+8)/(c(3c+4)) = -7.719796 (Chat 9).
C5 perturbed parameter: c -> c(1 +- 1e-3) on the original shell must shift mu2 smoothly (finite derivative),
   and on the +1 branch (background re-solved) must not create a root.
Output: CONTROLS_RESULTS.json
"""
import sys, json, time, math, hashlib
sys.dont_write_bytecode = True
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.special import eval_gegenbauer
from scipy.integrate import solve_ivp
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gi_core import Background, C_REG, growth_exponent

t0 = time.time()
B = json.loads((HERE/'BACKGROUNDS.json').read_text())
P = {r['delta']: r['gi_core'] for r in B['plus']}; O = {r['delta']: r['gi_core'] for r in B['original']}
out = {}

def M_variant(bg, mu2, variant):
    f = bg.scalar_M(mu2, full=True); sv = bg._sv
    X, Z, Bv = f['X_b'], f['Z_b'], f['B']; lam = mu2 + 4
    g = sv['g_b']; s2 = sv['sigma2_b']
    if variant == 'flip_sigma2':
        st = bg.scalar_XZ_renorm(mu2); X, Z = st[3], st[4]
        M = (g - s2/2)*X + 3*lam*Z/sv['rho_b']**2; n = abs((g - s2/2)*X) + abs(3*lam*Z/sv['rho_b']**2)
    elif variant == 'no_grav':
        M = Bv*X; n = abs(Bv*X)
    return M/n

def scan_roots(fun, grid):
    vals = [fun(m) for m in grid]; roots = []
    for i in range(len(grid) - 1):
        if vals[i]*vals[i + 1] < 0:
            roots.append(brentq(fun, grid[i], grid[i + 1], xtol=1e-12))
    return roots, vals

# C1 Gegenbauer known limit.  The deviation from the decoupled linear AdS scalar is expected at O(eta_b) = O(delta)
# (nonlinear U''(eta) = 28/9 + (100/9) eta + ... and backreaction), so the test is that
# (relative difference)/delta is ~constant, i.e. the difference extrapolates to zero as delta -> 0.
c1 = {}
for d in (0.0003, 0.001, 0.003):
    bg = Background(+1, d, e_h=P[d]['e_h'], y_b=P[d]['y_b'], polish=False); sv = bg.shell_values()
    f = bg.scalar_M(0.0, full=True)
    Xp = sv['g_b']*f['X_b'] + (12/sv['rho_b']**2 - 2*sv['phi_y_b']**2)*f['Z_b']
    x = math.cosh(bg.y_b/9)
    geg = (1/9)*math.sinh(bg.y_b/9)*4*eval_gegenbauer(13, 3, x)/eval_gegenbauer(14, 2, x)
    rel = Xp/f['X_b']/geg - 1
    c1[str(d)] = dict(Xprime_over_X=Xp/f['X_b'], gegenbauer_logderivative=geg, relative_difference=rel,
                      relative_difference_over_delta=rel/d, eta_b=sv['eta_b'])
ratios = [c1[k]['relative_difference_over_delta'] for k in c1]
out['C1_gegenbauer_limit'] = dict(rows=c1, spread_of_ratio=(max(ratios) - min(ratios))/abs(np.mean(ratios)),
                                  pass_=bool((max(ratios) - min(ratios))/abs(np.mean(ratios)) < 0.05))
print(out['C1_gegenbauer_limit'], flush=True)

# C2 flipped shell curvature on +1 branch (detection power)
grid = np.concatenate([-np.logspace(math.log10(400), math.log10(5), 30)[:-1], np.linspace(-5, 2.2499, 120)])
c2 = {}
for d in (0.001, 0.1):
    bg = Background(+1, d, e_h=P[d]['e_h'], y_b=P[d]['y_b'], polish=False); bg._sv = bg.shell_values()
    wide = np.concatenate([-np.logspace(6, math.log10(5), 80)[:-1], np.linspace(-5, 2.2499, 60)])
    roots, _ = scan_roots(lambda m: M_variant(bg, m, 'flip_sigma2'), wide)
    c2[str(d)] = dict(wrong_B=bg._sv['g_b'] - bg._sv['sigma2_b']/2, roots=roots,
                      growth=[growth_exponent(r) for r in roots])
out['C2_flipped_shell_curvature_plus_branch'] = c2
out['C2_detects_unstable_modes'] = bool(all(len(v['roots']) > 0 for v in c2.values()))
print(c2, flush=True)

# consistency of the renormalised integrator with the direct one (original shell, true operator)
bgc = Background(-1, 0.001, e_h=O[0.001]['e_h'], y_b=O[0.001]['y_b'], polish=False); svc = bgc.shell_values()
def M_renorm(m):
    st = bgc.scalar_XZ_renorm(m); return svc['B']*st[3] + 3*(m + 4)*st[4]/svc['rho_b']**2
out['renorm_integrator_check'] = dict(root_renorm=brentq(M_renorm, -7.8, -7.6, xtol=1e-12), root_direct=-7.7178716262)
print(out['renorm_integrator_check'], flush=True)

# C3 drop gravitational term on original shell
bg = Background(-1, 0.001, e_h=O[0.001]['e_h'], y_b=O[0.001]['y_b'], polish=False); bg._sv = bg.shell_values()
roots, _ = scan_roots(lambda m: M_variant(bg, m, 'no_grav'), grid)
out['C3_no_grav_term_original'] = dict(roots=roots, correct=-7.7178716252, moved=bool(all(abs(r + 7.7178716) > 1e-3 for r in roots)))
print(out['C3_no_grav_term_original'], flush=True)

# C4 small-detuning limit of the calibration
r = {}
for d in (0.001, 0.003):
    bgd = Background(-1, d, e_h=O[d]['e_h'], y_b=O[d]['y_b'], polish=False)
    r[d] = brentq(lambda m: bgd.scalar_M(m, rtol=1e-13, y0=1e-5), -7.8, -7.6, xtol=1e-13)
slope = (r[0.003] - r[0.001])/0.002; icpt = r[0.001] - slope*0.001
c = C_REG; closed = -4*(3*c*c - 4*c + 8)/(c*(3*c + 4))
out['C4_small_detuning_limit'] = dict(roots=r, slope=slope, intercept=icpt, closed_form_4D=closed,
                                      diff=icpt - closed, pass_=bool(abs(icpt - closed) < 1e-5))
print(out['C4_small_detuning_limit'], flush=True)

# C5 perturbed c
c5 = {}
for fac in (1 - 1e-3, 1 + 1e-3):
    cc = C_REG*fac
    bgo = Background(-1, 0.001, c=cc, e_h=O[0.001]['e_h'], y_b=O[0.001]['y_b'], polish=True)
    ro = brentq(lambda m: bgo.scalar_M(m), -7.9, -7.5, xtol=1e-12)
    bgp = Background(+1, 0.001, c=cc, e_h=P[0.001]['e_h'], y_b=P[0.001]['y_b'], polish=True)
    rp, vals = scan_roots(lambda m: bgp.scalar_M(m), grid)
    c5[str(fac)] = dict(original_root=ro, original_polish=bgo.polish_info['residual_after'],
                        plus_roots=rp, plus_min_Mhat=float(min(vals)), plus_polish=bgp.polish_info['residual_after'],
                        plus_phi_b=bgp.shell_values()['phi_b'])
c5['dmu2_dc_original'] = (c5[str(1 + 1e-3)]['original_root'] - c5[str(1 - 1e-3)]['original_root'])/(2e-3*C_REG)
out['C5_perturbed_c'] = c5
print(c5, flush=True)
out['runtime_s'] = time.time() - t0
out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(HERE/'CONTROLS_RESULTS.json').write_text(json.dumps(out, indent=2, default=float) + '\n')
print('done', out['runtime_s'])
