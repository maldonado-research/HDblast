"""Steps 3-4: calibration on the ORIGINAL unstable shell, then the boundary eigenvalue problem on the +1 branch.

For every background:
  * scalar sector : M_hat(mu2) on a grid of 360 points in [-400, 2.2499] (bound states must have mu2 < 9/4),
                    sign changes -> brentq roots; independent cross-check with the UNREDUCED longitudinal-gauge
                    equations (psi, chi, chi') and their Gauss-constraint residual
  * tensor sector : h'(y_b) mismatch on the same grid; exact zero mode at mu2 = 0 expected
  * convergence   : integrator rtol in {1e-10,1e-12,1e-13} x cone start y0 in {1e-3,1e-4,1e-5}
  * analytic flags: B = phi''/phi' + sigma''/2 at the shell (B>0 => spectrum real and mu2 > -4; see README)
  * l=1 harmonic  : physical only if B phi'_b = 0
  * l=0 static    : singular values of the static junction Jacobian (non-singular => no static zero mode)
Growth exponent of a mode: p = -3/2 + sqrt(9/4 - mu2) (e^{p H tau} on the shell).
Uses at most 2 worker processes.  Output: SPECTRUM_RESULTS.json
"""
import sys, json, time, math, hashlib, os
sys.dont_write_bytecode = True
os.environ.setdefault('OMP_NUM_THREADS', '1')
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from multiprocessing import Pool
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gi_core import Background, growth_exponent

GRID = np.concatenate([-np.logspace(math.log10(400), math.log10(5), 60)[:-1], np.linspace(-5, 2.2499, 301)])
CONV_MU2 = [-400, -10, -7.7, -4.0, -2.0, -1.0, 0.0, 1.0, 2.0, 2.2499]
VARIANTS = [(rt, y0) for rt in (1e-10, 1e-12, 1e-13) for y0 in (1e-3, 1e-4, 1e-5)]
RECORDED = dict(mu2_chat9_GN=-7.717871625176294, mu2_chat9_orchestrator=-7.7178716187365355,
                rate_prior=1.657193631245, rate_controls_spectral=1.657193663767,
                tscan_orchestrator={0.01: -7.700561458, 0.003: -7.7140236})

def roots_of(f, grid, vals, xtol=1e-13):
    out, exact = [], []
    for i in range(len(grid) - 1):
        a, b = vals[i], vals[i + 1]
        if a == 0.0:
            exact.append(float(grid[i]))
        elif a*b < 0:
            out.append(float(brentq(f, grid[i], grid[i + 1], xtol=xtol, rtol=1e-15)))
    if vals[-1] == 0.0: exact.append(float(grid[-1]))
    return out, exact

def analyse(task):
    branch, delta, e_h, y_b = task
    t0 = time.time()
    v = +1 if branch == 'plus' else -1
    bg = Background(v, delta, e_h=e_h, y_b=y_b, polish=False)
    sv = bg.shell_values()
    res = dict(branch=branch, delta=delta, e_h=e_h, y_b=y_b, shell={k: float(x) for k, x in sv.items()})
    res['H_brane'] = 1/sv['rho_b']
    # scalar + tensor + longitudinal scan
    Ms = np.array([bg.scalar_M(m) for m in GRID])
    Ts = np.array([bg.tensor_mismatch(m) for m in GRID])
    LGRID = GRID[::4]
    Ls = [bg.longitudinal_M(m) for m in LGRID]
    Lm = np.array([a for a, _ in Ls]); Lg = np.array([abs(b) for _, b in Ls])
    sroots, _ = roots_of(lambda m: bg.scalar_M(m), GRID, Ms)
    lroots, _ = roots_of(lambda m: bg.longitudinal_M(m)[0], LGRID, Lm)
    troots, texact = roots_of(lambda m: bg.tensor_mismatch(m), GRID, Ts)
    res['scalar'] = dict(grid_points=len(GRID), grid_range=[float(GRID[0]), float(GRID[-1])],
                         min_Mhat=float(Ms.min()), max_Mhat=float(Ms.max()), Mhat_at_0=float(Ms[np.argmin(abs(GRID))]),
                         Mhat_at_minus4=float(bg.scalar_M(-4.0)),
                         roots=sroots, growth_exponents=[growth_exponent(r) for r in sroots],
                         longitudinal_crosscheck=dict(roots=lroots, sign_agreement_fraction=float(np.mean(np.sign(Lm) == np.sign(Ms[::4]))), grid_points=len(LGRID),
                                                      max_gauss_residual=float(Lg.max())))
    # details of the M decomposition at a few points
    det = {}
    for m in (-4.0, -1.0, 0.0, 2.0):
        f = bg.scalar_M(m, full=True)
        det[str(m)] = dict(B_Xb=float(f['B']*f['X_b']), grav_term=float(3*(m + 4)*f['Z_b']/sv['rho_b']**2),
                           ratio=float(abs(3*(m + 4)*f['Z_b']/sv['rho_b']**2)/abs(f['B']*f['X_b'])))
    res['scalar']['M_decomposition'] = det
    res['tensor'] = dict(min=float(Ts.min()), max=float(Ts.max()), sign_change_roots=troots, exact_zeros_on_grid=texact,
                         zero_mode_mu2=0.0,
                         zero_mode_norm_int_rho2=float(np.trapezoid(bg.integrate(bg.e_h, bg.y_b, dense=True).sol(np.linspace(bg.y0, bg.y_b, 20001))[0]**2,
                                                                     np.linspace(bg.y0, bg.y_b, 20001))),
                         bound_states_between_0_and_9over4=[r for r in troots if 1e-9 < r < 2.25])
    # convergence
    conv = {}
    for m in CONV_MU2:
        vals = [bg.scalar_M(m, y0=y0, rtol=rt) for rt, y0 in VARIANTS]
        tv = [bg.tensor_mismatch(m, y0=y0, rtol=rt) for rt, y0 in VARIANTS]
        conv[str(m)] = dict(scalar_min=float(min(vals)), scalar_max=float(max(vals)), scalar_spread=float(max(vals) - min(vals)),
                            tensor_min=float(min(tv)), tensor_max=float(max(tv)), tensor_spread=float(max(tv) - min(tv)))
    res['convergence'] = conv
    # root convergence (calibration branch)
    if sroots:
        rc = {}
        for r in sroots:
            rc[repr(r)] = {f'rtol={rt:g},y0={y0:g}': float(brentq(lambda m: bg.scalar_M(m, y0=y0, rtol=rt), r - 0.05, r + 0.05, xtol=1e-13))
                           for rt, y0 in VARIANTS}
        res['root_convergence'] = rc
    # analytic flags
    B = sv['B']
    res['analytic'] = dict(B=float(B), B_positive=bool(B > 0),
                           implies=('spectrum real and mu2 > -4 (energy + Wronskian identities)' if B > 0 else
                                    'no bound from the identities; mu2 < -4 allowed'),
                           l1_physical_condition_B_phi_b=float(B*sv['phi_y_b']))
    # static (l=0, p=0) Jacobian
    sgn = math.copysign(1, bg.e_h)
    def F(x): return bg.junction_residual(sgn*10**x[0], x[1])
    x0 = np.array([math.log10(abs(bg.e_h)), bg.y_b]); J = np.zeros((2, 2)); h = 1e-6
    for k in range(2):
        dx = np.zeros(2); dx[k] = h; J[:, k] = (F(x0 + dx) - F(x0 - dx))/(2*h)
    res['static_jacobian'] = dict(J=J.tolist(), singular_values=np.linalg.svd(J, compute_uv=False).tolist(),
                                  cond=float(np.linalg.cond(J)))
    res['runtime_s'] = time.time() - t0
    print(branch, delta, 'done', res['runtime_s'], 'roots', sroots, 'troots', troots, flush=True)
    return res

if __name__ == '__main__':
    t0 = time.time()
    B = json.loads((HERE/'BACKGROUNDS.json').read_text())
    tasks = [('original', r['delta'], r['gi_core']['e_h'], r['gi_core']['y_b']) for r in B['original']]
    tasks += [('plus', r['delta'], r['gi_core']['e_h'], r['gi_core']['y_b']) for r in B['plus']]
    with Pool(2) as pool:
        results = pool.map(analyse, tasks, chunksize=1)
    orig = {r['delta']: r for r in results if r['branch'] == 'original'}
    plus = [r for r in results if r['branch'] == 'plus']
    o = orig[0.001]
    cal_mu2 = o['scalar']['roots'][0] if len(o['scalar']['roots']) == 1 else None
    finest = o['root_convergence'][repr(cal_mu2)]['rtol=1e-13,y0=1e-05'] if cal_mu2 is not None else None
    calibration = dict(n_bound_states=len(o['scalar']['roots']), mu2=finest, growth=growth_exponent(finest) if finest else None,
                       diff_vs_chat9_GN=finest - RECORDED['mu2_chat9_GN'] if finest else None,
                       diff_vs_chat9_orchestrator=finest - RECORDED['mu2_chat9_orchestrator'] if finest else None,
                       growth_diff_vs_prior=growth_exponent(finest) - RECORDED['rate_prior'] if finest else None,
                       growth_diff_vs_controls_spectral=growth_exponent(finest) - RECORDED['rate_controls_spectral'] if finest else None,
                       tscan={str(d): dict(ours=orig[d]['scalar']['roots'], recorded=RECORDED['tscan_orchestrator'][d]) for d in (0.003, 0.01)},
                       B_original=o['analytic']['B'])
    calibration['pass'] = bool(finest is not None and abs(calibration['diff_vs_chat9_GN']) < 1e-7 and
                               all(len(orig[d]['scalar']['roots']) == 1 and abs(orig[d]['scalar']['roots'][0] - RECORDED['tscan_orchestrator'][d]) < 1e-6
                                   for d in (0.003, 0.01)))
    summary = []
    for r in plus:
        summary.append(dict(delta=r['delta'], H_brane=r['H_brane'], B=r['analytic']['B'], scalar_roots=r['scalar']['roots'],
                            min_Mhat=r['scalar']['min_Mhat'], Mhat_0=r['scalar']['Mhat_at_0'], Mhat_m4=r['scalar']['Mhat_at_minus4'],
                            longitudinal_roots=r['scalar']['longitudinal_crosscheck']['roots'],
                            tensor_roots=r['tensor']['sign_change_roots'], tensor_exact_zero=r['tensor']['exact_zeros_on_grid'],
                            max_conv_spread=max(v['scalar_spread'] for v in r['convergence'].values()),
                            static_cond=r['static_jacobian']['cond'], l1_B_phi=r['analytic']['l1_physical_condition_B_phi_b']))
    out = dict(status='NUMERICAL', calibration=calibration, plus_summary=summary, details=results, recorded_inputs=RECORDED,
               grid=GRID.tolist(), variants=VARIANTS, runtime_s=time.time() - t0,
               script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               core_sha256=hashlib.sha256((HERE/'gi_core.py').read_bytes()).hexdigest())
    (HERE/'SPECTRUM_RESULTS.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(calibration, indent=1)); print(json.dumps(summary, indent=1)); print('total', out['runtime_s'])
