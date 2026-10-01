#!/usr/bin/env python3
"""B4 step 1: exact residual vacuum Lambda_res(d) and frozen scalar phi_s(d) of the tuned model at delta = 0.1.

Exact static '+1' branch (dS-sliced; copy of the M8 / 22 Sept plus-branch solver, sigma = 2W + delta(1 + c phi + d phi^2/2))
for d/d* in {0.90 ... 0.995}, where H^2 > 0, then cubic (and quadratic) polynomial fits in d, extrapolated through
H^2 = 0 (to d*_exact) and to d/d* = 1, 1.01, 1.05.  H0^2 = 1/rho_b^2 of the tuned initial shell (M8: rho_b = 7.838285538...).
Control C2: d = 0 reproduces M2 (H+^2/H0^2 = 0.406168, phi = 0.991435).  Output: B4_STATIC_LAMBDA.json.  Numerical."""
import math, json, warnings, hashlib
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
warnings.filterwarnings('ignore')
HERE = Path(__file__).resolve().parent
C = 2/1.0357712571566784 - 4/3
DELTA = 0.1
DSTAR_M8 = -3.106933495673783            # M8 series root, delta = 0.1 (used by A1)
RHO_B_DSTAR = 7.838285538074403          # M8 initial shell at d*, delta = 0.1
RHO_B_D0 = 7.836278852194686             # M8 initial shell at d = 0, delta = 0.1
ESC = 1e-3

def series(delta, c, d):
    A = 1 + c + d/2; B = c + d
    return delta*A/27 + delta**2*(A*A/36 - B*B/384)

def pot(e):
    w = 1/3 + e*e + e**3/3; wp = e*(2+e); p = 1+e
    u = .5*wp*wp - (2/3)*w*w; up = wp*(2*p - (4/3)*w); upp = 4*p*p + 2*wp - (4/3)*(wp*wp + 2*p*w)
    return u, up, upp

def plus_branch(delta, d, guess, rtol=2e-12, y0=1e-4):
    A = 1 + C + d/2; B = C + d
    sig = lambda e: 2*(1/3 + e*e + e**3/3) + delta*(1 + C*(1+e) + d*(1+e)**2/2)
    sig1 = lambda e: 2*e*(2+e) + delta*(C + d*(1+e))
    def integrate(params):
        eh = params[0]*ESC; yb = params[1]; u, up, upp = pot(eh)    # signed vertex offset phi_h - 1 (eh crosses 0 at d = -c)
        aa = -u/36; bb = u*u/4320 - up*up/750; cc = up/10; dd = up*(upp/280 + u/630)
        ini = [y0 + aa*y0**3 + bb*y0**5, eh + cc*y0*y0 + dd*y0**4, 2*cc*y0 + 4*dd*y0**3]
        def rhs(y, v):
            rho, eta, py = v; u, up, _ = pot(eta)
            rad = 1 + rho*rho*(py*py/12 - u/6)
            return [np.sqrt(max(rad, 1e-300)), py, up - 4*np.sqrt(max(rad, 1e-300))*py/rho]
        return solve_ivp(rhs, (y0, yb), ini, method='DOP853', rtol=rtol, atol=[1e-14, 1e-42, 1e-42], max_step=.1)
    def residual(params):
        if not (np.all(np.isfinite(params)) and params[1] > y0 and abs(params[0]*ESC) < 0.5): return np.array([1e3, 1e3])
        sol = integrate(params)
        if not sol.success or not np.all(np.isfinite(sol.y[:, -1])): return np.array([1e3, 1e3])
        rho, e, py = sol.y[:, -1]; u, _, _ = pot(e)
        ry = np.sqrt(1 + rho*rho*(py*py/12 - u/6))
        return np.array([ry/rho - sig(e)/6, py + sig1(e)/2])/delta
    fit = root(residual, guess, tol=1e-11, options={'eps': 1e-8})
    res = residual(fit.x)*delta
    rho, e, py = integrate(fit.x).y[:, -1]
    return dict(d=d, H2=float(1/rho**2), rho_b=float(rho), phi_s=float(1+e), sigma_s=float(sig(e)),
                junction_residual=float(np.max(np.abs(res))), ok=bool(np.max(np.abs(res)) < 1e-10)), fit.x

def main():
    out = dict(status='numerical (exact static branch + polynomial extrapolation through H^2 = 0)', delta=DELTA, c=C, d_star_M8=DSTAR_M8)
    # continuation from d = 0 toward d*
    from scipy.special import eval_gegenbauer
    sig0 = 2/3 + DELTA*(1 + C); hm2 = sig0*sig0/36 - 1/81
    yg = 9*np.arcsinh(1/(9*np.sqrt(hm2))); x = np.cosh(yg/9)
    f = eval_gegenbauer(14, 2, x)/eval_gegenbauer(14, 2, 1.)
    g = (1/9)*np.sinh(yg/9)*4*eval_gegenbauer(13, 3, x)/eval_gegenbauer(14, 2, x)
    eg = -DELTA*C/(2*(g+2)); guess = np.array([-abs(eg/f)/ESC, yg])
    r0, guess0 = plus_branch(DELTA, 0.0, guess)
    if not r0['ok']:   # try the other sign convention of the guess
        pass
    out['C2_d0'] = dict(r0, H2_over_H0sq=r0['H2']*RHO_B_D0**2, target_H2_over_H0sq=0.40616846475555957, target_phi=0.9914346595842017,
                        PASS=bool(abs(r0['H2']*RHO_B_D0**2 - 0.40616846475555957) < 1e-5 and abs(r0['phi_s'] - 0.9914346595842017) < 1e-5))
    print('C2', out['C2_d0'])
    path = list(np.linspace(0, 0.9, 37)) + [0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.975, 0.98, 0.985, 0.99, 0.9925, 0.995, 0.9975]
    rows = []; hist = [(0.0, guess0)]
    for q in path[1:]:
        if len(hist) >= 2:   # linear extrapolation of (eh, y_b) in q
            (q1, g1), (q2, g2) = hist[-2], hist[-1]; guess = g2 + (g2 - g1)*(q - q2)/(q2 - q1)
        else: guess = hist[-1][1]
        r, gnew = plus_branch(DELTA, q*DSTAR_M8, guess)
        if not r['ok']:
            r, gnew = plus_branch(DELTA, q*DSTAR_M8, hist[-1][1])
        if r['ok']: hist.append((q, gnew))
        r['params'] = [float(gnew[0]*ESC), float(gnew[1])]
        r['q'] = q; r['H2_over_H0sq'] = r['H2']*RHO_B_DSTAR**2; r['series_over_H0sq'] = series(DELTA, C, q*DSTAR_M8)*RHO_B_DSTAR**2
        rows.append(r); print(r['params'], 'q=%.4f d=%.5f H2/H0^2=%.6e series=%.6e phi_s=%.6f ok=%s res=%.1e' % (q, r['d'], r['H2_over_H0sq'], r['series_over_H0sq'], r['phi_s'], r['ok'], r['junction_residual']), flush=True)
    out['branch'] = rows
    fitq = [r for r in rows if r['q'] >= 0.899 and r['ok']]
    dd = np.array([r['d'] for r in fitq]); L = np.array([r['H2_over_H0sq'] for r in fitq]); P = np.array([r['phi_s'] for r in fitq])
    S = np.array([r['sigma_s'] for r in fitq])
    fits = {}
    for deg in (2, 3):
        cL = np.polyfit(dd, L, deg); cP = np.polyfit(dd, P, deg); cS = np.polyfit(dd, S, deg)
        roots = [x.real for x in np.roots(cL) if abs(x.imag) < 1e-12 and DSTAR_M8*1.1 < x.real < DSTAR_M8*0.9]
        dex = float(min(roots, key=lambda x: abs(x - DSTAR_M8)))
        ev = {('%.2f' % q): dict(d=q*DSTAR_M8, Lambda_over_H0sq=float(np.polyval(cL, q*DSTAR_M8)), phi_s=float(np.polyval(cP, q*DSTAR_M8)),
                                 sigma_s=float(np.polyval(cS, q*DSTAR_M8))) for q in (0.95, 0.99, 1.0, 1.01, 1.05)}
        fits['deg%d' % deg] = dict(coef_Lambda=cL.tolist(), coef_phi=cP.tolist(), coef_sigma=cS.tolist(), d_star_exact=dex,
                                   slope_at_dstar=float(np.polyval(np.polyder(cL), DSTAR_M8)), at=ev,
                                   max_fit_residual=float(np.max(np.abs(np.polyval(cL, dd) - L))))
    out['fits'] = fits
    out['extrapolation_error'] = {k: abs(fits['deg3']['at'][k]['Lambda_over_H0sq'] - fits['deg2']['at'][k]['Lambda_over_H0sq']) for k in fits['deg3']['at']}
    A = 1 + C + DSTAR_M8/2; B = C + DSTAR_M8
    out['series_slope_over_H0sq'] = (DELTA/54 + DELTA**2*(A/36 - B/192))*RHO_B_DSTAR**2
    out['H0sq_5D'] = 1/RHO_B_DSTAR**2
    out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (HERE/'B4_STATIC_LAMBDA.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps({k: out[k] for k in ('extrapolation_error', 'series_slope_over_H0sq')}, indent=1))
    for k in fits: print(k, fits[k]['d_star_exact'], fits[k]['slope_at_dstar'], fits[k]['at'])

if __name__ == '__main__':
    main()
