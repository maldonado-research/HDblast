"""Mechanism screen, part 8: removing the residual vacuum with a quadratic tension term (MODEL CHANGE).

sigma = 2W + delta(1 + c phi + d phi^2/2).  At phi = 1 the detuning becomes A = 1 + c + d/2 and the scalar
junction slope B = c + d.  Replacing (1+c) -> A and c -> B in the 22 Sept series gives the candidate
   H_vac^2 = delta A/27 + delta^2 [A^2/36 - B^2/384] + O(delta^3)          (checked numerically below)
so the '+1' endpoint is vacuum-free at d = d*(delta) (root of the series).
 1. d*(delta) for delta = 0.1, 0.01, 0.001.
 2. Numerical check of the generalised series with a '+1'-branch solver (a copy of the package solver
    with the generalised tension), for delta = 0.01 and several d on the path toward d*.
 3. Continuation of the INITIAL unstable shell (Chat13/14 solver) from d = 0 to d* at delta = 0.1 and 0.01:
    does the blast's initial state survive the tuning?
"""
import math, sys, time, warnings
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root, brentq
from scipy.special import eval_gegenbauer
from common import *
sys.path.insert(0, str(HERE/'pilot_5d'))
import rolloff5d_v1 as R_
warnings.filterwarnings('ignore')

def series(delta, c, d):
    A = 1 + c + d/2; B = c + d
    return delta*A/27 + delta**2*(A*A/36 - B*B/384)

def dstar(delta, c=C):
    return brentq(lambda d: series(delta, c, d), -6, 0)

def pot(e):
    w = 1/3 + e*e + e**3/3; wp = e*(2+e); p = 1+e
    u = .5*wp*wp - (2/3)*w*w; up = wp*(2*p - (4/3)*w); upp = 4*p*p + 2*wp - (4/3)*(wp*wp + 2*p*w)
    return u, up, upp

def plus_branch(delta, d, rtol=2e-12, y0=1e-4, guess=None):
    """Copy of the 22 Sept plus-branch solver with sigma generalised to include d."""
    A = 1 + C + d/2; B = C + d
    sig = lambda e: 2*(1/3 + e*e + e**3/3) + delta*(1 + C*(1+e) + d*(1+e)**2/2)
    sig1 = lambda e: 2*e*(2+e) + delta*(C + d*(1+e))
    sig0 = 2/3 + delta*A; hm2 = sig0*sig0/36 - 1/81
    if guess is None:
        yg = 9*np.arcsinh(1/(9*np.sqrt(hm2))); x = np.cosh(yg/9)
        f = eval_gegenbauer(14, 2, x)/eval_gegenbauer(14, 2, 1.)
        g = (1/9)*np.sinh(yg/9)*4*eval_gegenbauer(13, 3, x)/eval_gegenbauer(14, 2, x)
        eg = -delta*B/(2*(g+2)); guess = [math.log10(abs(eg/f)), yg]
    sgn = -1.0 if B > 0 else 1.0
    def integrate(params):
        eh = sgn*10**params[0]; yb = params[1]; u, up, upp = pot(eh)
        aa = -u/36; bb = u*u/4320 - up*up/750; cc = up/10; dd = up*(upp/280 + u/630)
        ini = [y0 + aa*y0**3 + bb*y0**5, eh + cc*y0*y0 + dd*y0**4, 2*cc*y0 + 4*dd*y0**3]
        def rhs(y, v):
            rho, eta, py = v; u, up, _ = pot(eta)
            rad = 1 + rho*rho*(py*py/12 - u/6)
            return [np.sqrt(max(rad, 1e-300)), py, up - 4*np.sqrt(max(rad, 1e-300))*py/rho]
        return solve_ivp(rhs, (y0, yb), ini, method='DOP853', rtol=rtol, atol=[1e-14, 1e-42, 1e-42], max_step=.1)
    def residual(params):
        rho, e, py = integrate(params).y[:, -1]; u, _, _ = pot(e)
        ry = np.sqrt(1 + rho*rho*(py*py/12 - u/6))
        return np.array([ry/rho - sig(e)/6, py + sig1(e)/2])/delta
    fit = root(residual, guess, tol=1e-9, options={'eps': 1e-8})
    res = residual(fit.x)*delta
    rho, e, py = integrate(fit.x).y[:, -1]
    return dict(delta=delta, d=d, H2=float(1/rho**2), rho_b=float(rho), phi_b=float(1+e), junction_residual=res.tolist(),
                series=float(series(delta, C, d)), ok=bool(np.max(np.abs(res)) < 1e-10)), fit.x

def initial_family(delta, n=12):
    g = (math.log10(0.2207*delta**1.8), 8.2282 + 0.9*math.log(1e-3/delta))
    ph, yb, r = R_.solve_shell(R_.Tension(delta, C, 0), guess=g); rows = []
    ds = dstar(delta)
    for d in np.linspace(0, ds, n+1):
        ph, yb, r = R_.solve_shell(R_.Tension(delta, C, d), guess=(math.log10(ph+1), yb))
        rho, phi, s = R_.integrate_y(ph, yb)
        rows.append(dict(d=float(d), phi_b=float(phi), rho_b=float(rho), H0sq_over_delta=float(1/rho**2/delta),
                         cone_phi_plus_1=float(ph+1), max_junction_residual=float(np.max(np.abs(r)))))
    return rows

def main():
    t0 = time.time()
    out = dict(status='series-based tuning + floating-point checks; the quadratic term is a CHANGE of the registered model',
               d_star={str(dl): dstar(dl) for dl in [0.1, 0.01, 0.001]})
    # 2. check the generalised series on the path toward d* at delta = 0.01
    chk = []; guess = None
    ds = out['d_star']['0.01']
    for frac in [0.0, 0.25, 0.5, 0.75, 0.9]:
        d = frac*ds
        try:
            r, guess = plus_branch(0.01, d, guess=guess); chk.append(r)
            print('plus-branch delta=0.01 d=%.4f H2=%.6e series=%.6e rel=%.2e ok=%s' % (d, r['H2'], r['series'], (r['H2']-r['series'])/r['H2'], r['ok']), flush=True)
        except Exception as e:
            chk.append(dict(d=d, error=str(e))); print('fail', d, e)
    out['plus_branch_series_check_delta_0p01'] = chk
    good = [r for r in chk if r.get('ok')]
    out['series_check_max_abs_H2_error_over_delta_cubed'] = max(abs(r['H2']-r['series'])/0.01**3 for r in good) if good else None
    out['control_d0_matches_package_H2'] = abs(chk[0]['H2'] - 0.0005986962494778341) < 1e-12 if chk and chk[0].get('ok') else False
    # 3. initial shell survives?
    for dl in [0.1, 0.01]:
        rows = initial_family(dl); out['initial_shell_family_delta_%s' % dl] = rows
        print('initial shell delta=%g: d=0 phi_b=%.5f H0^2/delta=%.5f ; d=d*=%.4f phi_b=%.5f H0^2/delta=%.5f res=%.1e' % (
            dl, rows[0]['phi_b'], rows[0]['H0sq_over_delta'], rows[-1]['d'], rows[-1]['phi_b'], rows[-1]['H0sq_over_delta'],
            rows[-1]['max_junction_residual']), flush=True)
    out['runtime_s'] = time.time() - t0
    dump('M8_QUADRATIC_TENSION_TUNING.json', out)
    print(out['d_star'], out['series_check_max_abs_H2_error_over_delta_cubed'], out['control_d0_matches_package_H2'])

if __name__ == '__main__':
    main()
