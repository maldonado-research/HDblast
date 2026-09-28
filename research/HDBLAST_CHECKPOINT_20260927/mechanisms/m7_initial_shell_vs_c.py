"""Mechanism screen, part 7: can the residual vacuum be tuned away by the tension slope c while keeping
the blast's initial state?

The residual RS vacuum of the '+1' endpoint is H_vac^2 = delta(1+c)/27 + O(delta^2), so it vanishes near
c* ~ -1 + 27 delta/384.  But c also fixes where the initial unstable shell sits.  This script follows the
registered initial-shell family (regular cone near phi=-1, shell near phi~0) by continuation in c from the
registered value downward, at delta = 0.1, 0.01, 0.001, using the Chat13/14 shell solver (copied module,
pilot_5d/rolloff5d_v1.py), and records phi_b(c), H0^2(c), junction residuals, and where the family ends.
Control: at the registered c the continuation start reproduces the archived rho_b.
"""
import math, sys, time, warnings
import numpy as np
from common import *
sys.path.insert(0, str(HERE/'pilot_5d'))
import rolloff5d_v1 as R_
warnings.filterwarnings('ignore')

def family(td, n_coarse=12, tmax=240):
    t0 = time.time()
    g = (math.log10(0.2207*td**1.8), 8.2282 + 0.9*math.log(1e-3/td))
    ph, yb, r = R_.solve_shell(R_.Tension(td, C, 0), guess=g)
    rows = []; c = C; step = (C - (-1.0))/n_coarse; last_ok = None
    def record(c, ph, yb, r):
        rho, phi, s = R_.integrate_y(ph, yb)
        rows.append(dict(c=float(c), phi_h_plus_1=float(ph+1), y_b=float(yb), rho_b=float(rho), phi_b=float(phi),
                         H0sq=float(1/rho**2), max_junction_residual=float(np.max(np.abs(r))),
                         Hvac2_series=float(td*(1+c)/27 + td**2*((1+c)**2/36 - c*c/384))))
    record(c, ph, yb, r)
    while step > 1e-4*abs(C+1) and time.time()-t0 < tmax:
        cn = c - step
        try:
            ph2, yb2, r2 = R_.solve_shell(R_.Tension(td, cn, 0), guess=(math.log10(ph+1), yb))
            if not np.all(np.isfinite(r2)) or np.max(np.abs(r2)) > 1e-9 or not (ph2 + 1 > 0):
                raise RuntimeError('bad root')
            c, ph, yb = cn, ph2, yb2; record(c, ph, yb, r2)
        except Exception:
            step /= 2
    return rows

def main():
    out = {}
    for td in [0.1, 0.01, 0.001]:
        rows = family(td)
        last = rows[-1]
        out[str(td)] = dict(rows=rows, c_registered=C, c_last_converged=last['c'], phi_b_at_last=last['phi_b'],
                            cone_displacement_at_last=last['phi_h_plus_1'],
                            c_star_series=float(-1 + 27*td/384),
                            family_reaches_c_star=bool(last['c'] <= -1 + 27*td/384),
                            fH_at_last=float(last['Hvac2_series']/last['H0sq']),
                            fH_at_registered=float(rows[0]['Hvac2_series']/rows[0]['H0sq']))
        print(td, 'rho_b(reg)=%.6f  last c=%.5f phi_b=%.5f cone=%.3e H0^2/td=%.4f fH_last=%.4f reaches c*: %s' % (
            rows[0]['rho_b'], last['c'], last['phi_b'], last['phi_h_plus_1'], last['H0sq']/td,
            out[str(td)]['fH_at_last'], out[str(td)]['family_reaches_c_star']), flush=True)
    checks = dict(registered_delta_1em3_rho_b_matches_archive=abs(out['0.001']['rows'][0]['rho_b'] - 78.82817714224423) < 1e-5,
                  registered_delta_0p1_rho_b_matches_archive=abs(out['0.1']['rows'][0]['rho_b'] - 7.836278852192727) < 1e-5)
    dump('M7_INITIAL_SHELL_VS_C.json', dict(status='numerical (floating-point continuation; end of family = last converged root, '
                                             'not a proven existence boundary)', families=out, checks=checks))
    print(checks)

if __name__ == '__main__':
    main()
