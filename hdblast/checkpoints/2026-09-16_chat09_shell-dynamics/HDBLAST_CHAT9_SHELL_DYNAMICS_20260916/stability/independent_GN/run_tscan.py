"""T2 / t-dependence: lowest mu2 and physical m2 = mu2*h as the detuning t -> 0 (flat BPS limit)."""
import json, math, numpy as np
from gn_lib import *
from run_roots import find_root, norm_data
Ip = I_plus(); cst = 2/Ip - 4/3
out = []
prev = None
for t in [1e-3, 3e-3, 1e-2, 3e-4, 1e-4, 3e-5, 1e-5]:
    al = 0.2206821 * t**1.8; yb = 8.2282 + 0.9*math.log(1e-3/t)
    sol, bg = polish_shell(t, cst, al, yb, n=1500, itmax=25); C = bg_coeffs(bg)
    mu2 = np.linspace(-40, 2.2, 212); X, _ = integrate_pert(mu2, C, sol['phi_h']); M = mismatch(X, C, sol)[0]
    idx = np.where(np.sign(M[1:]) != np.sign(M[:-1]))[0]
    roots = [find_root(C, sol, mu2[i], mu2[i+1], rounds=7) for i in idx]
    nd = norm_data(C, sol, roots[0]) if roots else None
    row = dict(t=t, alpha=sol['alpha'], y_b=float(sol['y_b']), h=float(sol['h']), h_over_hlin=float(sol['h']/(t/(6*Ip))), residual=sol['residual'],
               roots=[float(r) for r in roots], m2_phys=[float(r*sol['h']) for r in roots], N_scalar=float(nd[0]) if nd else None, N_metric=float(nd[1]) if nd else None,
               chi_s=float(nd[2]), psi_s=float(nd[3]), eps=float(nd[4]), Z_canonical=float(nd[0]/(sol['rho_b']**2*nd[2]**2)))
    print(row, flush=True); out.append(row)
json.dump(out, open('tscan.json', 'w'), indent=1)
