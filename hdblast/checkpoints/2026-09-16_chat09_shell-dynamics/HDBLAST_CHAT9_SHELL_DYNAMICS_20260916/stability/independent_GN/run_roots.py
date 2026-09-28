import json, math, numpy as np, time
from gn_lib import *
def find_root(C, sol, lo, hi, rounds=9, npts=11):
    for r in range(rounds):
        mu2 = np.linspace(lo, hi, npts); X, _ = integrate_pert(mu2, C, sol['phi_h']); M = mismatch(X, C, sol)[0]
        w = np.where(np.sign(M[1:]) != np.sign(M[:-1]))[0]
        if len(w) == 0: return None
        i = w[0]; lo, hi = mu2[i], mu2[i+1]
    # final linear interpolation
    return lo - M[i]*(hi-lo)/(M[i+1]-M[i])
def norm_data(C, sol, r, d=1e-4):
    mu2 = np.array([r-d, r, r+d]); X, _ = integrate_pert(mu2, C, sol['phi_h'])
    M, chis, eps = mismatch(X, C, sol); Es, psis, chis2, eps2 = metric_mismatch(X, C, sol)
    rb4 = C['rho'][-1]**4
    N_scalar = -2*rb4*chis[1]*(M[2]-M[0])/(2*d)
    N_metric = 6*(r+4)*rb4*psis[1]*(Es[2]-Es[0])/(2*d)
    return N_scalar, N_metric, chis[1], psis[1], eps[1], X[:,1]
if __name__ == "__main__":
    Ip = I_plus(); c = 2/Ip - 4/3; t = 1e-3; out = {}
    for n, y0 in [(1000, 2e-3), (2000, 2e-3), (4000, 2e-3), (4000, 5e-4), (4000, 8e-3)]:
        sol, bg = polish_shell(t, c, 8.785514694954969e-07, 8.228201550199884, n=n, y0=y0); C = bg_coeffs(bg)
        r1 = find_root(C, sol, -9, -6); nd1 = norm_data(C, sol, r1)
        print('n=%d y0=%.0e  mu2_0=%.10f  N_scalar=%.8e N_metric=%.8e chi_s=%.5e psi_s=%.5e eps=%.5e' % (n, y0, r1, nd1[0], nd1[1], nd1[2], nd1[3], nd1[4]))
        out['n%d_y0%g' % (n, y0)] = dict(mu2_0=r1, N_scalar_0=nd1[0], N_metric_0=nd1[1], chi_s=nd1[2], psi_s=nd1[3], eps=nd1[4], rho_b=sol['rho_b'], h=sol['h'])
    json.dump(out, open('roots_convergence.json', 'w'), indent=1)
