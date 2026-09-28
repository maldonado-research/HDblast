"""Gauge test: the residual GN gauge mode (xi^y = eps Y):  chi = phi' e, chi' = phi'' e, Q = E' = -e/rho^2, psi = H e
must (i) solve the bulk flow for every mu2 and (ii) leave the mismatch M unchanged.  We contaminate the cone data with a large
admixture of it and compare shell data with the analytic gauge mode and M with the uncontaminated M."""
import json, numpy as np, gn_lib
from gn_lib import *
Ip = I_plus(); c = 2/Ip - 4/3; t = 1e-3
sol, bg = polish_shell(t, c, 8.785514640315845e-07, 8.228201550106348, n=2000); C = bg_coeffs(bg)
mu2 = np.array([-20., -7.7178716, -4.3, -1.0, 0.7, 2.0])
X0, _ = integrate_pert(mu2, C, sol['phi_h']); M0, chis0, _ = mismatch(X0, C, sol)
orig = gn_lib.cone_start
def contaminated(mu2_, C_, phi_h):
    X, sp = orig(mu2_, C_, phi_h); e = 1e3*X[0]*C_['y'][0]      # psi_gauge(y0) = H e ~ 1000 chi(y0)
    g = np.array([C_['s'][0]*e, C_['spp'][0]*e, -e*C_['ir2'][0], C_['H'][0]*e]); contaminated.e = e
    return X + g, sp
gn_lib.cone_start = contaminated
X1, _ = integrate_pert(mu2, C, sol['phi_h']); M1, chis1, _ = mismatch(X1, C, sol)
e = contaminated.e
gb = np.array([C['s'][-1]*e, C['spp'][-1]*e, -e*C['ir2'][-1], C['H'][-1]*e])
flow_err = np.max(np.abs((X1 - X0) - gb)/np.abs(gb), axis=0)
out = dict(mu2=mu2.tolist(), gauge_mode_flow_relerr=flow_err.tolist(), M_clean=M0.tolist(), M_contaminated=M1.tolist(),
           rel_change_M=(np.abs(M1 - M0)/np.abs(M0).max()).tolist(), rel_change_chi_s=(np.abs(chis1 - chis0)/np.abs(chis0)).tolist(),
           size_of_gauge_part_vs_clean=(np.abs(gb[:, None]).max(axis=0)/np.abs(X0).max(axis=0)).tolist())
print(json.dumps(out, indent=1)); json.dump(out, open('gauge_test.json', 'w'), indent=1)
# (iii) pure gauge data alone: must reproduce the analytic gauge mode at the shell for every mu2 (tests MC + EV-E + scalar eq. for mu2 != 0)
def pure(mu2_, C_, phi_h):
    X, sp = orig(mu2_, C_, phi_h); e = np.ones_like(X[0])
    return np.array([C_['s'][0]*e, C_['spp'][0]*e, -e*C_['ir2'][0], C_['H'][0]*e]), sp
gn_lib.cone_start = pure
Xg, _ = integrate_pert(mu2, C, sol['phi_h'])
gb1 = np.array([C['s'][-1], C['spp'][-1], -C['ir2'][-1], C['H'][-1]])
pure_err = np.max(np.abs(Xg - gb1[:, None])/np.abs(gb1[:, None]), axis=0)
Mg = mismatch(Xg, C, sol)[0]
out['pure_gauge_flow_relerr'] = pure_err.tolist(); out['pure_gauge_M'] = Mg.tolist()
print('pure gauge flow rel err', pure_err, ' M on pure gauge', Mg); json.dump(out, open('gauge_test.json', 'w'), indent=1)
