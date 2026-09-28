"""Check of the identities proved in GN_DERIVATION.md sec. 7 on the numerical eigenmode:
   Q_FK = lam * N_FK,   N_action = 12 Q_FK,  with
   Q_FK = int rho^2 X^2/6 + int rho^2 Psi^2,  N_FK = int 3 Psi^2/(2 phi'^2) + 3 Psi_b^2/(2 phi_b'^2 B),  B = phi''/phi' + sigma''/2."""
import json, numpy as np
from gn_lib import *
from run_roots import find_root, norm_data
Ip = I_plus(); c = 2/Ip - 4/3; t = 1e-3; out = {}
for n in (2000, 4000):
    sol, bg = polish_shell(t, c, 8.785514640315845e-07, 8.228201550106348, n=n); C = bg_coeffs(bg)
    r = find_root(C, sol, -9, -6); nd = norm_data(C, sol, r); lam = r + 4
    X4, rec = integrate_pert(np.array([r]), C, sol['phi_h'], record_every=1)
    idx = np.array([0] + [i for i, _ in rec]); st, sp = cone_start(np.array([r]), C, sol['phi_h'])
    arr = np.array([st[:, 0]] + [x[:, 0] for _, x in rec]); chi, dchi, Q, psi = arr.T
    y = C['y'][idx]; rho = C['rho'][idx]; s = C['s'][idx]; H = C['H'][idx]
    X = chi + s*rho**2*Q; Psi = psi + H*rho**2*Q
    trap = lambda f: float(np.sum(0.5*(f[1:] + f[:-1])*np.diff(y)))
    f3 = 1.5*Psi**2/s**2
    cone_piece = f3[0]*y[0]/(2*sp[0] + 3)            # integrand ~ y^(2s+2) below y0
    B = C['spp'][-1]/s[-1] + 0.5*shell_data(sol)['ddsig']
    QFK = trap(rho**2*X**2)/6 + trap(rho**2*Psi**2)
    NFK_bulk = trap(f3) + cone_piece; NFK_bdy = 1.5*Psi[-1]**2/(s[-1]**2*B); NFK = NFK_bulk + NFK_bdy
    Yb2 = (rho[-1]**2*Psi[-1])**2
    out[n] = dict(mu2=float(r), lam=float(lam), B=float(B), B_rho_b2=float(B*rho[-1]**2), Q_FK=QFK, N_FK_bulk=NFK_bulk, N_FK_boundary=float(NFK_bdy),
                  Q_over_lamN=float(QFK/(lam*NFK)), N_action=float(nd[0]), N_action_over_12Q=float(nd[0]/(12*QFK)),
                  BC_residual_rel=float(abs(B*X[-1] + 3*lam*Psi[-1]/(rho[-1]**2*s[-1]))/abs(B*X[-1])))
    print(out[n])
json.dump(out, open('fk_identity.json', 'w'), indent=1)
