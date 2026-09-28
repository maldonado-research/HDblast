"""Third route to the norm: manifestly positive bulk integral in gauge-invariant variables
   X = chi + phi' rho^2 Q,  Psi = psi + H rho^2 Q   (invariant under the GN residual gauge eps),
   N_bulk = 2 int rho^2 X^2 dy + 12 int rho^2 Psi^2 dy      (two Z2 copies),
compared with the on-shell-action norms N_scalar, N_metric (d/dmu2 of the reduced boundary action).
Also evaluates the Frolov-Kofman-type quantities Q_FK = int rho^2 X^2/6 + int rho^2 Psi^2 and lam = mu2+4."""
import json, numpy as np
from gn_lib import *
from run_roots import find_root, norm_data
Ip = I_plus(); c = 2/Ip - 4/3; t = 1e-3
out = {}
for n in (2000, 4000):
    sol, bg = polish_shell(t, c, 8.785514640315845e-07, 8.228201550106348, n=n); C = bg_coeffs(bg)
    r = find_root(C, sol, -9, -6); nd = norm_data(C, sol, r)
    X4, rec = integrate_pert(np.array([r]), C, sol['phi_h'], record_every=1)
    idx = np.array([0] + [i for i, _ in rec]); st, _ = cone_start(np.array([r]), C, sol['phi_h'])
    arr = np.array([st[:, 0]] + [x[:, 0] for _, x in rec])
    chi, dchi, Q, psi = arr.T; y = C['y'][idx]; rho = C['rho'][idx]; s = C['s'][idx]; H = C['H'][idx]
    X = chi + s*rho**2*Q; Psi = psi + H*rho**2*Q
    trap = lambda f: float(np.sum(0.5*(f[1:] + f[:-1])*np.diff(y)))
    IX, IP = trap(rho**2*X**2), trap(rho**2*Psi**2)
    Nb = 2*IX + 12*IP
    # consistency of the gauge-invariant constraint X = -3 (rho^2 Psi)'/(rho^2 phi') via finite differences (interior)
    Yv = rho**2*Psi; dY = np.gradient(Yv, y); k = slice(len(y)//2, len(y)-5)
    cons = float(np.max(np.abs(X[k] + 3*dY[k]/(rho[k]**2*s[k]))/np.abs(X[k])))
    out[n] = dict(mu2=float(r), N_scalar=float(nd[0]), N_metric=float(nd[1]), N_bulk=Nb, ratio_bulk_over_action=Nb/float(nd[0]),
                  int_rho2_X2=IX, int_rho2_Psi2=IP, Q_FK=IX/6 + IP, lam=float(r) + 4, N_FK_implied=(IX/6 + IP)/(float(r) + 4),
                  Z_canonical=float(nd[0])/(sol['rho_b']**2*nd[2]**2), Z_EFT=c*Ip*(1 + 0.75*c), constraint_check=cons,
                  X_nodes=int(np.sum(np.sign(X[1:]) != np.sign(X[:-1]))), Psi_nodes=int(np.sum(np.sign(Psi[1:]) != np.sign(Psi[:-1]))))
    print(out[n])
json.dump(out, open('norm_bulk.json', 'w'), indent=1)
