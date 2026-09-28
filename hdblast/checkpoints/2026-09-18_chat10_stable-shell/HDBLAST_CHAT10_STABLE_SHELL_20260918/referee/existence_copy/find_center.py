#!/usr/bin/env python3
"""NON-RIGOROUS high-precision Newton iteration for the shell root (p = phi_h, y_b) of S_8/5.  Its only output is a
GUESS (centre + approximate inverse Jacobian) written to CENTER_GUESS.json; nothing proved here.  The rigorous
statement is checked afterwards by certify_existence.py.   Usage: python3 find_center.py [n_iter]"""
import sys, json, os, time
from fractions import Fraction as Fr
import cert_core as cc
from cert_core import TM, ONE
import exist_lib as L

HERE = os.path.dirname(os.path.abspath(__file__))
cc.set_degree(2)
T0 = time.time()
p = Fr(-0.9999991215629986)          # float guidance (stable_shell_recheck.solve(1e-3, c, 1.6, ...))
yb = Fr(8.228191724874382)
n_iter = int(sys.argv[1]) if len(sys.argv) > 1 else 4
RAD_BITS = 100
for it in range(n_iter):
    c_mid = cc.ffloor(p); c_rad = ONE >> RAD_BITS
    yb = Fr(int(yb * 2 ** 160), 2 ** 160)
    P_TM = TM([c_mid, c_rad])
    L.GATES.clear()
    st, _ = L.cone_start(P_TM, Fr(c_mid - c_rad, ONE), Fr(c_mid + c_rad, ONE))
    I = L.Integrator(st, L.Y0)
    st = I.advance_to(yb)
    g1, g2 = L.residuals(st)
    phi, s, H, w = [Fr(x.c[0], ONE) for x in st]
    G = [Fr(g1.c[0], ONE), Fr(g2.c[0], ONE)]
    Jp = [Fr(g1.c[1], c_rad), Fr(g2.c[1], c_rad)]
    dsig = 2 * (phi * phi - 1) + L.T_PAR * (L.C_PAR + L.D_PAR * phi)
    ddsig = 4 * phi + L.T_PAR * L.D_PAR
    U1 = (phi * phi - 1) * (Fr(-4, 3) + Fr(10, 3) * phi - Fr(4, 9) * phi ** 3)
    Jy = [(-w - s * s / 3) - dsig * s / 6, (U1 - 4 * H * s) + ddsig * s / 2]
    det = Jp[0] * Jy[1] - Jp[1] * Jy[0]
    Cinv = [[Jy[1] / det, -Jy[0] / det], [-Jp[1] / det, Jp[0] / det]]          # J^-1, rows: (dp, dy)
    dp = -(Cinv[0][0] * G[0] + Cinv[0][1] * G[1])
    dy = -(Cinv[1][0] * G[0] + Cinv[1][1] * G[1])
    print("iter %d  G=(%.3e, %.3e)  dp=%.3e dy=%.3e  remwidth=%.2e  steps=%d  %.0fs" % (
        it, float(G[0]), float(G[1]), float(dp), float(dy), max(x.rhi - x.rlo for x in st) / ONE, I.nstep, time.time() - T0), flush=True)
    print("        J = [[%.6e, %.6e],[%.6e, %.6e]]" % (float(Jp[0]), float(Jy[0]), float(Jp[1]), float(Jy[1])), flush=True)
    p_used, y_used = Fr(c_mid, ONE), yb
    p = Fr(c_mid, ONE) + dp; yb = yb + dy

rnd = lambda x, b: Fr(int(x * 2 ** b), 2 ** b)
out = {"note": "NON-RIGOROUS guess produced by find_center.py; used only to place the box and the preconditioner",
       "p_center": str(rnd(p, 300)), "y_center": str(rnd(yb, 160)),
       "C": [[str(rnd(Cinv[i][j], 80 + 0)) for j in range(2)] for i in range(2)],
       "last_residual_float": [float(G[0]), float(G[1])], "last_update_float": [float(dp), float(dy)],
       "J_float": [[float(Jp[0]), float(Jy[0])], [float(Jp[1]), float(Jy[1])]]}
json.dump(out, open(os.path.join(HERE, "CENTER_GUESS.json"), "w"), indent=1)
print("wrote CENTER_GUESS.json  p=%.20f  y=%.20f" % (float(p), float(yb)))
