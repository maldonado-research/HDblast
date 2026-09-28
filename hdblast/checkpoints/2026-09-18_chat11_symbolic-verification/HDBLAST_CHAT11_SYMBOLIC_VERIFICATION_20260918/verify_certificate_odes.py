#!/usr/bin/env python3
"""HDBLAST Chat 11 - SYMBOLIC check that the ODEs actually integrated by the interval certificates
(Chat 9 certify_tachyon.py, Chat 10 certify_existence.py / certify_stability.py) are the verified master system.
Certificate variables:  s = phi',  H = rho'/rho,  w = 1/rho^2,  G = U1/s,  X = s chi/3 = -(psi' + 2H psi),  R = X/psi,  lam = mu2 + 4.
Run:  ../../.hdblast_venv/bin/python verify_certificate_odes.py"""
import json, sys
import sympy as sp
y, mu2, U1, U2, s2 = sp.symbols('y mu2 U1 U2 sigma2')
rho, ph, psi = [sp.Function(n)(y) for n in ('rho', 'phi0', 'psi')]
lam = mu2 + 4
r1, p1 = sp.diff(rho, y), sp.diff(ph, y)
rho2 = r1**2/rho - 1/rho - rho*p1**2/3            # background rho''
ph2 = U1 - 4*r1/rho*p1                            # background phi0''
H, w, s, G = r1/rho, 1/rho**2, p1, U1/p1
def red(e):                                       # apply background ODEs (dU1/dy = U2 phi0') and the master equation
    master = -2*(H - ph2/p1)*sp.diff(psi, y) - (-sp.Rational(4, 3)*p1**2 - 4*H*ph2/p1 + (2 + mu2)/rho**2)*psi
    e = e.subs(sp.diff(psi, y, 2), master)
    for _ in range(2): e = e.subs({sp.diff(rho, y, 2): rho2, sp.diff(ph, y, 2): ph2})
    return sp.simplify(e)
def D(e): return sp.diff(e, y) + sp.diff(e, U1)*U2*p1     # total y-derivative including U1(phi0(y))
res = {}
def record(n, ok): res[n] = bool(ok); print("[%s] %s" % ("PASS" if ok else "FAIL", n), flush=True)
record("B1_background_H'=-w-s^2/3", red(D(H) - (-w - s**2/3)) == 0)
record("B2_background_w'=-2Hw", red(D(w) + 2*H*w) == 0)
record("B3_background_s'=U1-4Hs", red(D(s) - (U1 - 4*H*s)) == 0)
record("B4_identity_G'=U2-G^2+4HG", red(D(G) - (U2 - G**2 + 4*H*G)) == 0)
X = -(sp.diff(psi, y) + 2*H*psi)
record("C1_X'=[lam w-(2/3)s^2]psi+(2G-8H)X  (equivalent to the master equation)", red(D(X) - ((lam*w - sp.Rational(2, 3)*s**2)*psi + (2*G - 8*H)*X)) == 0)
R = X/psi
record("C2_Riccati_R'=R^2+(2G-6H)R+lam w-(2/3)s^2", red(D(R) - (R**2 + (2*G - 6*H)*R + lam*w - sp.Rational(2, 3)*s**2)) == 0)
r = y*R; c1 = 1 + 2*y*G - 6*y*H; c0 = lam*y**2*w - sp.Rational(2, 3)*y**2*s**2
record("C3_cone_form_y r'=r^2+c1 r+c0", red(y*D(r) - (r**2 + c1*r + c0)) == 0)
nu = sp.sqrt(sp.Rational(9, 4) - mu2); rp = -sp.Rational(5, 2) - nu
record("C4_indicial_root_r+=-5/2-nu_solves_r^2+5r+lam=0_(cone limit c1->5, c0->lam)", sp.simplify(rp**2 + 5*rp + lam) == 0)
B = ph2/p1 + s2/2
record("C5_B=G-4H+sigma''/2", sp.simplify(B - (G - 4*H + s2/2)) == 0)
chi = 3*X/s
mism = red(D(chi) + 2*s*psi + (s2/2)*chi)
record("C6_scalar_junction_=(3 psi/s)(lam w+B R)", sp.simplify(mism - (3*psi/s)*(lam*w + B*R)) == 0)
ok = all(res.values())
json.dump(dict(sympy=sp.__version__, all_pass=ok, checks=res), open("CERTIFICATE_ODES_RESULT.json", "w"), indent=1)
print("ALL PASS" if ok else "SOME CHECKS FAILED"); sys.exit(0 if ok else 1)
