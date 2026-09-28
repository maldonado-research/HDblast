#!/usr/bin/env python3
"""FLOAT check of Remark 1.2: with the junction conditions, B and 1/rho_b^2 are explicit functions of phi_b alone:
   B = -2 U_phi/sigma_t' - (2/3) sigma_t + sigma_t''/2,      1/rho_b^2 = sigma_t^2/36 - sigma_t'^2/48 + U/6   (at phi_b)."""
import json
t = 1e-3; c = 0.5975949350280
W = lambda p: 1 - p + p**3/3; W1 = lambda p: p*p - 1; W2 = lambda p: 2*p
U = lambda p: 0.5*W1(p)**2 - (2/3)*W(p)**2; U1 = lambda p: W1(p)*W2(p) - (4/3)*W(p)*W1(p)
for sh in json.load(open("shells.json")):
    d, p = sh["d"], sh["phi_b"]
    sig = 2*W(p) + t*(1 + c*p + d*p*p/2); s1 = 2*W1(p) + t*(c + d*p); s2 = 2*W2(p) + t*d
    B = -2*U1(p)/s1 - (2/3)*sig + s2/2; h = sig**2/36 - s1**2/48 + U(p)/6
    print("d=%.2f  B_closed=%+.8e  B_ode=%+.8e   rho_b_closed=%.6f rho_b_ode=%.6f" % (d, B, sh["B"], h**-0.5, sh["rho_b"]))
