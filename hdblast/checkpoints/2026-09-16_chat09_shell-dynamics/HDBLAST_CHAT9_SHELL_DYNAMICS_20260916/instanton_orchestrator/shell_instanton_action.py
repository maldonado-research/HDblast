#!/usr/bin/env python3
"""Euclidean on-shell action of the O(5)-symmetric continuation of the registered shell solution
(two mirrored 5-balls ds^2 = dy^2 + rho(y)^2 dOmega_4^2 glued at the shell): the 'creation of the shell universe' instanton.
  S_E = -int sqrt(g)[R/2 - (d phi)^2/2 - U] + int sqrt(h) sigma_t + GHY  ->  on shell:
  S_E = -(4/3) Omega_4 int_0^{y_b} rho^4 U dy - (1/3) Omega_4 rho_b^4 sigma_t(phi_b),   Omega_4 = 8 pi^2/3.
Compared with the 4D effective-theory (de Sitter / Hawking-Moss) value  S_E = -24 pi^2 f^2/V = -24 pi^2/V_E."""
import sys, json, math
import numpy as np
sys.path.insert(0, "../background")
import hdblast_background as bg
Ip = bg.I_plus(); c = 2/Ip - 4/3; Om4 = 8*math.pi**2/3
rows = []
for t in [1e-2, 3e-3, 1e-3, 3e-4, 1e-4]:
    guess = (math.log10(0.2207*t**1.8), 8.2282 + 0.9*math.log(1e-3/t))
    sol = bg.solve_shell(t, c, guess, dv=4e-4)
    yend, rec = bg.integrate(sol["phi_h"], sol["v_b"], dv=2e-4, record=True)
    y, rho, phi = rec[:, 0], rec[:, 1], rec[:, 2]
    integrand = rho**4*bg.U(phi)
    bulk = np.trapz(integrand, y) + (rec[0, 0]**5/5)*bg.U(sol["phi_h"])      # + analytic piece on [0, y0] (rho ~ y)
    sig = 2*bg.W(sol["phi_b"]) + t*(1 + c*sol["phi_b"])
    SE = -(4/3)*Om4*bulk - (1/3)*Om4*sol["rho_b"]**4*sig
    # 4D estimates
    SE_eft = -24*math.pi**2*(2*Ip)**2/t                       # leading order: f = 2 I_plus, V = t
    Ikept = np.trapz((rho/sol["rho_b"])**2, y)
    SE_dS = -24*math.pi**2*(2*Ikept)/(3*sol["h"])             # -24 pi^2 M4^2/(3 H^2) with M4^2 = 2 I_kept: -8 pi^2 M4^2/H^2
    rows.append(dict(t=t, S_E_5D=SE, S_E_EFT_leading=SE_eft, ratio_5D_over_EFT=SE/SE_eft, S_E_dS_with_Ikept=SE_dS, ratio_5D_over_dS=SE/SE_dS,
                     bulk_term=-(4/3)*Om4*bulk, brane_term=-(1/3)*Om4*sol["rho_b"]**4*sig))
    print("t=%.0e  S_E(5D) = %.8e   EFT -24pi^2 f^2/V = %.8e  ratio = %.8f   | -8pi^2 (2 I_kept)/h = %.8e ratio = %.10f" % (t, SE, SE_eft, SE/SE_eft, SE_dS, SE/SE_dS))
json.dump(rows, open("SHELL_INSTANTON_ACTION.json", "w"), indent=1)
