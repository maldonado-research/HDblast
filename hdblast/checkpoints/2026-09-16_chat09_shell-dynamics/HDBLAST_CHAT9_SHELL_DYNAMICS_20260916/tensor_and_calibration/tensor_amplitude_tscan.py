#!/usr/bin/env python3
"""HDBLAST Chat 9 - zero-mode normalisation vs detuning t: F_shell^2 = I_plus/I_kept compared with Langlois-Maartens-Wands
F^-2 = sqrt(1+x^2) - x^2 asinh(1/x), x = H/k (hep-th/0006007), evaluated with the kept-side throat curvature k_- = 5/9."""
import json, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import tensor_spectrum as ts
bg = ts.bg
Ip = ts.sol["I_plus"]; c = ts.sol["c"]; km = 5/9
LMW = lambda x: 1/(math.sqrt(1 + x*x) - x*x*math.asinh(1/x))
rows = []
for t in (1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5):   # t = 1e-5: float Newton fails (phi_h + 1 ~ 2e-10, singular Jacobian)
    s = bg.solve_shell(t, c, (math.log10(0.2207*t**1.8), 8.2282 + 0.9*math.log(1e-3/t)), dv=4e-4)
    B = ts.background(2e-4, phi_h=s["phi_h"], v_b=s["v_b"]); rb = B["rho"][-1]
    Ik = (ts.simpson(B["rho"]**2, B["y"]) + B["y"][0]**3/3)/rb**2; H = 1/rb
    F2 = Ip/Ik; x = H/km
    Cfit = (F2 - 1)*2*km**3*Ip/H**2 - math.log(1/H)
    rows.append(dict(t=t, H=H, x=x, I_kept=Ik, F2_shell=F2, F2_LMW=LMW(x), C=Cfit, residual=s["residual"]))
    print("t=%.0e  H=%.6e  H/k_-=%.5f  I_kept=%.9f  F^2_shell-1=%.6e  F^2_LMW(H/k_-)-1=%.6e  ratio=%.4f   C=[(F^2-1) 2k^3 I+/H^2 - ln(1/H)]=%.4f" % (t, H, x, Ik, F2 - 1, LMW(x) - 1, (F2 - 1)/(LMW(x) - 1), Cfit), flush=True)
json.dump(rows, open(os.path.join(HERE, "TENSOR_AMPLITUDE_TSCAN.json"), "w"), indent=1)
