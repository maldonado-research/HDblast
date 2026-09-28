#!/usr/bin/env python3
"""HDBLAST Chat 9 - part B: doubly tuned inflection-point variant  sigma_t - 2W = t(1 + c phi + d phi^2/2), c = c_star + dc (dc<0), d ~ d0.
Near phi_b = 0: V'/V = lam1 + gamma Theta^2/2, roll toward phi -> +1. Analytic: n_s - 1 = -4a cot(a N), a = sqrt(lam1 gamma/2).
N is counted back from the point of maximal eps_H (inflation does NOT end in this family: max eps_H ~ 0.50)."""
import json, math, os
import numpy as np
import shell_inflation as si
HERE = si.HERE; d0 = 1.1134966; out = {}

def run(dc, d, Nstars=(50, 55, 60)):
    tr = si.evolve(si.c_star + dc, d, 0.03, dN=5e-3, Nmax=3000.0)
    imax = int(np.argmax(tr["e1"]))
    return tr, imax, [si.observables(tr, n, iend=imax) for n in Nstars]

def ns55(dc, d=d0):
    o = run(dc, d, (55,))[2][0]
    return o["n_s"] if o else float("nan")

target = 0.9649
lo, hi = -2.0e-4, -2.0e-5            # n_s increases with |dc|
flo, fhi = ns55(lo), ns55(hi); print("bracket: n_s55(dc=%.1e)=%.5f  n_s55(dc=%.1e)=%.5f" % (lo, flo, hi, fhi), flush=True)
for _ in range(22):
    mid = 0.5*(lo + hi); fm = ns55(mid)
    if fm > target: lo = mid
    else: hi = mid
dc = 0.5*(lo + hi)
tr, imax, obs = run(dc, d0)
print("tuned dc = c - c_star = %.6e (relative %.3e), d = d0 = %.7f ; max eps_H = %.4f ; e-folds from Theta=0.03 to max-eps point: %.1f ; ends? %s" % (dc, dc/si.c_star, d0, tr["e1"][imax], tr["N"][imax], tr["ended"]))
lnV, G, GT = si.tables(si.c_star + dc, d0); GTT = np.gradient(GT, si.TG); i0 = int(np.argmin(np.abs(si.TG)))
jm = int(np.argmin(np.abs(G[i0-3000:i0+3000]))) + i0 - 3000
j = int(np.argmin(np.abs(GT[i0-3000:i0+3000]))) + i0 - 3000
lam1, gam = G[j], GTT[j]; a = math.sqrt(lam1*gam/2)
print("inflection at Theta=%.5f: lam1 = V'/V = %.5e, gamma = %.4f, a = %.5f, N_tot = pi/a = %.1f ; analytic n_s(55) = %.5f" % (si.TG[j], lam1, gam, a, math.pi/a, 1 - 4*a/math.tan(55*a)))
res = []
for o in obs:
    print(" N*=%d: n_s=%.5f  r=%.3e  alpha_s=%+.3e  Theta*=%+.5f phi_b*=%+.5f f*=%.5f  eps1=%.3e  H/M_Pl=%.3e  t/(M5 L0)^3=%.4e" % (o["Nstar"], o["n_s"], o["r"], o["alpha_s"], o["Theta"], o["phi_b"], o["f"], o["e1"], o["H_over_MPl"], o["t_over_M5L0_cubed"]))
    cal = [si.calibrate(o, t) for t in (1e-2, 1e-3, 1e-4, 1e-6)]
    for q in cal: print("      t=%.0e: M5 L0=%.4g  M5=%.4g GeV  L0=%.4g GeV^-1 = %.4g m  H=%.4g GeV  k_-=%.4g GeV  H/k_-=%.4g" % (q["t"], q["M5L0"], q["M5_GeV"], q["L0_invGeV"], q["L0_m"], q["H_GeV"], q["k_minus_GeV"], q["H_over_kminus"]))
    res.append(dict(obs=o, calibration=cal))
out["tuned"] = dict(dc=dc, d=d0, lam1=float(lam1), gamma=float(gam), results=res, inflation_ends=bool(tr["ended"]), max_epsH=float(tr["e1"][imax]))
# tuning widths
print("\ntuning sensitivity at N*=55 (Planck 1 sigma = 0.0042):")
rows = []
for fdc in (0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4):
    v = ns55(dc*fdc); rows.append(("dc", fdc, v)); print("   dc x %.2f (d=d0): n_s = %.5f" % (fdc, v), flush=True)
for dd in (-3e-3, -1e-3, -3e-4, 3e-4, 1e-3, 3e-3):
    v = ns55(dc, d0 + dd); rows.append(("d", dd, v)); print("   d = d0 %+.0e (dc fixed): n_s = %s" % (dd, v), flush=True)
out["sensitivity"] = rows
json.dump(out, open(os.path.join(HERE, "SHELL_INFLATION_B.json"), "w"), indent=1)
