#!/usr/bin/env python3
"""FLOAT evaluation (guidance for the interval certificate, NOT a proof) of the criterion of Theorem 2 for the d = 8/5 shell:
at mu*^2 in {0, 1, 2, 2.2, 9/4}: min over the grid of psi/y^alpha-normalised psi (positivity), F, lam/beta, G = F - lam/beta, m/psi_b."""
import json, numpy as np
import osc_validate as ov, sys
sys.path.insert(0, ov.BG); import hdblast_background as bg
sh = [s for s in json.load(open("shells.json")) if abs(s["d"] - 1.6) < 1e-12][0]
mu2 = np.array([-4.0, 0.0, 1.0, 2.0, 2.2, 2.25]); r = ov.shoot(mu2, sh["phi_h"], sh["y_b"])
B, beta = sh["B"], sh["beta"]; lam = mu2 + 4; G = r["F"] - lam/beta
m = B*r["chi"] + 3*lam*r["psi"]/(r["rho"]**2*r["s"])
rows = [dict(mu2=float(mu2[k]), zeros=int(r["Z"][k]), psi_b_over_cone_norm=float(r["psi"][k]), F=float(r["F"][k]), lam_over_beta=float(lam[k]/beta),
             G=float(G[k]), G_over_F=float(G[k]/r["F"][k]), m_over_psi_b=float(m[k]/r["psi"][k]), criterion_no_state_below=bool(r["Z"][k] == 0 and G[k] > 0)) for k in range(len(mu2))]
for x in rows: print(x)
json.dump(dict(shell=sh, rows=rows), open("osc_criterion_S85_output.json", "w"), indent=1)
