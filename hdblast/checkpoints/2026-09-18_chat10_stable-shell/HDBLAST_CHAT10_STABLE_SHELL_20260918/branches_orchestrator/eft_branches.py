#!/usr/bin/env python3
"""Chat 10 (orchestrator): static de Sitter shells of the quadratic-detuning family in the 4D effective theory (t -> 0).
V_E(phi; d) = (1 + c phi + d phi^2/2)/f(phi)^2.  For each d list ALL stationary points phi_b in (-1,1), their mu2 = 3 (ln V_E)''/Z_E,
and V_E.  Shows how the unstable hilltop and the stable shell are connected (bifurcation structure near d0 = 1.1135)."""
import json, numpy as np
P = "/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/effective_theory_orchestrator/eft_landscape.npz"
D = np.load(P); y, phi, I = D["y"], D["phi"], D["I"]
m = (y < 7) & (y > -9); y, phi, I = y[m][::10], phi[m][::10], I[m][::10]
W = 1 - phi + phi**3/3; W1 = phi**2 - 1; W2 = 2*phi
I1 = (-1 + 2*W*I/3)/W1
I2 = ((2*W1*I/3 + 2*W*I1/3)*W1 - (-1 + 2*W*I/3)*W2)/W1**2
f, f1, f2 = 2*I, 2*I1, 2*I2
Z = -W*f1/W1; ZE = Z/f + 1.5*(f1/f)**2
c = 2/1.0357712571566782 - 4/3
rows = []
for d in [0.0, 0.5, 0.9, 1.0, 1.05, 1.08, 1.10, 1.11, 1.1135, 1.12, 1.15, 1.2, 1.3, 1.6, 2.0, 3.0]:
    V = 1 + c*phi + d*phi**2/2; Vp = c + d*phi
    g = Vp/V - 2*f1/f                                   # (ln V_E)'
    ok = V > 0
    idx = np.where((g[:-1]*g[1:] < 0) & ok[:-1] & ok[1:])[0]
    pts = []
    for i in idx:
        w = g[i]/(g[i] - g[i+1]); lin = lambda a: a[i] + w*(a[i+1] - a[i])
        p = lin(phi); lnV2 = d/lin(V) - (lin(Vp)/lin(V))**2 - 2*(lin(f2)/lin(f) - (lin(f1)/lin(f))**2)
        pts.append(dict(phi_b=float(p), mu2=float(3*lnV2/lin(ZE)), VE_over_t=float(lin(V)/lin(f)**2)))
    rows.append(dict(d=d, stationary=pts))
    print("d=%.4f : " % d + " | ".join("phi_b=%+.5f mu2=%+.4f V_E/t=%.5f" % (q["phi_b"], q["mu2"], q["VE_over_t"]) for q in pts))
json.dump(dict(c=c, rows=rows), open("EFT_BRANCHES.json", "w"), indent=1)
