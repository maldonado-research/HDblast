#!/usr/bin/env python3
"""Quantum metastability of the classically stable shells (4D effective theory, t -> 0, M_Pl = 1 in the Einstein frame).
Barrier top = the unstable branch at phi_b > 0.  Hawking-Moss exponent  B_HM = 24 pi^2 (1/V_min - 1/V_top) = (b_HM)/t.
Coleman-De Luccia bounces exist only if |V''|/H^2 > 4 at the top, i.e. mu2_top < -4."""
import json, math
rows = json.load(open("EFT_BRANCHES.json"))["rows"]
out = []
print("   d     phi_top    mu2_top   (V_top-V_min)/V_min    b_HM = t*B_HM     CDL possible?   V_E(+1)/t = (1+c+d/2)/81")
c = 2/1.0357712571566782 - 4/3
for r in rows:
    if r["d"] <= 1.1135 or len(r["stationary"]) < 2: continue
    mn = min(r["stationary"], key=lambda q: abs(q["phi_b"])); top = max(r["stationary"], key=lambda q: q["phi_b"])
    b = 24*math.pi**2*(1/mn["VE_over_t"] - 1/top["VE_over_t"])
    print(" %.2f   %+.5f   %+.4f     %.5f              %8.3f          %s            %.5f" % (
        r["d"], top["phi_b"], top["mu2"], top["VE_over_t"]/mn["VE_over_t"] - 1, b, "yes" if top["mu2"] < -4 else "no (HM only)", (1 + c + r["d"]/2)/81))
    out.append(dict(d=r["d"], phi_top=top["phi_b"], mu2_top=top["mu2"], barrier_fraction=top["VE_over_t"]/mn["VE_over_t"] - 1, b_HM=b, CDL=top["mu2"] < -4))
json.dump(out, open("METASTABILITY_EFT.json", "w"), indent=1)
