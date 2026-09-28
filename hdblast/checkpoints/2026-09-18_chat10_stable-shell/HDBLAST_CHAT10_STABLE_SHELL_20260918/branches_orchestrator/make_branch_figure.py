#!/usr/bin/env python3
"""Figure: static-shell branches (4D effective theory curves, 5D points).  Hand-written SVG (no matplotlib)."""
import json, numpy as np
P = "/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/effective_theory_orchestrator/eft_landscape.npz"
D = np.load(P); y, phi, I = D["y"], D["phi"], D["I"]
m = (y < 7) & (y > -9); phi, I = phi[m][::10], I[m][::10]
W = 1 - phi + phi**3/3; W1 = phi**2 - 1
f = 2*I; f1 = 2*(-1 + 2*W*I/3)/W1
c = 2/1.0357712571566782 - 4/3
# for each phi_b != 0 the curvature d that makes it stationary:  (c + d p) = g (1 + c p + d p^2/2),  g = 2 f'/f
g = 2*f1/f
with np.errstate(divide='ignore', invalid='ignore'):
    d_of_phi = (g*(1 + c*phi) - c)/(phi - g*phi**2/2)
sel = (np.abs(phi) > 2e-3) & (phi > -0.32) & (phi < 0.46) & np.isfinite(d_of_phi) & (d_of_phi > 0.8) & (d_of_phi < 2.08)
bp, bd = phi[sel], d_of_phi[sel]; o = np.argsort(bd); bp, bd = bp[o], bd[o]
d0 = 1.1134966
pts = json.load(open("BRANCHES_5D_CHECK.json"))
SURF, INK, INK2, GRID, BLUE, ORANGE = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3de", "#2a78d6", "#eb6834"
Wd, Ht = 760, 480; x0, y0, w, h = 80, 70, 620, 330; xl, yl = (0.85, 2.08), (-0.30, 0.46)
sx = lambda x: x0 + (x - xl[0])/(xl[1] - xl[0])*w; sy = lambda v: y0 + h - (v - yl[0])/(yl[1] - yl[0])*h
out = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" font-family="-apple-system,Helvetica,Arial,sans-serif">' % (Wd, Ht), '<rect width="%d" height="%d" fill="%s"/>' % (Wd, Ht, SURF)]
T = lambda x, y_, s, size=12, fill=INK2, anchor="start", weight="400": out.append('<text x="%.1f" y="%.1f" font-size="%d" fill="%s" text-anchor="%s" font-weight="%s">%s</text>' % (x, y_, size, fill, anchor, weight, s))
T(x0, 30, "Stable and unstable shells are two branches that swap at d₀ = 1.1135", 16, INK, weight="600")
T(x0, 50, "Shell position φ_b versus tension curvature d  (lines: 4D effective theory; dots: full 5D solutions, t = 0.001)", 12)
for gy in (-0.2, 0.0, 0.2, 0.4):
    out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (x0, sy(gy), x0 + w, sy(gy), GRID)); T(x0 - 8, sy(gy) + 4, "%+.1f" % gy if gy else "0", 11, INK2, "end")
for gx in (1.0, 1.2, 1.4, 1.6, 1.8, 2.0): T(sx(gx), y0 + h + 18, "%.1f" % gx, 11, INK2, "middle")
T(x0 + w/2, y0 + h + 40, "tension curvature d", 12, INK2, "middle"); T(x0 - 52, y0 - 12, "φ_b", 12, INK2)
seg = lambda X, Y: " ".join("%s%.1f,%.1f" % ("M" if i == 0 else "L", sx(a), sy(b)) for i, (a, b) in enumerate(zip(X, Y)))
# stable pieces: phi_b<0 branch for d<d0 ; phi_b=0 for d>d0.  unstable: phi_b=0 for d<d0 ; phi_b>0 branch for d>d0
neg = bp < 0; pos = bp > 0
out.append('<path d="%s" fill="none" stroke="%s" stroke-width="2"><title>stable shell (throat side), d &lt; d0</title></path>' % (seg(bd[neg], bp[neg]), BLUE))
out.append('<path d="M%.1f,%.1f L%.1f,%.1f" fill="none" stroke="%s" stroke-width="2"><title>stable shell at the wall centre, d &gt; d0</title></path>' % (sx(d0), sy(0), sx(xl[1]), sy(0), BLUE))
out.append('<path d="%s" fill="none" stroke="%s" stroke-width="2" stroke-dasharray="6 4"><title>unstable shell (barrier top), d &gt; d0</title></path>' % (seg(bd[pos], bp[pos]), ORANGE))
out.append('<path d="M%.1f,%.1f L%.1f,%.1f" fill="none" stroke="%s" stroke-width="2" stroke-dasharray="6 4"><title>unstable hilltop at the wall centre, d &lt; d0 (registered shell is d = 0)</title></path>' % (sx(xl[0]), sy(0), sx(d0), sy(0), ORANGE))
for q in pts:
    stable = (len(q["mu2_5D"]) == 0) or (q["mu2_5D"][0] > 0)
    lab = "5D: d = %.2f, φ_b = %+.5f, m²/H² = %s" % (q["d"], q["phi_b_5D"], ("%+.4f" % q["mu2_5D"][0]) if q["mu2_5D"] else "no bound state")
    out.append('<circle cx="%.1f" cy="%.1f" r="5" fill="%s" stroke="%s" stroke-width="2"><title>%s</title></circle>' % (sx(q["d"]), sy(q["phi_b_5D"]), BLUE if stable else ORANGE, SURF, lab))
T(sx(1.62), sy(0) - 10, "stable shell S₈⁄₅ (no scalar bound state)", 12, INK, weight="600")
T(sx(1.62), sy(0.293) - 12, "barrier top, m² = −2.54 H²", 12, INK, "end")
T(sx(0.93), sy(-0.2477) + 22, "stable, m² = +1.71 H²", 12, INK)
lx, ly = x0 + 16, y0 + 18
out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2"/>' % (lx, ly, lx + 28, ly, BLUE)); T(lx + 36, ly + 4, "stable (m² &gt; 0 or no bound state)", 12)
out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2" stroke-dasharray="6 4"/>' % (lx, ly + 20, lx + 28, ly + 20, ORANGE)); T(lx + 36, ly + 24, "unstable (one growing mode)", 12)
T(x0, Ht - 12, "Floating-point calculations, not an interval certificate. The registered shell is d = 0 (off-scale, on the unstable branch).", 11, INK2)
out.append("</svg>"); open("SHELL_BRANCHES.svg", "w").write("\n".join(out)); print("wrote SHELL_BRANCHES.svg", len(bd))
