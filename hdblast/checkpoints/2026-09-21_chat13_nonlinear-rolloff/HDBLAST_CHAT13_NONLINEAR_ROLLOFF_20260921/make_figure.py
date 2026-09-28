#!/usr/bin/env python3
"""Figure: the two fates of the shell (brane-frame Hubble rate and shell position versus brane proper time), t_det = 0.1, dz = 1e-3."""
import numpy as np, json
Dp = np.load("runs/t01_plus_v2_timeseries.npz")["rec"]; Dm = np.load("runs/t01_minus_ext_timeseries.npz")["rec"]
SURF, INK, INK2, GRID, BLUE, ORANGE = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3de", "#2a78d6", "#eb6834"
Wd, Ht = 1040, 470
out = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" font-family="-apple-system,Helvetica,Arial,sans-serif">' % (Wd, Ht), '<rect width="%d" height="%d" fill="%s"/>' % (Wd, Ht, SURF)]
T = lambda x, y, s, size=12, fill=INK2, anchor="start", weight="400": out.append('<text x="%.1f" y="%.1f" font-size="%d" fill="%s" text-anchor="%s" font-weight="%s">%s</text>' % (x, y, size, fill, anchor, weight, s))
def panel(x0, y0, w, h, xl, yl):
    return (lambda x: x0 + (x - xl[0])/(xl[1] - xl[0])*w), (lambda v: y0 + h - (v - yl[0])/(yl[1] - yl[0])*h)
def path(sx, sy, X, Y, xl, yl):
    pts = [(sx(a), sy(b)) for a, b in zip(X, Y) if xl[0] <= a <= xl[1] and yl[0] <= b <= yl[1]]
    return " ".join("%s%.1f,%.1f" % ("M" if i == 0 else "L", a, b) for i, (a, b) in enumerate(pts))
# panel A: shell position phi_b vs H0 tau
x0, y0, w, h = 70, 70, 420, 320; xl, yl = (0, 7.5), (-2.2, 1.2)
sx, sy = panel(x0, y0, w, h, xl, yl)
T(x0, 30, "A. Where the shell goes", 16, INK, weight="600"); T(x0, 50, "Shell position φ_b versus brane proper time H₀τ  (5D non-linear, t = 0.1)", 12)
for gy in (-2, -1, 0, 1):
    out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (x0, sy(gy), x0 + w, sy(gy), GRID)); T(x0 - 8, sy(gy) + 4, "%+d" % gy if gy else "0", 11, INK2, "end")
for gx in range(0, 8): T(sx(gx), y0 + h + 18, str(gx), 11, INK2, "middle")
T(x0 + w/2, y0 + h + 40, "H₀τ  (brane proper time, Hubble units of the static shell)", 12, INK2, "middle")
out.append('<path d="%s" fill="none" stroke="%s" stroke-width="2"><title>seed toward φ = +1: φ_b = 0.991 when the chart freezes</title></path>' % (path(sx, sy, Dp[:, 9], Dp[:, 1], xl, yl), BLUE))
out.append('<path d="%s" fill="none" stroke="%s" stroke-width="2"><title>seed toward the throat: crosses φ = −1 at H₀τ = 5.70, turnaround at 5.99, collapse</title></path>' % (path(sx, sy, Dm[:, 9], Dm[:, 1], xl, yl), ORANGE))
out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-dasharray="4 4"/>' % (sx(0), sy(-1), sx(7.5), sy(-1), INK2)); T(sx(0.1), sy(-1) - 6, "φ = −1 (AdS₋ throat value)", 11, INK2)
T(sx(4.6), sy(0.75), "→ toward a new de Sitter brane, φ_b → 0.99", 12, INK, weight="600")
T(sx(6.45), sy(-1.6), "collapse", 12, INK, "end", "600"); T(sx(6.45), sy(-1.85), "singularity ≈ H₀τ 6.43 (extrapolated)", 11, INK2, "end")
# panel B: H_J / H0
x0, y0, w, h = 590, 70, 410, 320; xl, yl = (0, 7.5), (-3.0, 1.3)
sx, sy = panel(x0, y0, w, h, xl, yl)
T(x0, 30, "B. What the brane universe sees", 16, INK, weight="600"); T(x0, 50, "Brane-frame Hubble rate H_J / H₀ versus H₀τ", 12)
for gy in (-3, -2, -1, 0, 1):
    out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (x0, sy(gy), x0 + w, sy(gy), GRID)); T(x0 - 8, sy(gy) + 4, "%+d" % gy if gy else "0", 11, INK2, "end")
for gx in range(0, 8): T(sx(gx), y0 + h + 18, str(gx), 11, INK2, "middle")
T(x0 + w/2, y0 + h + 40, "H₀τ", 12, INK2, "middle")
out.append('<path d="%s" fill="none" stroke="%s" stroke-width="2"><title>plus side: H_J/H₀ = 0.649 at H₀τ = 7.3, still relaxing toward the exact RS value 0.638</title></path>' % (path(sx, sy, Dp[:, 9], Dp[:, 2], xl, yl), BLUE))
out.append('<path d="%s" fill="none" stroke="%s" stroke-width="2"><title>throat side: expansion reverses at H₀τ = 5.99; collapse; singularity near 6.43 by extrapolation</title></path>' % (path(sx, sy, Dm[:, 9], Dm[:, 2], xl, yl), ORANGE))
T(sx(4.0), sy(0.52), "still inflating: H_J/H₀ ≈ 0.65, relaxing toward 0.638", 12, INK, weight="600"); T(sx(5.85), sy(-2.2), "collapse", 12, INK, "end", "600")
lx, ly = x0 + 12, y0 + 250
out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2"/>' % (lx, ly, lx + 28, ly, BLUE)); T(lx + 36, ly + 4, "seed toward φ = +1", 12)
out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2"/>' % (lx, ly + 20, lx + 28, ly + 20, ORANGE)); T(lx + 36, ly + 24, "seed toward the throat", 12)
T(70, Ht - 14, "Floating-point 1+1 numerical relativity in the registered 5D model (κ₅ = 1), detuning t = 0.1 (registered value 0.001). Not observations; not a certificate.", 11, INK2)
out.append("</svg>"); open("TWO_FATES.svg", "w").write("\n".join(out)); print("wrote TWO_FATES.svg")
