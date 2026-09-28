#!/usr/bin/env python3
"""Build figures/SHELL_HILLTOP_AND_SPECTRUM.svg (no matplotlib needed: hand-written SVG).
Panel A: Einstein-frame modulus potential V_E/t versus canonical field Theta (registered linear detuning).
Panel B: 5D tachyon mass mu^2(t) versus detuning t, with the Hamilton-Jacobi closed-form limit.
Colours: slots 1-2 of the documented validated default categorical palette (blue #2a78d6, orange #eb6834);
text uses ink tokens, never series colours."""
import json, os
import numpy as np
here = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(here, "..", "effective_theory_orchestrator", "eft_landscape.npz"))
ys, Th, V = d["y"], d["Theta"], d["V"]
m = (ys <= 20) & (ys >= -26); Th, V = Th[m][::400], V[m][::400]
scan_path = os.path.join(here, "..", "stability", "orchestrator_derivation", "T_SCAN_ORCHESTRATOR.json")
mu0 = -7.719795918397119
pts = [(1e-2, -7.700561458338), (3e-3, -7.714023601543), (1e-3, -7.717871618737), (3e-4, -7.719218614232), (1e-4, -7.719603430014)]
if os.path.exists(scan_path):
    try:
        S = json.load(open(scan_path)); mu0 = S["hj_closed_form_limit"]
        pts = [(r["t"], r["mu2_n12000"][0]) for r in S["scan"] if r["mu2_n12000"]]
    except Exception: pass
slope = np.polyfit([p[0] for p in pts], [p[1] - mu0 for p in pts], 1)[0]

SURF, INK, INK2, GRID, BLUE, ORANGE = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3de", "#2a78d6", "#eb6834"
Wd, Ht = 1040, 470
def panel(x0, y0, w, h, xlim, ylim):
    sx = lambda x: x0 + (x - xlim[0])/(xlim[1] - xlim[0])*w
    sy = lambda y: y0 + h - (y - ylim[0])/(ylim[1] - ylim[0])*h
    return sx, sy
out = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" font-family="-apple-system,Helvetica,Arial,sans-serif">' % (Wd, Ht),
       '<rect width="%d" height="%d" fill="%s"/>' % (Wd, Ht, SURF)]
T = lambda x, y, s, size=13, fill=INK2, anchor="start", weight="400": out.append(
    '<text x="%.1f" y="%.1f" font-size="%d" fill="%s" text-anchor="%s" font-weight="%s">%s</text>' % (x, y, size, fill, anchor, weight, s))

# ---------------- Panel A
x0, y0, w, h = 70, 70, 420, 320; xlim, ylim = (-3.4, 0.7), (0.0, 0.26)
sx, sy = panel(x0, y0, w, h, xlim, ylim)
T(x0, 30, "A. The registered shell sits on a hilltop", 16, INK, weight="600")
T(x0, 50, "Modulus potential V_E / t versus canonical field Θ (Planck units)", 12)
for gy in (0.0, 0.05, 0.10, 0.15, 0.20, 0.25):
    out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1"/>' % (x0, sy(gy), x0 + w, sy(gy), GRID))
    T(x0 - 8, sy(gy) + 4, "%.2f" % gy, 11, INK2, "end")
for gx in (-3, -2, -1, 0):
    T(sx(gx), y0 + h + 18, "%d" % gx, 11, INK2, "middle")
T(x0, y0 + h + 40, "← shell moves toward φ = +1", 12, INK2)
T(x0 + w/2 + 40, y0 + h + 40, "Θ", 12, INK2, "middle")
T(x0 + w + 30, y0 + h + 40, "toward the AdS₋ throat →", 12, INK2, "end")
path = " ".join("%s%.1f,%.1f" % ("M" if i == 0 else "L", sx(a), sy(b)) for i, (a, b) in enumerate(zip(Th, V)) if xlim[0] <= a <= xlim[1])
out.append('<path d="%s" fill="none" stroke="%s" stroke-width="2" stroke-linejoin="round"><title>V_E/t along the registered wall</title></path>' % (path, BLUE))
out.append('<circle cx="%.1f" cy="%.1f" r="6" fill="%s" stroke="%s" stroke-width="2"><title>registered shell: Θ = 0, V_E/t = 0.23303, m²/H² = −7.72</title></circle>' % (sx(0), sy(0.233030), ORANGE, SURF))
T(sx(0) - 12, sy(0.233030) - 12, "registered shell:  m² = −7.72 H²", 12, INK, "end", "600")
out.append('<circle cx="%.1f" cy="%.1f" r="4" fill="%s" stroke="%s" stroke-width="2"><title>φ_b = −1 at Θ = 0.548, V_E/t = 0.1242; the effective theory continues beyond (V_E turns negative); 5D fate unknown</title></circle>' % (sx(0.5479), sy(0.124199), BLUE, SURF))
T(sx(0.5479) + 4, sy(0.124199) + 22, "φ_b = −1 (throat);", 11, INK2, "middle"); T(sx(0.5479) + 4, sy(0.124199) + 36, "continues, 5D fate unknown", 11, INK2, "middle")
T(sx(-3.3), sy(0.0198) - 10, "ends in a stable dS brane (0.0197)", 11, INK2)

# ---------------- Panel B
x0, y0, w, h = 590, 70, 410, 320; xlim, ylim = (0.0, 0.0108), (-7.7225, -7.6985)
sx, sy = panel(x0, y0, w, h, xlim, ylim)
T(x0, 30, "B. Five-dimensional spectrum meets the closed form", 16, INK, weight="600")
T(x0, 50, "Tachyon mass m² / H² versus tension detuning t", 12)
for gy in (-7.720, -7.715, -7.710, -7.705, -7.700):
    out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1"/>' % (x0, sy(gy), x0 + w, sy(gy), GRID))
    T(x0 - 8, sy(gy) + 4, "%.3f" % gy, 11, INK2, "end")
for gx in (0.0, 0.002, 0.004, 0.006, 0.008, 0.010):
    T(sx(gx), y0 + h + 18, ("%.3f" % gx) if gx else "0", 11, INK2, "middle")
T(x0 + w/2, y0 + h + 40, "detuning t   (registered value 0.001)", 12, INK2, "middle")
out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2" stroke-dasharray="5 4"><title>linear fit: μ² = %.6f + %.4f t</title></line>' % (
    sx(0), sy(mu0), sx(0.0105), sy(mu0 + slope*0.0105), INK2, mu0, slope))
for t, v in pts:
    out.append('<circle cx="%.1f" cy="%.1f" r="5" fill="%s" stroke="%s" stroke-width="2"><title>5D shooting: t = %g, μ² = %.7f</title></circle>' % (sx(t), sy(v), BLUE, SURF, t, v))
out.append('<rect x="%.1f" y="%.1f" width="11" height="11" transform="rotate(45 %.1f %.1f)" fill="%s" stroke="%s" stroke-width="2"><title>Hamilton–Jacobi closed form at t → 0: μ² = %.9f</title></rect>' % (
    sx(0) - 5.5, sy(mu0) - 5.5, sx(0), sy(mu0), ORANGE, SURF, mu0))
T(sx(0) + 14, sy(mu0) + 16, "closed form  −4(3c²−4c+8)/(c(3c+4)) = %.6f" % mu0, 12, INK, weight="600")
T(sx(0.0062), sy(mu0 + slope*0.0062) - 12, "slope ≈ +1.92  (small-t value 1.9243)", 12, INK2, "end")
# legend
lx, ly = x0 + 12, y0 + 14
out.append('<circle cx="%.1f" cy="%.1f" r="5" fill="%s"/>' % (lx, ly, BLUE)); T(lx + 12, ly + 4, "5D linear perturbation theory (shooting)", 12, INK2)
out.append('<rect x="%.1f" y="%.1f" width="9" height="9" transform="rotate(45 %.1f %.1f)" fill="%s"/>' % (lx - 4.5, ly + 15.5, lx, ly + 20, ORANGE)); T(lx + 12, ly + 24, "4D Hamilton–Jacobi effective theory (t → 0)", 12, INK2)
T(70, Ht - 14, "Floating-point calculations in the registered dimensionless model (κ₅ = 1). Not observations; not an interval certificate.", 11, INK2)
out.append("</svg>")
open(os.path.join(here, "SHELL_HILLTOP_AND_SPECTRUM.svg"), "w").write("\n".join(out))
print("wrote SVG; slope =", slope, " mu0 =", mu0, " n points =", len(pts))
