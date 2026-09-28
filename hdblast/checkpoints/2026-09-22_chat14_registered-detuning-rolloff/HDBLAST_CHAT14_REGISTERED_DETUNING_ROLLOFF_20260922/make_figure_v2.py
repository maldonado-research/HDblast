#!/usr/bin/env python3
"""Figure: the two fates at the REGISTERED detuning t = 0.001 (solid), with the Chat 13 curves at t = 0.1 overlaid (faint) for comparison."""
import sys, numpy as np
P150 = "/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 150/HDBLAST_CHAT13_NONLINEAR_ROLLOFF_20260921/runs/"
plus_tag, minus_tag = (sys.argv[1], sys.argv[2]) if len(sys.argv) > 2 else ("runs/v2_t1e3_plus_c4", "runs/v2_t1e3_minus_c4")
Dp = np.load(plus_tag + "_timeseries.npz")["rec"]; Dm = np.load(minus_tag + "_timeseries.npz")["rec"]
T_CLEAN = float(sys.argv[3]) if len(sys.argv) > 3 else 8.0   # +1 side: plot only the window where the bulk constraint stays small (see 00_READ_FIRST.md)
Dp = Dp[Dp[:, 0] <= T_CLEAN]
Rp = np.load(P150 + "t01_plus_v2_timeseries.npz")["rec"]; Rm = np.load(P150 + "t01_minus_ext_timeseries.npz")["rec"]
SURF, INK, INK2, GRID, BLUE, ORANGE, BLUE_L, ORANGE_L = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3de", "#2a78d6", "#eb6834", "#b7d3f6", "#f6c9b7"
Wd, Ht = 1040, 470
out = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" font-family="-apple-system,Helvetica,Arial,sans-serif">' % (Wd, Ht), '<rect width="%d" height="%d" fill="%s"/>' % (Wd, Ht, SURF)]
T = lambda x, y, s, size=12, fill=INK2, anchor="start", weight="400": out.append('<text x="%.1f" y="%.1f" font-size="%d" fill="%s" text-anchor="%s" font-weight="%s">%s</text>' % (x, y, size, fill, anchor, weight, s))
def panel(x0, y0, w, h, xl, yl): return (lambda x: x0 + (x - xl[0])/(xl[1] - xl[0])*w), (lambda v: y0 + h - (v - yl[0])/(yl[1] - yl[0])*h)
def path(sx, sy, X, Y, xl, yl):
    pts = [(sx(a), sy(b)) for a, b in zip(X, Y) if xl[0] <= a <= xl[1] and yl[0] <= b <= yl[1] and np.isfinite(a) and np.isfinite(b)]
    return " ".join("%s%.1f,%.1f" % ("M" if i == 0 else "L", a, b) for i, (a, b) in enumerate(pts))
def draw(x0, y0, w, h, xl, yl, col, title, sub, ylab, yticks, xoff=0.0):
    sx, sy = panel(x0, y0, w, h, xl, yl)
    T(x0, 30, title, 16, INK, weight="600"); T(x0, 50, sub, 12)
    for gy in yticks:
        out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s"/>' % (x0, sy(gy), x0 + w, sy(gy), GRID)); T(x0 - 8, sy(gy) + 4, ("%+g" % gy) if gy else "0", 11, INK2, "end")
    for gx in range(int(xl[0]), int(xl[1]) + 1): T(sx(gx), y0 + h + 18, str(gx), 11, INK2, "middle")
    T(x0 + w/2, y0 + h + 40, "H₀τ  (brane proper time, Hubble units of the static shell)", 12, INK2, "middle")
    # faint t = 0.1 reference, shifted so that the roll-offs line up (seeds differ in amplitude); shift = difference of the times at which phi_b crosses 0.05 / -0.05
    for R, D, cl, cf in ((Rp, Dp, BLUE_L, BLUE), (Rm, Dm, ORANGE_L, ORANGE)):
        i1 = np.argmax(np.abs(R[:, 1]) > 0.05); i2 = np.argmax(np.abs(D[:, 1]) > 0.05); sh = D[i2, 9] - R[i1, 9]
        out.append('<path d="%s" fill="none" stroke="%s" stroke-width="2"><title>t = 0.1 (Chat 13), time-shifted by %.2f to align the roll-off</title></path>' % (path(sx, sy, R[:, 9] + sh, R[:, col], xl, yl), cl, sh))
    out.append('<path d="%s" fill="none" stroke="%s" stroke-width="2.5"><title>t = 0.001, seed toward φ = +1</title></path>' % (path(sx, sy, Dp[:, 9], Dp[:, col], xl, yl), BLUE))
    out.append('<path d="%s" fill="none" stroke="%s" stroke-width="2.5"><title>t = 0.001, seed toward the throat</title></path>' % (path(sx, sy, Dm[:, 9], Dm[:, col], xl, yl), ORANGE))
    return sx, sy
sx, sy = draw(70, 70, 420, 320, (0, 8), (-2.4, 1.2), 1, "A. Where the shell goes — registered detuning", "Shell position φ_b versus H₀τ  (t = 0.001 solid; t = 0.1 faint, time-aligned)", "φ_b", (-2, -1, 0, 1))
out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-dasharray="4 4"/>' % (sx(0), sy(-1), sx(8), sy(-1), INK2)); T(sx(0.1), sy(-1) - 6, "φ = −1 (AdS₋ throat value)", 11, INK2)
sx, sy = draw(590, 70, 410, 320, (0, 8), (-3.0, 1.3), 2, "B. What the brane universe sees", "Brane-frame Hubble rate H_J / H₀ versus H₀τ", "H_J/H₀", (-3, -2, -1, 0, 1))
out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-dasharray="4 4"/>' % (sx(0), sy(0.6067), sx(8), sy(0.6067), INK2)); T(sx(0.15), sy(0.6067) - 6, "exact Randall–Sundrum brane in AdS₊: 0.607", 11, INK2)
T(sx(0.15), sy(0.26), "blue curve: 0.63 at H₀τ = 7.05, still falling slowly toward 0.607;", 10, INK, weight="600"); T(sx(0.15), sy(0.08), "the chart's lapse collapses there (brane time stops advancing)", 10, INK2)
lx, ly = 602, 320
for cl, lab in ((BLUE, "t = 0.001, seed toward φ = +1"), (ORANGE, "t = 0.001, seed toward the throat"), (BLUE_L, "t = 0.1 reference (Chat 13), time-aligned")):
    out.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="2.5"/>' % (lx, ly, lx + 28, ly, cl)); T(lx + 36, ly + 4, lab, 11); ly += 18
T(70, Ht - 14, "Floating-point 1+1 numerical relativity in the registered 5D model (κ₅ = 1) at the registered tension detuning t = 0.001. Not observations; not a certificate.", 11, INK2)
out.append("</svg>"); open("TWO_FATES_REGISTERED.svg", "w").write("\n".join(out)); print("wrote TWO_FATES_REGISTERED.svg")
