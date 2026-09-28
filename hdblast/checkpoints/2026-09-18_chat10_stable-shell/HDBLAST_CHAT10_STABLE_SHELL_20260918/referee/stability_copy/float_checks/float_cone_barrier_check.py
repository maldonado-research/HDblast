#!/usr/bin/env python3
"""FLOAT cross-check (NOT a proof) of the certified cone-barrier coefficients k (r = y R = r_+ + k y^2 + ...):
independent RK4 integration in t = ln y of the background (phi, s, beta = H - 1/y) and of r from y = 1e-5 (r = r_+ there)
to y00 = 1/320; prints (r(y00) - r_+)/y00^2 next to the certified bracket [k_lo, k_hi]."""
import json, math, os
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
cert = json.load(open(os.path.join(HERE, "..", "STABILITY_CERTIFICATE_S85.json")))
A = 8.784370018889562540e-7       # alpha = phi_h + 1 (Task A enclosure midpoint); the state variable is v = phi + 1 so that phi^2 - 1 = v (v - 2) has no cancellation
p = A
U = lambda v: 0.5 * (v * (v - 2)) ** 2 - (2 / 3) * (1 - (v - 1) + (v - 1) ** 3 / 3) ** 2
U1 = lambda v: (v * (v - 2)) * (-4 / 3 + 10 * (v - 1) / 3 - 4 * (v - 1) ** 3 / 9)
def rhs(t, v, lam):
    y = math.exp(t); phi, s, be, r = v
    H = 1 / y + be
    y2w = (1 + y * be) ** 2 - y * y * s * s / 12 + y * y * U(phi) / 6
    c1 = 1 + 2 * y * U1(phi) / s - 6 * y * H
    c0 = lam * y2w - (2 / 3) * y * y * s * s
    return [y * s, y * (U1(phi) - 4 * H * s), y * (-2 * be / y - be * be - s * s / 4 - U(phi) / 6), r * r + c1 * r + c0]
for key, bar in cert["cone_barrier"]["per_mu2"].items():
    mu2 = float(Fr(key)); lam = mu2 + 4; nu = math.sqrt(max(2.25 - mu2, 0.0)); rp = -2.5 - nu
    yi = 1e-5; y00 = 1 / 320; n = 6000; h = math.log(y00 / yi) / n; t = math.log(yi)
    v = [p + U1(p) * yi * yi / 10, U1(p) * yi / 5, -U(p) * yi / 18, rp]
    for _ in range(n):
        k1 = rhs(t, v, lam); k2 = rhs(t + h / 2, [a + h / 2 * b for a, b in zip(v, k1)], lam)
        k3 = rhs(t + h / 2, [a + h / 2 * b for a, b in zip(v, k2)], lam); k4 = rhs(t + h, [a + h * b for a, b in zip(v, k3)], lam)
        v = [a + h / 6 * (b + 2 * c + 2 * d + e) for a, b, c, d, e in zip(v, k1, k2, k3, k4)]; t += h
    kf = (v[3] - rp) / y00 ** 2
    print("mu2=%-5s float k = %.5f   certified bracket [%.5f, %.5f]   inside: %s" % (key, kf, float(Fr(bar["k_lo"])), float(Fr(bar["k_hi"])), float(Fr(bar["k_lo"])) < kf < float(Fr(bar["k_hi"]))))
