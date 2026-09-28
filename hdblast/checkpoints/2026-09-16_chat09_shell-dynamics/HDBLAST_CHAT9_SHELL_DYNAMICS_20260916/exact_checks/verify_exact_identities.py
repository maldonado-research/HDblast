#!/usr/bin/env python3
"""HDBLAST Chat 9 - exact rational-arithmetic checks (python fractions) of the analytic identities used by the orchestrator.
No floating point is used in the PASS/FAIL decisions.  Run: python3 verify_exact_identities.py"""
import json, random
from fractions import Fraction as Fr
random.seed(20260916)
def rnd(lo=-3, hi=3, nz=True):
    while True:
        x = Fr(random.randint(lo*997, hi*997), random.randint(1, 997))
        if not nz or x != 0: return x
results = {}
def record(name, n_fail, n): results[name] = dict(samples=n, failures=n_fail, status="PASS" if n_fail == 0 else "FAIL")

# ---- model polynomials (exact)
W  = lambda p: 1 - p + p**3/3
W1 = lambda p: p*p - 1
W2 = lambda p: 2*p
U  = lambda p: W1(p)**2/2 - Fr(2, 3)*W(p)**2
U_registered = lambda p: -Fr(1,6) + Fr(4,3)*p - Fr(5,3)*p**2 - Fr(4,9)*p**3 + Fr(17,18)*p**4 - Fr(2,27)*p**6
N = 400
# E1: registered sextic equals the superpotential form
f = sum(1 for _ in range(N) if (lambda p: U(p) != U_registered(p))(rnd()))
record("E1_U_equals_superpotential_form", f, N)

# E2: Hamilton-Jacobi order-2 consistency.  Unknown functions at a point: Phi (free), W, W', W'' (free rationals).
#     transport:  W' Phi' - (2/3) W Phi = -1/2 ;  M = 2 W Phi'/W' ;  identity: -2 W Phi'' + (1/3) W M + (1/2) W' M' + 1/2 = 0
f = 0
for _ in range(N):
    w, w1, w2, w3, Ph = rnd(), rnd(), rnd(), rnd(), rnd()
    Ph1 = (-Fr(1,2) + Fr(2,3)*w*Ph)/w1
    Ph2 = (Fr(2,3)*w1*Ph + Fr(2,3)*w*Ph1 - w2*Ph1)/w1            # derivative of the transport equation
    M   = 2*w*Ph1/w1
    M1  = 2*(w1*Ph1 + w*Ph2)/w1 - 2*w*Ph1*w2/w1**2
    if -2*w*Ph2 + Fr(1,3)*w*M + Fr(1,2)*w1*M1 + Fr(1,2) != 0: f += 1
record("E2_HJ_second_order_consistency_identity", f, N)

# E3: HJ order-0: with pi^{mu nu} = (1/2) sqrt(g) W g^{mu nu}, pi_phi = sqrt(g) W':  (1/2)[(4/3)(2W)^2 - 4 W^2] - (1/2) W'^2 + U = 0
f = sum(1 for _ in range(N) if (lambda p: Fr(1,2)*(Fr(4,3)*(2*W(p))**2 - 4*W(p)**2) - Fr(1,2)*W1(p)**2 + U(p) != 0)(rnd()))
record("E3_HJ_zeroth_order_gives_registered_U", f, N)

# E4: modulus-mass formulae.  General point: W, W', W'' and I free; f = 2I, I' = (-1 + 2WI/3)/W', I'' from differentiation,
#     Z = -W f'/W', Z_E = Z/f + (3/2)(f'/f)^2, linear detuning with stationarity V'/V = 2f'/f:
#     mu2 = 3[-(V'/V)^2 - 2 f''/f + 2 (f'/f)^2]/Z_E   must equal   -4 - Delta,
#     Delta = (8/3) f [ f - 2Z - (3/2) f' W''/W' ] / [ f'^2 + (2/3) f Z ].
f = 0
for _ in range(N):
    w, w1, w2, I = rnd(), rnd(), rnd(), rnd()
    I1 = (-1 + Fr(2,3)*w*I)/w1
    I2 = ((Fr(2,3)*w1*I + Fr(2,3)*w*I1)*w1 - (-1 + Fr(2,3)*w*I)*w2)/w1**2
    F, F1, F2 = 2*I, 2*I1, 2*I2
    Z = -w*F1/w1; ZE = Z/F + Fr(3,2)*(F1/F)**2
    den = F1**2 + Fr(2,3)*F*Z
    if ZE == 0 or den == 0: continue
    g = 2*F1/F
    mu2 = 3*(-(g**2) - 2*F2/F + 2*(F1/F)**2)/ZE
    Delta = Fr(8,3)*F*(F - 2*Z - Fr(3,2)*F1*w2/w1)/den
    if mu2 != -4 - Delta: f += 1
record("E4_mu2_equals_minus4_minus_Delta", f, N)

# E5: registered point phi_b = 0 (W=1, W'=-1, W''=0): stationarity gives c = 2/I - 4/3 and mu2 = -4(3c^2-4c+8)/(c(3c+4)), Z = f' = c I
f = 0
for _ in range(N):
    I = Fr(random.randint(1, 1400), 1000)                                   # 0 < I < 1.4 (< 3/2 so Z > 0)
    w, w1, w2 = Fr(1), Fr(-1), Fr(0)
    I1 = (-1 + Fr(2,3)*w*I)/w1; I2 = ((Fr(2,3)*w1*I + Fr(2,3)*w*I1)*w1)/w1**2
    F, F1, F2 = 2*I, 2*I1, 2*I2; c = 2*F1/F
    Z = -w*F1/w1; ZE = Z/F + Fr(3,2)*(F1/F)**2
    mu2 = 3*(-(c**2) - 2*F2/F + 2*(F1/F)**2)/ZE
    ok = (c == 2/I - Fr(4,3)) and (Z == c*I) and (ZE == c/2 + 3*c*c/8) and (mu2 == -4*(3*c*c - 4*c + 8)/(c*(3*c + 4)))
    if not ok: f += 1
record("E5_registered_closed_forms_c_star_Z_ZE_mu2", f, N)

# E6: instanton total derivative.  rho'^2 := 1 + rho^2 (s^2/12 - U/6);  rho'' := rho[ -1/rho^2 - s^2/3 + rho'^2/rho^2 ]
#     6 rho^2 rho'^2 + 2 rho^3 rho''  ==  6 rho^2 - (4/3) rho^4 U
f = 0
for _ in range(N):
    r, s, u = rnd(1, 9), rnd(), rnd()
    rp2 = 1 + r*r*(s*s/12 - u/6); rpp = r*(-1/(r*r) - s*s/3 + rp2/(r*r))
    if 6*r*r*rp2 + 2*r**3*rpp != 6*r*r - Fr(4,3)*r**4*u: f += 1
record("E6_instanton_total_derivative_identity", f, N)

# E7: master-equation bookkeeping: conformal (z) form -> proper-distance (y) form.
#     coefficient of psi in y-form times rho^2 must equal  4 H_z - 4 H g + 6 + mu2  with H = rho', H_z = rho rho'', g = rho phi''/phi' + rho'
f = 0
for _ in range(N):
    r, s, u, spp, mu2 = rnd(1, 9), rnd(), rnd(), rnd(), rnd()
    rp2 = 1 + r*r*(s*s/12 - u/6); rpp = r*(-1/(r*r) - s*s/3 + rp2/(r*r))
    # use rp symbolically through rp2 only where it appears squared; terms linear in rho' cancel pairwise:
    # z-form zeroth-order coefficient: 4 r rpp - 4 rp (r spp/s + rp) + 6 + mu2 ; y-form * r^2: -(4/3) s^2 r^2 - 4 r rp spp/s + 2 + mu2
    # difference = 4 r rpp - 4 rp2 + 4 + (4/3) s^2 r^2  must vanish
    if 4*r*rpp - 4*rp2 + 4 + Fr(4,3)*s*s*r*r != 0: f += 1
record("E7_master_equation_z_to_y_conversion", f, N)

# E8: shell identity M462 (7.1) from the static Friedmann form h = sigma^2/36 + U/6 - sigma'^2/48
f = 0
for _ in range(N):
    p, t, c = rnd(), rnd(), rnd()
    sig = 2*W(p) + t*(1 + c*p); dsig = 2*W1(p) + t*c
    lhs = sig**2/36 + U(p)/6 - dsig**2/48
    rhs = t*(W(p)*(1 + c*p)/9 - c*W1(p)/12) + t*t*((1 + c*p)**2/36 - c*c/48)
    if lhs != rhs: f += 1
record("E8_static_shell_identity_M462_7_1", f, N)

ok = all(v["status"] == "PASS" for v in results.values())
json.dump(dict(all_pass=ok, checks=results), open("EXACT_IDENTITY_CHECKS.json", "w"), indent=1)
for k, v in results.items(): print("%-52s %s  (%d samples, %d failures)" % (k, v["status"], v["samples"], v["failures"]))
print("ALL PASS" if ok else "SOME CHECKS FAILED")
