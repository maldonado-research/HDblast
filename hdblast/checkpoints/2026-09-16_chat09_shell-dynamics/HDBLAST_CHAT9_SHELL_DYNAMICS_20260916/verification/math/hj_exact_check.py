#!/usr/bin/env python3
"""VERIFICATION (math agent) - exact rational-arithmetic checks of the Hamilton-Jacobi algebra (orchestrator R3).

My own derivation of the order-2 HJ equations (see VERIFICATION_REPORT.md section D) gives, for
S_1 = int sqrt(-g)[ W + Phi R + (1/2) M (d phi)^2 ],
  total order-2 constraint = (2/3) W [Phi R - 3 Box Phi + (1/2) M (dphi)^2] - W'(Phi' R - (1/2) M' (dphi)^2 - M Box phi) - R/2 + (dphi)^2/2
  (R):        (2/3) W Phi - W' Phi' - 1/2 = 0
  (Box phi):  -2 W Phi' + W' M = 0
  (dphi^2):   -2 W Phi'' + (1/3) W M + (1/2) W' M' + 1/2 = 0
Checks done here with python Fractions at random rational points, for random rational cubic superpotentials AND the registered W:
  1. the (dphi^2) equation is an identity once (R) and (Box phi) hold [Phi(phi0) is a free rational number = integration constant];
  2. Z = -2M*... : with f = 4 Phi (two copies, f R/2) and -Z/2 = 2*(M/2):  Z = -2M = -W f'/W';
  3. f'' = (2/3) f - (2/3) Z - f' W''/W'   (used in the Delta formula);
  4. mu2 = 3 (ln V_E)''/Z_E  ==  3 (V''/V)/Z_E - 4 - Delta  with Delta = (8/3) f [f - 2Z - (3/2) f' W''/W']/[f'^2 + (2/3) f Z]
     at a stationary point V'/V = 2 f'/f (V'' random);
  5. registered point phi_b = 0: mu2 = -4(3c^2-4c+8)/(c(3c+4)) with c = 2/I - 4/3, as an identity in I (random rational I);
     also Z_E = 3/(2 I^2) - 1/I  (the Brax-van de Bruck-Davis-Rhodes moduli-space metric in M_pl = 1 units, as transcribed by the literature agent)
     and order-0: U = W'^2/2 - (2/3) W^2 is what the code uses.
"""
from fractions import Fraction as Fr
import random
random.seed(9)
rr = lambda a=-3, b=3: Fr(random.randint(a*1000, b*1000), 1000)

def polyval(c, x, d=0):
    # c = [c0, c1, c2, c3]; d-th derivative
    c = list(c)
    for _ in range(d): c = [i*ci for i, ci in enumerate(c)][1:]
    return sum(ci*x**i for i, ci in enumerate(c))

def check(Wc, p, Phi):
    W, W1, W2, W3 = [polyval(Wc, p, d) for d in range(4)]
    assert W1 != 0
    # (R):  W1 Phi' = (2/3) W Phi - 1/2 ; differentiate twice exactly
    P1 = (Fr(2, 3)*W*Phi - Fr(1, 2))/W1
    # W2 P1 + W1 P2 = (2/3)(W1 Phi + W P1)
    P2 = (Fr(2, 3)*(W1*Phi + W*P1) - W2*P1)/W1
    # W3 P1 + 2 W2 P2 + W1 P3 = (2/3)(W2 Phi + 2 W1 P1 + W P2)
    P3 = (Fr(2, 3)*(W2*Phi + 2*W1*P1 + W*P2) - W3*P1 - 2*W2*P2)/W1
    M = 2*W*P1/W1
    M1 = 2*((W1*P1 + W*P2)*W1 - W*P1*W2)/W1**2
    id3 = -2*W*P2 + Fr(1, 3)*W*M + Fr(1, 2)*W1*M1 + Fr(1, 2)
    f, f1, f2 = 4*Phi, 4*P1, 4*P2
    Z = -2*M
    idZ = Z - (-W*f1/W1)
    idf2 = f2 - (Fr(2, 3)*f - Fr(2, 3)*Z - f1*W2/W1)
    # mu2 forms
    V = Fr(1); V1 = 2*f1/f*V; V2 = rr()
    ZE = Z/f + Fr(3, 2)*(f1/f)**2
    lnVE2 = V2/V - (V1/V)**2 - 2*(f2/f - (f1/f)**2)
    mu2_a = 3*lnVE2/ZE
    Delta = Fr(8, 3)*f*(f - 2*Z - Fr(3, 2)*f1*W2/W1)/(f1**2 + Fr(2, 3)*f*Z)
    mu2_b = 3*(V2/V)/ZE - 4 - Delta
    return id3, idZ, idf2, mu2_a - mu2_b

bad = 0
for trial in range(200):
    Wc = [Fr(1), Fr(-1), Fr(0), Fr(1, 3)] if trial % 2 == 0 else [rr(), rr(), rr(), rr()]
    p = rr(-2, 2); Phi = rr(1, 3)
    try:
        res = check(Wc, p, Phi)
    except (AssertionError, ZeroDivisionError):
        continue
    if any(r != 0 for r in res): bad += 1; print("FAIL", trial, res)
print("200 random exact checks of [consistency identity, Z=-W f'/W', f'' relation, mu2 == 3V''/(V Z_E) - 4 - Delta]: failures =", bad)

# registered point phi_b = 0 as identities in I
bad = 0
for trial in range(100):
    I = Fr(random.randint(500, 1490), 1000)       # 0 < I < 3/2 (Z > 0 region)
    W, W1, W2 = Fr(1), Fr(-1), Fr(0)
    Phi = I/2
    P1 = (Fr(2, 3)*W*Phi - Fr(1, 2))/W1; P2 = (Fr(2, 3)*(W1*Phi + W*P1) - W2*P1)/W1
    f, f1, f2 = 4*Phi, 4*P1, 4*P2; Z = -W*f1/W1
    c = 2*f1/f                                     # stationarity for V = t(1 + c phi) at phi = 0
    ZE = Z/f + Fr(3, 2)*(f1/f)**2
    mu2 = 3*(-c*c - 2*(f2/f - (f1/f)**2))/ZE
    ok = (c == 2/I - Fr(4, 3)) and (mu2 == -4*(3*c*c - 4*c + 8)/(c*(3*c + 4))) and (ZE == Fr(3, 2)/I**2 - 1/I) \
         and (ZE == c/2 + 3*c*c/8) and (Z == c*I) and (f2/f == (2 - c)/3)
    lnV_zz_BBDR = -2/I**2 + 4/I - Fr(28, 9)
    ok = ok and (3*lnV_zz_BBDR/(Fr(3, 2)/I**2 - 1/I) == mu2)
    if not ok: bad += 1
print("100 random-I exact checks at phi_b=0 [c=2/I-4/3, closed-form mu2, Z_E=3/(2I^2)-1/I=c/2+3c^2/8, Z=cI, I''/I=(2-c)/3, BBDR form]: failures =", bad)
# discriminant statement
print("discriminant of 3c^2-4c+8 =", 16 - 96, "(<0: numerator never vanishes); pure AdS check: I = 3/(2W) gives Z =",
      -Fr(1)*(4*((Fr(2, 3)*Fr(3, 4)) - Fr(1, 2))/Fr(-1))/Fr(-1))
c = 0.5975949350280132
print("closed form at c_star: %.12f" % (-4*(3*c*c - 4*c + 8)/(c*(3*c + 4))))
