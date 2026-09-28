#!/usr/bin/env python3
"""HDBLAST Chat 10 - shared rigorous machinery for the background (cone -> shell) problem with QUADRATIC detuning
        sigma_t(phi) = 2 W(phi) + t (1 + c phi + d phi^2 / 2),     t, c, d exact rationals.
Python 3.9 standard library only.  Re-uses cert_core.py (byte-identical copy of the Chat 9 arithmetic core:
outward-rounded 2^-384 fixed point, one-parameter Taylor models, validated Taylor integrator with verified a-priori
boxes).  The cone start is the Chat 9 Volterra/disc-algebra contraction argument (certify_tachyon.py, step 1),
transcribed as a function; it does not depend on the detuning.

State vector: (phi, s = phi', H = rho'/rho, w = 1/rho^2),   phi' = s, s' = U1 - 4Hs, H' = -w - s^2/3, w' = -2Hw.
Junction residuals at a point y:   G1 = H - sigma_t(phi)/6,    G2 = s + sigma_t'(phi)/2.
"""
import math
from fractions import Fraction as Fr
import cert_core as cc
from cert_core import IV, TM, ONE, P

T_PAR = Fr(1, 1000)
C_PAR = Fr(5975949350280, 10 ** 13)
D_PAR = Fr(8, 5)

Y0 = Fr(1, 20)               # end of the cone chart
NC = 46                      # cone series degree
KAPPA = Fr(1, 8)             # step <= KAPPA * y
HMAX = Fr(1, 10)
NORD_MAX = 64
TOL = ONE >> 140             # accepted scaled Lagrange remainder per step

GATES = []


def gate(name, ok, detail=""):
    GATES.append({"gate": name, "pass": bool(ok), "detail": detail})
    if not ok:
        raise SystemExit("GATE FAILED: %s  %s" % (name, detail))


# ------------------------------------------------------------------ polynomial helpers (exact)
def poly_shift(c, a):
    out = [Fr(0)] * len(c)
    for k, ck in enumerate(c):
        for j in range(k + 1):
            out[j] += ck * math.comb(k, j) * a ** (k - j)
    return out


def poly_der(c):
    return [k * c[k] for k in range(1, len(c))]


def majorant(c, x):
    return sum(abs(ck) * x ** k for k, ck in enumerate(c))


def pmul(a, b):
    out = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def U_of(x):       # generic (TM or IV)
    x2 = x * x; x3 = x2 * x
    Wv = (x3.scale(1, 3) - x).addc(1); W1v = x2.addc(-1)
    return (W1v * W1v).scale(1, 2) - (Wv * Wv).scale(2, 3)


def U1_of(x):
    x2 = x * x; x3 = x2 * x
    return x2.addc(-1) * (x.scale(10, 3) - x3.scale(4, 9)).addc(Fr(-4, 3))


def sigma_of(x):
    """sigma_t = 2(1 - phi + phi^3/3) + t(1 + c phi + d phi^2/2)"""
    x2 = x * x; x3 = x2 * x
    Wv = (x3.scale(1, 3) - x).addc(1)
    tc = T_PAR * C_PAR; td2 = T_PAR * D_PAR / 2
    return (Wv.scale(2) + x.scale(tc.numerator, tc.denominator) + x2.scale(td2.numerator, td2.denominator)).addc(T_PAR)


def dsigma_of(x):
    """sigma_t' = 2(phi^2 - 1) + t(c + d phi)"""
    td = T_PAR * D_PAR
    return ((x * x).addc(-1).scale(2) + x.scale(td.numerator, td.denominator)).addc(T_PAR * C_PAR)


def residuals(st):
    phi, s, H, w = st[:4]
    return H - sigma_of(phi).scale(1, 6), s + dsigma_of(phi).scale(1, 2)


def B_of(st):
    """B = phi''/phi' + sigma_t''/2 = U1/s - 4H + 2 phi + t d / 2"""
    phi, s, H, w = st[:4]
    return (U1_of(phi) * s.recip() - H.scale(4) + phi.scale(2)).addc(T_PAR * D_PAR / 2)


# ------------------------------------------------------------------ cone start (Chat 9 step 1, as a function)
def cone_start(P_TM, p_lo, p_hi):
    """P_TM: Taylor model of p = phi_h over the parameter range; p_lo, p_hi exact Fractions bounding that range.
    Returns ([phi, s, H, w] at Y0 as TMs, info dict).  All inequalities are exact rational gates."""
    Ucoef = [Fr(-1, 6), Fr(4, 3), Fr(-5, 3), Fr(-4, 9), Fr(17, 18), Fr(0), Fr(-2, 27)]
    Wc = [Fr(1), Fr(-1), Fr(0), Fr(1, 3)]
    W1c = [Fr(-1), Fr(0), Fr(1)]
    Uchk = [Fr(1, 2) * a for a in pmul(W1c, W1c)] + [Fr(0)] * 2
    for i, a in enumerate(pmul(Wc, Wc)):
        Uchk[i] -= Fr(2, 3) * a
    gate("U polynomial equals W_phi^2/2 - (2/3)W^2", Uchk == Ucoef)
    gate("U_phi = (phi^2-1)(-4/3 + 10 phi/3 - 4 phi^3/9)", pmul(W1c, [Fr(-4, 3), Fr(10, 3), Fr(0), Fr(-4, 9)]) == poly_der(Ucoef))
    U2chk = [a + b for a, b in zip(pmul([Fr(0), Fr(2)], [Fr(-4, 3), Fr(10, 3), Fr(0), Fr(-4, 9)]), pmul(W1c, [Fr(10, 3), Fr(0), Fr(-4, 3)]))]
    gate("U_phiphi = 2 phi Bq + (phi^2-1)(10/3 - 4 phi^2/3)  (used in the a-priori G equation)", U2chk == poly_der(poly_der(Ucoef)))
    A0 = poly_shift(Ucoef, Fr(-1)); A1 = poly_der(A0); A2 = poly_der(A1)
    RC = Fr(1, 2); V0 = Fr(1, 10 ** 6)
    r_be, r_s, r_ph = Fr(6, 100), Fr(2, 10 ** 6), Fr(11, 10 ** 7)
    gate("|p+1| <= v0 on the parameter range", abs(p_lo + 1) <= V0 and abs(p_hi + 1) <= V0)
    vs = V0 + r_ph
    selfmap = [RC / 3 * (r_be ** 2 + r_s ** 2 / 4 + majorant(A0, vs) / 6), RC / 5 * (majorant(A1, vs) + 4 * r_be * r_s), RC * r_s]
    gate("cone map sends the ball into itself", selfmap[0] < r_be and selfmap[1] < r_s and selfmap[2] < r_ph,
         str([float(x) for x in selfmap]))
    lip = [RC / 3 * (2 * r_be * r_be + r_s * r_s / 2 + majorant(A1, vs) * r_ph / 6) / r_be,
           RC / 5 * (4 * r_s * r_be + 4 * r_be * r_s + majorant(A2, vs) * r_ph) / r_s,
           RC * r_s / r_ph]
    gate("cone map is a contraction (weighted max norm)", max(lip) < 1, str([float(x) for x in lip]))
    cau = TM.cauchy
    zero = TM([0])
    be = [zero]; s_ = [zero]; ph = [P_TM]
    p2 = []; p3 = []; Ws = []; W1s = []; Bq = []; Us = []; U1s = []; be2 = []; s2 = []; bes = []
    for n in range(1, NC + 1):
        k = n - 1
        p2.append(cau(ph, ph, k)); p3.append(cau(p2, ph, k))
        Wk = p3[k].scale(1, 3) - ph[k]
        Ws.append(Wk.addc(1) if k == 0 else Wk)
        W1s.append(p2[k].addc(-1) if k == 0 else p2[k])
        b = ph[k].scale(10, 3) - p3[k].scale(4, 9)
        Bq.append(b.addc(Fr(-4, 3)) if k == 0 else b)
        Us.append(cau(W1s, W1s, k).scale(1, 2) - cau(Ws, Ws, k).scale(2, 3))
        U1s.append(cau(W1s, Bq, k))
        be2.append(cau(be, be, k)); s2.append(cau(s_, s_, k)); bes.append(cau(be, s_, k))
        f = -be2[k] - s2[k].scale(1, 4) - Us[k].scale(1, 6)
        g = U1s[k] - bes[k].scale(4)
        be.append(f.scale(Y0.numerator, Y0.denominator * (n + 2)))
        s_.append(g.scale(Y0.numerator, Y0.denominator * (n + 4)))
        ph.append(s_[k].scale(Y0.numerator, Y0.denominator * n))
    xq = Y0 / RC

    def tm_sum(lst, r):
        acc = lst[0]
        for a in lst[1:]:
            acc = acc + a
        tl = cc.fceil(r * xq ** (NC + 1) / (1 - xq))
        return TM(acc.c, acc.rlo - tl, acc.rhi + tl)

    be0 = tm_sum(be, r_be); s0 = tm_sum(s_, r_s); ph0 = tm_sum(ph, r_ph)
    H0 = be0.addc(1 / Y0)
    w0 = H0 * H0 - (s0 * s0).scale(1, 12) + U_of(ph0).scale(1, 6)
    gate("s > 0 at y0", s0.bound().lo > 0)
    # sigma = s/y > 0 on the whole cone chart (0, y0]  (Horner over tau = y/y0 in [0,1], tail bound divided by y0)
    tauI = IV(0, ONE); acc = IV(0, 0)
    for a in reversed(s_[1:]):
        acc = acc * tauI + a.bound()
    tls = cc.fceil(r_s * xq ** (NC + 1) / (1 - xq))
    sgI = (acc + IV(-tls, tls)).scale(Y0.denominator, Y0.numerator)
    gate("sigma = s/y > 0 on (0,y0]  (phi strictly increasing on the cone chart)", sgI.lo > 0)
    info = {"y0": str(Y0), "series_degree": NC, "disc_radius": str(RC), "ball_radii_beta_s_phi": [str(r_be), str(r_s), str(r_ph)],
            "self_map_bounds": [str(x) for x in selfmap], "contraction_constant_upper": "%.6f" % float(max(lip))}
    return [ph0, s0, H0, w0], info


# ------------------------------------------------------------------ validated integration
class Integrator:
    def __init__(self, state, y, verbose=False):
        self.state = state; self.y = y; self.N = 24; self.h_prev = Y0 / 64
        self.nstep = 0; self.max_rem = 0; self.verbose = verbose

    def _adaptive(self, state, h, tube=False, allow_halving=True):
        while True:
            try:
                out = cc.step(state, [], h, self.N, tau_interval=tube)
            except ArithmeticError:
                out = None
            if out is not None and out[2] <= TOL:
                if out[2] < (TOL >> 40) and self.N > 12:
                    self.N -= 2
                return out, h
            if out is not None and self.N + 4 <= NORD_MAX:
                self.N += 4
                continue
            if not allow_halving:
                raise SystemExit("fixed-size step failed")
            h = h / 2
            if h < Fr(1, 2 ** 30):
                raise SystemExit("step size underflow")

    def advance_to(self, y_target):
        while self.y < y_target:
            h = min(KAPPA * self.y, HMAX, 2 * self.h_prev)
            h = Fr(int(h * 2 ** 16), 2 ** 16)
            truncated = self.y + h >= y_target
            if truncated:
                h = y_target - self.y
            (self.state, _, rem), h2 = self._adaptive(self.state, h)
            if not (truncated and h2 == h):
                self.h_prev = h2                 # a step cut short by the target must not shrink the next proposals
            h = h2; self.max_rem = max(self.max_rem, rem)
            self.y += h; self.nstep += 1
            if self.verbose and self.nstep % 40 == 0:
                print("   step %3d y=%.5f N=%d phi+1=%.6e remwidth=%.2e" % (self.nstep, float(self.y), self.N, self.state[0].c[0] / ONE + 1,
                      max(x.rhi - x.rlo for x in self.state) / ONE), flush=True)
        return self.state

    def fixed_step(self, state, h, tube=False):
        (new, tb, rem), _ = self._adaptive(state, h, tube=tube, allow_halving=False)
        self.max_rem = max(self.max_rem, rem)
        return new, tb


def eval_on_segment(state, h, N):
    """state: TMs at y_a whose polynomial parts are CONSTANT (c[k] = 0 for k >= 1; all uncertainty in the remainder).
    Returns TM enclosures (in the free parameter u) of the solution at  y = y_a + h (1+u)/2,  u in [-1,1], i.e. on the
    whole segment [y_a, y_a + h], keeping the y-dependence as a polynomial in u (needed for the preconditioned test)."""
    for x in state:
        assert not any(x.c[1:]), "eval_on_segment needs parameter-free input"
    G_tm = U1_of(state[0]) * state[1].recip()
    x0 = [x.bound() for x in state] + [G_tm.bound()]
    Bh = apriori_box_combined(x0, [], h)
    if Bh is None:
        raise SystemExit("no a-priori enclosure on the segment")
    GB = Bh[-1]; Bh = Bh[:-1]
    X = cc.taylor(list(state), [], h, N - 1)
    XB = cc.taylor(Bh, [], h, N, G0=GB)
    tau = TM([ONE // 2, ONE // 2])             # tau = (1+u)/2 exactly
    unit = IV(0, ONE)
    out = []
    for i in range(4):
        acc = X[i][N - 1]
        for k in range(N - 2, -1, -1):
            acc = acc * tau + X[i][k]
        out.append(acc.add_iv(XB[i][N] * unit))
    return out, max(XB[i][N].mag() for i in range(4))


def tm_at(x, sign):
    """enclosure (IV) of the Taylor model x at u = +1 or u = -1"""
    v = 0
    for k, ck in enumerate(x.c):
        v += ck if (sign > 0 or k % 2 == 0) else -ck
    return IV(v + x.rlo, v + x.rhi)


def iv_sqrt(x):
    if x.lo <= 0:
        raise ArithmeticError("sqrt of non-positive interval")
    return IV(math.isqrt(x.lo << P), math.isqrt(x.hi << P) + 1)


def fr_of(n):
    return Fr(n, ONE)


# ------------------------------------------------------------------ a-priori box with a more robust inflation heuristic
def apriori_box_v2(x0, mu2s, h, maxit=60):
    """Same rigorous acceptance test as cert_core.apriori_box  ( x0 + [0,h] f(Bi) subset of Bi  =>  every solution
    starting in x0 exists on [0,h] and stays in the returned box ), but the trial boxes are inflated by
    width/8 + (largest component width)/64.  The Chat 9 heuristic inflates each component only relative to its own
    width and stalls where s' = s (G - 4H) changes sign, which happens right at this shell (phi'' ~ 0 at y_b).
    Only the search heuristic differs; the inclusion that is finally verified is identical."""
    hI = IV(0, cc.fceil(h))
    f0 = cc.vector_field(x0, mu2s)
    B = [x + hI * f for x, f in zip(x0, f0)]
    for it in range(maxit):
        wmax = max(b.width() for b in B)
        Bi = [IV(b.lo - (b.width() // 8 + wmax // 64 + 1), b.hi + (b.width() // 8 + wmax // 64 + 1)) for b in B]
        try:
            f = cc.vector_field(Bi, mu2s)
        except ArithmeticError:
            return None
        Bn = [x + hI * g for x, g in zip(x0, f)]
        if all(bn.subset_of(bi) for bn, bi in zip(Bn, Bi)):
            return Bn
        B = [bn.hull(bi) for bn, bi in zip(Bn, Bi)]
        if any(b.mag() > (ONE << 24) for b in B):
            return None
    return None


_apriori_v1 = cc.apriori_box


def apriori_box_combined(x0, mu2s, h, maxit=60):
    """Chat 9 heuristic first; if it finds no box, the v2 inflation.  Either result has passed the same inclusion test."""
    B = _apriori_v1(x0, mu2s, h)
    return B if B is not None else apriori_box_v2(x0, mu2s, h, maxit)


cc.apriori_box = apriori_box_combined
