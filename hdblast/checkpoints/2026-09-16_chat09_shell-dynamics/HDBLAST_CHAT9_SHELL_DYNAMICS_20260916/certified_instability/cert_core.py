#!/usr/bin/env python3
"""Rigorous arithmetic core for the HDBLAST tachyon certificate.  Python 3.9 standard library only.

All real numbers are enclosed with outward-rounded fixed-point integers (unit 2^-P).

IV  : interval [lo, hi]  (integers, unit 2^-P).
TM  : Taylor model in ONE parameter u in [-1, 1]:  f(u) in  sum_k c[k] u^k + [rlo, rhi]   for every u in [-1,1]
      (c[k] exact fixed-point numbers, remainder an interval).  The parameter is the horizon value phi_h.
Both types offer the same interface (+, -, scale by an exact rational, add an exact rational, reciprocal,
Cauchy product of coefficient lists), so the Taylor-coefficient recurrences are written once.
"""
from fractions import Fraction

P = 384
ONE = 1 << P
D = 6            # parameter degree of the Taylor models (set by the driver before use)


def set_degree(d):
    global D
    D = d


def ffloor(fr):
    fr = Fraction(fr)
    return (fr.numerator << P) // fr.denominator


def fceil(fr):
    fr = Fraction(fr)
    return -((-fr.numerator << P) // fr.denominator)


def cdiv(a, b):          # ceil(a/b), b > 0
    return -((-a) // b)


class IV:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi):
        if lo > hi:
            raise ArithmeticError("empty interval")
        self.lo = lo
        self.hi = hi

    @staticmethod
    def frac(a, b=None):
        if b is None:
            b = a
        return IV(ffloor(a), fceil(b))

    def __add__(self, o):
        return IV(self.lo + o.lo, self.hi + o.hi)

    def __sub__(self, o):
        return IV(self.lo - o.hi, self.hi - o.lo)

    def __neg__(self):
        return IV(-self.hi, -self.lo)

    def __mul__(self, o):
        a = self.lo * o.lo; b = self.lo * o.hi; c = self.hi * o.lo; d = self.hi * o.hi
        return IV(min(a, b, c, d) >> P, -((-max(a, b, c, d)) >> P))

    def scale(self, num, den=1):
        """multiply by the exact rational num/den (den > 0)"""
        if num >= 0:
            return IV((self.lo * num) // den, cdiv(self.hi * num, den))
        return IV((self.hi * num) // den, cdiv(self.lo * num, den))

    def addc(self, fr):
        return IV(self.lo + ffloor(fr), self.hi + fceil(fr))

    def recip(self):
        if self.lo > 0 or self.hi < 0:
            return IV((ONE * ONE) // self.hi if self.hi > 0 else -cdiv(ONE * ONE, -self.hi),
                      cdiv(ONE * ONE, self.lo) if self.lo > 0 else -((ONE * ONE) // (-self.lo)))
        raise ArithmeticError("division by interval containing zero")

    def bound(self):
        return self

    def hull(self, o):
        return IV(min(self.lo, o.lo), max(self.hi, o.hi))

    def subset_of(self, o):
        return o.lo <= self.lo and self.hi <= o.hi

    def mag(self):
        return max(abs(self.lo), abs(self.hi))

    def width(self):
        return self.hi - self.lo

    @staticmethod
    def cauchy(a, b, k):
        lo = 0; hi = 0
        for j in range(k + 1):
            x = a[j]; y = b[k - j]
            p1 = x.lo * y.lo; p2 = x.lo * y.hi; p3 = x.hi * y.lo; p4 = x.hi * y.hi
            lo += min(p1, p2, p3, p4); hi += max(p1, p2, p3, p4)
        return IV(lo >> P, -((-hi) >> P))

    def fl(self):
        return (self.lo / ONE, self.hi / ONE)


def _imul_raw(alo, ahi, blo, bhi):
    p1 = alo * blo; p2 = alo * bhi; p3 = ahi * blo; p4 = ahi * bhi
    return min(p1, p2, p3, p4), max(p1, p2, p3, p4)


class TM:
    __slots__ = ("c", "rlo", "rhi", "plo", "phi_")

    def __init__(self, c, rlo=0, rhi=0):
        c = list(c) + [0] * (D + 1 - len(c))
        if rlo > rhi:
            raise ArithmeticError("empty remainder")
        self.c = c
        self.rlo = rlo
        self.rhi = rhi
        m = 0
        for x in c[1:]:
            m += abs(x)
        self.plo = c[0] - m         # range of the polynomial part over u in [-1,1]
        self.phi_ = c[0] + m

    @staticmethod
    def const_iv(iv):
        return TM([iv.lo], 0, iv.hi - iv.lo)

    def bound(self):
        return IV(self.plo + self.rlo, self.phi_ + self.rhi)

    def __add__(self, o):
        return TM([x + y for x, y in zip(self.c, o.c)], self.rlo + o.rlo, self.rhi + o.rhi)

    def __sub__(self, o):
        return TM([x - y for x, y in zip(self.c, o.c)], self.rlo - o.rhi, self.rhi - o.rlo)

    def __neg__(self):
        return TM([-x for x in self.c], -self.rhi, -self.rlo)

    def add_iv(self, iv):
        return TM(self.c, self.rlo + iv.lo, self.rhi + iv.hi)

    def scale(self, num, den=1):
        c = [(x * num) // den for x in self.c]           # each floor error in [0,1) ulp, times u^k in [-1,1]
        if num >= 0:
            rlo = (self.rlo * num) // den; rhi = cdiv(self.rhi * num, den)
        else:
            rlo = (self.rhi * num) // den; rhi = cdiv(self.rlo * num, den)
        e = 0 if den == 1 else D + 1
        return TM(c, rlo - e, rhi + e)

    def addc(self, fr):
        c = list(self.c)
        c[0] += ffloor(fr)
        return TM(c, self.rlo, self.rhi + 1)

    @staticmethod
    def cauchy(a, b, k):
        """enclosure of sum_{j<=k} a[j]*b[k-j]"""
        full = [0] * (2 * D + 1)
        rlo = 0; rhi = 0
        for j in range(k + 1):
            A = a[j]; B = b[k - j]
            Bc = B.c
            for i, ai in enumerate(A.c):
                if ai:
                    for l, bl in enumerate(Bc):
                        if bl:
                            full[i + l] += ai * bl
            # (pa + Ia)(pb + Ib) = pa pb + (pa + Ia) Ib + pb Ia
            if B.rlo or B.rhi:
                lo, hi = _imul_raw(A.plo + A.rlo, A.phi_ + A.rhi, B.rlo, B.rhi)
                rlo += lo; rhi += hi
            if A.rlo or A.rhi:
                lo, hi = _imul_raw(B.plo, B.phi_, A.rlo, A.rhi)
                rlo += lo; rhi += hi
        hi_part = 0
        for x in full[D + 1:]:
            hi_part += abs(x)
        rlo -= hi_part; rhi += hi_part
        c = [x >> P for x in full[:D + 1]]
        return TM(c, (rlo >> P) - (D + 1), -((-rhi) >> P) + (D + 1))

    def __mul__(self, o):
        return TM.cauchy([self], [o], 0)

    def recip(self):
        c0 = self.c[0]
        if c0 == 0:
            raise ArithmeticError("TM reciprocal: zero centre")
        a0 = (ONE * ONE) // c0                      # exact number a0/2^P, approx 1/c0
        q = self.scale(a0, ONE).addc(-1)            # 1/self = (a0/2^P) / (1+q)
        th = q.bound().mag()
        if 2 * th >= ONE:
            raise ArithmeticError("TM reciprocal: relative spread too large")
        S = TM([ONE])
        mq = -q
        for _ in range(D):
            S = (mq * S).addc(1)
        t = ONE
        for _ in range(D + 1):
            t = cdiv(t * th, ONE)
        t = cdiv(t * ONE, ONE - th)                 # |q|^(D+1)/(1-|q|)
        S = TM(S.c, S.rlo - t, S.rhi + t)
        return S.scale(a0, ONE)


# ----------------------------------------------------------------------------------------------
# Taylor-coefficient recurrences (scaled by the step: X_k = x_k h^k) for
#   phi' = s,  s' = U1(phi) - 4 H s,  H' = -w - s^2/3,  w' = -2 H w,
#   R_j' = R_j^2 + (2G - 6H) R_j + (mu2_j + 4) w - (2/3) s^2,     G = U1(phi)/s,
#   U1 = (phi^2 - 1)(-4/3 + (10/3) phi - (4/9) phi^3).
# ----------------------------------------------------------------------------------------------
def taylor(state, mu2s, h, N, G0=None):
    """state = [phi, s, H, w, R_0, R_1, ...]; mu2s[j] (Fraction) belongs to R_j.  Returns the list of
    coefficient lists X[i][k], k = 0..N, of the solution in the scaled time tau = (y - y_k)/h.
    G0 (optional): any valid enclosure of G = U1(phi)/s at the expansion point (used for a-priori boxes, where the
    quotient of the two boxes would be needlessly wide)."""
    T = type(state[0])
    cau = T.cauchy
    hn, hd = h.numerator, h.denominator
    nR = len(state) - 4
    phi = [state[0]]; s = [state[1]]; H = [state[2]]; w = [state[3]]
    R = [[state[4 + j]] for j in range(nR)]
    p2 = []; p3 = []; A = []; Bq = []; U1 = []; Hs = []; s2 = []; Hw = []; G = []; C = []
    inv_s0 = state[1].recip()
    mu4 = [Fraction(m) + 4 for m in mu2s]
    for k in range(N):
        p2.append(cau(phi, phi, k))
        p3.append(cau(p2, phi, k))
        A.append(p2[k].addc(-1) if k == 0 else p2[k])
        b = phi[k].scale(10, 3) - p3[k].scale(4, 9)
        Bq.append(b.addc(Fraction(-4, 3)) if k == 0 else b)
        U1.append(cau(A, Bq, k))
        Hs.append(cau(H, s, k))
        s2.append(cau(s, s, k))
        Hw.append(cau(H, w, k))
        # G_k = (U1_k - sum_{j<k} G_j s_{k-j}) / s_0
        if k == 0:
            num = U1[0]
        else:
            num = U1[k] - cau(G + [ZERO(T)], s, k)      # sum_{j<k} G_j s_{k-j}
        G.append(G0 if (k == 0 and G0 is not None) else num * inv_s0)
        C.append(G[k].scale(2) - H[k].scale(6))
        sc_n, sc_d = hn, hd * (k + 1)
        phi.append(s[k].scale(sc_n, sc_d))
        s.append((U1[k] - Hs[k].scale(4)).scale(sc_n, sc_d))
        H.append((-(w[k]) - s2[k].scale(1, 3)).scale(sc_n, sc_d))
        w.append(Hw[k].scale(-2 * sc_n, sc_d))
        for j in range(nR):
            Rj = R[j]
            f = cau(Rj, Rj, k) + cau(C, Rj, k) + w[k].scale(mu4[j].numerator, mu4[j].denominator) - s2[k].scale(2, 3)
            Rj.append(f.scale(sc_n, sc_d))
    return [phi, s, H, w] + R


def ZERO(T):
    return IV(0, 0) if T is IV else TM([0])


def _mono(lam, x, f):
    """range of f over x in the box when df/dx = lam has a definite sign (then the extremes sit at the x-endpoints)"""
    if (lam.hi < 0 or lam.lo > 0) and x.lo < x.hi:
        return f(IV(x.lo, x.lo)).hull(f(IV(x.hi, x.hi)))
    return f(x)


def vector_field(box, mu2s):
    """IV evaluation of the right-hand side on a box  [phi, s, H, w, R_0.., G].
    G = U1(phi)/s is carried as an extra component with its own equation  G' = U2(phi) - G^2 + 4 H G  (an identity
    along solutions), and s' = U1 - 4Hs is written as s (G - 4H): this removes the phi/s dependency problem."""
    phi, s, H, w = box[:4]
    G = box[-1]
    p2 = phi * phi
    Bq = (phi.scale(10, 3) - (p2 * phi).scale(4, 9)).addc(Fraction(-4, 3))
    U2 = (phi * Bq).scale(2) + p2.addc(-1) * p2.scale(-4, 3).addc(Fraction(10, 3))
    C = G.scale(2) - H.scale(6)
    s2 = s * s
    out = [s, s * (G - H.scale(4)), -w - s2.scale(1, 3), (H * w).scale(-2)]
    for j, m in enumerate(mu2s):
        R = box[4 + j]
        m4 = Fraction(m) + 4
        rest = w.scale(m4.numerator, m4.denominator) - s2.scale(2, 3)
        out.append(_mono(R.scale(2) + C, R, lambda r: r * r + C * r + rest))
    H4 = H.scale(4)
    out.append(_mono(H4 - G.scale(2), G, lambda g: U2 - g * (g - H4)))
    return out


def apriori_box(x0, mu2s, h, maxit=40):
    """Find a box Bh with  x0 + [0,h] f(Bh)  subset of  Bh  (then every solution starting in x0 exists on [0,h]
    and stays in Bh).  x0: list of IV."""
    hI = IV(0, fceil(h))
    f0 = vector_field(x0, mu2s)
    B = [x + hI * f for x, f in zip(x0, f0)]
    for it in range(maxit):
        # inflate
        Bi = []
        for b, x in zip(B, x0):
            wd = b.width()
            e = wd // 8 + 1
            Bi.append(IV(b.lo - e, b.hi + e))
        try:
            f = vector_field(Bi, mu2s)
        except ArithmeticError:
            return None
        Bn = [x + hI * g for x, g in zip(x0, f)]
        if all(bn.subset_of(bi) for bn, bi in zip(Bn, Bi)):
            return Bn
        B = [bn.hull(bi) for bn, bi in zip(Bn, Bi)]
        if any(b.mag() > (ONE << 24) for b in B):      # diverging Picard iteration: give up, caller halves the step
            return None
    return None


def step(state, mu2s, h, N, tau_interval=False):
    """One validated Taylor step of size h (Fraction).  state: list of TM at y_k  (4 background components, then Riccati variables).
    Returns (new_state at y_k + h, tube, max scaled Lagrange remainder); tube (if tau_interval) is a TM enclosure valid
    for ALL y in [y_k, y_k+h].

    Background components: direct Taylor-model evaluation of the Taylor step (remainders propagate by interval arithmetic).
    Riccati components (scalar, do not feed back): mean-value form in the R-direction.  With R_k(u) = P(u) + e, e in I_k, 0 in I_k,
        R_{k+1} = Phi_h(x, P) + dPhi/dR * e,   dPhi/dR = exp( int_0^h (2R + 2G - 6H) dy ) in (0, exp(h * Lam_hi)],
    where Lam_hi = sup of 2R + 2G - 6H over the verified a-priori box.  Phi_h(x, P) is evaluated with a thin R input."""
    ph0 = state[0]; p2 = ph0 * ph0
    G_tm = (p2.addc(-1) * (ph0.scale(10, 3) - (p2 * ph0).scale(4, 9)).addc(Fraction(-4, 3))) * state[1].recip()
    x0 = [x.bound() for x in state] + [G_tm.bound()]
    Bh = apriori_box(x0, mu2s, h)
    if Bh is None:
        raise ArithmeticError("no a-priori enclosure; reduce the step")
    GB = Bh[-1]; Bh = Bh[:-1]
    thin = list(state[:4]) + [TM(x.c, 0, 0) for x in state[4:]]
    X = taylor(thin, mu2s, h, N - 1)             # coefficients 0..N-1 as Taylor models
    XB = taylor(Bh, mu2s, h, N, G0=GB)           # coefficient N over the a-priori box: Lagrange remainder
    # contraction factors for the Riccati components
    CB = GB.scale(2) - Bh[2].scale(6)
    hc = fceil(h); hf = ffloor(h)
    new = []; tube = []
    for i in range(len(state)):
        acc = X[i][0]
        for k in range(1, N):
            acc = acc + X[i][k]
        rem = XB[i][N]
        res = acc.add_iv(rem)
        extra = IV(0, 0)
        if i >= 4:
            I = IV(min(state[i].rlo, 0), max(state[i].rhi, 0))
            lam_hi = (Bh[i].scale(2) + CB).hi
            if lam_hi < 0:
                d_hi = cdiv(ONE * ONE, ONE + (((-lam_hi) * hf) >> P))        # exp(-x) <= 1/(1+x), x = h|lam| rounded down
            else:
                xx = cdiv(lam_hi * hc, ONE)
                if 2 * xx >= ONE:
                    raise ArithmeticError("Riccati expansion bound too large; reduce the step")
                d_hi = cdiv(ONE * ONE, ONE - xx)                              # exp(x) <= 1/(1-x)
            extra = IV(0, d_hi) * I                                          # tube: factor in (0, max(1,d_hi)]
            res = res.add_iv(IV(d_hi, d_hi) * I)
            res = TM(res.c, min(res.rlo, 0), max(res.rhi, 0))      # keep 0 in the remainder: P(u) itself stays in the box
        new.append(res)
        if tau_interval:
            # x(y_k + tau h) in X_0 + sum_{k>=1} X_k [0,1] + rem [0,1]   (+ [0, max(1,d_hi)] * I for Riccati components)
            lo = 0; hi = 0
            for k in range(1, N):
                b = X[i][k].bound()
                lo += min(b.lo, 0); hi += max(b.hi, 0)
            lo += min(rem.lo, 0); hi += max(rem.hi, 0)
            tb = X[i][0].add_iv(IV(lo, hi))
            if i >= 4:
                tb = tb.add_iv(IV(0, max(d_hi, ONE)) * I)
            tube.append(tb)
    return new, tube, max(XB[i][N].mag() for i in range(len(state)))


def shell_point(state, mu2s, hnom, N, HBiv):
    """Rigorous evaluation of the state AT the shell point  y_b(u, h) defined by  w(y_b) = h,  h in HBiv,  assuming
    y_b in [y_a, y_a + hnom] (verified here).  With tau = (y - y_a)/hnom, the Taylor step polynomial of w is solved for
    tau by a non-rigorous Newton iteration in Taylor-model arithmetic; the result tau_c(u) is then VERIFIED:
        w(tau_c(u) - eps) > sup HB   and   w(tau_c(u) + eps) < inf HB   for all u,
    so (w strictly decreasing) tau*(u,h) = tau_c(u) + [-eps, eps].  All state components are then evaluated at this
    Taylor-model-valued tau.  Returns (list of TM at the shell, eps as Fraction of hnom, info)."""
    ph0 = state[0]; p2 = ph0 * ph0
    G_tm = (p2.addc(-1) * (ph0.scale(10, 3) - (p2 * ph0).scale(4, 9)).addc(Fraction(-4, 3))) * state[1].recip()
    x0 = [x.bound() for x in state] + [G_tm.bound()]
    Bh = apriori_box(x0, mu2s, hnom)
    if Bh is None:
        raise ArithmeticError("no a-priori enclosure for the shell step")
    GB = Bh[-1]; Bh = Bh[:-1]
    if not (Bh[2].lo > 0 and Bh[3].lo > 0):
        raise ArithmeticError("H, w not positive on the shell step")       # w' = -2Hw < 0: w strictly decreasing
    thin = list(state[:4]) + [TM(x.c, 0, 0) for x in state[4:]]
    X = taylor(thin, mu2s, hnom, N - 1)
    XB = taylor(Bh, mu2s, hnom, N, G0=GB)
    unit = IV(0, ONE)

    def ev(i, tau):
        acc = X[i][N - 1]
        for k in range(N - 2, -1, -1):
            acc = acc * tau + X[i][k]
        return acc.add_iv(XB[i][N] * unit)            # Lagrange remainder * tau^N, tau in [0,1]

    # non-rigorous Newton for tau_c(u)
    hmid = TM([(HBiv.lo + HBiv.hi) // 2])
    tau = TM([ONE // 2])
    for _ in range(8):
        F = ev(3, tau) - hmid
        dF = X[3][N - 1].scale(N - 1)
        for k in range(N - 2, 0, -1):
            dF = dF * tau + X[3][k].scale(k)
        tau = tau - F * dF.recip()
        tau = TM(tau.c)
    eps = None
    for e in range(60, 8, -2):
        ee = ONE >> e
        tm_lo = TM(tau.c, -ee, -ee); tm_hi = TM(tau.c, ee, ee)
        if tm_lo.bound().lo < 0 or tm_hi.bound().hi > ONE:
            continue
        if ev(3, tm_lo).bound().lo > HBiv.hi and ev(3, tm_hi).bound().hi < HBiv.lo:
            eps = ee; break
    if eps is None:
        raise ArithmeticError("shell point not verified")
    tau_tm = TM(tau.c, -eps, eps)
    out = []
    lamB = None
    for i in range(len(state)):
        v = ev(i, tau_tm)
        if i >= 4:
            I = IV(min(state[i].rlo, 0), max(state[i].rhi, 0))
            lam_hi = (Bh[i].scale(2) + GB.scale(2) - Bh[2].scale(6)).hi
            if lam_hi >= 0:
                raise ArithmeticError("Riccati not contracting on the shell step")
            v = v.add_iv(unit * I)                     # dPhi/dR in (0,1]
        out.append(v)
    return out, Fraction(eps, ONE), {"tau_center_coeffs": [dec(c, 30) for c in tau.c[:3]], "H_box_lo_positive": Bh[2].lo > 0}


def dec(n, digits=40):
    """decimal string of the fixed-point integer n (truncated toward -inf at the given digits; display only)"""
    fr = Fraction(n, ONE)
    sgn = "-" if fr < 0 else ""
    fr = abs(fr)
    ip = fr.numerator // fr.denominator
    fp = fr - ip
    ds = str((fp.numerator * 10 ** digits) // fp.denominator).rjust(digits, "0")
    return "%s%d.%s" % (sgn, ip, ds)


def iv_str(iv, digits=40):
    """outward-rounded decimal strings [lo, hi]"""
    from decimal import Decimal, getcontext, ROUND_FLOOR, ROUND_CEILING
    getcontext().prec = digits + 10
    def conv(n, mode):
        getcontext().rounding = mode
        v = Decimal(n) / Decimal(ONE)
        q = Decimal(1).scaleb(v.adjusted() - digits + 1) if v != 0 else Decimal(1)
        return str(v.quantize(q, rounding=mode))
    return [conv(iv.lo, ROUND_FLOOR), conv(iv.hi, ROUND_CEILING)]
