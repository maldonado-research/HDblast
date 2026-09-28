#!/usr/bin/env python3
"""HDBLAST Chat 10 - computer-assisted proof that the stable shell S_8/5 has NO scalar bound state (mu^2 < 9/4).

Python 3.9 standard library only.  All rigorous steps use cert_core.py (byte-identical Chat 9 core: outward-rounded 2^-384
fixed point, Taylor models, validated Taylor integrator) and exist_lib.py (byte-identical Task A library).  Floats are
used only to GUESS the barrier coefficients k (every guess is verified by an interval gate) and in print-outs.
Run:  python3 certify_stability.py        (any failed gate -> SystemExit, non-zero exit status)

SETTING (inherited, see PREMISES in the certificate)
  Background  phi' = s, s' = U1(phi) - 4Hs, H' = -w - s^2/3, w' = -2Hw   (H = rho'/rho, w = 1/rho^2), regular cone at y = 0,
  phi(0) = p.  Shell S_8/5: sigma_t = 2W + t(1 + c phi + d phi^2/2), t = 1/1000, c = 5975949350280/10^13, d = 8/5 (exact).
  Task A (EXISTENCE_CERTIFICATE_S85.json): a root (p*, y*) of both junction conditions lies in BOX = [p_lo,p_hi] x [ya,yb].
  Scalar perturbations (longitudinal gauge, X := s chi/3 = -(psi' + 2H psi), lam := mu^2 + 4, G := U1/s):
      psi' = -2H psi - X,    X' = [lam w - (2/3) s^2] psi + (2G - 8H) X,
      Riccati variable R = X/psi:   R' = R^2 + (2G - 6H) R + lam w - (2/3) s^2,
  regular (Frobenius, normalisable) solution psi ~ y^(1/2+nu), nu = sqrt(9/4 - mu^2):  r := y R -> r_+ := -5/2 - nu  (y -> 0).
  Shell condition  B chi + 3 lam psi/(rho^2 phi') = 0   <=>   (3 psi/s) * b = 0,   b := lam w + B R,
      B = phi''/phi' + sigma_t''/2 = G - 4H + 2 phi + t d/2.

THEOREM CHECKED.  For EVERY (p, y_s) in BOX (hence for the shell (p*, y*)) and EVERY real mu^2 < 9/4, the regular solution
  satisfies  (i) R finite on (0, y_s]  (so psi has no zero on (0, y_s]),   (ii) b(mu^2) < 0 at y_s,   and (iii) B > 0.
  Hence the shell condition is violated for every mu^2 < 9/4: no scalar bound state, in particular no tachyon and no
  zero mode.  The same statement is proved independently (without the threshold barrier) for mu^2 <= mu*^2 with
  mu*^2 = 0, 2, 11/5.

PROOF STEPS CHECKED HERE
  1. Cone chart: Chat 9/Task A Volterra contraction gates (exist_lib.cone_start).  PARITY LEMMA (Lemma P below): beta = H - 1/y and s are odd in y, phi is even; the Taylor recursion is re-run with the
     vanishing coefficients set to exact zeros.  The series is evaluated at y00 = 1/320 and enclosed on (0, y00].
  2. Cone barrier (new; works at the threshold nu = 0).  With  y r' = (r - r_+)(r - r_-) + delta1 r + delta0,
         delta1 = c1 - 5 = 2 y sigma'/sigma + 2 y beta          (sigma = s/y; identity U1 = 5 sigma + y sigma' + 4 y beta sigma),
         delta0 = c0 - lam = lam (2 y beta + y^2 beta^2 - y^2 s^2/12 + y^2 U/6) - (2/3) y^2 s^2,
     D1 = delta1/y^2 and D0 = delta0/y^2 are enclosed on (0, y00].  For rbar(y) = r_+ + k y^2 the gate
         M(k) := 2k - [ k (k y^2 - 2 nu) + D1 (r_+ + k y^2) + D0 ]   > 0  on (0,y00]   (k = k_hi)      resp.  < 0  (k = k_lo)
     says y rbar' > Q(rbar) (resp. <).  Gronwall argument (Lemma C below): the regular solution obeys
     r_+ + k_lo y^2 < r(y) < r_+ + k_hi y^2 on (0, y00]; in particular R is finite there.
  3. Validated integration from y00 to the shell tube [ya, yb] of the background (Taylor model in p) together with, for each
     listed mu*^2, the Riccati solutions started at the lower / upper barrier value at y00 (mean-value form, cert_core.step).
     Scalar comparison: R_lower(y) <= R(y; mu*^2) <= R_upper(y); and R(y; mu^2) <= R(y; mu*^2) for mu^2 <= mu*^2 (Lemma M).
  4. Shell: B > 0 and w > 0 on the tube, and  b_up := lam* w + B R_upper < 0  on the tube.  Then for mu^2 <= mu*^2 (resp. < 9/4
     for the threshold value): b(mu^2) = lam w + B R <= lam* w + B R_upper < 0.

ANALYTIC LEMMAS (pencil-and-paper; NOT machine-checked; they are the only non-computational steps besides the premises)
  Lemma P (parity).  The Volterra map of the cone chart, (beta, s, phi-p) -> (y int_0^1 u^2 f(uy) du, y int_0^1 u^4 g(uy) du,
     y int_0^1 s(uy) du), f = -beta^2 - s^2/4 - U(phi)/6, g = U1(phi) - 4 beta s, maps the closed non-empty subset
     {beta odd, s odd, phi - p even} of the contraction ball into itself (f, g are then even).  The unique fixed point lies
     in it: beta_n = s_n = 0 for even n, phi_n = 0 for odd n.
  Lemma C (cone barrier).  Let lam <= 25/4, nu = sqrt(25/4 - lam) >= 0, psi = y^(1/2+nu) h(y), h analytic, h(0) = 1, the regular
     solution (the first Frobenius solution; it exists also at the double root nu = 0), r = yR = -y psi'/psi - 2yH, so
     r - r_+ = O(y).  Exactly, y r' = Q(y,r) := (r - r_+)(r - r_-) + delta1 r + delta0, r_- = -5/2 + nu.  Let
     rbar = r_+ + k y^2 with y rbar' - Q(y, rbar) = y^2 M(y), M >= M0 > 0 on (0,y00] (gate).  On the maximal interval
     J = (0,y1) of (0,y00] where psi > 0, e := r - rbar satisfies  y e' = A e - y^2 M,  A = (r - r_+) + (rbar - r_+) - 2nu + delta1
     = -2nu + O(y).  With E(y) = exp(int_y^c A(t)/t dt) (c in J fixed):  (eE)' = -y M E < 0;  E is bounded near 0 (nu >= 0,
     A + 2nu = O(t)), e -> 0, hence e(y)E(y) = -int_0^y t M E dt < 0:  r < rbar on J.  If y1 <= y00 were a zero of psi, then
     X(y1) = -psi'(y1) > 0 (X(y1) != 0 for a nontrivial solution of the linear first-order system) and R -> +infinity
     as y -> y1-, contradicting r < rbar.  So psi > 0 and r < rbar on (0,y00].  Reversed inequalities give r > r_+ + k_lo y^2.
  Lemma M (monotone comparison).  (a) Same lam: if R_a, R_b solve the Riccati equation on [y00, Y], R_a(y00) <= R_b(y00)
     and R_b is finite on [y00,Y], then R_a <= R_b and R_a is finite on [y00, Y]:  D = R_b - R_a obeys the linear homogeneous
     equation D' = (R_a + R_b + 2G - 6H) D, so D keeps its sign while both are finite; R_a cannot reach +infinity below R_b,
     and no solution of R' = R^2 + ... reaches -infinity at the right end of an interval (v = 1/R has v' = -1 at v = 0).
     (b) lam_1 < lam_2 <= 25/4, regular solutions R_1, R_2, R_2 finite on (0,Y]:  r_1(0+) = -5/2 - nu_1 < -5/2 - nu_2 = r_2(0+),
     so D = R_2 - R_1 > 0 near 0;  D' = (R_1 + R_2 + 2G - 6H) D + (lam_2 - lam_1) w  with w > 0 forbids a first zero of D.
     Hence R_1 < R_2 and R_1 finite on (0,Y].
  ASSEMBLY for lam* = 25/4 (mu^2 < 9/4):  R(.;25/4) is finite on (0,y00] and <= rbar/y (Lemma C), hence <= R_upper on
     [y00, y_s] (Lemma M a; R_upper is finite because the validated integrator encloses it in bounded boxes); for lam < 25/4,
     R(.;lam) < R(.;25/4) (Lemma M b) is finite on (0,y_s], so psi > 0 there, and with B > 0, w > 0:
     b(lam) = lam w + B R(lam) < (25/4) w + B R_upper < 0  at y_s.   (For mu*^2 = 0, 2, 11/5 the same with nu* > 0.)
"""
LEMMAS = __doc__[__doc__.index("ANALYTIC LEMMAS"):]
import sys, os, json, time, hashlib, math
from fractions import Fraction as Fr
import cert_core as cc
from cert_core import IV, TM, ONE, P
import exist_lib as L            # also installs the combined a-priori-box search into cert_core (same verified inclusion test)
from exist_lib import gate

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
DEG = 2
cc.set_degree(DEG)
Y0 = L.Y0; NC = L.NC
T00 = Fr(1, 16)                  # y00 = T00 * Y0 = 1/320 : end of the barrier chart
Y00 = T00 * Y0
RC = Fr(1, 2); r_be, r_s, r_ph = Fr(6, 100), Fr(2, 10 ** 6), Fr(11, 10 ** 7)      # as in exist_lib.cone_start (gated there)
MU2S = [Fr(0), Fr(2), Fr(11, 5), Fr(9, 4)]
K_MARGIN = float(os.environ.get("HDB_K_MARGIN", "1e-3"))
NEG_D = os.environ.get("HDB_NEG_CONTROL_D")          # negative control: replace d (e.g. 13/10: a bound state exists, gates must FAIL)
D_EFF = Fr(NEG_D) if NEG_D else L.D_PAR

# ------------------------------------------------------------------ inherited existence certificate (Task A)
exist_path = os.path.join(HERE, "inputs", "EXISTENCE_CERTIFICATE_S85.json")
EX = json.load(open(exist_path))
gate("existence certificate is for the exact S_8/5 parameters", EX["parameters_exact"]["t"] == str(L.T_PAR) and
     EX["parameters_exact"]["c"] == str(L.C_PAR) and EX["parameters_exact"]["d"] == str(L.D_PAR) and
     EX["status"] == "PASS_EXISTENCE_CERTIFIED_BY_PRECONDITIONED_MIRANDA")
c_mid = int(EX["BOX_exact"]["p_mid_fixed_point_numerator_over_2^384"]); c_rad = ONE >> 100
p_lo = Fr(c_mid - c_rad, ONE); p_hi = Fr(c_mid + c_rad, ONE)
ya = Fr(EX["BOX_exact"]["ya"]); yb = Fr(EX["BOX_exact"]["yb"])
gate("BOX of the existence certificate reproduced exactly", str(p_lo) == EX["BOX_exact"]["p_lo"] and str(p_hi) == EX["BOX_exact"]["p_hi"] and ya < yb)
P_TM = TM([c_mid, c_rad])

# ------------------------------------------------------------------ 1. cone: Task A gates, then parity-enforced recursion
state_Y0_ref, cone_info = L.cone_start(P_TM, p_lo, p_hi)
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
    b_ = ph[k].scale(10, 3) - p3[k].scale(4, 9)
    Bq.append(b_.addc(Fr(-4, 3)) if k == 0 else b_)
    Us.append(cau(W1s, W1s, k).scale(1, 2) - cau(Ws, Ws, k).scale(2, 3))
    U1s.append(cau(W1s, Bq, k))
    be2.append(cau(be, be, k)); s2.append(cau(s_, s_, k)); bes.append(cau(be, s_, k))
    f = -be2[k] - s2[k].scale(1, 4) - Us[k].scale(1, 6)
    g = U1s[k] - bes[k].scale(4)
    # PARITY LEMMA: beta_n = s_n = 0 for even n, (phi)_n = 0 for odd n  (exact zeros)
    be.append(f.scale(Y0.numerator, Y0.denominator * (n + 2)) if n % 2 else zero)
    s_.append(g.scale(Y0.numerator, Y0.denominator * (n + 4)) if n % 2 else zero)
    ph.append(s_[k].scale(Y0.numerator, Y0.denominator * n) if n % 2 == 0 else zero)


def tail_frac(r, y):
    x = y / RC
    return r * x ** (NC + 1) / (1 - x)


def tm_eval(lst, t, r, y):
    """sum_n lst[n] t^n (t exact rational) + Cauchy tail at y"""
    acc = lst[0]
    tn = Fr(1)
    for a in lst[1:]:
        tn *= t
        acc = acc + a.scale(tn.numerator, tn.denominator)
    tl = cc.fceil(tail_frac(r, y))
    return TM(acc.c, acc.rlo - tl, acc.rhi + tl)


# consistency of the parity-enforced series with Task A's cone state at Y0
for nm, lst, r, ref in (("phi", ph, r_ph, state_Y0_ref[0]), ("s", s_, r_s, state_Y0_ref[1])):
    a = tm_eval(lst, Fr(1), r, Y0).bound(); b_ = ref.bound()
    gate("parity-enforced cone series agrees with Task A cone state at Y0 (%s)" % nm, a.lo <= b_.hi and b_.lo <= a.hi)
be00 = tm_eval(be, T00, r_be, Y00); s00 = tm_eval(s_, T00, r_s, Y00); ph00 = tm_eval(ph, T00, r_ph, Y00)
H00 = be00.addc(1 / Y00)
w00 = H00 * H00 - (s00 * s00).scale(1, 12) + L.U_of(ph00).scale(1, 6)
gate("s > 0 at y00", s00.bound().lo > 0)

# enclosures on (0, y00]  (tau = y/Y0 in [0, T00]); all p in the box at once
tauI = IV(0, cc.fceil(T00))


def horner(lst, start, weights=None):
    acc = IV(0, 0)
    for n in range(len(lst) - 1, start - 1, -1):
        a = lst[n].bound()
        if weights:
            a = a.scale(weights(n))
        acc = acc * tauI + a
    return acc


def sym(fr):
    v = cc.fceil(fr)
    return IV(-v, v)


x00 = Y00 / RC
beI = horner(be, 0) + sym(tail_frac(r_be, Y00))                                        # beta
sI = horner(s_, 0) + sym(tail_frac(r_s, Y00))                                          # s
phI = horner(ph, 0) + sym(tail_frac(r_ph, Y00))                                        # phi
be_y = (horner(be, 1) + sym(tail_frac(r_be, Y00) / T00)).scale(Y0.denominator, Y0.numerator)      # beta/y   (tail/y00: powers n-1 >= 0)
sgI = (horner(s_, 1) + sym(tail_frac(r_s, Y00) / T00)).scale(Y0.denominator, Y0.numerator)        # sigma = s/y
gate("sigma = s/y > 0 on (0, y00]", sgI.lo > 0)
tail_ds = r_s / Y00 ** 3 * x00 ** (NC + 1) * (Fr(NC) / (1 - x00) + 1 / (1 - x00) ** 2)            # sum_{n>NC} (n-1)|s_n| y^(n-3)
dsg_y = horner(s_, 3, weights=lambda n: n - 1).scale(Y0.denominator ** 3, Y0.numerator ** 3) + sym(tail_ds)   # sigma'/y (s_2 = 0 exactly)
gate("s_2 coefficient is an exact zero (parity), so sigma'/y has no 1/y term", s_[2].c == zero.c and s_[2].rlo == 0 and s_[2].rhi == 0)
D1 = (dsg_y * sgI.recip()).scale(2) + be_y.scale(2)
y2I = IV(0, cc.fceil(Y00 * Y00))
E0 = be_y.scale(2) + beI * beI - (sI * sI).scale(1, 12) + L.U_of(phI).scale(1, 6)      # (y^2 w - 1)/y^2
S2 = (sI * sI).scale(2, 3)

BARRIER = {}
R_START = []
for m in MU2S:
    lam = m + 4
    nu2 = Fr(9, 4) - m
    nuI = IV(0, 0) if nu2 == 0 else L.iv_sqrt(IV.frac(nu2))
    rpI = (-nuI).addc(Fr(-5, 2))                       # r_+ = -5/2 - nu
    D0 = E0.scale(lam.numerator, lam.denominator) - S2
    kI = (D0 + rpI * D1) * (nuI.scale(2).addc(2)).recip()
    kmid = (kI.lo + kI.hi) / 2 / ONE                    # float guess (guidance only)
    k_hi = Fr(math.ceil((kmid + K_MARGIN) * 2 ** 24), 2 ** 24)
    k_lo = Fr(math.floor((kmid - K_MARGIN) * 2 ** 24), 2 ** 24)

    def M(k):
        kyy = y2I.scale(k.numerator, k.denominator)
        inner = (kyy - nuI.scale(2)).scale(k.numerator, k.denominator) + D1 * (rpI + kyy) + D0
        return (-inner).addc(2 * k)

    Mh = M(k_hi); Ml = M(k_lo)
    gate("cone barrier mu2=%s: y rbar' > Q(rbar) for rbar = r_+ + k_hi y^2 on (0,y00]" % m, Mh.lo > 0, str(cc.iv_str(Mh, 8)))
    gate("cone barrier mu2=%s: y rbar' < Q(rbar) for rbar = r_+ + k_lo y^2 on (0,y00]" % m, Ml.hi < 0, str(cc.iv_str(Ml, 8)))
    y00I = IV.frac(Y00 * Y00)
    st_lo = (rpI + y00I.scale(k_lo.numerator, k_lo.denominator)).scale(Y00.denominator, Y00.numerator)
    st_hi = (rpI + y00I.scale(k_hi.numerator, k_hi.denominator)).scale(Y00.denominator, Y00.numerator)
    gate("barrier start values ordered mu2=%s" % m, st_lo.lo < st_hi.hi)
    R_START.append(st_lo.lo); R_START.append(st_hi.hi)       # lower solution starts at/below the lower barrier, upper at/above the upper
    BARRIER[str(m)] = {"nu": cc.iv_str(nuI, 25) if nu2 else ["0", "0"], "k_enclosure_(guidance)": cc.iv_str(kI, 10), "k_lo": str(k_lo), "k_hi": str(k_hi),
                       "M(k_hi) on (0,y00]": cc.iv_str(Mh, 8), "M(k_lo) on (0,y00]": cc.iv_str(Ml, 8),
                       "R_lower_start_at_y00": cc.dec(st_lo.lo, 25), "R_upper_start_at_y00": cc.dec(st_hi.hi, 25)}
MU_LIST = [m for m in MU2S for _ in (0, 1)]
print("cone + barriers done %.1fs   D1 in %s" % (time.time() - T0, cc.iv_str(D1, 8)), flush=True)

# ------------------------------------------------------------------ 3. validated integration  y00 -> ya -> tube [ya, yb]
if NEG_D:
    L.D_PAR = D_EFF          # only affects B_of / sigma (negative control)
state = [ph00, s00, H00, w00] + [TM([r], 0, 0) for r in R_START]
TOL = L.TOL; NORD_MAX = L.NORD_MAX
NORD = 24; y = Y00; h_prev = Y00 / 64; nstep = 0; max_rem = 0
rw = lambda st: max(x.rhi - x.rlo for x in st)


def adaptive(st, h, tube=False, allow_halving=True):
    global NORD
    while True:
        try:
            out = cc.step(st, MU_LIST, h, NORD, tau_interval=tube)
        except ArithmeticError:
            out = None
        if out is not None and out[2] <= TOL:
            if out[2] < (TOL >> 40) and NORD > 12:
                NORD -= 2
            return out, h
        if out is not None and NORD + 4 <= NORD_MAX:
            NORD += 4
            continue
        if not allow_halving:
            raise SystemExit("fixed-size step failed")
        h = h / 2
        if h < Fr(1, 2 ** 40):
            raise SystemExit("step size underflow")


CHECK = []
targets = [Y0, Fr(1), Fr(2), Fr(4), Fr(6), Fr(8), ya]
for yt in targets:
    while y < yt:
        h = min(L.KAPPA * y, L.HMAX, 2 * h_prev)
        h = Fr(int(h * 2 ** 20), 2 ** 20)
        trunc = y + h >= yt
        if trunc:
            h = yt - y
        (state, _, rem), h2 = adaptive(state, h)
        if not (trunc and h2 == h):
            h_prev = h2
        y += h2; nstep += 1; max_rem = max(max_rem, rem)
        if nstep % 40 == 0:
            print("   step %3d y=%.5f N=%d R(9/4,up)=%.6f remwidth=%.2e  %.0fs" % (nstep, float(y), NORD, state[-1].c[0] / ONE, rw(state) / ONE, time.time() - T0), flush=True)
    if yt == Y0:
        for i, nm in enumerate(("phi", "s", "H", "w")):
            a = state[i].bound(); b_ = state_Y0_ref[i].bound()
            gate("integrated state at Y0 agrees with Task A cone series (%s)" % nm, a.lo <= b_.hi and b_.lo <= a.hi)
    CHECK.append({"y": str(yt), "phi": cc.iv_str(state[0].bound(), 25), "s": cc.iv_str(state[1].bound(), 25), "H": cc.iv_str(state[2].bound(), 25),
                  "w": cc.iv_str(state[3].bound(), 25), "R_lower_upper_per_mu2": [cc.iv_str(x.bound(), 18) for x in state[4:]]})
state_a = state
(state_b, tube, rem2), _ = adaptive(state_a, yb - ya, tube=True, allow_halving=False)
max_rem = max(max_rem, rem2)
print("integration done %.1fs steps=%d" % (time.time() - T0, nstep), flush=True)

# ------------------------------------------------------------------ 4. shell gates on the whole tube (all p in the box, all y in [ya,yb])
tb = [x.bound() for x in tube]
phi_t, s_t, H_t, w_t = tb[:4]
gate("s > 0 on the tube", s_t.lo > 0)
gate("H > 0 on the tube (H' < 0, hence H > 0 and w' = -2Hw < 0 on (0,yb])", H_t.lo > 0)
gate("w > 0 on the tube", w_t.lo > 0)
for i, key in enumerate(("phi_b", "phi'_b = s_b", "H_b = rho'_b/rho_b", "h = w_b = 1/rho_b^2")):
    lo, hi = [Fr(v) for v in EX["enclosures"][key]]
    gate("tube intersects the Task A enclosure (%s)" % key, Fr(tb[i].lo, ONE) <= hi and lo <= Fr(tb[i].hi, ONE))
G1_t, G2_t = L.residuals(tb[:4])
if not NEG_D:
    gate("both junction residual enclosures over the tube contain 0", G1_t.lo <= 0 <= G1_t.hi and G2_t.lo <= 0 <= G2_t.hi)
B_t = L.B_of(tb[:4])
gate("B = phi''/phi' + sigma_t''/2 > 0 on the tube", B_t.lo > 0, str(cc.iv_str(B_t, 20)))
RESULT = {}
for j, m in enumerate(MU2S):
    lam = m + 4
    Rl = tb[4 + 2 * j]; Rh = tb[5 + 2 * j]
    gate("Riccati sandwich ordered at the shell mu2=%s" % m, Rl.lo <= Rh.hi)
    lw = w_t.scale(lam.numerator, lam.denominator)
    b_up = lw + B_t * IV(Rh.hi, Rh.hi)
    b_dn = lw + B_t * IV(Rl.lo, Rl.lo)
    bI = IV(b_dn.lo, b_up.hi)
    RESULT[str(m)] = {"lam": str(lam), "R_lower_solution_on_tube": cc.iv_str(Rl, 20), "R_upper_solution_on_tube": cc.iv_str(Rh, 20),
                      "R_regular_solution_in": [cc.iv_str(Rl, 20)[0], cc.iv_str(Rh, 20)[1]],
                      "b_enclosure = lam w + B R": cc.iv_str(bI, 20), "relative_margin |b|/(lam w) >=": "%.4f" % (-(bI.hi / ONE) / (lw.hi / ONE))}
    gate("junction mismatch b = lam w + B R_upper < 0 on the tube, mu*2=%s  (=> no bound state with mu^2 %s %s)" % (m, "<" if m == Fr(9, 4) else "<=", m),
         bI.hi < 0, str(cc.iv_str(bI, 16)))
    print("mu*^2 = %-5s  R in [%s, %s]   b in %s" % (m, cc.iv_str(Rl, 12)[0], cc.iv_str(Rh, 12)[1], cc.iv_str(bI, 12)), flush=True)

cert = {
    "artifact": "HDBLAST Chat 10: certified absence of scalar bound states on the stable de Sitter shell S_8/5",
    "status": "PASS_NO_SCALAR_BOUND_STATE_BELOW_9/4_CERTIFIED",
    "parameters_exact": {"t": str(L.T_PAR), "c": str(L.C_PAR), "d": str(D_EFF), "sigma_t": "2W + t(1 + c phi + d phi^2/2)", "W": "1 - phi + phi^3/3"},
    "statement": "For every (p, y_s) in BOX (Task A box; contains the certified junction root (p*, y*)) and every real mu^2 < 9/4: the Frobenius-regular "
                 "solution of the longitudinal-gauge scalar system has R = X/psi finite on (0, y_s] (psi has no zero), B > 0, and the junction mismatch "
                 "b = (mu^2+4) w + B R is < 0 at y_s.  Hence no scalar bound state with mu^2 < 9/4 exists (no tachyon, no zero mode, no massive bound state). "
                 "Independently of the threshold barrier the same holds for mu^2 <= 0, <= 2 and <= 11/5.",
    "BOX_exact": EX["BOX_exact"],
    "B_on_BOX": cc.iv_str(B_t, 30),
    "result_per_mu_star2": RESULT,
    "cone_barrier": {"y00": str(Y00), "D1 = delta1/y^2 on (0,y00]": cc.iv_str(D1, 12), "(y^2 w - 1)/y^2 on (0,y00]": cc.iv_str(E0, 12), "per_mu2": BARRIER,
                     "k_margin_float": K_MARGIN},
    "cone": cone_info,
    "integration": {"tm_degree": DEG, "fixed_point_bits": P, "steps": nstep, "taylor_order_final": NORD, "remainder_tolerance": "2^-140",
                    "max_scaled_Lagrange_remainder": "%.3e" % (max_rem / ONE), "TM_remainder_width_at_ya": "%.3e" % (rw(state_a) / ONE)},
    "checkpoints_all_p": CHECK,
    "premises_inherited": [
        "Task A existence certificate (inputs/EXISTENCE_CERTIFICATE_S85.json): a root of both junction conditions lies in BOX",
        "the longitudinal-gauge scalar system and the scalar junction condition B chi + 3 lam psi/(rho^2 phi') = 0 (Chat 9 derivation in two gauges, "
        "checked against the full linearised Einstein equations to 1e-9 in floating point; identical to Frolov-Kofman hep-th/0309002)",
        "Frobenius theory at the regular-singular cone point: bound state <=> mu^2 real < 9/4 and psi = y^(1/2+nu) h(y), h analytic, h(0) = 1 "
        "(Task B Lemmas 1-2; reality of mu^2 for B > 0: Task B Theorem 1).  Analyticity of the background at the cone is supplied by the gated "
        "disc-algebra contraction (the fixed point is holomorphic in |y| < 1/2).",
        "analytic lemmas P (parity), C (cone barrier, Gronwall), M (monotone comparison in mu^2) in the docstring of certify_stability.py (reproduced under analytic_lemmas): pencil-and-paper proofs, not machine-checked",
        "identification of the contraction fixed point with the regular cone solution (as in Chat 9 / Task A)"],
    "not_proved_here": [
        "uniqueness of the shell root in BOX (Task A caveat); the theorem holds for every root in BOX",
        "tensor and vector sectors; nonlinear stability; completeness of the scalar mode expansion (self-adjointness is cited in Task B, not proved)",
        "the derivation of the perturbation equations themselves"],
    "analytic_lemmas": LEMMAS.split("\n"),
    "source_sha256": {f: hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()
                      for f in ("cert_core.py", "exist_lib.py", "certify_stability.py", os.path.join("inputs", "EXISTENCE_CERTIFICATE_S85.json"))},
    "gates": L.GATES,
    "runtime_seconds": round(time.time() - T0, 1),
}
name = "STABILITY_CERTIFICATE_S85" + ("_NEGCONTROL_SHOULD_NOT_EXIST" if NEG_D else "") + os.environ.get("HDB_CERT_SUFFIX", "") + ".json"
json.dump(cert, open(os.path.join(HERE, name), "w"), indent=1)
print("ALL %d GATES PASSED.  wrote %s  (%.0fs)" % (len(L.GATES), name, time.time() - T0))
