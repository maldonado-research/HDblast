#!/usr/bin/env python3
"""HDBLAST Chat 9 - computer-assisted proof that the registered de Sitter shell has a tachyonic scalar bound state.

Python 3.9 standard library only (integers / fractions; floats are used ONLY to guess barrier points and the shell
position, every guess is verified afterwards).  Run:
    python3 certify_tachyon.py                                   -> TACHYON_CERTIFICATE.json        (mu^2 = -8, -15/2;  ~2 min)
    python3 certify_tachyon.py -7.7180 -7.7177 -7.71788 -7.71786 -> TACHYON_CERTIFICATE_TIGHT.json  (extra values;  ~4 min)
Any failed gate raises SystemExit (non-zero exit status).

THEOREM CHECKED (conditional only on the inherited premises listed under PREMISES and on the analytic derivation of the
perturbation equations, see "not_proved_here" in the certificate):
  Let p = phi_h range over the whole archived M462 root interval PH and let h range over the whole interval HB
  (M462 Krawczyk image for eta, times the archived h_lin enclosure).  For each such (p, h) let (phi, s, H, w)(y) be the
  regular-cone bulk solution (s = phi', H = rho'/rho, w = 1/rho^2) and let y_b be the unique point with w(y_b) = h. Then
     (i)   ya < y_b < yb for the explicit rationals ya, yb;
     (ii)  for every listed mu^2 the Frobenius-regular solution (psi ~ y^alpha, alpha = 1/2 + sqrt(9/4 - mu^2)) of the
           longitudinal-gauge scalar system has psi > 0 on (0, y_b];
     (iii) the scalar-junction mismatch  b(mu^2) = (mu^2+4) w + (G - 4H + 2 phi) R  at y_b is < 0 for mu^2 = -8 and > 0 for
           mu^2 = -15/2 (and has the certified signs listed in the certificate for the extra values).
  The true registered shell has (p, h) in PH x HB.  The Frobenius-normalised solution is continuous in mu^2, so the
  boundary-value problem has an eigenvalue in the certified bracket, which lies entirely in mu^2 < 0:
  a tachyonic scalar bound state, m^2 = mu^2 H^2 < 0.

EQUATIONS (kappa_5 = 1; ' = d/dy; y = proper distance from the first cone):
  phi' = s,  s' = U1(phi) - 4 H s,  H' = -w - s^2/3,  w' = -2 H w        (polynomial; the constraint
  H^2 = w + s^2/12 - U/6 is propagated by the flow and imposed at the cone; it is re-checked at the shell).
  Perturbation (xi = -2 psi, X := s chi/3 = -(psi' + 2 H psi), from the orchestrator's equations (C), (H), (BC)):
      psi' = -2 H psi - X,     X' = [(mu^2+4) w - (2/3) s^2] psi + (2G - 8H) X,     G := U1(phi)/s,
  Riccati variable R = X/psi:   R' = R^2 + (2G - 6H) R + (mu^2+4) w - (2/3) s^2.
  Scalar junction (linear detuning, sigma_t'' = 2 W''):  chi' + 2 s psi + W''(phi) chi = (3 psi/s) * b,
      b = (mu^2+4) w + (G - 4H + 2 phi) R.      With psi > 0, s > 0: sign(junction mismatch) = sign(b).

PROOF STEPS CHECKED BY THIS SCRIPT
  1. Cone: beta = H - 1/y, s, phi - p solve the Volterra system  beta = y int_0^1 u^2 f(uy) du,  s = y int_0^1 u^4 g(uy) du,
     phi - p = y int_0^1 s(uy) du,  f = -beta^2 - s^2/4 - U/6,  g = U1 - 4 beta s.  Exact rational inequalities show that
     the map is a contraction of a closed polydisc ball in the disc algebra |y| <= 1/2, uniformly for p in PH; Cauchy
     estimates bound the Taylor tails.  The Taylor coefficients are enclosed with Taylor models in p.  w(y0) follows from
     the constraint (which the solution of the Volterra system satisfies, with w rho^2 = 1).
  2. Riccati start: r = y R obeys  y r' = r^2 + c1 r + c0,  c1 = 1 + 2yG - 6yH -> 5,  c0 = (mu^2+4) y^2 w - (2/3) y^2 s^2.
     With enclosures of c1, c0 on (0, y0] a forward-invariant interval [r_lo, r_hi] containing the regular exponent
     r_+ = -5/2 - sqrt(9/4 - mu^2) traps the Frobenius solution.  Scalar ODE comparison then sandwiches R(y) between the
     solutions started at r_lo/y0 and r_hi/y0 (both are integrated; both stay finite, hence psi has no zero).
  3. Validated Taylor integration (adaptive order, Lagrange remainder on a verified a-priori box) of the polynomial system
     with Taylor models in p (degree D, interval remainder), fixed-point outward rounding with 2^-384 resolution.
     The Riccati components are propagated in mean-value form (their flow contracts; plain interval evaluation would
     turn the contraction into an expansion).  For a-priori boxes G is carried with its own equation
     G' = U2 - G^2 + 4HG (an identity along solutions) to avoid the phi/s dependency problem.
  4. Shell localisation: w(ya) > sup HB, w(yb) < inf HB for all p, w strictly decreasing (H > 0).  Sharp version: a
     Taylor-model-valued shell position tau_c(u) +- eps is verified and all fields are evaluated there.
  5. Sign of b at the shell for all (p, h) in PH x HB.
"""
import sys, os, json, time, hashlib
from fractions import Fraction as Fr
import cert_core as cc
from cert_core import IV, TM, ONE, P

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------ run parameters
DEG = 6                      # Taylor-model degree in the parameter u
NORD = 24                    # initial Taylor order in y (adapted)
NORD_MAX = 64
Y0 = Fr(1, 20)               # end of the cone chart
NC = 46                      # cone series degree
KAPPA = Fr(1, 8)            # step <= KAPPA * y
HMAX = Fr(1, 10)
Y1 = Fr(41, 5)               # 8.2: last regular grid point before the shell
MARGIN = Fr(1, 2 ** 25)
MU2S = [Fr(-8), Fr(-15, 2)]      # main pair first; further values (command line) tighten the bracket
if len(sys.argv) > 1:
    MU2S += [Fr(a) for a in sys.argv[1:]]
cc.set_degree(DEG)

# ------------------------------------------------------------------ inherited premises (archived project certificates)
PREMISES = {
    "phi_h_exact_enclosure": [
        "-14757382293836628409055945955638376056067244952601439453213988003982854014198008114112139564980144014991928783004985291113591199160241975371/14757395258967641292800000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",
        "-2361181167013806746700395108892597146300344827675010661093596338499171546322394485578886830451294183515168773799015215320850233940877724247773/2361183241434822606848000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"],
    "eta_center": "0.12006349132995513",
    "eta_K_minus_center": ["-1.72733124632899741302044101481753812787900415059993974864482879638671875E-8",
                           "1.70627498578617397617549907256999375260164697465370409190654754638671875E-8"],
    "h_lin": ["0.1609106890301122079438195943542259003570638392159557779696258416417065709726635459398300926784087278",
              "0.1609106890301122079438195943542259027318480501841361260022540269380224005584525375018980713522020316"],
    "c_star": ["0.597594935028013161992501798917377470951432737258136002302176766367145518338629217944627778807571400",
               "0.597594935028013161992501798917377499448843268876300178693714989922935473368097116689443522893091046"],
    "t": "1/1000",
    "sources": ["untitled folder 120/M462_CERTIFICATE_RESEND.json  (K_minus_center[0] = eta Krawczyk image; existence+uniqueness of the shooting root)",
                "untitled folder 139/CONTINUUM_CERTIFICATION_STEP_20260904/background/inputs/BACKGROUND_ROOT_CONTRACT.json  (phi_h_exact_enclosure, h_lin, c_star)"],
    "meaning": "h = h_lin * t * (1 + t*eta) is the brane Hubble^2 and equals w = 1/rho^2 at the shell (M462 normalisation q = sqrt(h) rho = 1); "
               "phi_h = -1 + t^(9/5) delta.  These enclosures are NOT re-proved here.  c_star is used only for a consistency diagnostic."}


def check_premise_files():
    """If the archived files are present, assert that the embedded numbers are byte-identical to them."""
    out = {}
    base = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
    f1 = os.path.join(base, "untitled folder 139", "CONTINUUM_CERTIFICATION_STEP_20260904", "background", "inputs", "BACKGROUND_ROOT_CONTRACT.json")
    f2 = os.path.join(base, "untitled folder 120", "M462_CERTIFICATE_RESEND.json")
    if os.path.exists(f1):
        d = json.load(open(f1))
        assert d["phi_h_exact_enclosure"] == PREMISES["phi_h_exact_enclosure"]
        assert [d["constant_enclosures"]["h_lin"]["lo"], d["constant_enclosures"]["h_lin"]["hi"]] == PREMISES["h_lin"]
        assert [d["constant_enclosures"]["c_star"]["lo"], d["constant_enclosures"]["c_star"]["hi"]] == PREMISES["c_star"]
        out["contract_sha256"] = hashlib.sha256(open(f1, "rb").read()).hexdigest()
    if os.path.exists(f2):
        d = json.load(open(f2))
        assert d["box"]["center"][0] == PREMISES["eta_center"]
        assert [d["K_minus_center"][0]["lo"], d["K_minus_center"][0]["hi"]] == PREMISES["eta_K_minus_center"]
        out["m462_certificate_sha256"] = hashlib.sha256(open(f2, "rb").read()).hexdigest()
    return out


GATES = []


def gate(name, ok, detail=""):
    GATES.append({"gate": name, "pass": bool(ok), "detail": detail})
    if not ok:
        raise SystemExit("GATE FAILED: %s  %s" % (name, detail))


def ivs(x, digits=30):
    return cc.iv_str(x.bound(), digits)


# ------------------------------------------------------------------ parameter box
p_lo, p_hi = [Fr(x) for x in PREMISES["phi_h_exact_enclosure"]]
c_mid = cc.ffloor((p_lo + p_hi) / 2)
c_rad = cc.fceil((p_hi - p_lo) / 2) + 2
gate("PH covered by the Taylor-model parameter range", Fr(c_mid - c_rad, ONE) <= p_lo and p_hi <= Fr(c_mid + c_rad, ONE))
P_TM = TM([c_mid, c_rad])                     # p(u) = c_mid + c_rad u  (exact), u in [-1,1]
t_par = Fr(PREMISES["t"])
eta_lo = Fr(PREMISES["eta_center"]) + Fr(PREMISES["eta_K_minus_center"][0])
eta_hi = Fr(PREMISES["eta_center"]) + Fr(PREMISES["eta_K_minus_center"][1])
h_lo = Fr(PREMISES["h_lin"][0]) * t_par * (1 + t_par * eta_lo)
h_hi = Fr(PREMISES["h_lin"][1]) * t_par * (1 + t_par * eta_hi)
HB = IV.frac(h_lo, h_hi)

# ------------------------------------------------------------------ 1. cone contraction (exact rationals)
Ucoef = [Fr(-1, 6), Fr(4, 3), Fr(-5, 3), Fr(-4, 9), Fr(17, 18), Fr(0), Fr(-2, 27)]     # U(phi) ascending


def poly_shift(c, a):
    """coefficients of c(a + v) in v"""
    out = [Fr(0)] * len(c)
    from math import comb
    for k, ck in enumerate(c):
        for j in range(k + 1):
            out[j] += ck * comb(k, j) * a ** (k - j)
    return out


def poly_der(c):
    return [k * c[k] for k in range(1, len(c))]


def majorant(c, x):
    return sum(abs(ck) * x ** k for k, ck in enumerate(c))


# U = W_phi^2/2 - (2/3) W^2 check
def pmul(a, b):
    out = [Fr(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


Wc = [Fr(1), Fr(-1), Fr(0), Fr(1, 3)]
W1c = [Fr(-1), Fr(0), Fr(1)]
Uchk = [Fr(1, 2) * a for a in pmul(W1c, W1c)] + [Fr(0)] * 2
for i, a in enumerate(pmul(Wc, Wc)):
    Uchk[i] -= Fr(2, 3) * a
gate("U polynomial equals W_phi^2/2 - (2/3)W^2", Uchk == Ucoef)
U1chk = pmul(W1c, [Fr(-4, 3), Fr(10, 3), Fr(0), Fr(-4, 9)])
gate("U_phi = (phi^2-1)(-4/3 + 10 phi/3 - 4 phi^3/9)", U1chk == poly_der(Ucoef))

U2chk = [a + b for a, b in zip(pmul([Fr(0), Fr(2)], [Fr(-4, 3), Fr(10, 3), Fr(0), Fr(-4, 9)]) , pmul(W1c, [Fr(10, 3), Fr(0), Fr(-4, 3)]))]
gate("U_phiphi = 2 phi Bq + (phi^2-1)(10/3 - 4 phi^2/3)  (used in the a-priori G equation)", U2chk == poly_der(poly_der(Ucoef)))
A0 = poly_shift(Ucoef, Fr(-1)); A1 = poly_der(A0); A2 = poly_der(A1)
RC = Fr(1, 2); V0 = Fr(1, 10 ** 6)
r_be, r_s, r_ph = Fr(6, 100), Fr(2, 10 ** 6), Fr(11, 10 ** 7)
gate("|p+1| <= v0 on the parameter range", abs(Fr(c_mid - c_rad, ONE) + 1) <= V0 and abs(Fr(c_mid + c_rad, ONE) + 1) <= V0)
vs = V0 + r_ph
selfmap = [RC / 3 * (r_be ** 2 + r_s ** 2 / 4 + majorant(A0, vs) / 6), RC / 5 * (majorant(A1, vs) + 4 * r_be * r_s), RC * r_s]
gate("cone map sends the ball into itself", selfmap[0] < r_be and selfmap[1] < r_s and selfmap[2] < r_ph,
     str([float(x) for x in selfmap]))
lip = [RC / 3 * (2 * r_be * r_be + r_s * r_s / 2 + majorant(A1, vs) * r_ph / 6) / r_be,
       RC / 5 * (4 * r_s * r_be + 4 * r_be * r_s + majorant(A2, vs) * r_ph) / r_s,
       RC * r_s / r_ph]
gate("cone map is a contraction (weighted max norm)", max(lip) < 1, str([float(x) for x in lip]))
CONTRACTION = max(lip)

# ------------------------------------------------------------------ cone Taylor coefficients as Taylor models (scaled by Y0^n)
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
tail = lambda r: cc.fceil(r * xq ** (NC + 1) / (1 - xq))
def tm_sum(lst, r):
    acc = lst[0]
    for a in lst[1:]:
        acc = acc + a
    tl = tail(r)
    return TM(acc.c, acc.rlo - tl, acc.rhi + tl)

be0 = tm_sum(be, r_be); s0 = tm_sum(s_, r_s); ph0 = tm_sum(ph, r_ph)
H0 = be0.addc(1 / Y0)
def U_of(x):       # generic (TM or IV)
    x2 = x * x; x3 = x2 * x
    Wv = (x3.scale(1, 3) - x).addc(1); W1v = x2.addc(-1)
    return (W1v * W1v).scale(1, 2) - (Wv * Wv).scale(2, 3)
def U1_of(x):
    x2 = x * x; x3 = x2 * x
    return x2.addc(-1) * (x.scale(10, 3) - x3.scale(4, 9)).addc(Fr(-4, 3))
w0 = H0 * H0 - (s0 * s0).scale(1, 12) + U_of(ph0).scale(1, 6)
gate("s > 0 at y0", s0.bound().lo > 0)

# ------------------------------------------------------------------ 2. Riccati barrier on (0, y0]
tauI = IV(0, ONE)
def horner(lst, start=0):
    acc = IV(0, 0)
    for a in reversed(lst[start:]):
        acc = acc * tauI + a.bound()
    return acc
tl = lambda r: IV(-tail(r), tail(r))
beI = horner(be) + tl(r_be)                                   # beta(y), y in [0,y0]
phI = horner(ph) + tl(r_ph)
sgI = (horner(s_, 1) + tl(r_s)).scale(Y0.denominator, Y0.numerator)   # sigma = s/y  (tail bound divided by y0 as well)
gate("sigma = s/y > 0 on (0,y0]", sgI.lo > 0)
yI = IV(0, cc.fceil(Y0))
ybe = yI * beI
y2 = yI * yI
sI = yI * sgI
y2s2 = y2 * (sI * sI)
yG = U1_of(phI) * sgI.recip()
c1I = (yG.scale(2) - ybe.scale(6)).addc(-5)
y2w = (ybe.addc(1) * ybe.addc(1)) - y2s2.scale(1, 12) + (y2 * U_of(phI)).scale(1, 6)
BARRIER = {}
R_START = []
for m in MU2S:
    m4 = m + 4
    c0I = y2w.scale(m4.numerator, m4.denominator) - y2s2.scale(2, 3)
    # floating guidance for the barrier points
    import math
    c1l, c1h = c1I.fl(); c0l, c0h = c0I.fl()
    rl = (-c1h - math.sqrt(c1h * c1h - 4 * c0l)) / 2
    rh = (-c1l - math.sqrt(c1l * c1l - 4 * c0h)) / 2
    r_lo = Fr(math.floor((rl - 1e-3) * 2 ** 20), 2 ** 20)
    r_hi = Fr(math.ceil((rh + 1e-3) * 2 ** 20), 2 ** 20)
    q_lo = IV.frac(r_lo * r_lo) + c1I.scale(r_lo.numerator, r_lo.denominator) + c0I
    q_hi = IV.frac(r_hi * r_hi) + c1I.scale(r_hi.numerator, r_hi.denominator) + c0I
    gate("barrier mu2=%s: y r' > 0 at r_lo on (0,y0]" % m, q_lo.lo > 0)
    gate("barrier mu2=%s: y r' < 0 at r_hi on (0,y0]" % m, q_hi.hi < 0)
    Q0 = lambda r: r * r + 5 * r + m + 4
    gate("barrier mu2=%s: r_lo < r_+ < r_hi < -5/2" % m, Q0(r_lo) > 0 and Q0(r_hi) < 0 and r_lo < r_hi < Fr(-5, 2))
    BARRIER[str(m)] = {"r_lo": str(r_lo), "r_hi": str(r_hi), "c1_on_(0,y0]": cc.iv_str(c1I, 12), "c0_on_(0,y0]": cc.iv_str(c0I, 12),
                       "r_plus_float": -2.5 - math.sqrt(2.25 - float(m))}
    R_START.append(r_lo / Y0); R_START.append(r_hi / Y0)
MU_LIST = [m for m in MU2S for _ in (0, 1)]           # R_{2j} = lower solution, R_{2j+1} = upper solution for MU2S[j]

state = [ph0, s0, H0, w0] + [TM([cc.ffloor(r)], 0, 1) for r in R_START]
CHECKPOINTS = [{"y": str(Y0), "phi": ivs(ph0), "s": ivs(s0), "H": ivs(H0), "w": ivs(w0)}]
print("cone done  %.1fs   contraction %.4f   s(y0) in %s" % (time.time() - T0, float(CONTRACTION), ivs(s0, 20)), flush=True)

# ------------------------------------------------------------------ 3. validated integration to Y1
y = Y0; nstep = 0; max_rem = 0
def remainder_width(st):
    return max((x.rhi - x.rlo) for x in st)
TOL = ONE >> 140                       # accepted scaled Lagrange remainder per step (about 7e-43)
h_prev = Y0 / 64
def adaptive_step(state, h, tube=False):
    """validated step with order control; halves h when no a-priori box exists or the remainder is too large"""
    global NORD
    while True:
        try:
            out = cc.step(state, MU_LIST, h, NORD, tau_interval=tube)
        except ArithmeticError:
            out = None
        if out is not None and out[2] <= TOL:
            if out[2] < (TOL >> 40) and NORD > 12:
                NORD -= 2
            return out, h
        if out is not None and NORD + 4 <= NORD_MAX:
            NORD += 4
            continue
        if tube:
            raise SystemExit("tube step failed")
        h = h / 2                          # exact dyadic rational
        if h < Fr(1, 2 ** 30):
            raise SystemExit("step size underflow")
while y < Y1:
    h = min(KAPPA * y, HMAX, 2 * h_prev)
    h = Fr(int(h * 2 ** 16), 2 ** 16)
    if y + h > Y1:
        h = Y1 - y
    (state, _, rem), h = adaptive_step(state, h)
    h_prev = h
    max_rem = max(max_rem, rem)
    y += h; nstep += 1
    if nstep % 10 == 0 or y == Y1:
        CHECKPOINTS.append({"y": str(y), "phi": ivs(state[0]), "s": ivs(state[1]), "H": ivs(state[2]), "w": ivs(state[3]),
                            "R": [ivs(x, 16) for x in state[4:]],
                            "max_TM_remainder_width": "%.3e" % (remainder_width(state) / ONE)})
        print("step %3d  y=%.5f h=%.4f N=%d  phi+1=%.6e  R=%.4f remainder width %.2e  local trunc %.1e   %.0fs" %
              (nstep, float(y), float(h), NORD, state[0].c[0] / ONE + 1, state[4].c[0] / ONE, remainder_width(state) / ONE, rem / ONE, time.time() - T0), flush=True)

# ------------------------------------------------------------------ 4. locate the shell (floating guidance, then rigorous gates)
cen = [IV(x.c[0], x.c[0]) for x in state]
hh = Fr(1, 16)
ser = cc.taylor(cen, MU_LIST, hh, 30)[3]
wc = [a.lo / ONE for a in ser]
target = float((h_lo + h_hi) / 2)
tau = 0.4
for _ in range(60):
    f = sum(c * tau ** k for k, c in enumerate(wc)) - target
    df = sum(k * c * tau ** (k - 1) for k, c in enumerate(wc) if k)
    tau -= f / df
yb_guess = Y1 + Fr(tau) * hh
ya = Fr(int((yb_guess - MARGIN) * 2 ** 40), 2 ** 40)
yb = ya + 2 * MARGIN
ycur = Y1; state_a = state
while ycur < ya:
    h = min(HMAX, 2 * h_prev, ya - ycur)
    (state_a, _, rem), h = adaptive_step(state_a, h)
    h_prev = h; ycur += h; nstep += 1; max_rem = max(max_rem, rem)
(state_b, tube, rem2), _ = adaptive_step(state_a, yb - ya, tube=True)
gate("w(ya) > sup HB for all p in PH", state_a[3].bound().lo > HB.hi)
gate("w(yb) < inf HB for all p in PH", state_b[3].bound().hi < HB.lo)
gate("H > 0 and w > 0 on the tube (w strictly decreasing)", tube[2].bound().lo > 0 and tube[3].bound().lo > 0)

# ------------------------------------------------------------------ 5. sign of the junction mismatch
gate("H > 0 at yb (H' = -w - s^2/3 < 0, so H > 0 and w' = -2Hw < 0 on all of (0, yb])", state_b[2].bound().lo > 0)
# 5a. coarse version: enclosure on the whole tube y in [ya, yb] (no use of the relation between y_b and phi_h)
phi_t, s_t, H_t, w_t = tube[:4]
gate("s > 0 on the tube", s_t.bound().lo > 0)
G_t = U1_of(phi_t) * s_t.recip()
coefI = (G_t - H_t.scale(4) + phi_t.scale(2)).bound()
COARSE = {}
for j, m in enumerate(MU2S):
    Rl = tube[4 + 2 * j].bound(); Rh = tube[5 + 2 * j].bound()
    RI = IV(min(Rl.lo, Rh.lo), max(Rl.hi, Rh.hi))
    m4 = m + 4
    bI = w_t.bound().scale(m4.numerator, m4.denominator) + coefI * RI
    COARSE[str(m)] = {"b_enclosure_on_tube": cc.iv_str(bI, 20),
                      "sign": "negative" if bI.hi < 0 else ("positive" if bI.lo > 0 else "undetermined")}
# 5b. sharp version: evaluation AT the shell point w(y_b) = h (Taylor-model-valued y_b(phi_h), verified)
shell, eps_tau, shinfo = cc.shell_point(state_a, MU_LIST, yb - ya, NORD, HB)
phi_s, s_s, H_s, w_s = shell[:4]
gate("s > 0 at the shell", s_s.bound().lo > 0)
gate("sanity: enclosure of w at the shell point meets HB and exceeds it by < 1e-6 relative",
     w_s.bound().lo <= HB.hi and w_s.bound().hi >= HB.lo and w_s.bound().lo > HB.lo - (HB.lo >> 20) and w_s.bound().hi < HB.hi + (HB.hi >> 20))
G_s = U1_of(phi_s) * s_s.recip()
coef_s = G_s - H_s.scale(4) + phi_s.scale(2)
gate("junction coefficient G - 4H + 2 phi < 0 at the shell (b decreasing in R)", coef_s.bound().hi < 0)
RESULT = {}
for j, m in enumerate(MU2S):
    Rl = shell[4 + 2 * j]; Rh = shell[5 + 2 * j]
    gate("ordering of the Riccati sandwich mu2=%s" % m, Rl.bound().lo <= Rh.bound().hi)
    m4 = m + 4
    b_up = (w_s.scale(m4.numerator, m4.denominator) + coef_s * Rl).bound().hi      # R >= R_lower-solution, coef < 0
    b_dn = (w_s.scale(m4.numerator, m4.denominator) + coef_s * Rh).bound().lo
    bI = IV(b_dn, b_up)
    sign = "negative" if bI.hi < 0 else ("positive" if bI.lo > 0 else "undetermined")
    RESULT[str(m)] = {"b_enclosure_at_shell": cc.iv_str(bI, 20), "sign": sign, "mu2_float": float(m),
                      "R_lower_solution_at_shell": ivs(Rl, 20), "R_upper_solution_at_shell": ivs(Rh, 20),
                      "coarse_tube_version": COARSE[str(m)]}
    print("mu2 = %s :  b in %s   (%s)   [coarse tube: %s]" % (m, cc.iv_str(bI, 12), sign, COARSE[str(m)]["sign"]), flush=True)
signs = [RESULT[str(m)]["sign"] for m in MU2S]
gate("main pair: signs determined and opposite", {signs[0], signs[1]} == {"negative", "positive"})
gate("all mu^2 values are negative", all(m < 0 for m in MU2S))
neg = [m for m in MU2S if RESULT[str(m)]["sign"] == "negative"]
pos = [m for m in MU2S if RESULT[str(m)]["sign"] == "positive"]
gate("negative-sign values lie below positive-sign values", max(neg) < min(pos))
BRACKET = [max(neg), min(pos)]
print("certified eigenvalue bracket: %s < mu^2 < %s   (%.7f, %.7f)" % (BRACKET[0], BRACKET[1], float(BRACKET[0]), float(BRACKET[1])))
G_t = G_s; phi_t, s_t, H_t, w_t = phi_s, s_s, H_s, w_s          # diagnostic below is evaluated at the shell point

# consistency diagnostic (not needed for the proof): the background junction residuals must contain zero on the tube
cs = TM.const_iv(IV.frac(Fr(PREMISES["c_star"][0]), Fr(PREMISES["c_star"][1])))
tT = TM.const_iv(IV.frac(t_par))
Wt = ((phi_t * phi_t * phi_t).scale(1, 3) - phi_t).addc(1)
sig = Wt.scale(2) + tT * (cs * phi_t).addc(1)
dsig = (phi_t * phi_t).addc(-1).scale(2) + tT * cs
F1 = (H_t - sig.scale(1, 6)).bound(); F2 = (s_t + dsig.scale(1, 2)).bound()
cons = (H_t * H_t - w_t - (s_t * s_t).scale(1, 12) + U_of(phi_t).scale(1, 6)).bound()
gate("integrator sanity: the constraint H^2 - w - s^2/12 + U/6 = 0 is inside its enclosure at the shell", cons.lo <= 0 <= cons.hi)
DIAG = {"constraint_residual_at_shell_point": cc.iv_str(cons, 6), "H - sigma_t/6 at shell point": cc.iv_str(F1, 12), "s + sigma_t'/2 at shell point": cc.iv_str(F2, 12),
        "both_contain_zero": F1.lo <= 0 <= F1.hi and F2.lo <= 0 <= F2.hi}

cert = {
    "artifact": "HDBLAST Chat 9 certified tachyonic scalar bound state of the registered de Sitter shell",
    "status": "PASS_SIGN_CHANGE_CERTIFIED_CONDITIONAL_ON_INHERITED_ROOT_ENCLOSURES",
    "statement": "For every phi_h in PH and every h in HB (boxes containing the M462-certified shell root), with y_b defined by w(y_b)=1/rho_b^2=h: "
                 "y_b in [ya,yb]; the Frobenius-regular scalar mode has psi>0 on (0,y_b]; the scalar-junction mismatch b has the listed strict signs "
                 "at the listed mu^2 values. Hence an eigenvalue mu^2 (m^2 = mu^2 H^2 < 0) exists with %s < mu^2 < %s." % (BRACKET[0], BRACKET[1]),
    "mu2_values": [str(m) for m in MU2S],
    "result": RESULT,
    "PH_covered": [cc.dec(c_mid - c_rad, 40), cc.dec(c_mid + c_rad, 40)],
    "HB": cc.iv_str(HB, 30),
    "shell_interval_Y": {"ya": str(ya), "yb": str(yb), "ya_float": float(ya), "yb_float": float(yb)},
    "certified_bracket": {"lower": str(BRACKET[0]), "upper": str(BRACKET[1]), "float": [float(BRACKET[0]), float(BRACKET[1])]},
    "shell_point_enclosures": {"phi_b": ivs(phi_s), "s_b": ivs(s_s), "H_b": ivs(H_s), "w_b": ivs(w_s), "G_b": ivs(G_s),
                               "G-4H+2phi": ivs(coef_s, 20), "y_b_in": [str(ya), str(yb)],
                               "tau_halfwidth_eps": str(eps_tau), "tau_center": shinfo["tau_center_coeffs"]},
    "tube_coefficient_G-4H+2phi": cc.iv_str(coefI, 20),
    "cone": {"y0": str(Y0), "series_degree": NC, "disc_radius": str(RC), "ball_radii_beta_s_phi": [str(r_be), str(r_s), str(r_ph)],
             "self_map_bounds": [str(x) for x in selfmap], "contraction_constant_upper": "%.6f" % float(CONTRACTION), "barrier": BARRIER},
    "integration": {"taylor_order_final": NORD, "taylor_order_max": NORD_MAX, "remainder_tolerance": "2^-140", "tm_degree": DEG, "fixed_point_bits": P, "steps_to_Y1": nstep, "Y1": str(Y1),
                    "max_scaled_Lagrange_remainder": "%.3e" % (max(max_rem, rem, rem2) / ONE),
                    "final_TM_remainder_width": "%.3e" % (remainder_width(state_b) / ONE)},
    "junction_consistency_diagnostic": DIAG,
    "checkpoints": CHECKPOINTS,
    "premises": PREMISES,
    "premise_file_binding": check_premise_files(),
    "gates": GATES,
    "not_proved_here": [
        "the inherited M462/M463 root enclosures and the I_plus-derived constants (premises), and the identification w(y_b) = 1/rho_b^2 = h of the M462 normalisation q = sqrt(h) rho",
        "identification of the contraction fixed point with 'the' regular cone solution uses the standard uniqueness of bounded solutions of the Volterra system (argued in the text, not machine-checked)",
        "uniqueness of the bound state / absence of other eigenvalues (floating argument-principle evidence only)",
        "continuity of the Frobenius-normalised solution in mu^2 (standard regular-singular-point theory, used for the intermediate value step)",
        "the derivation of the master system and of the junction condition (analytic, in stability/orchestrator_derivation; cross-checked in floating point only)",
        "positivity of the mode norm (tachyon rather than ghost)"],
    "source_sha256": {f: hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest() for f in ("cert_core.py", "certify_tachyon.py")},
    "runtime_seconds": round(time.time() - T0, 1),
}
name = ("TACHYON_CERTIFICATE" if len(sys.argv) == 1 else "TACHYON_CERTIFICATE_TIGHT") + os.environ.get("HDB_CERT_SUFFIX", "") + ".json"
json.dump(cert, open(os.path.join(HERE, name), "w"), indent=1)
print("ALL %d GATES PASSED.  wrote %s   (%.0fs)" % (len(GATES), name, time.time() - T0))
