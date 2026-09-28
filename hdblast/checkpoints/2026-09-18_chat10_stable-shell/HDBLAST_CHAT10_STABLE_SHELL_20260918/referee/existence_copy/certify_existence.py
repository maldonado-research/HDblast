#!/usr/bin/env python3
"""HDBLAST Chat 10 - computer-assisted EXISTENCE proof and enclosure of the stable de Sitter shell S_8/5.

Python 3.9 standard library only; integers / fractions with outward rounding (cert_core.py, 2^-384 fixed point).
Floats appear only in print-outs.  Run:   python3 certify_existence.py      (~10 s; any failed gate -> non-zero exit)

MODEL   phi' = s, s' = U1(phi) - 4Hs, H' = -w - s^2/3, w' = -2Hw  (H = rho'/rho, w = 1/rho^2), W = 1 - phi + phi^3/3,
        U = W'^2/2 - (2/3)W^2, regular cone at y = 0 (rho ~ y, phi(0) = p, phi'(0) = 0), y = proper distance.
SHELL   sigma_t = 2W + t(1 + c phi + d phi^2/2),  t = 1/1000, c = 5975949350280/10^13, d = 8/5   (EXACT rationals).
        Junction residuals at y:  G1 = H - sigma_t(phi)/6,   G2 = s + sigma_t'(phi)/2.

THEOREM CHECKED.  Let BOX = [p_lo, p_hi] x [y_a, y_b] (exact rationals in the certificate).  There is (p*, y*) in BOX
such that the regular-cone solution with phi(0) = p* exists on (0, y*], has s > 0, H > 0, w > 0 there, and satisfies
G1 = G2 = 0 at y = y*.  All shell quantities listed under "enclosures" are rigorous bounds valid for EVERY point of
BOX, hence for (p*, y*).

PROOF STEPS
  1. Cone chart (0, 1/20]: Chat 9 Volterra contraction in the disc algebra, uniform in p (exact rational gates);
     Taylor coefficients enclosed as Taylor models in u, p = p_mid + p_rad u, u in [-1,1].
  2. Validated Taylor integration to y_a for all p in [p_lo, p_hi] at once (Taylor model in u + interval remainder;
     Lagrange remainder on verified a-priori boxes).  One further validated step gives the state at y_b.
  3. p-faces: the Taylor model at y_a is evaluated at u = -1 / u = +1 (rigorous interval), and the solution on the whole
     segment y = y_a + (y_b - y_a)(1+v)/2, v in [-1,1], is enclosed as a Taylor model in v (Taylor step polynomial
     composed with tau = (1+v)/2, Lagrange remainder over a verified a-priori box).
  4. Preconditioned Miranda test.  With the exact rational matrix C (approximately J^-1, det C != 0 gated) put
     Gt = C (G1, G2)^T.  Gates:  Gt_1 < 0 on {p = p_lo} x [y_a,y_b],  Gt_1 > 0 on {p = p_hi} x [y_a,y_b],
     Gt_2 < 0 on [p_lo,p_hi] x {y_a},  Gt_2 > 0 on [p_lo,p_hi] x {y_b}.   (p,y) -> G is continuous on BOX (continuous
     dependence of the contraction fixed point and of the ODE flow on p; the flow exists up to y_b for all p in the
     range by step 2).  Poincare-Miranda gives a zero of Gt in BOX; C invertible gives G1 = G2 = 0 there.
  5. Enclosures: validated tube over [y_a, y_b] for all p in the range.
NOT PROVED HERE: uniqueness of the root (globally or inside BOX); see "not_proved_here".
"""
import sys, os, json, time, hashlib
from fractions import Fraction as Fr
import cert_core as cc
from cert_core import IV, TM, ONE, P
import exist_lib as L
from exist_lib import gate

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
DEG = 2
RP_BITS = 100                 # p-radius 2^-100  ~ 7.9e-31
RY = Fr(1, 2 ** 84)           # y-radius 2^-84   ~ 5.2e-26
cc.set_degree(DEG)

guess = json.load(open(os.path.join(HERE, "CENTER_GUESS.json")))          # NON-RIGOROUS guidance only
p_c = Fr(guess["p_center"]); y_c = Fr(guess["y_center"])
C = [[Fr(x) for x in row] for row in guess["C"]]
detC = C[0][0] * C[1][1] - C[0][1] * C[1][0]
gate("preconditioner C is an exact rational matrix with det C != 0", detC != 0, "det C = %.6e" % float(detC))

c_mid = cc.ffloor(p_c); c_rad = ONE >> RP_BITS
NEG_SHIFT = int(os.environ.get("HDB_NEG_CONTROL_SHIFT", "0"))      # negative control: move the box off the root (gates must FAIL)
c_mid += NEG_SHIFT * c_rad
p_lo = Fr(c_mid - c_rad, ONE); p_hi = Fr(c_mid + c_rad, ONE)
P_TM = TM([c_mid, c_rad])                      # p(u) = (c_mid + c_rad u) / 2^P exactly
ya = Fr(int(y_c * 2 ** 160), 2 ** 160) - RY; yb = ya + 2 * RY

# 1. cone
state0, cone_info = L.cone_start(P_TM, p_lo, p_hi)
print("cone done %.1fs" % (time.time() - T0), flush=True)

# 2. integration to ya, then one step to yb (with tube)
I = L.Integrator(state0, L.Y0, verbose=True)
CHECK = []
for yt in [Fr(1), Fr(2), Fr(4), Fr(6), Fr(7), Fr(8), Fr(41, 5)]:
    st = I.advance_to(yt)
    CHECK.append({"y": str(yt), "phi": cc.iv_str(st[0].bound(), 30), "s": cc.iv_str(st[1].bound(), 30),
                  "H": cc.iv_str(st[2].bound(), 30), "w": cc.iv_str(st[3].bound(), 30)})
state_a = I.advance_to(ya)
state_b, tube = I.fixed_step(state_a, yb - ya, tube=True)
N_last = I.N
rw = lambda st: max(x.rhi - x.rlo for x in st)
print("integration done %.1fs  steps=%d  remainder width at ya %.2e" % (time.time() - T0, I.nstep, rw(state_a) / ONE), flush=True)
gate("s > 0 on every a-priori box of the integration (implicit: the recurrences divide by s) and on the tube", tube[1].bound().lo > 0)
gate("H > 0 on the tube [ya,yb] for all p (H' = -w - s^2/3 < 0, hence H > 0 on (0,yb])", tube[2].bound().lo > 0)
gate("w > 0 on the tube", tube[3].bound().lo > 0)


def Gt(st):
    g1, g2 = L.residuals(st)
    return [g1.scale(C[i][0].numerator, C[i][0].denominator) + g2.scale(C[i][1].numerator, C[i][1].denominator) for i in (0, 1)]


# 3./4. Miranda faces
FACES = {}
# y-faces: all p at once
gt_a = Gt(state_a)[1].bound(); gt_b = Gt(state_b)[1].bound()
FACES["Gt2 on y=ya (all p)"] = cc.iv_str(gt_a, 12); FACES["Gt2 on y=yb (all p)"] = cc.iv_str(gt_b, 12)
gate("Gt_2 < 0 on the face y = ya", gt_a.hi < 0, str(FACES["Gt2 on y=ya (all p)"]))
gate("Gt_2 > 0 on the face y = yb", gt_b.lo > 0, str(FACES["Gt2 on y=yb (all p)"]))
# p-faces: all y at once
seg_rem = 0
for sign, nm in ((-1, "p=p_lo"), (+1, "p=p_hi")):
    st_face = [TM.const_iv(L.tm_at(x, sign)) for x in state_a]
    seg, rem = L.eval_on_segment(st_face, yb - ya, max(N_last, 12))
    seg_rem = max(seg_rem, rem)
    g = Gt(seg)[0].bound()
    FACES["Gt1 on %s (all y)" % nm] = cc.iv_str(g, 12)
    if sign < 0:
        gate("Gt_1 < 0 on the face p = p_lo", g.hi < 0, str(cc.iv_str(g, 12)))
    else:
        gate("Gt_1 > 0 on the face p = p_hi", g.lo > 0, str(cc.iv_str(g, 12)))
    # consistency: the segment enclosure at v = +1 must intersect the independent validated step to yb
    for i in range(4):
        e1 = L.tm_at(seg[i], +1); e2 = L.tm_at(state_b[i], sign)
        gate("consistency segment vs step (%s, comp %d)" % (nm, i), e1.lo <= e2.hi and e2.lo <= e1.hi)
gate("segment Lagrange remainder below tolerance", seg_rem <= L.TOL)

# 5. enclosures over the whole BOX
phi_t, s_t, H_t, w_t = [x.bound() for x in tube]
tubeIV = [phi_t, s_t, H_t, w_t]
G1_t, G2_t = L.residuals(tubeIV)
gate("sanity: both residual enclosures over BOX contain 0", G1_t.lo <= 0 <= G1_t.hi and G2_t.lo <= 0 <= G2_t.hi)
cons = H_t * H_t - w_t - (s_t * s_t).scale(1, 12) + L.U_of(phi_t).scale(1, 6)
gate("sanity: constraint H^2 - w - s^2/12 + U/6 contains 0 over BOX", cons.lo <= 0 <= cons.hi)
G_t = L.U1_of(phi_t) * s_t.recip()
B_t = L.B_of(tubeIV)
gate("B = phi''/phi' + sigma_t''/2 > 0 on BOX", B_t.lo > 0, str(cc.iv_str(B_t, 20)))
rho_t = L.iv_sqrt(w_t).recip()
rhop_t = H_t * rho_t
sig_t = L.sigma_of(phi_t); dsig_t = L.dsigma_of(phi_t)
ddsig_t = phi_t.scale(4).addc(L.T_PAR * L.D_PAR)
dG1dy = (-w_t - (s_t * s_t).scale(1, 3)) - (dsig_t * s_t).scale(1, 6)
gate("dG1/dy = H' - sigma_t' s/6 < 0 on BOX (for each p at most one zero of G1 in [ya,yb])", dG1dy.hi < 0)
alpha_iv = IV(c_mid - c_rad + ONE, c_mid + c_rad + ONE)
ENC = {
    "alpha = phi_h + 1": cc.iv_str(alpha_iv, 40), "phi_h": cc.iv_str(IV(c_mid - c_rad, c_mid + c_rad), 40),
    "y_b": cc.iv_str(IV.frac(ya, yb), 40),
    "phi_b": cc.iv_str(phi_t, 40), "phi'_b = s_b": cc.iv_str(s_t, 40), "H_b = rho'_b/rho_b": cc.iv_str(H_t, 40),
    "h = w_b = 1/rho_b^2": cc.iv_str(w_t, 40), "rho_b": cc.iv_str(rho_t, 40), "rho'_b": cc.iv_str(rhop_t, 40),
    "G_b = phi''_b/phi'_b": cc.iv_str(G_t - H_t.scale(4), 30), "U_phi(phi_b)/phi'_b": cc.iv_str(G_t, 30),
    "B = phi''/phi' + sigma_t''/2": cc.iv_str(B_t, 30),
    "sigma_t(phi_b)": cc.iv_str(sig_t, 30), "sigma_t'(phi_b)": cc.iv_str(dsig_t, 30), "sigma_t''(phi_b)": cc.iv_str(ddsig_t, 30),
    "dG1/dy": cc.iv_str(dG1dy, 20),
}
for k, v in ENC.items():
    print("  %-32s %s" % (k, v))

cert = {
    "artifact": "HDBLAST Chat 10: certified existence and enclosure of the stable de Sitter shell S_8/5",
    "status": "PASS_EXISTENCE_CERTIFIED_BY_PRECONDITIONED_MIRANDA",
    "parameters_exact": {"t": str(L.T_PAR), "c": str(L.C_PAR), "d": str(L.D_PAR), "sigma_t": "2W + t(1 + c phi + d phi^2/2)", "W": "1 - phi + phi^3/3"},
    "statement": "There is (p*, y*) in BOX = [p_lo,p_hi] x [ya,yb] such that the regular-cone bulk solution with phi(0) = p* exists on (0,y*], "
                 "has phi' > 0, rho'/rho > 0 there, and satisfies both junction conditions rho'/rho = sigma_t/6, phi' = -sigma_t'/2 at y*. "
                 "Every enclosure below holds on all of BOX.",
    "BOX_exact": {"p_lo": str(p_lo), "p_hi": str(p_hi), "ya": str(ya), "yb": str(yb),
                  "p_mid_fixed_point_numerator_over_2^384": str(c_mid), "p_radius": "2^-%d" % RP_BITS, "y_radius": str(RY)},
    "enclosures": ENC,
    "miranda_faces": FACES,
    "preconditioner_C": guess["C"], "det_C_float": float(detC),
    "cone": cone_info,
    "integration": {"tm_degree": DEG, "fixed_point_bits": P, "steps": I.nstep, "taylor_order_final": I.N, "remainder_tolerance": "2^-140",
                    "max_scaled_Lagrange_remainder": "%.3e" % (max(I.max_rem, seg_rem) / ONE),
                    "TM_remainder_width_at_ya": "%.3e" % (rw(state_a) / ONE)},
    "checkpoints_all_p": CHECK,
    "guess_file_is_guidance_only": guess,
    "gates": L.GATES,
    "not_proved_here": [
        "uniqueness of the root, inside BOX or globally (floating Jacobian at the centre: det J = %.4e, non-degenerate; not interval-certified)" % (
            guess["J_float"][0][0] * guess["J_float"][1][1] - guess["J_float"][0][1] * guess["J_float"][1][0]),
        "identification of the contraction fixed point with 'the' regular cone solution uses the standard uniqueness of bounded solutions of the Volterra system (as in Chat 9)",
        "continuity of (p,y) -> (phi,s,H,w) on BOX is the standard continuous dependence theorem (uniform contraction on the cone chart + smooth ODE flow); used for Poincare-Miranda, not machine-checked",
        "that y* is the FIRST point along the flow where the junction conditions hold (not needed for existence)",
        "the a-priori-box search uses a modified inflation heuristic (exist_lib.apriori_box_v2) when the Chat 9 heuristic stalls; the verified inclusion test is unchanged"],
    "source_sha256": {f: hashlib.sha256(open(os.path.join(HERE, f), "rb").read()).hexdigest()
                      for f in ("cert_core.py", "exist_lib.py", "certify_existence.py", "CENTER_GUESS.json")},
    "runtime_seconds": round(time.time() - T0, 1),
}
name = "EXISTENCE_CERTIFICATE_S85" + ("_NEGCONTROL_SHOULD_NOT_EXIST" if NEG_SHIFT else "") + os.environ.get("HDB_CERT_SUFFIX", "") + ".json"
json.dump(cert, open(os.path.join(HERE, name), "w"), indent=1)
print("ALL %d GATES PASSED.  wrote %s  (%.0fs)" % (len(L.GATES), name, time.time() - T0))
