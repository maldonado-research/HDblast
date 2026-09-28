#!/usr/bin/env python3
"""Sweep 3 (methods), theme E: a validated-numerics pilot for certifying the static "+1 branch".

Question. What would a computer-assisted existence proof of the +1 branch start from, and how well
conditioned is it? This script does NOT prove existence of the nonlinear solution. It

  1. (EXACT, SymPy) verifies the explicit regular linear cone solution on the pure-AdS (phi=+1)
     background: eta = a * C_14^(2)(cosh u)/C_14^(2)(1), u = y/9, plus the growth exponents 14 and -18,
     and the leading-order junction coefficient eta_b = -(9c/64) delta;
  2. (INTERVAL, mpmath.iv, 60 digits) encloses the leading-order ("LO") truncated model
     [thin-brane metric junction + linear scalar germ + linearised scalar junction] at the six
     detunings of the 22 Sept checkpoint, and encloses the germ amplification G(u_b) at the
     checkpoint's own shell position;
  3. (NUMERICAL) compares the LO enclosures with the checkpoint's floating-point nonlinear solutions
     (PLUS_BRANCH_RESULTS.json, read-only, SHA-256 checked) and fits how the relative defect scales
     with delta;
  4. (CONDITIONAL planning numbers) digits needed by a direct phi-formulation, and the ratio by
     which the unwanted (cone-singular) linear solution is suppressed when integrating outward.

Pre-declared expectations (written before the first run):
  - E1: every exact identity holds and every wrong-formula control fails;
  - E2: every 100-digit point evaluation lies inside the corresponding interval enclosure;
  - E3: LO metric-only H^2 reproduces the checkpoint's `metric_only_H2` to relative 1e-12;
  - E4: the relative defects of eta_b, eta_h, H^2 and of the ratio eta_b/eta_h vs G(u_b) are O(delta),
        i.e. log-log slope in [0.8, 1.2] over the three smallest detunings (the first neglected terms
        are O(delta) relative).
Controls (must fail / be detected): wrong Gegenbauer degree, wrong index, wrong dimension coefficient,
1% perturbed c, wrong-sign scalar junction, and a containment self-test on a deliberately shifted point.

Runtime: a few seconds, one core. Output: certification_pilot.json next to this script.
"""
import hashlib
import json
import math
import os
import sys
import time

import mpmath
import sympy as sp
from mpmath import iv, mp

HERE = os.path.dirname(os.path.abspath(__file__))
CKPT = ("/home/user/unified-theory-maldonado/new-files/latest-work/"
        "HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/"
        "static_branch/PLUS_BRANCH_RESULTS.json")
CKPT_SHA = "b03e1c5569afb3567873fa990d16af9321b06e12a80ed89f9088a1bc00cd4fc7"
I_PLUS = "1.0357712571566784"          # c = 2/I_+ - 4/3 (registered definition)
N_DEG, ALPHA = 14, 2                   # Gegenbauer degree/index derived below
IV_DPS, MP_DPS = 60, 100
U0_GERM = mpmath.mpf("0.25")           # interface used by the analytic_structure Frobenius series

t_start = time.time()
out = {"script": os.path.basename(__file__), "status_labels": {}, "controls": {}, "expectations": {}}
with open(__file__, "rb") as fh:
    out["script_sha256"] = hashlib.sha256(fh.read()).hexdigest()

# ---------------------------------------------------------------- input
raw = open(CKPT, "rb").read()
sha = hashlib.sha256(raw).hexdigest()
if sha != CKPT_SHA:
    sys.exit(f"checkpoint JSON hash mismatch: {sha}")
rows = json.loads(raw)["rows"]
out["input"] = {"path": CKPT, "sha256": sha, "n_rows": len(rows)}

# ---------------------------------------------------------------- 1. exact checks (SymPy)
ph, x, t, s = sp.symbols("phi x t s")
W = 1 - ph + ph**3 / 3
Wp = sp.diff(W, ph)
U = Wp**2 / 2 - sp.Rational(2, 3) * W**2
k = W.subs(ph, 1) / 3                                   # A' = -W/3 -> AdS scale k = 1/9
m2_over_k2 = sp.simplify(sp.diff(U, ph, 2).subs(ph, 1) / k**2)
exact = {}
exact["k_equals_1_over_9"] = bool(k == sp.Rational(1, 9))
exact["Upp1_over_k2_equals_252"] = bool(m2_over_k2 == 252)


def ode_residual(P, coef=5, rhs=252):
    """(x^2-1)P'' + coef*x P' - rhs*P, i.e. eta_uu + 4 coth(u) eta_u - 252 eta with x=cosh u (coef=5)."""
    return sp.expand((x**2 - 1) * sp.diff(P, x, 2) + coef * x * sp.diff(P, x) - rhs * P)


C14 = sp.expand(sp.gegenbauer(N_DEG, ALPHA, x))
exact["C14_2_solves_linear_eq"] = bool(ode_residual(C14) == 0)
exact["C14_2_at_1_equals_680"] = bool(C14.subs(x, 1) == 680)
lhs = sp.expand(sp.cancel(C14.subs(x, (t + 1 / t) / 2) * t**14))
rhs_exp = sp.expand(sum((j + 1) * (15 - j) * t**(28 - 2 * j) for j in range(15)))
exact["exponential_form_identity"] = bool(sp.expand(lhs - rhs_exp) == 0)
roots = sorted(sp.solve(sp.Eq(s * (s - 1) + 5 * s - 252, 0), s))
exact["growth_exponents_at_infinity"] = [int(r) for r in roots]
exact["growth_exponents_are_minus18_and_14"] = roots == [-18, 14]
# leading-order scalar junction coefficient: eta_b (k*14 + 2) = -delta c / 2  -> -9 c delta / 64
cc, dd = sp.symbols("c delta")
eta_lo = sp.simplify(-dd * cc / 2 / (k * 14 + 2))
exact["eta_b_leading_equals_minus_9c_over_64_delta"] = bool(sp.simplify(eta_lo + sp.Rational(9, 64) * cc * dd) == 0)
# Where the evolution chart ends (methods observation for P3): with dz = dy/rho on the phi=+1 AdS
# background rho = sinh(k y)/k, z = ln tanh(k y / 2) + const -> -infinity at the regular cone y = 0,
# and the 2D part e^{2B}(-dt^2+dz^2) = -rho^2 dt^2 + dy^2 with rho/y -> 1: a Rindler-type horizon that
# the conformal chart reaches only as z -> -infinity (so any finite far boundary z = -L is artificial).
yy = sp.symbols("y", positive=True)
rho_ads = sp.sinh(k * yy) / k
z_of_y = sp.log(sp.tanh(k * yy / 2))
exact["dz_dy_equals_1_over_rho"] = bool(sp.simplify(sp.diff(z_of_y, yy) - 1 / rho_ads) == 0)
exact["z_to_minus_infinity_at_cone"] = bool(sp.limit(z_of_y, yy, 0, "+") == -sp.oo)
exact["rho_over_y_to_1_at_cone"] = bool(sp.limit(rho_ads / yy, yy, 0, "+") == 1)
# wrong-formula controls (each must give a NONZERO residual)
ctrl = {}
ctrl["wrong_degree_13_fails"] = ode_residual(sp.gegenbauer(13, ALPHA, x)) != 0
ctrl["wrong_degree_15_fails"] = ode_residual(sp.gegenbauer(15, ALPHA, x)) != 0
ctrl["wrong_index_3half_fails"] = ode_residual(sp.gegenbauer(N_DEG, sp.Rational(3, 2), x)) != 0
ctrl["wrong_dimension_coef4_fails"] = ode_residual(C14, coef=4) != 0     # 4D-like 3*coth(u) damping
out["exact"] = {kk: (bool(v) if isinstance(v, (bool,)) else v) for kk, v in exact.items()}
out["status_labels"]["exact"] = "EXACT (SymPy %s)" % sp.__version__
exact_ok = all(v for kk, v in exact.items() if isinstance(v, bool))

# ---------------------------------------------------------------- 2. interval enclosures (mpmath.iv)
iv.dps = IV_DPS
mp.dps = MP_DPS
W_J = [(j + 1) * (15 - j) for j in range(15)]          # C_14^(2)(cosh u) = sum W_J[j] e^{(14-2j)u}
kI = iv.mpf(1) / 9
cI = 2 / iv.mpf(I_PLUS) - iv.mpf(4) / 3
kP = mpmath.mpf(1) / 9
cP = 2 / mpmath.mpf(I_PLUS) - mpmath.mpf(4) / 3


def G_iv(u, deg_weights=W_J, deg=14, norm=680):
    return sum(iv.mpf(w) * iv.exp((deg - 2 * j) * u) for j, w in enumerate(deg_weights)) / norm


def Gu_iv(u):
    return sum(iv.mpf(w * (14 - 2 * j)) * iv.exp((14 - 2 * j) * u) for j, w in enumerate(W_J)) / 680


def G_mp(u):
    return sum(mpmath.mpf(w) * mpmath.exp((14 - 2 * j) * u) for j, w in enumerate(W_J)) / 680


def Gu_mp(u):
    return sum(mpmath.mpf(w * (14 - 2 * j)) * mpmath.exp((14 - 2 * j) * u) for j, w in enumerate(W_J)) / 680


def lo_model_iv(delta_str, c_val, junction_sign=+1):
    d = iv.mpf(delta_str)
    q = 1 + iv.mpf(3) / 2 * d * (1 + c_val)              # coth(u_b) from k coth u_b = sigma(1)/6
    u_b = iv.log((q + 1) / (q - 1)) / 2
    H2 = kI**2 * (q**2 - 1)
    G, Gu = G_iv(u_b), Gu_iv(u_b)
    # linearised scalar junction  eta_y = -sigma'/2 = -(2 eta_b) - delta c/2   (junction_sign=-1: wrong sign)
    a = -junction_sign * d * c_val / 2 / (kI * Gu + junction_sign * 2 * G)
    return {"u_b": u_b, "H2": H2, "eta_h": a, "eta_b": a * G, "G": G}


def lo_model_mp(delta_str):
    d = mpmath.mpf(delta_str)
    q = 1 + mpmath.mpf(3) / 2 * d * (1 + cP)
    u_b = mpmath.log((q + 1) / (q - 1)) / 2
    G, Gu = G_mp(u_b), Gu_mp(u_b)
    a = -d * cP / 2 / (kP * Gu + 2 * G)
    return {"u_b": u_b, "H2": kP**2 * (q**2 - 1), "eta_h": a, "eta_b": a * G, "G": G}


def lo(I):
    return mp.mpf(I._mpi_[0])          # exact endpoint (mp precision > iv precision)


def hi(I):
    return mp.mpf(I._mpi_[1])


def contains(I, p):
    return bool(lo(I) <= p <= hi(I))


def mid(I):
    return (lo(I) + hi(I)) / 2


def width(I):
    return float(hi(I) - lo(I))


def rel(a, b):
    return float((mpmath.mpf(a) - mpmath.mpf(b)) / mpmath.mpf(b))


table = []
containment = []
for r in rows:
    dstr = repr(r["delta"])
    L = lo_model_iv(dstr, cI)
    P = lo_model_mp(dstr)
    for key in ("u_b", "H2", "eta_h", "eta_b", "G"):
        containment.append(contains(L[key], P[key]))
    u_ck = mpmath.mpf(r["y_b"]) / 9                        # checkpoint's own shell position (point input)
    Gck = G_iv(iv.mpf(u_ck))
    ratio_ck = mpmath.mpf(r["eta_b"]) / mpmath.mpf(r["eta_h"])
    row = {
        "delta": r["delta"],
        "LO_u_b": float(mid(L["u_b"])), "checkpoint_u_b": float(u_ck),
        "LO_H2": float(mid(L["H2"])), "LO_H2_width": width(L["H2"]),
        "LO_eta_b": float(mid(L["eta_b"])), "LO_eta_b_width": width(L["eta_b"]),
        "LO_eta_h": float(mid(L["eta_h"])), "LO_eta_h_width": width(L["eta_h"]),
        "G_at_checkpoint_u_b": float(mid(Gck)), "G_at_checkpoint_u_b_relwidth": width(Gck) / float(mid(Gck)),
        "checkpoint_ratio_eta_b_over_eta_h": float(ratio_ck),
        "reldev_ratio_vs_G": float(ratio_ck / mid(Gck) - 1),
        "reldev_eta_b_LO_vs_checkpoint": rel(mid(L["eta_b"]), r["eta_b"]),
        "reldev_eta_h_LO_vs_checkpoint": rel(mid(L["eta_h"]), r["eta_h"]),
        "reldev_H2_LO_vs_checkpoint_H2": rel(mid(L["H2"]), r["H2"]),
        "reldev_H2_LO_vs_checkpoint_metric_only_H2": rel(mid(L["H2"]), r["metric_only_H2"]),
        "reldev_eta_b_LO_vs_checkpoint_eta_finite_curvature_linear": rel(mid(L["eta_b"]), r["eta_finite_curvature_linear"]),
        "reldev_eta_b_refined_vs_default_checkpoint": rel(r["refined"]["eta_b"], r["eta_b"]),
        "reldev_eta_h_refined_vs_default_checkpoint": rel(r["refined"]["eta_h"], r["eta_h"]),
    }
    table.append(row)
out["rows"] = table
out["status_labels"]["enclosures"] = ("INTERVAL (mpmath.iv %s digits): encloses the closed-form LO truncated model only, "
                                     "not the nonlinear solution" % IV_DPS)
out["status_labels"]["comparisons"] = "NUMERICAL (checkpoint values are floating-point, rtol 2e-12 and refined 8e-14)"


def slope(keys):
    small = sorted(table, key=lambda z: z["delta"])[:3]
    xs = [math.log(z["delta"]) for z in small]
    ys = [math.log(abs(z[keys])) for z in small]
    n = len(xs); mx = sum(xs) / n; my = sum(ys) / n
    return sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / sum((a - mx) ** 2 for a in xs)


slopes = {kk: slope(kk) for kk in ("reldev_ratio_vs_G", "reldev_eta_b_LO_vs_checkpoint",
                                    "reldev_eta_h_LO_vs_checkpoint", "reldev_H2_LO_vs_checkpoint_H2")}
coeffs = {kk: [z[kk] / z["delta"] for z in sorted(table, key=lambda z: z["delta"])] for kk in slopes}
out["scaling"] = {"loglog_slope_three_smallest_delta": slopes, "reldev_over_delta_by_increasing_delta": coeffs}

# ---------------------------------------------------------------- 3. controls on the comparison
reg = [z for z in rows if abs(z["delta"] - 1e-3) < 1e-15][0]
reg_row = [z for z in table if abs(z["delta"] - 1e-3) < 1e-15][0]
# C_a: wrong degree 13 in the amplification factor (normalised to 1 at u=0)
W13 = [int(v) for v in sp.Poly(sp.expand(sp.gegenbauer(13, ALPHA, (t + 1 / t) / 2) * t**13), t).all_coeffs()]
W13 = [w for w in W13]  # coefficients of t^26 ... t^0 (odd degree -> only even powers t^(26-2j))
W13 = W13[0::2]
G13 = sum(iv.mpf(w) * iv.exp((13 - 2 * j) * iv.mpf(mpmath.mpf(reg["y_b"]) / 9)) for j, w in enumerate(W13)) / int(sp.gegenbauer(13, ALPHA, 1))
ctrl["wrong_degree_13_amplification_misses_by_over_90pct"] = abs(float(reg_row["checkpoint_ratio_eta_b_over_eta_h"] / mid(G13)) - 1) > 0.9
# C_b: c perturbed by +1%
small = sorted(rows, key=lambda z: z["delta"])[:3]
pert_ok = []
for r in small:
    base = abs([z for z in table if z["delta"] == r["delta"]][0]["reldev_eta_b_LO_vs_checkpoint"])
    Lp = lo_model_iv(repr(r["delta"]), cI * iv.mpf("1.01"))
    pert_ok.append(abs(rel(mid(Lp["eta_b"]), r["eta_b"])) > 10 * base)
ctrl["c_perturbed_1pct_detected_at_three_smallest_delta"] = all(pert_ok)
# C_c: wrong-sign scalar junction (eta_y = +sigma'/2).  DESIGN NOTE: the first version of this control
# expected a sign flip of eta_b.  That expectation was wrong: with k*G_u/G ~ 14/9 < 2 the wrong-sign
# denominator (k G_u - 2 G) is also negative, so eta_b keeps its sign and is ~8x too large instead
# (-(9/8) c delta vs -(9/64) c delta).  The first run recorded the sign criterion as not met; the
# criterion is now "relative deviation from the checkpoint exceeds 100%", and both are reported.
Lw = lo_model_iv(repr(reg["delta"]), cI, junction_sign=-1)
wrong_sign_reldev = rel(mid(Lw["eta_b"]), reg["eta_b"])
ctrl["wrong_sign_scalar_junction_detected_reldev_over_100pct"] = abs(wrong_sign_reldev) > 1.0
out["control_design_note"] = {
    "wrong_sign_junction_original_sign_criterion_met": bool((mid(Lw["eta_b"]) > 0) != (reg["eta_b"] > 0)),
    "wrong_sign_junction_eta_b": float(mid(Lw["eta_b"])),
    "wrong_sign_junction_reldev_vs_checkpoint": wrong_sign_reldev,
    "note": ("First run: sign-flip criterion was a design error of the control (the sign does not flip); "
             "replaced by a magnitude criterion. Reported, not hidden."),
}
# C_d: containment self-test with a deliberately shifted point
Lreg = lo_model_iv(repr(reg["delta"]), cI)
shifted = mid(Lreg["eta_b"]) + 1000 * (hi(Lreg["eta_b"]) - lo(Lreg["eta_b"])) + mpmath.mpf("1e-50")
ctrl["containment_detects_shifted_point"] = not contains(Lreg["eta_b"], shifted)
out["controls"] = {kk: bool(v) for kk, v in ctrl.items()}
controls_ok = all(ctrl.values())

# ---------------------------------------------------------------- 4. conditioning / planning numbers
eta_h_reg = reg["eta_h"]
iv.dps = 15
phi_h_15 = iv.mpf(1) + iv.mpf(repr(eta_h_reg))
iv.dps = IV_DPS
u_b_reg = mpmath.mpf(reg["y_b"]) / 9
plan = {
    "binary64_1_plus_eta_h_equals_1": (1.0 + eta_h_reg) == 1.0,
    "interval_15_digits_phi_h_contains_1": bool(lo(phi_h_15) <= 1 <= hi(phi_h_15)),
    "digits_for_16_significant_digits_of_eta_h_in_direct_phi_form": int(math.ceil(-math.log10(abs(eta_h_reg)))) + 16,
    "amplification_G_at_registered_u_b": float(reg_row["G_at_checkpoint_u_b"]),
    "log10_amplification": math.log10(reg_row["G_at_checkpoint_u_b"]),
    "unwanted_mode_suppression_log10_exp32_times_ub_minus_u0": float(32 * (u_b_reg - U0_GERM) / mpmath.log(10)),
    "u0_germ_interface": float(U0_GERM),
    "registered_u_b": float(u_b_reg),
    "a_hat_eta_h_times_G_over_delta_registered": float(mpmath.mpf(eta_h_reg) * reg_row["G_at_checkpoint_u_b"] / mpmath.mpf(reg["delta"])),
    "minus_9c_over_64": float(-9 * cP / 64),
}
# known O(delta) relative H^2 defect of the thin-brane LO model: +27 c^2 / (384 (1+c)) (from -c^2 delta^2/384)
out["H2_defect_coefficient_expected"] = float(27 * cP**2 / (384 * (1 + cP)))
out["H2_defect_coefficient_observed_smallest_delta"] = coeffs["reldev_H2_LO_vs_checkpoint_H2"][0]
out["planning"] = plan
out["status_labels"]["planning"] = ("CONDITIONAL: linear-regime conditioning estimates; they assume the nonlinear "
                                    "correction stays small (|eta| < 1e-4 at delta=1e-3) and say nothing rigorous about "
                                    "the nonlinear problem")

# ---------------------------------------------------------------- expectations and summary
exp = {}
exp["E1_exact_and_controls"] = exact_ok and controls_ok
exp["E2_containment_all"] = all(containment)
exp["E2_n_containment_checks"] = len(containment)
exp["E3_metric_only_H2_reproduced_1e-12"] = all(abs(z["reldev_H2_LO_vs_checkpoint_metric_only_H2"]) < 1e-12 for z in table)
exp["E4_slopes_in_0.8_1.2"] = {kk: (0.8 <= v <= 1.2) for kk, v in slopes.items()}
out["expectations"] = exp
out["controls_ok"] = controls_ok
out["exact_ok"] = exact_ok
out["all_expectations_met"] = (exp["E1_exact_and_controls"] and exp["E2_containment_all"]
                               and exp["E3_metric_only_H2_reproduced_1e-12"] and all(exp["E4_slopes_in_0.8_1.2"].values()))
out["precision"] = {"iv_dps": IV_DPS, "mp_dps": MP_DPS, "sympy": sp.__version__, "mpmath": mpmath.__version__}
out["scope"] = ("Pilot for a computer-assisted proof. Encloses only the closed-form leading-order truncated model and "
                "the linear germ amplification. It does not enclose the nonlinear +1 branch, does not prove existence or "
                "uniqueness, and adds no stability or dynamical statement.")
out["runtime_s"] = round(time.time() - t_start, 2)
with open(os.path.join(HERE, "certification_pilot.json"), "w") as fh:
    json.dump(out, fh, indent=1)
print(json.dumps({"exact_ok": exact_ok, "controls": out["controls"], "expectations": exp, "slopes": slopes,
                  "planning": plan, "runtime_s": out["runtime_s"]}, indent=1))
