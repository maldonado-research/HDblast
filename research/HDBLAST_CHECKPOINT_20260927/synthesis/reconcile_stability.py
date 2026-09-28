#!/usr/bin/env python3
"""Cross-method reconciliation for the 27 Sept 2026 HDBLAST checkpoint.

Reads (read-only) the saved JSON outputs of three workstreams
  - stability_gauge_invariant  (Method A: gauge-invariant shooting, all dS4 harmonics)
  - stability_time_domain      (Method B: discrete operator + time domain, homogeneous sector)
  - analytic_structure         (exact series + leading-order spectrum)
and checks, with stated tolerances, whether they agree where they overlap.

Also includes deliberate controls (a wrong mu^2 <-> growth-rate map, a wrong
analytic limit) that must FAIL.  Writes synthesis/RECONCILIATION.json.
No computation here re-solves any ODE/PDE; it only compares saved results and
evaluates exact expressions with sympy.
"""
import hashlib
import json
import os
import sys

import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def load(rel):
    with open(os.path.join(ROOT, rel)) as f:
        return json.load(f)


def sha(rel):
    with open(os.path.join(ROOT, rel), "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def mu2_from_growth(p):
    """dS4 harmonic: e^{p H tau}, p^2 + 3p + mu^2 = 0  =>  mu^2 = -p(p+3)."""
    return -p * (p + 3.0)


def mu2_wrong(p):
    """Deliberately wrong map (flat-space intuition mu^2 = -p^2). Must fail."""
    return -p * p


inputs = {
    "A_spectrum": "stability_gauge_invariant/SPECTRUM_RESULTS.json",
    "A_controls": "stability_gauge_invariant/CONTROLS_RESULTS.json",
    "A_backgrounds": "stability_gauge_invariant/BACKGROUNDS.json",
    "B_key": "stability_time_domain/results/KEY_RESULTS.json",
    "S_series": "analytic_structure/SERIES_COEFFICIENTS.json",
    "S_spectrum": "analytic_structure/SPECTRUM_LEADING_ORDER.json",
    "P_summary": "preheating/SUMMARY.json",
}
A = load(inputs["A_spectrum"])
AC = load(inputs["A_controls"])
B = load(inputs["B_key"])["key_results"]
S = load(inputs["S_series"])
P = load(inputs["P_summary"])

checks = []


def record(name, value, target, tol, expect_pass=True, note=""):
    diff = abs(value - target)
    rel = diff / abs(target) if target != 0 else diff
    ok = rel <= tol
    checks.append({
        "name": name, "value": value, "target": target, "rel_diff": rel,
        "tolerance_rel": tol, "agrees": ok, "expected_to_agree": expect_pass,
        "passes": ok == expect_pass, "note": note,
    })
    return ok


# ---------------------------------------------------------------- 1. calibration
pA = A["calibration"]["growth"]
muA = A["calibration"]["mu2"]
cal_B = {r["run"]: r for r in B["calibration_time_domain"]}
pB_pencil = cal_B["original_lin_h5e-5"]["rate_pencil"]
pB_nl = cal_B["original_nl_eps1e-8_h5e-5"]["rate_pencil"]
pB_eig = B["calibration_semi_discrete_eigenvalue"]["5e-5"]
record("calibration growth: A (shooting) vs B (time-domain pencil, h=5e-5)", pB_pencil, pA, 1e-6)
record("calibration growth: A vs B (nonlinear RHS, eps=1e-8, h=5e-5)", pB_nl, pA, 1e-6)
record("calibration growth: A vs B (semi-discrete eigenvalue, h=5e-5)", pB_eig, pA, 1e-7)
record("calibration: mu^2 from A's growth via p(p+3) reproduces A's mu^2", mu2_from_growth(pA), muA, 1e-12)
record("CONTROL (must fail): wrong map mu^2=-p^2 applied to A's growth", mu2_wrong(pA), muA, 1e-3,
       expect_pass=False, note="flat-space map ignores the 3H friction of dS4 harmonics")

# ------------------------------------------------ 2. wrong-sign sigma'' control
pA_flip = AC["C2_flipped_shell_curvature_plus_branch"]["0.001"]["growth"][0]
muA_flip = AC["C2_flipped_shell_curvature_plus_branch"]["0.001"]["roots"][0]
pB_flip = B["control_sigma2_flipped_growth_rate"]["shift_h5e-5"]
record("flipped-sigma'' control growth (delta=0.001): A vs B shift-invert h=5e-5", pB_flip, pA_flip, 1e-7,
       note="two independent codes, one detects the same artificial tachyon")
record("flipped-sigma'' control: mu^2 from B's rate vs A's root", mu2_from_growth(pB_flip), muA_flip, 1e-7)
record("CONTROL (must fail): wrong map on B's flipped rate", mu2_wrong(pB_flip), muA_flip, 1e-3, expect_pass=False)

# B's residual gauge mode lambda->1 corresponds to the l=1 harmonic mu^2=-4, which Method A shows is pure gauge
lam_g = B["gauge_lambda1_L6_vs_L8_h1e-4"]["L8"]
record("B's gauge eigenvalue (L=8) mapped to mu^2 vs A's pure-gauge l=1 value -4", mu2_from_growth(lam_g), -4.0, 1e-6,
       note="consistency of interpretation: dS4 translation harmonic p=1 <-> mu^2=-4")

# ----------------------------------------- 3. leading physical rate on +1 branch
rich = B["plus_td_rate_richardson"]["extrapolated"]
rich_b = B["plus_td_rate_richardson_bumpB"]["extrapolated"]
record("+1 branch leading amplitude rate (B, Richardson) vs continuum edge -3/2", rich, -1.5, 1e-3,
       note="verifier: only -1.500 +/- 1e-3 is supported")
record("+1 branch leading rate, second profile (B) vs -3/2", rich_b, -1.5, 1e-3)
h_vac = P["constants"]["h_vac"]  # H_+/H_0 from the +1 branch (preheating constants block)
edge_H0 = -1.5 * h_vac
e_fold_H0tau = 1.0 / (1.5 * h_vac)

# A: no scalar roots below 9/4 at any delta; B: no shell-supported mode above the line
A_no_roots = all(len(r["scalar_roots"]) == 0 and len(r["longitudinal_roots"]) == 0 for r in A["plus_summary"])
B_no_shell_modes = (len(B["plus_dense_shell_supported_modes_above_line_excluding_grid_scale"]) == 0)
checks.append({"name": "A: zero scalar/longitudinal roots below 9/4 at all six delta", "agrees": A_no_roots,
               "expected_to_agree": True, "passes": A_no_roots})
checks.append({"name": "B: no shell-supported mode above Re(lambda)=-3/2 (dense spectra, excl. grid-scale)",
               "agrees": B_no_shell_modes, "expected_to_agree": True, "passes": B_no_shell_modes})

# ----------------------------------- 4. analytic delta->0 limit of B = g_b + sigma''/2
phi = sp.symbols("phi")
W = 1 - phi + phi**3 / 3
U = sp.Rational(1, 2) * sp.diff(W, phi)**2 - sp.Rational(2, 3) * W**2
U2 = sp.simplify(sp.diff(U, phi, 2).subs(phi, 1))          # bulk mass^2 at phi=+1
ell = sp.simplify(3 / W.subs(phi, 1))                      # AdS radius: A' = -W/3 => 1/ell = W/3
Delta = sp.simplify(2 + sp.sqrt(4 + U2 * ell**2))          # conformal dimension in AdS5
n_gegen = sp.simplify(Delta - 4)                           # Gegenbauer degree
g_limit = sp.simplify(n_gegen / ell)                       # regular-solution growth rate d ln X/dy
sigma2_half = sp.simplify(sp.diff(2 * W, phi, 2).subs(phi, 1) / 2)
B_limit = sp.simplify(g_limit + sigma2_half)
rows = sorted(A["plus_summary"], key=lambda r: r["delta"])[:3]
# linear fit B(delta) = B0 + B1 delta through the three smallest delta
import numpy as np  # noqa: E402
d_arr = np.array([r["delta"] for r in rows])
b_arr = np.array([r["B"] for r in rows])
B1, B0 = np.polyfit(d_arr, b_arr, 1)
fit_resid = float(np.max(np.abs(b_arr - (B0 + B1 * d_arr))))
record("B(delta->0) from A's numerics (linear fit, 3 smallest delta) vs exact 32/9", float(B0), float(B_limit), 1e-4)
record("CONTROL (must fail): same fit vs wrong limit 14/9 (bulk term only)", float(B0), float(g_limit), 1e-2,
       expect_pass=False)

# --------------------------------------------------- 5. tensor sector agreement
max_tensor = max(abs(t) for r in A["plus_summary"] for t in r["tensor_roots"])
n_tensor = [len(r["tensor_roots"]) for r in A["plus_summary"]]
tensor_ok = max_tensor < 1e-12 and all(n == 1 for n in n_tensor)
checks.append({"name": "A tensor roots: exactly one per delta, |mu^2|<1e-12 (graviton); analytic_structure: spectrum {0} U [9/4, inf)",
               "value": max_tensor, "agrees": tensor_ok, "expected_to_agree": True, "passes": tensor_ok})

# ----------------------------------------- 6. small-delta 4D limit of calibration
c4 = AC["C4_small_detuning_limit"]
record("calibration delta->0 extrapolation (A) vs closed-form 4D value", c4["intercept"], c4["closed_form_4D"], 1e-6)

# ------------------------------------------------ 7. series vs numerical branch
c = sp.Rational(0)  # placeholder to keep sympy namespace clear
csym = sp.symbols("c")
h1 = sp.sympify(S["H2"][1], locals={"c": csym})
h2 = sp.sympify(S["H2"][2], locals={"c": csym})
h2_pkg = (1 + csym)**2 / 36 - csym**2 / 384
record("series h2 (analytic_structure) minus 22 Sept package form, symbolic", float(sp.simplify(h2 - h2_pkg).subs(csym, 1)) + 1.0,
       1.0, 1e-15, note="value-1 is the symbolic difference evaluated at c=1 (identically 0)")
creg = 0.5975949350280132
H2_series = sum(float(sp.sympify(S["H2"][k], locals={"c": csym}).subs(csym, creg)) * 0.001**k for k in range(1, len(S["H2"])))
bg = [r for r in load(inputs["A_backgrounds"])["plus"] if abs(r["delta"] - 0.001) < 1e-12][0]
record("H^2 at delta=0.001: exact series through delta^8 vs package float solution", H2_series, bg["package"]["H2"], 1e-13)
H2_two_terms = float((h1 * 0.001 + h2 * 0.001**2).subs(csym, creg))
record("CONTROL (must fail at 1e-13): two-term series vs package H^2", H2_two_terms, bg["package"]["H2"], 1e-13,
       expect_pass=False, note="truncation error O(delta^3) ~ 4e-8 relative")

out = {
    "status": "numerical comparison of saved outputs + exact sympy limits; no new ODE/PDE solves",
    "inputs": {k: {"path": v, "sha256": sha(v)} for k, v in inputs.items()},
    "checks": checks,
    "all_pass": all(ch["passes"] for ch in checks),
    "n_checks": len(checks),
    "exact": {"U''(1)": str(U2), "ell_plus": str(ell), "Delta": str(Delta), "gegenbauer_degree": str(n_gegen),
              "g_limit": str(g_limit), "sigma''/2 at phi=1": str(sigma2_half), "B_limit": str(B_limit)},
    "fit_B_vs_delta": {"B0": float(B0), "B1": float(B1), "max_fit_residual": fit_resid,
                       "deltas": d_arr.tolist()},
    "derived_for_next_tests": {
        "status": "conditional (linear theory; assumes the far-region evolution approaches the +1 branch)",
        "H_plus_over_H0": h_vac,
        "slowest_linear_decay_rate_in_H0_units": edge_H0,
        "e_folding_time_in_H0tau": e_fold_H0tau,
        "note": "If Chat 14's positive-side run relaxes onto the +1 branch linearly, deviations should decay no slower "
                "than exp(-1.5 H_+ tau) (continuum edge, possibly with power-law prefactors from the branch cut).",
    },
}
with open(os.path.join(HERE, "RECONCILIATION.json"), "w") as f:
    json.dump(out, f, indent=1)
for ch in checks:
    print(("PASS " if ch["passes"] else "FAIL ") + ch["name"])
print("all_pass:", out["all_pass"], "B0 fit:", B0, "B_limit:", B_limit, "edge rate (H0):", edge_H0)
sys.exit(0 if out["all_pass"] else 1)
