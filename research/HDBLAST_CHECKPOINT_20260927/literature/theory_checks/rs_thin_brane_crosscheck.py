#!/usr/bin/env python3
"""Literature anchoring check for the HDBLAST static "+1 branch".

Question: does the registered model's static +1 branch reduce, at leading orders in the
detuning delta, to the standard thin-brane de Sitter relation for a single Z2 brane in AdS5,
    H^2 = (sigma/6)^2 - 1/ell^2      (kappa_5 = 1, junction K = sigma/6),
found in the braneworld literature (Kaloper 1999; Nihei 1999; Kim & Kim 2000; Garriga & Sasaki 2000)?
And is the registered "balanced" tension 2W exactly the Randall-Sundrum/Karch-Randall critical
tension 6/ell at each AdS extremum of W?

Inputs (read-only): the checkpoint's PLUS_BRANCH_RESULTS.json (22 Sept 2026 package).
Output: rs_thin_brane_crosscheck.json in this folder.
Precision: exact rational/symbolic algebra (SymPy) for identities; mpmath at 30 digits for numbers.
Controls: three deliberately wrong variants must FAIL.
"""
import json, hashlib, os, sys
import sympy as sp
import mpmath as mp

mp.mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
PKG = ("/home/user/unified-theory-maldonado/new-files/latest-work/"
       "HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/"
       "HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/static_branch/PLUS_BRANCH_RESULTS.json")

phi, d, c = sp.symbols('phi delta c', real=True)
W = 1 - phi + phi**3/3
Wp = sp.diff(W, phi)
U = sp.Rational(1, 2)*Wp**2 - sp.Rational(2, 3)*W**2
sigma = 2*W + d*(1 + c*phi)

checks = {}

# 1. AdS radii at the extrema: Einstein eq R_AB = (2/3) U g_AB and AdS5 R_AB = -(4/ell^2) g_AB  =>  U = -6/ell^2
ell = {}
for p0 in (1, -1):
    Uval = sp.nsimplify(U.subs(phi, p0))
    ell2 = sp.nsimplify(-6/Uval)
    ell[p0] = sp.sqrt(ell2)
    crit = sp.nsimplify(6/ell[p0])           # RS / Karch-Randall critical (flat-brane) tension, junction K = sigma/6
    bal = sp.nsimplify(2*W.subs(phi, p0))    # registered balanced tension
    checks[f"phi={p0}"] = dict(U=str(Uval), ell=str(ell[p0]), critical_tension_6_over_ell=str(crit),
                               balanced_tension_2W=str(bal), equal=bool(sp.simplify(crit - bal) == 0))

# 2. Thin-brane dS relation evaluated with the tension at phi=1 (constant-phi, metric-only junction)
sig1 = sp.expand(sigma.subs(phi, 1))
H2_thin = sp.expand((sig1/6)**2 - 1/ell[1]**2)
# checkpoint expansion: H^2 = delta(1+c)/27 + delta^2[(1+c)^2/36 - c^2/384] + O(delta^3)
H2_ckpt_2nd = d*(1 + c)/27 + d**2*((1 + c)**2/36 - c**2/384)
diff = sp.simplify(sp.expand(H2_ckpt_2nd - H2_thin))
checks["thin_brane_H2_exact_in_delta"] = str(H2_thin)
checks["checkpoint_minus_thin_brane_through_delta2"] = str(diff)
checks["leading_order_agrees"] = bool(sp.simplify(sp.diff(H2_thin, d).subs(d, 0) - sp.diff(H2_ckpt_2nd, d).subs(d, 0)) == 0)
checks["difference_is_pure_scalar_profile_term"] = bool(sp.simplify(diff + d**2*c**2/384) == 0)

# 3. Numbers at the registered point, compared against the checkpoint JSON (read-only)
cR = mp.mpf('0.5975949350280132'); dR = mp.mpf('0.001')
rho0 = mp.mpf('78.82817714224423')    # original registered shell radius, quoted in the checkpoint review
H2_thin_num = (((mp.mpf(2)/3 + dR*(1 + cR))/6)**2 - mp.mpf(1)/81)
with open(PKG, 'rb') as fh:
    raw = fh.read()
pkg = json.loads(raw)
row = [r for r in pkg['rows'] if abs(r['delta'] - 0.001) < 1e-15][0]
num = dict(
    H2_thin_brane=mp.nstr(H2_thin_num, 20),
    checkpoint_metric_only_H2=row['metric_only_H2'],
    rel_diff_thin_vs_checkpoint_metric_only=mp.nstr(abs(H2_thin_num - row['metric_only_H2'])/H2_thin_num, 5),
    checkpoint_full_H2=row['H2'],
    H_over_H0_thin_brane=mp.nstr(mp.sqrt(H2_thin_num)*rho0, 12),
    H_over_H0_full_branch=mp.nstr(mp.sqrt(mp.mpf(row['H2']))*rho0, 12),
    ppm_shift_in_H_full_vs_thin=mp.nstr((mp.sqrt(mp.mpf(row['H2'])/H2_thin_num) - 1)*1e6, 6),
    predicted_ppm_from_c2_over_384=mp.nstr(-(dR**2*cR**2/384)/(2*H2_thin_num)*1e6, 6),
)
checks["numbers_at_registered_point"] = num
checks["input_sha256"] = hashlib.sha256(raw).hexdigest()

# 4. dS-brane KK continuum threshold: m^2 = (3H/2)^2 = 9H^2/4 (Garriga-Sasaki; 2608.01762 snippet)
checks["KK_threshold_over_H2"] = str(sp.Rational(3, 2)**2)


# 5. Geometry identities behind the literature comparison (exact, SymPy Ricci tensor from the metric)
y, t, x1, x2, x3, L = sp.symbols('y t x1 x2 x3 ell', positive=True)
def ricci(g, X):
    n = len(X); ginv = g.inv()
    Gam = [[[sp.simplify(sum(ginv[a, e]*(sp.diff(g[e, b], X[cc]) + sp.diff(g[e, cc], X[b]) - sp.diff(g[b, cc], X[e]))
             for e in range(n))/2) for cc in range(n)] for b in range(n)] for a in range(n)]
    R = sp.zeros(n)
    for b in range(n):
        for cc in range(n):
            R[b, cc] = sp.simplify(sum(sp.diff(Gam[a][b][cc], X[a]) for a in range(n))
                                   - sum(sp.diff(Gam[a][b][a], X[cc]) for a in range(n))
                                   + sum(Gam[a][a][e]*Gam[e][b][cc] for a in range(n) for e in range(n))
                                   - sum(Gam[a][cc][e]*Gam[e][b][a] for a in range(n) for e in range(n)))
    return R
X = [y, t, x1, x2, x3]
# (a) dS-sliced AdS5 with a regular cone: rho = ell*sinh(y/ell), dS4 in flat slicing, signature (-++++) on the slice
rho = L*sp.sinh(y/L)
g = sp.diag(1, -rho**2, rho**2*sp.exp(2*t), rho**2*sp.exp(2*t), rho**2*sp.exp(2*t))
R = ricci(g, X)
checks["dS_sliced_AdS5_is_Einstein_R_eq_minus4_over_ell2_g"] = bool(sp.simplify(R + 4/L**2*g) == sp.zeros(5))
checks["regular_cone_rho_prime_at_0"] = str(sp.diff(rho, y).subs(y, 0))
# (b) flat wall with the registered superpotential: R_AB = phi_A phi_B + (2/3) U g_AB with phi' = W_phi, A' = -W/3
A = sp.Function('A')(y); ph = sp.Function('ph')(y)
g2 = sp.diag(1, -sp.exp(2*A), sp.exp(2*A), sp.exp(2*A), sp.exp(2*A))
R2 = ricci(g2, X)
Wy = W.subs(phi, ph); Wpy = Wp.subs(phi, ph); Uy = U.subs(phi, ph)
subsd = {sp.Derivative(A, (y, 2)): sp.diff(-Wy/3, y), sp.Derivative(A, y): -Wy/3}
eqs_ok = True
for i in range(5):
    lhs = R2[i, i] - (sp.diff(ph, y)**2 if i == 0 else 0) - sp.Rational(2, 3)*Uy*g2[i, i]
    lhs = lhs.subs(subsd).subs(sp.Derivative(ph, y), Wpy).doit().subs(sp.Derivative(ph, y), Wpy)
    eqs_ok = eqs_ok and bool(sp.simplify(sp.expand(lhs)) == 0)
checks["first_order_flow_solves_Einstein_eqs_(DeWolfe_et_al_structure)"] = eqs_ok
kink = -sp.tanh(y)
checks["flat_kink_phi_eq_minus_tanh_solves_phi_prime_eq_W_phi"] = bool(sp.simplify(sp.diff(kink, y) - Wp.subs(phi, kink)) == 0)
checks["ell_equals_3_over_W_at_extrema"] = bool(all(sp.simplify(ell[p0] - 3/W.subs(phi, p0)) == 0 for p0 in (1, -1)))
lead = sp.limit(U/phi**6, phi, sp.oo)
checks["U_large_phi_leading_coefficient_of_phi6"] = str(lead)
checks["U_unbounded_below"] = bool(lead < 0)

# Controls (must fail)
ctrl = {}
# (a) one-sided junction K = sigma/3 instead of sigma/6
H2_bad = sp.expand((sig1/3)**2 - 1/ell[1]**2)
ctrl["one_sided_junction_sigma_over_3_matches"] = bool(sp.simplify(sp.diff(H2_bad, d).subs(d, 0) - sp.Rational(1, 27)*(1 + c)) == 0)
# (b) wrong AdS radius (phi=-1 vacuum, ell=9/5) with the phi=+1 tension
H2_bad2 = sp.expand((sig1/6)**2 - 1/ell[-1]**2)
ctrl["wrong_vacuum_radius_gives_positive_H2_at_delta0"] = bool(sp.nsimplify(H2_bad2.subs(d, 0)) > 0)
ctrl["wrong_vacuum_radius_matches"] = bool(sp.simplify(H2_bad2 - H2_thin) == 0)
# (c) perturbed c: numbers must visibly disagree with the checkpoint
H2_pert = (((mp.mpf(2)/3 + dR*(1 + cR*mp.mpf('1.01')))/6)**2 - mp.mpf(1)/81)
ctrl["perturbed_c_rel_diff_vs_checkpoint_metric_only"] = mp.nstr(abs(H2_pert - row['metric_only_H2'])/H2_pert, 5)

subs_bad = {sp.Derivative(A, (y, 2)): sp.diff(+Wy/3, y), sp.Derivative(A, y): +Wy/3}
bad_ok = True
for i in range(5):
    lhs = R2[i, i] - (sp.diff(ph, y)**2 if i == 0 else 0) - sp.Rational(2, 3)*Uy*g2[i, i]
    lhs = lhs.subs(subs_bad).subs(sp.Derivative(ph, y), Wpy).doit().subs(sp.Derivative(ph, y), Wpy)
    bad_ok = bad_ok and bool(sp.simplify(sp.expand(lhs)) == 0)
ctrl["wrong_sign_flow_A_prime_eq_plus_W_over_3_solves"] = bad_ok
checks["controls"] = ctrl

ok = (checks["dS_sliced_AdS5_is_Einstein_R_eq_minus4_over_ell2_g"] and checks["first_order_flow_solves_Einstein_eqs_(DeWolfe_et_al_structure)"] and checks["flat_kink_phi_eq_minus_tanh_solves_phi_prime_eq_W_phi"] and checks["ell_equals_3_over_W_at_extrema"] and checks["U_unbounded_below"]
      and checks["phi=1"]["equal"] and checks["phi=-1"]["equal"] and checks["leading_order_agrees"]
      and checks["difference_is_pure_scalar_profile_term"]
      and float(num["rel_diff_thin_vs_checkpoint_metric_only"]) < 1e-12
      and not ctrl["one_sided_junction_sigma_over_3_matches"]
      and not ctrl["wrong_vacuum_radius_matches"]
      and not ctrl["wrong_sign_flow_A_prime_eq_plus_W_over_3_solves"]
      and float(ctrl["perturbed_c_rel_diff_vs_checkpoint_metric_only"]) > 1e-4)
checks["all_pass"] = bool(ok)
checks["label"] = ("EXACT (identities) + NUMERICAL (comparison with checkpoint floats); "
                   "a consistency check against standard literature formulas, not a new result")
out = os.path.join(HERE, 'rs_thin_brane_crosscheck.json')
with open(out, 'w') as fh:
    json.dump(checks, fh, indent=1)
print(json.dumps(checks, indent=1))
sys.exit(0 if ok else 1)
