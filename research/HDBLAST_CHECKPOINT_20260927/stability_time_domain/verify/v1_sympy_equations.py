"""Audit V1 (symbolic): independent derivation of the 1+1 Einstein-scalar system used by the
frozen solver, its linearisation, the junction data and the potential expansions.
Writes verify/v1_sympy_equations.json.  kappa_5^2 = 1, canonical scalar, metric
e^{2B}(-dt^2+dz^2)+e^{2A}dx_3^2.  Everything is derived here from the metric (Christoffels ->
Ricci), not copied from the audited code; the audited code's expressions are then compared."""
import json, sympy as sp
from pathlib import Path

t, z = sp.symbols('t z', real=True)
x1, x2, x3 = sp.symbols('x1 x2 x3', real=True)
X = [t, z, x1, x2, x3]
A = sp.Function('A')(t, z); B = sp.Function('B')(t, z); ph = sp.Function('phi')(t, z)
Uf = sp.Function('U')
g = sp.diag(-sp.exp(2*B), sp.exp(2*B), sp.exp(2*A), sp.exp(2*A), sp.exp(2*A))
gi = g.inv()
n = 5
Gam = [[[sp.simplify(sum(gi[a, d]*(sp.diff(g[d, b], X[c]) + sp.diff(g[d, c], X[b]) - sp.diff(g[b, c], X[d]))
                          for d in range(n))/2) for c in range(n)] for b in range(n)] for a in range(n)]
def ric(b, c):
    return sp.simplify(sum(sp.diff(Gam[a][b][c], X[a]) for a in range(n)) - sum(sp.diff(Gam[a][b][a], X[c]) for a in range(n))
                       + sum(Gam[a][a][d]*Gam[d][b][c] for a in range(n) for d in range(n))
                       - sum(Gam[a][c][d]*Gam[d][b][a] for a in range(n) for d in range(n)))
R = sp.Matrix(n, n, lambda i, j: ric(i, j) if (i == j or (i, j) in [(0, 1), (1, 0)]) else 0)
Rs = sp.simplify(sum(gi[i, i]*R[i, i] for i in range(n)))
dph = [sp.diff(ph, x) for x in X]
kin = sum(gi[i, i]*dph[i]**2 for i in range(n))
T = sp.Matrix(n, n, lambda i, j: dph[i]*dph[j] - g[i, j]*(kin/2 + Uf(ph)))
G = R - Rs*g/2
E = sp.simplify(G - T)   # Einstein equations E=0 (kappa^2=1)

At, Az, Bt, Bz = [sp.diff(A, t), sp.diff(A, z), sp.diff(B, t), sp.diff(B, z)]
pt, pz = sp.diff(ph, t), sp.diff(ph, z)
Up = sp.Subs(sp.Derivative(Uf(sp.Symbol('s')), sp.Symbol('s')), sp.Symbol('s'), ph).doit()
# Evolution equations as implemented (full, nonlinear form reconstructed from rhs(); see VERIFICATION.md)
Att = sp.diff(A, z, 2) - 3*(At**2 - Az**2) + sp.Rational(2, 3)*sp.exp(2*B)*Uf(ph)
Btt = sp.diff(B, z, 2) + 3*(At**2 - Az**2) - pt**2/2 + pz**2/2 - sp.Rational(1, 3)*sp.exp(2*B)*Uf(ph)
checks = {}
# Scalar equation: box phi - U'(phi) = 0
box = sum(sp.diff(sp.sqrt(-g.det())*gi[i, i]*dph[i], X[i]) for i in range(n))/sp.sqrt(-g.det())
s_ = sp.Symbol('s'); dU = sp.diff(Uf(s_), s_).subs(s_, ph)
ptt_impl = sp.diff(ph, z, 2) - 3*At*pt + 3*Az*pz - sp.exp(2*B)*dU
scalar_eq = sp.simplify(box - dU)
# solve scalar_eq for phi_tt and compare
PTT = sp.Symbol('PTT'); ATT = sp.Symbol('ATT'); BTT = sp.Symbol('BTT')
rep = {sp.diff(ph, t, 2): PTT, sp.diff(A, t, 2): ATT, sp.diff(B, t, 2): BTT}
sol_ptt = sp.solve(sp.simplify(scalar_eq.subs(rep)), PTT)[0]
checks['scalar_evolution_matches_rhs'] = sp.simplify(sol_ptt - ptt_impl) == 0
# wrong-formula control: flip the sign of the 3 A_z phi_z friction term -> must be detected
checks['CONTROL_wrong_scalar_friction_detected'] = sp.simplify(sol_ptt - (ptt_impl - 6*Az*pz)) != 0
# Einstein: xx and zz/tt combinations -> solve for ATT, BTT
Exx = sp.simplify(E[2, 2].subs(rep)); Ezz = sp.simplify(E[1, 1].subs(rep)); Ett = sp.simplify(E[0, 0].subs(rep))
# The dynamical equations: xx component and (tt - zz) combination contain second time derivatives.
sol = sp.solve([Exx, sp.simplify(Ett - Ezz)], [ATT, BTT], dict=True)
if not sol:
    sol = sp.solve([Exx, Ezz], [ATT, BTT], dict=True)
sATT = sol[0][ATT]; sBTT = sol[0][BTT]
# The rhs equations may differ from G=T dynamical eqs by multiples of the constraints; test both.
Ham = sp.simplify(Ett)     # G_tt - T_tt: must contain no second time derivatives
Mom = sp.simplify(E[0, 1])
checks['Ham_has_no_second_time_derivs'] = not any(Ham.has(d) for d in [ATT, BTT, PTT])
dA = sp.simplify(sATT - Att); dB = sp.simplify(sBTT - Btt)
# Try to express the differences as combinations of the Hamiltonian constraint (Ett+Ezz) and momentum
c1, c2, c3, c4 = sp.symbols('c1 c2 c3 c4')
def as_constraint_combo(diff):
    for cand in [Ett, Ezz, Ett + Ezz, Ett - Ezz]:
        cand = sp.simplify(cand.subs(rep))
        if cand.has(ATT) or cand.has(BTT): continue
        r = sp.simplify(diff/cand)
        if r.free_symbols <= set() or sp.simplify(sp.diff(r, t)) == 0 and sp.simplify(sp.diff(r, z)) == 0:
            return str(r)
    return None
checks['A_evolution_matches_Einstein'] = (dA == 0)
checks['B_evolution_matches_Einstein'] = (dB == 0)
checks['CONTROL_wrong_B_coefficient_detected'] = sp.simplify(sBTT - (Btt + sp.exp(2*B)*Uf(ph)/6)) != 0
extra = {}
if dA != 0: extra['A_diff'] = str(dA); extra['A_diff_as_constraint_multiple'] = as_constraint_combo(dA)
if dB != 0: extra['B_diff'] = str(dB); extra['B_diff_as_constraint_multiple'] = as_constraint_combo(dB)

# Constraints as implemented in diagnostics() (full, nonlinear, written in A,B,phi):
pa = At - 1; pb = Bt
Mimpl = -3*sp.diff(A, t, z) - 3*At*(Az - Bz) + 3*Az*Bt - pt*pz
# Hamiltonian as implemented, rebuilt in full variables: H = -2*(e^{2B}U - rho^2 U_bg) + 12 pa + 6 pa^2 + 6(1+pa) pb
#   - 18 hc az + 6 hc bz - 12 az^2 + 6 az bz - 6 azz - pf^2 - 2 phz fz - fz^2 ;  the background pieces cancel
# only on-shell for the background, so compare against G-T up to a constant factor instead:
Hfull = -2*sp.exp(2*B)*Uf(ph) + 6*At**2 + 6*At*Bt - 12*Az**2 + 6*Az*Bz - 6*sp.diff(A, z, 2) - pt**2 - pz**2
Hfull = sp.expand(Hfull)
ratioM = sp.simplify(Mimpl/sp.simplify(E[0, 1]))
checks['momentum_constraint_is_multiple_of_G_tz'] = bool(ratioM.is_number)
extra['momentum_ratio'] = str(ratioM)
# find constant k so that Hfull = k*(Ett+Ezz) (both free of second time derivatives)
HE = sp.simplify(Ett)
ratioH = sp.simplify(Hfull/HE)
checks['hamiltonian_constraint_is_multiple_of_Gtt_minus_Ttt'] = bool(ratioH.is_number)
extra['hamiltonian_ratio'] = str(ratioH)

# ---- consistency of the audited Hfull with the implemented perturbative H (background pieces):
rho = sp.Function('rho')(z); a_, b_, f_ = [sp.Function(nm)(t, z) for nm in ('a', 'b', 'f')]
phb = sp.Function('phib')(z)
hc = sp.diff(rho, z)/rho; phz = sp.diff(phb, z)
subsd = {A: t + sp.log(rho) + a_, B: sp.log(rho) + b_, ph: phb + f_}

# ---- static background: y-equations (rho_yy=-rho(p^2/4+U/6), first integral) vs Einstein
y = sp.Symbol('y'); r = sp.Function('r')(y); P = sp.Function('P')(y)
ry = sp.diff(r, y); ryy = sp.diff(r, y, 2)
Hs = ry**2 - 1 - r**2*(sp.diff(P, y)**2/12 - Uf(P)/6)
dHs = sp.diff(Hs, y).subs(sp.diff(r, y, 2), -r*(sp.diff(P, y)**2/4 + Uf(P)/6)).subs(
    sp.diff(P, y, 2), sp.diff(Uf(s_), s_).subs(s_, P) - 4*ry*sp.diff(P, y)/r)
checks['static_first_integral_conserved_by_second_order_system'] = sp.simplify(dHs) == 0

# ---- superpotential / BPS / junction normalisation
p = sp.Symbol('p'); d, c = sp.symbols('delta c')
W = 1 - p + p**3/3; Wp = sp.diff(W, p); U = Wp**2/2 - sp.Rational(2, 3)*W**2
sig = 2*W + d*(1 + c*p)
checks['sigma1'] = sp.simplify(sp.diff(sig, p) - (2*(p**2 - 1) + d*c)) == 0
checks['sigma2_eq_4phi'] = sp.simplify(sp.diff(sig, p, 2) - 4*p) == 0
checks['sigma3_eq_4'] = sp.diff(sig, p, 3) == 4
# BPS flat brane: A_y = W/3, phi_y = -W_p solves flat-sliced Hamiltonian A_y^2 = phi_y^2/12 - U/6 and
# the junctions A_y = sigma/6, phi_y = -sigma'/2 at delta=0.
checks['BPS_hamiltonian'] = sp.simplify((W/3)**2 - (Wp**2/12 - U/6)) == 0
checks['BPS_junction_metric'] = sp.simplify(W/3 - sig.subs(d, 0)/6) == 0
checks['BPS_junction_scalar'] = sp.simplify(-Wp + sp.diff(sig, p).subs(d, 0)/2) == 0
# nonlinear junction expansions used in bc(): sigma(pb+f)-sigma(pb) and sigma'(pb+f)-sigma'(pb)
pb_, f = sp.symbols('pb f')
ds = sp.expand(sig.subs(p, pb_ + f) - sig.subs(p, pb_)); ds1 = sp.expand(sp.diff(sig, p).subs(p, pb_ + f) - sp.diff(sig, p).subs(p, pb_))
s10 = sp.diff(sig, p).subs(p, pb_); s20 = sp.diff(sig, p, 2).subs(p, pb_)
checks['bc_ds_matches'] = sp.simplify(ds - (s10*f + 2*pb_*f**2 + sp.Rational(2, 3)*f**3)) == 0
checks['bc_ds1_matches'] = sp.simplify(ds1 - (s20*f + 2*f**2)) == 0
# U coefficients: registered UC (in phi) and audited UETA (in eta=phi-1)
UC = [-sp.Rational(1, 6), sp.Rational(4, 3), -sp.Rational(5, 3), -sp.Rational(4, 9), sp.Rational(17, 18), 0, -sp.Rational(2, 27)]
checks['registered_UC_coefficients'] = sp.expand(U - sum(cc*p**k for k, cc in enumerate(UC))) == 0
e = sp.Symbol('eta')
UETA = [-sp.Rational(2, 27), 0, sp.Rational(14, 9), sp.Rational(50, 27), -sp.Rational(1, 6), -sp.Rational(4, 9), -sp.Rational(2, 27)]
checks['audited_UETA_coefficients'] = sp.expand(U.subs(p, 1 + e) - sum(cc*e**k for k, cc in enumerate(UETA))) == 0
checks['Upp_at_1_eq_28_over_9'] = sp.diff(U, p, 2).subs(p, 1) == sp.Rational(28, 9)
checks['U_at_1_eq_-2/27'] = U.subs(p, 1) == -sp.Rational(2, 27)
# pot_eta closed forms (td_core.pot_eta)
w = sp.Rational(1, 3) + e**2 + e**3/3; wp = e*(2 + e); pp = 1 + e
u = wp**2/2 - sp.Rational(2, 3)*w**2; up = wp*(2*pp - sp.Rational(4, 3)*w)
upp = 4*pp**2 + 2*wp - sp.Rational(4, 3)*(wp**2 + 2*pp*w)
checks['pot_eta_U'] = sp.expand(u - U.subs(p, 1 + e)) == 0
checks['pot_eta_Up'] = sp.expand(up - sp.diff(U, p).subs(p, 1 + e)) == 0
checks['pot_eta_Upp'] = sp.expand(upp - sp.diff(U, p, 2).subs(p, 1 + e)) == 0
# regular-cone series used in plus_background: rho = y + A y^3 + B y^5, eta = eh + C y^2 + D y^4
yy = sp.Symbol('y'); eh = sp.Symbol('eh')
u0 = U.subs(p, 1 + eh); u1 = sp.diff(U, p).subs(p, 1 + eh); u2 = sp.diff(U, p, 2).subs(p, 1 + eh)
Aa = -u0/36; Bb = u0**2/4320 - u1**2/750; Cc = u1/10; Dd = u1*(u2/280 + u0/630)
rs = yy + Aa*yy**3 + Bb*yy**5; es = eh + Cc*yy**2 + Dd*yy**4
Ue = lambda q: U.subs(p, 1 + q)
res_r = sp.series(sp.diff(rs, yy, 2) + rs*(sp.diff(es, yy)**2/4 + Ue(es)/6), yy, 0, 4).removeO()  # y^5 residual needs a y^7 term
res_e = sp.series(sp.diff(es, yy, 2) - sp.diff(U, p).subs(p, 1 + es) + 4*sp.diff(rs, yy)*sp.diff(es, yy)/rs, yy, 0, 4).removeO()
checks['cone_series_rho_residual_vanishes_through_y3'] = sp.simplify(sp.expand(res_r)) == 0
checks['cone_series_eta_through_y3'] = sp.simplify(sp.expand(res_e)) == 0
# ---- linearised junction data (linear_matrix ga, gf) from the nonlinear conditions
bsym, fsym, rb, s0 = sp.symbols('b f rb s0')
ga_nl = rb*(s0*(sp.exp(bsym) - 1) + sp.exp(bsym)*(s10*fsym + 2*pb_*fsym**2 + sp.Rational(2, 3)*fsym**3))/6
gf_nl = -rb*(s10*(sp.exp(bsym) - 1) + sp.exp(bsym)*(s20*fsym + 2*fsym**2))/2
lin = lambda ex: sp.expand(sp.series(sp.series(ex, bsym, 0, 2).removeO(), fsym, 0, 2).removeO())
lin_ga = lin(ga_nl); lin_gf = lin(gf_nl)
checks['linear_ga'] = sp.simplify(lin_ga.subs(bsym*fsym, 0) - rb*(s0*bsym + s10*fsym)/6) == 0
checks['linear_gf'] = sp.simplify(lin_gf.subs(bsym*fsym, 0) - (-rb*(s10*bsym + s20*fsym)/2)) == 0
checks['CONTROL_wrong_BPS_junction_sigma_over_3_detected'] = sp.simplify(W/3 - sig.subs(d, 0)/3) != 0
# ---- de Sitter separation: decoupled scalar f=e^{lam t}g(z), A=t+ln rho:  lam^2+3 lam = -mu^2
lam, mu2 = sp.symbols('lam mu2')
roots = sp.solve(lam**2 + 3*lam + mu2, lam)
checks['Re_lambda_is_-3/2_for_mu2_gt_9/4'] = all(sp.simplify(sp.re(r_.subs(mu2, sp.Rational(9, 4) + sp.Symbol('k', positive=True)**2)) + sp.Rational(3, 2)) == 0 for r_ in roots)
# Calibration consistency: GN mu2 -> lambda
mu2GN = -7.717871625176294
lamGN = -1.5 + (2.25 - mu2GN)**0.5
checks['GN_mu2_to_lambda'] = abs(lamGN - 1.6571936312453648) < 1e-12
# Flip control: lam=167.48925 -> partner
lamF = 167.4892453380038
checks['flip_partner_sum_-3'] = abs((-3 - lamF) - (-170.4892453380038)) < 1e-12

# ---- series expansions of the +1 branch quoted by the checkpoint (sanity of numbers only)
dd = 1e-3; cc = 0.5975949350280132
phi_series = 1 - 9*cc*dd/64
H2_series = dd*(1 + cc)/27 + dd**2*((1 + cc)**2/36 - cc**2/384)
out = dict(checks={k: bool(v) for k, v in checks.items()}, all_pass=all(bool(v) for v in checks.values()),
           extra=extra, lambda_from_GN_mu2=lamGN,
           plus_branch_series=dict(phi_b_first_order=phi_series, H2_second_order=H2_series,
                                   phi_b_archived=0.9999159473169134, H2_archived=5.92401479433e-5,
                                   phi_b_diff=phi_series - 0.9999159473169134, H2_rel_diff=H2_series/5.92401479433e-5 - 1))
Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=1, default=str) + '\n')
print(json.dumps(out, indent=1, default=str))
