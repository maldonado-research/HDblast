#!/usr/bin/env python3
"""Chat 13 - derive (SymPy, exact) the 1+1 evolution system used by the nonlinear code.
Metric (2D conformal gauge, flat 3-slices):   ds^2 = e^{2B(t,z)} (-dt^2 + dz^2) + e^{2A(t,z)} dx_3^2 ,   phi = phi(t,z),
field equations  E_AB := R_AB - d_A phi d_B phi - (2/3) U(phi) g_AB = 0,   Box phi = U'(phi).
Static registered solution in this chart:  B = ln rho(z),  A = ln rho(z) + t,  phi = phi0(z),  dz = dy/rho.
Z2 shell at z = 0 (bulk on z < 0), unit normal e^{-B} d_z:   K^x_x = e^{-B} A',  K^t_t = e^{-B} B'  =>  A' = B' = (sigma_t/6) e^{B},  phi' = -(sigma_t'/2) e^{B}."""
import sympy as sp, json
t, z, x1, x2, x3 = sp.symbols('t z x1 x2 x3'); X = [t, z, x1, x2, x3]
A, B, ph = [sp.Function(n)(t, z) for n in ('A', 'B', 'phi')]
U, U1 = sp.symbols('U U1')
g = sp.diag(-sp.exp(2*B), sp.exp(2*B), sp.exp(2*A), sp.exp(2*A), sp.exp(2*A)); gi = g.inv(); n = 5
Gam = [[[sum(sp.Rational(1, 2)*gi[a, d]*(sp.diff(g[d, c], X[b]) + sp.diff(g[d, b], X[c]) - sp.diff(g[b, c], X[d])) for d in range(n)) for c in range(n)] for b in range(n)] for a in range(n)]
def ric(b, c):
    return sp.simplify(sum(sp.diff(Gam[a][b][c], X[a]) - sp.diff(Gam[a][b][a], X[c]) + sum(Gam[a][a][d]*Gam[d][b][c] - Gam[a][c][d]*Gam[d][b][a] for d in range(n)) for a in range(n)))
dphi = [sp.diff(ph, v) for v in X]
E = {(b, c): sp.simplify(ric(b, c) - dphi[b]*dphi[c] - sp.Rational(2, 3)*U*g[b, c]) for b in range(n) for c in range(b, n)}
At, Az, Bt, Bz, pt, pz = [sp.diff(f, v) for f in (A, B, ph) for v in (t, z)]
Att, Azz, Btt, Bzz, ptt, pzz = [sp.diff(f, v, 2) for f in (A, B, ph) for v in (t, z)]
# evolution equations: solve for the wave operators
box = lambda f: sp.diff(f, t, 2) - sp.diff(f, z, 2)
eA = sp.solve(sp.Eq(E[(2, 2)], 0), Att)[0]                          # from the x-x equation
print("A_tt =", sp.simplify(eA))
sqrtg = sp.exp(2*B + 3*A)
KG = sp.simplify(sum(sp.diff(sqrtg*gi[a, a]*dphi[a], X[a]) for a in range(n))/sqrtg - U1)
ep = sp.solve(sp.Eq(KG, 0), ptt)[0]
print("phi_tt =", sp.simplify(ep))
comb = sp.simplify((E[(0, 0)] - E[(1, 1)]))                          # tt - zz combination contains B_tt - B_zz
eB = sp.solve(sp.Eq(comb.subs(Att, eA), 0), Btt)[0]
print("B_tt =", sp.simplify(eB))
C1 = sp.simplify(E[(0, 1)]); C2 = sp.simplify((E[(0, 0)] + E[(1, 1)]).subs(Att, eA))
print("momentum constraint  E_tz =", C1)
print("Hamiltonian-type constraint  E_tt + E_zz (A_tt eliminated) =", C2)
# checks against the closed forms used in the code
fA = Azz - 3*At**2 + 3*Az**2 + sp.Rational(2, 3)*sp.exp(2*B)*U
fp = pzz - 3*At*pt + 3*Az*pz - sp.exp(2*B)*U1
fB = Bzz + 3*At**2 - 3*Az**2 - pt**2/2 + pz**2/2 - sp.Rational(1, 3)*sp.exp(2*B)*U
ok = dict(A=sp.simplify(eA - fA) == 0, phi=sp.simplify(ep - fp) == 0, B=sp.simplify(eB - fB) == 0)
cm = -3*sp.diff(A, t, z) - 3*At*Az + 3*At*Bz + 3*Az*Bt - pt*pz
ok["momentum"] = sp.simplify(C1 - cm) == 0
print("closed forms verified:", ok)
# static solution check: B = ln rho(z), A = ln rho(z) + t, phi0(z) with the registered background ODEs in z
rho = sp.Function('rho')(z); p0 = sp.Function('phi0')(z)
sub = {A: sp.log(rho) + t, B: sp.log(rho), ph: p0}
Hc = sp.diff(rho, z)/rho
bgsub = {sp.diff(rho, z, 2): rho*(2*Hc**2 - 1 - sp.diff(p0, z)**2/3), sp.diff(p0, z, 2): rho**2*U1 - 3*Hc*sp.diff(p0, z)}     # background ODEs in z
Ucon = (sp.diff(p0, z)**2/2 - 6*Hc**2 + 6)/rho**2                                                                          # Hamiltonian constraint
fH = -2*U*sp.exp(2*B) + 6*At**2 + 6*At*Bt - 12*Az**2 + 6*Az*Bz - 6*Azz - pt**2 - pz**2
ok["hamiltonian"] = sp.simplify(C2 - fH) == 0
stat = [sp.simplify(e.subs(sub).doit().subs(bgsub).subs(U, Ucon)) for e in (fA - Att, fB - Btt, fp - ptt, cm, fH)]
ok["static_solution_solves_everything"] = all(s_ == 0 for s_ in stat)
print("static residuals on the background ODEs:", stat)
print("all checks:", ok)
json.dump({k: bool(v) for k, v in ok.items()}, open("EVOLUTION_EQUATIONS_CHECK.json", "w"), indent=1)
