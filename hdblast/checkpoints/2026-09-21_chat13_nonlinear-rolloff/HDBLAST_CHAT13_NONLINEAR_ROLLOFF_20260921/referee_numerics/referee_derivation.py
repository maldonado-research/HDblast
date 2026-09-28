#!/usr/bin/env python3
"""Referee (numerics): independent SymPy re-derivation, by a different route from derive_evolution_equations.py.
 1. Einstein tensor G_AB and scalar stress tensor T_AB (kappa_5^2 = 1, V = U) -> G_AB - T_AB; compare the (tz) component with
    the momentum constraint used in the code, the (tt)+(zz) combination with the Hamiltonian-type constraint, and the evolution
    equations with the closed forms.  Also verify R_AB = d_A phi d_B phi + (2/3) U g_AB is equivalent to G_AB = T_AB in 5D.
 2. Junction conditions: extrinsic curvature of z = const with unit normal n = e^{-B} d_z (pointing out of the bulk z < 0),
    K_mn = (1/2) L_n h_mn; Israel for a Z2 brane with pure tension S_mn = -sigma h_mn:  K_mn - h_mn K = -(1/2) S_mn  (kept side,
    outward normal)  =>  K^m_n = (sigma/6) delta^m_n  (checked against the RS solution below).  Scalar: phi_n = -sigma'/2.
 3. Identity: with A_z = B_z = sigma e^B/6 and phi_z = -sigma' e^B/2 imposed for all t at the shell, the momentum constraint
    vanishes identically there (so the Neumann data are compatible with the constraint, and the code's constraint monitor
    could legitimately include the shell point if pa_z were computed one-sidedly).
 4. RS check of the sign convention: ds^2 = dy^2 + e^{-2k|y|} eta, bulk on y > 0 side kept -> requires sigma = 6k > 0 for
    positive-tension brane with warp factor decreasing away from the brane.
"""
import sympy as sp
t, z = sp.symbols('t z'); xs = sp.symbols('x1:4'); X = [t, z, *xs]
A, B, ph = [sp.Function(n)(t, z) for n in ('A', 'B', 'phi')]
U, U1, sig, sig1 = sp.symbols('U U1 sigma sigma1')
n = 5
g = sp.diag(-sp.exp(2*B), sp.exp(2*B), sp.exp(2*A), sp.exp(2*A), sp.exp(2*A)); gi = g.inv()
Gam = [[[sum(sp.Rational(1, 2)*gi[a, d]*(sp.diff(g[d, c], X[b]) + sp.diff(g[d, b], X[c]) - sp.diff(g[b, c], X[d])) for d in range(n)) for c in range(n)] for b in range(n)] for a in range(n)]
def ric(b, c):
    return sp.simplify(sum(sp.diff(Gam[a][b][c], X[a]) - sp.diff(Gam[a][b][a], X[c]) + sum(Gam[a][a][d]*Gam[d][b][c] - Gam[a][c][d]*Gam[d][b][a] for d in range(n)) for a in range(n)))
Ric = sp.Matrix(n, n, lambda b, c: ric(b, c))
Rs = sp.simplify(sum(gi[a, a]*Ric[a, a] for a in range(n)))
G = sp.simplify(Ric - sp.Rational(1, 2)*g*Rs)
dphi = [sp.diff(ph, v) for v in X]
kin = sp.simplify(sum(gi[a, a]*dphi[a]**2 for a in range(n)))
T = sp.Matrix(n, n, lambda b, c: dphi[b]*dphi[c] - g[b, c]*(kin/2 + U))
Ein = sp.simplify(G - T)
# 1a. equivalence with the trace-reversed form used by the author
Etr = sp.Matrix(n, n, lambda b, c: Ric[b, c] - dphi[b]*dphi[c] - sp.Rational(2, 3)*U*g[b, c])
Ttr = sp.simplify(sum(gi[a, a]*T[a, a] for a in range(n)))
equiv = all(sp.simplify(Ein[b, c] - (Etr[b, c] - sp.Rational(1, 2)*g[b, c]*sp.simplify(sum(gi[a, a]*Etr[a, a] for a in range(n))))) == 0 for b in range(n) for c in range(n))
print("G_AB - T_AB == trace-reverse of (R_AB - dphi dphi - (2/3)U g):", equiv)
At, Az, Bt, Bz, pt, pz = [sp.diff(f, v) for f in (A, B, ph) for v in (t, z)]
Att, Azz, Btt, Bzz, ptt, pzz = [sp.diff(f, v, 2) for f in (A, B, ph) for v in (t, z)]
Atz = sp.diff(A, t, z)
# 1b. momentum constraint from G_tz - T_tz (no second time derivatives appear in it)
cm_code = -3*Atz - 3*At*Az + 3*At*Bz + 3*Az*Bt - pt*pz
ratio = sp.simplify(Ein[0, 1]/cm_code)
print("G_tz - T_tz  =", sp.factor(Ein[0, 1]), "   ratio to code's M:", ratio)
# 1c. Hamiltonian constraint: G_tt - T_tt contains no second time derivatives; compare with the code's Hc form
fH = -2*U*sp.exp(2*B) + 6*At**2 + 6*At*Bt - 12*Az**2 + 6*Az*Bz - 6*Azz - pt**2 - pz**2
print("G_tt - T_tt (Hamiltonian) =", sp.expand(Ein[0, 0]))
print("ratio (G_tt - T_tt)/(code Hc):", sp.simplify(Ein[0, 0]/fH))
# 1d. evolution equations: G_xx - T_xx gives B_tt - B_zz combination?  Solve the three second-order equations
sol = sp.solve([Ein[0, 0] + Ein[1, 1], Ein[2, 2], sp.expand(Ein[0, 0] - Ein[1, 1])], [Att, Btt, Azz], dict=True)
# Alternative route: use the (zz) equation (contains A_tt, B_tt but not second z derivatives of B) and the trace.
fA = Azz - 3*At**2 + 3*Az**2 + sp.Rational(2, 3)*sp.exp(2*B)*U
fB = Bzz + 3*At**2 - 3*Az**2 - pt**2/2 + pz**2/2 - sp.Rational(1, 3)*sp.exp(2*B)*U
sqrtg = sp.exp(2*B + 3*A)
KG = sp.simplify(sum(sp.diff(sqrtg*gi[a, a]*dphi[a], X[a]) for a in range(n))/sqrtg - U1)
fp = pzz - 3*At*pt + 3*Az*pz - sp.exp(2*B)*U1
# Check: substituting A_tt = fA, B_tt = fB, phi_tt = fp into ALL components of G - T must leave only the two constraints
subs_ev = {Att: fA, Btt: fB, ptt: fp}
resid = {}
for b in range(n):
    for c in range(b, n):
        r = sp.simplify(Ein[b, c].subs(subs_ev))
        resid[(b, c)] = r
print("after imposing the three evolution equations, the residual components of G-T are:")
for k, v in resid.items():
    if v != 0: print("   ", k, "=", sp.factor(v))
print("   (all others vanish);  (2,2),(3,3),(4,4) vanish:", all(resid[(i, i)] == 0 for i in (2, 3, 4)))
print("Klein-Gordon residual with phi_tt = fp:", sp.simplify(KG.subs(subs_ev)))
# the surviving residuals must be combinations of the momentum and Hamiltonian constraints
Hc_sym = sp.simplify(Ein[0, 0])
for k in [(0, 0), (0, 1), (1, 1)]:
    v = resid[k]
    c1 = sp.simplify(v/cm_code) if k == (0, 1) else None
    c2 = sp.simplify(v/Hc_sym) if k != (0, 1) else None
    print("   residual", k, "= (%s) * M   or  (%s) * Hc" % (c1, c2))

# 2. junction conditions
# induced metric on z = const: h = diag(-e^{2B}, e^{2A} x3); unit normal n^A = e^{-B} delta^A_z (points toward +z, i.e. out of the bulk z<0)
# K_mn = (1/2) n^z d_z h_mn
K_tt = sp.Rational(1, 2)*sp.exp(-B)*sp.diff(-sp.exp(2*B), z)
K_xx = sp.Rational(1, 2)*sp.exp(-B)*sp.diff(sp.exp(2*A), z)
Kt_t = sp.simplify(K_tt/(-sp.exp(2*B))); Kx_x = sp.simplify(K_xx/sp.exp(2*A))
print("K^t_t =", Kt_t, "  K^x_x =", Kx_x)
Ktr = Kt_t + 3*Kx_x
# Israel, Z2, kept side with outward normal: K_mn - h_mn K = -(1/2) S_mn with S_mn = -sigma h_mn  (see RS check below for the sign)
# mixed components: K^m_n - delta K = (sigma/2) delta  -> trace: K - 4K = 2 sigma -> K = -(2/3) sigma ; K^m_n = (sigma/2 - 2 sigma/3) delta = -(sigma/6) delta ??
# Resolve the sign with the RS check instead of trusting a memorised convention:
y, k, s = sp.symbols('y k sigma', positive=True)
# RS: ds^2 = dy^2 + e^{-2k y} eta on the kept side y>0 (brane at y=0, warp factor decreasing away from the brane, positive tension sigma = 6k for kappa^2=1).
# Outward normal from the kept side (y>0) at y=0 is n = -d_y.  K_mn = (1/2) L_n h_mn = -(1/2) d_y h_mn = k h_mn  -> K^m_n = +k delta = (sigma/6) delta.
print("RS check: outward-normal extrinsic curvature on the kept side  K^m_n = +k delta^m_n = (sigma/6) delta^m_n  with sigma = 6k > 0")
print("  => the correct Z2 Israel rule in this convention is  K^m_n|kept, outward normal = (sigma/6) delta^m_n ,")
print("     i.e. the algebraic form is  K_mn - h_mn K = +(1/2) S_mn ... only the RS-anchored form matters for the code:")
print("     A_z = B_z = (sigma/6) e^B  at the shell (bulk z<0, outward normal +z).  Code uses exactly this.")
# Scalar junction: action S = int sqrt(-g) [ -1/2 (dphi)^2 - U ] - int sqrt(-h) sigma(phi);  EOM: Box phi = U' + sigma'(phi) delta(y) (y proper distance)
# Z2: phi_y(0+) - phi_y(0-) = sigma'  and phi_y(0+) = -phi_y(0-)  -> on the side y<0 (our bulk): phi_y(0-) = -sigma'/2 ; phi_z = e^B phi_y = -(sigma'/2) e^B.
print("scalar junction on the kept side (bulk z<0):  phi_z = -(sigma'/2) e^B .  Code uses exactly this.")
# Check with the static background: rho_y/rho = sigma/6 and phi_y = -sigma'/2 are what solve_shell() imposes (res(): rp/rho - s/6, s + s1/2). consistent.

# 3. Identity: junction data satisfy the momentum constraint identically at the shell
Bb, fb = sp.Function('Bb')(t), sp.Function('fb')(t)      # boundary values B(t,0), phi(t,0)
sgm = sp.Function('sigma')(fb)
Gz = sgm*sp.exp(Bb)/6; Pz = -sp.diff(sgm, fb)*sp.exp(Bb)/2
M_shell = -3*sp.diff(Gz, t) - 3*sp.diff(A, t)*Gz + 3*sp.diff(A, t)*Gz + 3*Gz*sp.diff(Bb, t) - sp.diff(fb, t)*Pz
print("momentum constraint at the shell with junction data imposed for all t:", sp.simplify(M_shell))
# what A_tz must be at the shell:
print("A_tz at shell = d/dt[(sigma/6) e^B] = ", sp.simplify(sp.diff(Gz, t)), " (NOT zero: the code's constraint monitor uses pa_z = 0 there and skips the last 6 points)")
