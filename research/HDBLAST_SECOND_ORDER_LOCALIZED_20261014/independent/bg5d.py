"""Independent re-derivation of the background relations from the 5D field
equations (shares no code with engine.py).

g = Om(S)^2 eta (signature -++++), phi = Ph(S), S = -t^2 + x1^2+...+x4^2.
Field equations (kappa_5 = 1): R_AB = d_A phi d_B phi + (2/3) U g_AB,
Box phi = U'(phi).
We compute R_AB symbolically, write E_AB = R_AB - dphi dphi = alpha eta + beta X X,
and return U(S) = (3/2) alpha / Om^2, the ODE beta = 0, and U_phi = Box phi.
"""
import sympy as sp

t, x1, x2, x3, x4 = sp.symbols('t x1 x2 x3 x4', real=True)
C = [t, x1, x2, x3, x4]
eta = sp.diag(-1, 1, 1, 1, 1)
S = sp.Symbol('S', positive=True)
Om = sp.Function('Om'); Ph = sp.Function('Ph')

def derive():
    Sx = -t**2 + x1**2 + x2**2 + x3**2 + x4**2
    g = Om(Sx)**2 * eta
    gi = eta / Om(Sx)**2
    n = 5
    dg = [[[sp.diff(g[a, b], C[c]) for c in range(n)] for b in range(n)] for a in range(n)]
    Gam = [[[sum(gi[k, l] * (dg[l][j][i] + dg[l][i][j] - dg[i][j][l]) for l in range(n)) / 2
             for j in range(n)] for i in range(n)] for k in range(n)]
    def Ric(i, j):
        r = 0
        for k in range(n):
            r += sp.diff(Gam[k][i][j], C[k]) - sp.diff(Gam[k][i][k], C[j])
            for l in range(n):
                r += Gam[k][k][l] * Gam[l][i][j] - Gam[k][j][l] * Gam[l][i][k]
        return r
    phi = Ph(Sx)
    dphi = [sp.diff(phi, c) for c in C]
    # evaluate at a generic point on the t = 0 slice: components (x1,x1) and (x1,x2)
    def at(e):
        e = e.subs(t, 0).doit()
        return e
    E11 = at(Ric(1, 1) - dphi[1] * dphi[1])
    E12 = at(Ric(1, 2) - dphi[1] * dphi[2])
    # E_ij = alpha delta_ij + beta x_i x_j  on t=0 (eta_ij = delta_ij, X_i = x_i)
    rep = {x1: 0, x2: 0, x3: 0, x4: sp.sqrt(S)}
    alpha = sp.simplify(E11.subs(rep))
    beta = sp.simplify((E12 / (x1 * x2)).subs({x1: 1, x2: 1, x3: 0, x4: sp.sqrt(S - 2)}))
    # Box phi = (1/sqrt|g|) d_A (sqrt|g| g^AB d_B phi)
    sqg = Om(Sx)**5
    box = sum(sp.diff(sqg * gi[a, a] * dphi[a], C[a]) for a in range(n)) / sqg
    box = sp.simplify(at(box).subs(rep))
    return alpha, beta, box

if __name__ == '__main__':
    a, b, bx = derive()
    print('alpha =', a); print('beta =', b); print('box =', bx)
