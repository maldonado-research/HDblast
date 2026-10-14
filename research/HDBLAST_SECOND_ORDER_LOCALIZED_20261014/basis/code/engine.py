"""Exact jet engine for radial Einstein-scalar initial data on the T=0 slice.

Background (two identical regular interiors of the HDBLAST shell model):
    5D metric g = Omega(s)^2 eta,  phi = Phi(s),  s = eta_AB X^A X^B.
At T=0 the slice is time symmetric: gamma0 = Omega(r^2)^2 delta_ij on R^4,
phi0 = Phi(r^2), k = 0, d_T phi = 0.

Conventions (match the frozen HDBLAST capsules): kappa_5 = 1,
    S = int sqrt(-g) [R/2 - (grad phi)^2/2 - U(phi)],
    R_AB = phi_A phi_B + (2/3) U g_AB,   Box phi = U'.

Every scalar field used here is a polynomial in x1..x4 whose coefficients are
rational functions of "jet symbols": values at s=|x|^2 of the s-derivatives
of the radial functions.  d/dx_i acts by the chain rule
    d/dx_i J_k = 2 x_i J_{k+1}.
Background equations are imposed by substitution rules (see bg_rules), so an
identity that reduces to zero holds for every actual background solution.
"""
import sympy as sp

X = sp.symbols('x1:5')
EPS = sp.Symbol('eps')
NJ = 9  # number of jet orders kept

def jets(name, n=NJ, start=0):
    return [sp.Symbol(f'{name}{k}') for k in range(start, n)]

O = jets('O')          # Omega^(k)(s)
F = jets('F')          # Phi^(k)(s); F0 never used (only derivatives matter)
U = jets('U')          # U^(k)(phi) evaluated on the background
_families = {}

def register(name, n=NJ, kind='s'):
    """Register a radial function of s; returns its jet list."""
    J = jets(name, n)
    _families[name] = (J, kind)
    return J

for nm in ('O', 'F'):
    _families[nm] = (O if nm == 'O' else F, 's')
_families['U'] = (U, 'U')

def _succ():
    m = {}
    for name, (J, kind) in _families.items():
        for k in range(len(J) - 1):
            if kind == 's':
                m[J[k]] = J[k + 1]
            else:  # U^(k)(phi(s)) -> U^(k+1) * Phi'
                m[J[k]] = J[k + 1] * F[1]
    return m

def d(expr, i):
    """Exact x_i derivative using the radial chain rule."""
    expr = sp.sympify(expr)
    succ = _succ()
    out = sp.diff(expr, X[i])
    syms = expr.free_symbols
    for J, Jn in succ.items():
        if J in syms:
            out += 2 * X[i] * Jn * sp.diff(expr, J)
    return out

S2 = sum(xx**2 for xx in X)  # s as a polynomial

# ---------------------------------------------------------------- background
def bg_rules():
    """Substitution rules imposing the 5D background equations.

    For g = e^{2w} eta in five dimensions with w = ln Omega(s):
        Omega'' = 2 Omega'^2/Omega - Omega Phi'^2/3,
        U   = -(3/(2 Omega^2)) (16 w' + 4 s w'' + 12 s w'^2),
        U'  = Omega^-2 (10 Phi' + 4 s Phi'' + 12 s Phi' Omega'/Omega).
    Higher U-derivatives follow by differentiating along s and dividing by
    Phi' (U is a function of phi and Phi' != 0 on the domain).
    """
    s = sp.Symbol('s_')
    O0, O1 = O[0], O[1]
    rules = {}
    O2 = 2 * O1**2 / O0 - O0 * F[1]**2 / 3
    rules[O[2]] = O2
    # successive Omega derivatives
    cur = O2
    for k in range(3, NJ):
        cur = sp.together(dS(cur, rules))
        rules[O[k]] = cur
    w1 = O1 / O0
    w2 = O2 / O0 - O1**2 / O0**2
    U0 = -sp.Rational(3, 2) / O0**2 * (16 * w1 + 4 * s * w2 + 12 * s * w1**2)
    U1 = (10 * F[1] + 4 * s * F[2] + 12 * s * F[1] * O1 / O0) / O0**2
    rules[U[0]] = U0
    rules[U[1]] = U1
    cur = U1
    for k in range(2, 5):
        cur = sp.together(dS(cur, rules) / F[1])
        rules[U[k]] = cur
    return rules, s

def dS(expr, rules=None):
    """Total s-derivative of a jet expression (with explicit symbol s_)."""
    s = sp.Symbol('s_')
    succ = _succ()
    out = sp.diff(expr, s)
    for J, Jn in succ.items():
        if J in expr.free_symbols:
            out += Jn * sp.diff(expr, J)
    if rules:
        out = out.subs(rules)
    return out

_RULES = None

def reduce_bg(expr):
    """Impose background equations and s=|x|^2; return a canonical form."""
    global _RULES
    if _RULES is None:
        _RULES = bg_rules()
    rules, s = _RULES
    e = expr.subs(rules).subs(s, S2)
    return e

def is_zero(expr):
    e = sp.together(reduce_bg(expr))
    num = sp.numer(e)
    num = sp.expand(num)
    return num == 0, num

# ------------------------------------------------------- epsilon-jet algebra
class J3:
    """Second-order Taylor jet in eps: (value, first, second derivative).

    If J3.FIRST_ONLY is True the second-order slot is not propagated
    (exact for linearized operators, much cheaper)."""
    FIRST_ONLY = False
    __slots__ = ('a', 'b', 'c')
    def __init__(self, a, b=0, c=0):
        self.a, self.b, self.c = sp.sympify(a), sp.sympify(b), sp.sympify(c)
    def __add__(self, o):
        o = o if isinstance(o, J3) else J3(o)
        return J3(self.a + o.a, self.b + o.b, self.c + o.c)
    __radd__ = __add__
    def __neg__(self):
        return J3(-self.a, -self.b, -self.c)
    def __sub__(self, o):
        return self + (-(o if isinstance(o, J3) else J3(o)))
    def __rsub__(self, o):
        return J3(o) - self
    def __mul__(self, o):
        o = o if isinstance(o, J3) else J3(o)
        if J3.FIRST_ONLY:
            return J3(self.a * o.a, self.a * o.b + self.b * o.a, 0)
        return J3(self.a * o.a, self.a * o.b + self.b * o.a,
                  self.a * o.c + 2 * self.b * o.b + self.c * o.a)
    __rmul__ = __mul__
    def d(self, i):
        return J3(d(self.a, i), d(self.b, i), d(self.c, i))
    def map(self, f):
        return J3(f(self.a), f(self.b), f(self.c))

def j3zero():
    return J3(0, 0, 0)

def inv_metric_jet(g0, h):
    """eps-jet of (g0 + eps h)^{-1}; g0 diagonal-invertible sympy Matrix."""
    gi = g0.inv()
    first = -gi * h * gi
    second = 2 * gi * h * gi * h * gi
    return [[J3(gi[i, j], first[i, j], second[i, j]) for j in range(4)]
            for i in range(4)]

def christoffel(g, gi):
    """g, gi: 4x4 lists of J3. Returns Gam[k][i][j]."""
    dg = [[[g[i][j].d(k) for k in range(4)] for j in range(4)] for i in range(4)]
    Gam = [[[j3zero() for _ in range(4)] for _ in range(4)] for _ in range(4)]
    for k in range(4):
        for i in range(4):
            for j in range(i, 4):
                acc = j3zero()
                for l in range(4):
                    acc = acc + gi[k][l] * (dg[l][j][i] + dg[l][i][j] - dg[i][j][l])
                acc = acc * sp.Rational(1, 2)
                Gam[k][i][j] = acc
                Gam[k][j][i] = acc
    return Gam

def ricci(Gam):
    R = [[j3zero() for _ in range(4)] for _ in range(4)]
    for i in range(4):
        for j in range(i, 4):
            acc = j3zero()
            for k in range(4):
                acc = acc + Gam[k][i][j].d(k) - Gam[k][i][k].d(j)
                for l in range(4):
                    acc = acc + Gam[k][k][l] * Gam[l][i][j] - Gam[k][j][l] * Gam[l][i][k]
            R[i][j] = acc
            R[j][i] = acc
    return R
