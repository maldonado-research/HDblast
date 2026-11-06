"""Exact finite ADM polynomial identities, using only rational arithmetic.

Monomials are sorted tuples of variable names.  The epsilon coefficient
operation retains every other variable, so these are polynomial identities,
not identities inferred from a collection of numerical test points.
"""
from fractions import Fraction as F


class Polynomial:
    def __init__(self, terms=0):
        if isinstance(terms, Polynomial):
            self.terms = dict(terms.terms)
        elif isinstance(terms, dict):
            self.terms = {m: F(c) for m, c in terms.items() if c}
        else:
            self.terms = {(): F(terms)} if terms else {}

    def __add__(self, other):
        result = dict(self.terms)
        for monomial, coefficient in Polynomial(other).terms.items():
            result[monomial] = result.get(monomial, F(0)) + coefficient
        return Polynomial(result)

    __radd__ = __add__

    def __neg__(self):
        return Polynomial({m: -c for m, c in self.terms.items()})

    def __sub__(self, other):
        return self + -Polynomial(other)

    def __rsub__(self, other):
        return Polynomial(other) + -self

    def __mul__(self, other):
        result = {}
        for m, c in self.terms.items():
            for n, d in Polynomial(other).terms.items():
                monomial = tuple(sorted(m + n))
                result[monomial] = result.get(monomial, F(0)) + c*d
        return Polynomial(result)

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        return self * (1/F(scalar))

    def __pow__(self, degree):
        if type(degree) is not int or degree < 0:
            raise ValueError("Nonnegative polynomial power required")
        result = Polynomial(1)
        for _ in range(degree):
            result = result*self
        return result

    def __eq__(self, other):
        return self.terms == Polynomial(other).terms

    def coefficient(self, variable, degree):
        return Polynomial({tuple(v for v in m if v != variable): c
                           for m, c in self.terms.items()
                           if m.count(variable) == degree})


def variable(name):
    return Polynomial({(name,): F(1)})


def symbolic(check, data):
    """Recompute formal identities; check receives ``(name, boolean)``.

    ``data`` contains the four numeric claims returned below.  Each claim
    is substituted on the right of a recomputed polynomial equality.
    """
    def eq(actual, expected, name):
        check(name, actual == expected)

    def matrix():
        return [[Polynomial() for _ in range(4)] for _ in range(4)]

    def product(left, right):
        return [[sum(left[i][k]*right[k][j] for k in range(4))
                 for j in range(4)] for i in range(4)]

    def trace(value):
        return sum(value[i][i] for i in range(4))

    eps, A, a, f, S, p = [variable(name) for name in
                          ("eps", "A", "a", "f", "S", "p")]
    df = [variable(name) for name in ("fy", "f1", "f2", "f3")]
    base = [-4, 2, 2, 2]
    h, L, k, M = [matrix() for _ in range(4)]
    for i in range(4):
        h[i][i] = base[i]*f
        L[i][i] = (-2 if i == 0 else 1)*S
        for j in range(i, 4):
            k[i][j] = k[j][i] = variable("k" + str(i) + str(j))
            M[i][j] = M[j][i] = variable("M" + str(i) + str(j))
    connection = [matrix() for _ in range(4)]
    for b in range(1, 4):
        connection[b][0][b] = a
        connection[b][b][0] = -a
    dh = [matrix() for _ in range(4)]
    dk = [matrix() for _ in range(4)]
    for d in range(4):
        for i in range(4):
            for j in range(4):
                dh[d][i][j] = (df[d]*base[i] if i == j else Polynomial()) - sum(
                    connection[d][i][b]*h[b][j]
                    + connection[d][j][b]*h[i][b] for b in range(4))
            for j in range(i, 4):
                dk[d][i][j] = dk[d][j][i] = variable(
                    "D" + str(d) + "k" + str(i) + str(j))
    eh, Q, ek = [], [], []
    for i in range(4):
        eh.append(sum(dh[j][j][i] - dh[i][j][j]/2 for j in range(4)))
        ek.append(sum(dk[j][j][i] - dk[i][j][j]/2 for j in range(4)))
        Q.append(sum(h[j][l]*(dh[j][l][i] - dh[i][j][l]/2)
                     for j in range(4) for l in range(4)))
    eq(eh[0], -5*df[0]-18*a*f, "warped divergence radial")
    eq(Q[0], 2*f*df[0]-36*a*f*f, "warped quadratic contraction radial")
    for i in range(1, 4):
        eq(eh[i], df[i], "warped divergence angular " + str(i))
        eq(Q[i], -10*f*df[i], "warped quadratic angular " + str(i))
    hh = product(h, h)
    gamma = [[Polynomial(i == j)+eps*h[i][j]+eps**2*k[i][j]/2
              for j in range(4)] for i in range(4)]
    inverse = [[Polynomial(i == j)-eps*h[i][j]+eps**2*(hh[i][j]-k[i][j]/2)
                for j in range(4)] for i in range(4)]
    residual = product(gamma, inverse)
    for i in range(4):
        for j in range(4):
            eq(residual[i][j].coefficient("eps", 1), 0,
               "inverse first " + str(i) + str(j))
            eq(residual[i][j].coefficient("eps", 2), 0,
               "inverse second " + str(i) + str(j))
    radial = [F(x) for x in data["mixed_second_radial_coefficients"]]
    angular = F(data["mixed_second_angular_coefficient"])
    temporal_cross = F(data["temporal_second_cross_coefficient"])
    scalar_pair = [F(x) for x in data["coordinate_second_scalar_pair"]]
    for i in range(4):
        C = sum(inverse[j][l]*(eps*(dh[j][l][i]-dh[i][j][l]/2)
                    + eps**2*(dk[j][l][i]-dk[i][j][l]/2)/2)
                for j in range(4) for l in range(4))
        reference = gamma[i][0]*A*A*a
        lapse_gradient = A*A*(1+eps*f)*(
            a*(1+eps*f)*(1 if i == 0 else 0)+eps*df[i])
        J = A*A*(1+eps*f)**2*C+reference-lapse_gradient
        eq(J.coefficient("eps", 0), 0,
           "background shift cancellation " + str(i))
        first = -6*A*A*(df[0]+4*a*f) if i == 0 else Polynomial()
        eq(J.coefficient("eps", 1), first, "first mixed velocity " + str(i))
        quad = (radial[0]*f*df[0]+radial[1]*a*f*f if i == 0
                else angular*f*df[i])
        expected = A*A*(ek[i]+a*k[i][0]+quad)
        eq(2*J.coefficient("eps", 2), expected,
           "general second mixed velocity " + str(i))
        check("omitted-reference corruption rejected " + str(i),
              2*(J-reference).coefficient("eps", 2) != expected)
    physical_K = [[eps*L[i][j]+eps**2*M[i][j]/2
                   for j in range(4)] for i in range(4)]
    trace_K = trace(product(inverse, physical_K))
    eq(trace(L), S, "signed first trace K")
    eq(trace(product(h, L)), 14*f*S,
       "first metric extrinsic-curvature contraction")
    eq(2*trace_K.coefficient("eps", 2), trace(M)-28*f*S,
       "second physical trace")
    alpha_t = A*A*(1+eps*f)**2*trace_K
    eq(2*alpha_t.coefficient("eps", 2), A*A*(trace(M)-24*f*S),
       "second lapse velocity")
    gtt_t = -2*A*(1+eps*f)*alpha_t
    eq(2*gtt_t.coefficient("eps", 2),
       -2*A**3*trace(M)+temporal_cross*A**3*f*S,
       "second temporal metric velocity")
    eq(2*(-A*A*(1+eps*f)**2).coefficient("eps", 2), -2*A*A*f*f,
       "second temporal metric value")
    phi_t = A*(1+eps*f)*(eps*p-eps**2*f*p)
    eq(2*phi_t.coefficient("eps", 2), scalar_pair[1],
       "coordinate scalar velocity cancellation")
    phi_initial = variable("phi_background")+eps*variable("phi_first")
    eq(2*phi_initial.coefficient("eps", 2), scalar_pair[0],
       "coordinate scalar value cancellation")
    check("generic second lapse freedom changes wave data",
          -2*A*variable("free_alpha2") != 0)
    return {
        "mixed_second_radial_coefficients": [str(x) for x in radial],
        "mixed_second_angular_coefficient": str(angular),
        "temporal_second_cross_coefficient": str(temporal_cross),
        "coordinate_second_scalar_pair": [str(x) for x in scalar_pair],
        "method": "formal sparse rational polynomial identities",
    }
