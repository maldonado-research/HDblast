"""Independent exact S^3 angular checks; no producer/oracle imports."""
from fractions import Fraction
from math import prod
import sympy as sp

def sphere_monomial_moment(powers):
    """Uniform S^3 moment, independently via the Gaussian radial identity."""
    if any((p % 2 for p in powers)):
        return Fraction(0)
    halves = [p // 2 for p in powers]
    degree = sum(halves)
    numerator = prod((prod(range(1, 2 * p, 2)) for p in halves))
    denominator = prod((4 + 2 * k for k in range(degree)))
    return Fraction(numerator, denominator)

def need(condition, label):
    if not condition:
        raise ValueError(label)

def main():
    z = sp.symbols('z', positive=True)
    phi, chi = sp.symbols('phi chi', real=True)
    n = sp.Matrix([sp.sqrt(z) * sp.cos(phi), sp.sqrt(z) * sp.sin(phi), sp.sqrt(1 - z) * sp.cos(chi), sp.sqrt(1 - z) * sp.sin(chi)])
    dn = n.jacobian([z, phi, chi])
    metric = sp.simplify(sp.trigsimp(dn.T * dn))
    expected = sp.diag(1 / (4 * z * (1 - z)), z, 1 - z)
    need(sp.simplify(metric - expected) == sp.zeros(3), 'ANGULAR_IDENTITY_LINE_33')
    need(sp.simplify(metric.det()) == sp.Rational(1, 4), 'ANGULAR_IDENTITY_LINE_34')
    y_hopf = sp.trigsimp(n[0] ** 2 - n[1] ** 2)
    need(sp.trigsimp(y_hopf - z * sp.cos(2 * phi)) == 0, 'ANGULAR_IDENTITY_LINE_36')
    n0, n1, n2, n3 = sp.symbols('n0 n1 n2 n3', real=True)
    coords = sp.Matrix([n0, n1, n2, n3])
    y = n0 ** 2 - n1 ** 2
    zz = n0 ** 2 + n1 ** 2
    ambient = sp.Matrix([sp.diff(y, v) for v in coords])
    projected_norm = sp.expand(ambient.dot(ambient) - coords.dot(ambient) ** 2)
    need(sp.expand(projected_norm - 4 * (zz - y ** 2)) == 0, 'ANGULAR_IDENTITY_LINE_45')
    e_n04 = sphere_monomial_moment([4, 0, 0, 0])
    e_n02_n12 = sphere_monomial_moment([2, 2, 0, 0])
    e_y2 = 2 * e_n04 - 2 * e_n02_n12
    e_z = 2 * sphere_monomial_moment([2, 0, 0, 0])
    e_grad2 = 4 * (e_z - e_y2)
    need(e_n04 == Fraction(1, 8), 'ANGULAR_IDENTITY_LINE_52')
    need(e_n02_n12 == Fraction(1, 24), 'ANGULAR_IDENTITY_LINE_53')
    need(e_y2 == Fraction(1, 6), 'ANGULAR_IDENTITY_LINE_54')
    need(e_grad2 == Fraction(4, 3), 'ANGULAR_IDENTITY_LINE_55')
    r = sp.symbols('r', positive=True)
    theta = sp.symbols('theta', real=True)
    primitive = theta - r * sp.tan(theta)
    need(sp.trigsimp(sp.diff(primitive, theta) - (1 - r / sp.cos(theta) ** 2)) == 0, 'ANGULAR_IDENTITY_LINE_60')
    alpha = sp.acos(sp.sqrt(r))
    endpoint = sp.acos(sp.sqrt(r)) - sp.sqrt(r * (1 - r))
    need(sp.simplify(primitive.subs(theta, alpha).rewrite(sp.sqrt) - endpoint) == 0, 'ANGULAR_IDENTITY_LINE_64')
    f = 2 * endpoint / sp.pi
    need(sp.limit(f, r, 0, dir='+') == 1, 'ANGULAR_IDENTITY_LINE_68')
    need(sp.limit(f, r, 1, dir='-') == 0, 'ANGULAR_IDENTITY_LINE_69')
    x = sp.symbols('x', positive=True)
    tail_bracket = sp.asin(x) + x * sp.sqrt(1 - x * x)
    need(sp.simplify(sp.diff(tail_bracket, x) - 2 * sp.sqrt(1 - x * x)) == 0, 'ANGULAR_IDENTITY_LINE_73')
    need(tail_bracket.subs(x, 0) == 0, 'ANGULAR_IDENTITY_LINE_74')
    delta = Fraction(1, 1000)
    need(1 - Fraction(4, 3) * delta == Fraction(749, 750), 'ANGULAR_IDENTITY_LINE_76')
    a, b = sp.symbols('a b', real=True)
    g = sp.expand(-a * a * y * y + b * b * projected_norm)
    need(sp.expand(g - (4 * b * b * zz - (a * a + 4 * b * b) * y * y)) == 0, 'ANGULAR_IDENTITY_LINE_80')
    need(g.subs({n0: 1, n1: 0, n2: 0, n3: 0}) == -a * a, 'ANGULAR_IDENTITY_LINE_81')
    need(sp.simplify(g.subs({n0: sp.sqrt(2) / 2, n1: sp.sqrt(2) / 2, n2: 0, n3: 0})) == 4 * b * b, 'ANGULAR_IDENTITY_LINE_82')
    need(g.subs({a: 1, b: 1, n0: sp.Rational(3, 5), n1: sp.Rational(4, 5), n2: 0, n3: 0}) == sp.Rational(2255, 625), 'ANGULAR_IDENTITY_LINE_84')
    need(g.subs({a: 0, b: 0}) == 0, 'ANGULAR_IDENTITY_LINE_86')
    need(sp.expand(g.subs(a, 0) - b * b * projected_norm) == 0, 'ANGULAR_IDENTITY_LINE_87')
    print('PASS: Hopf metric, uniform Z law, Y and gradient identities')
    print('PASS: exact independent moments <Y²>=1/6 and <|DY|²>=4/3')
    print('PASS: angular-fraction antiderivative, endpoints, and rational tail bound')
    print('PASS: explicit mixed-sign and degenerate-case controls')
if __name__ == '__main__':
    main()
