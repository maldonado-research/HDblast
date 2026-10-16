#!/usr/bin/env python3
"""Independent exact algebra for the numerical-response initial-data review.

Standard-library only. Sparse Laurent polynomials retain every independent
point-jet symbol, so the comparisons below are polynomial identities, not
finite point samples. No producer code or scientific report is imported.
"""
import argparse
from fractions import Fraction as F
from math import comb
import json
from pathlib import Path
import subprocess
import sys


NAMES = ('e', 'rho', 'U', 'Uy', 'Ui', 'T', 'f', 'fy', 'fi', 'q', 'a',
         'p', 'py', 'Vp', 'fyy', 'zeta', 'zetay', 'ell', 'C', 'B', 'vy',
         'v', 'h', 'Vpp', 'vt', 'E')
ZERO = (0,) * len(NAMES)


class Failure(Exception):
    pass


def need(condition, message):
    if not condition:
        raise Failure(message)


class Poly:
    def __init__(self, value=0):
        if isinstance(value, dict):
            self.d = {k: F(v) for k, v in value.items() if v}
        else:
            self.d = {ZERO: F(value)} if value else {}

    def __add__(self, other):
        other = as_poly(other)
        d = self.d.copy()
        for key, value in other.d.items():
            d[key] = d.get(key, F(0)) + value
        return Poly(d)

    __radd__ = __add__

    def __neg__(self):
        return Poly({k: -v for k, v in self.d.items()})

    def __sub__(self, other):
        return self + -as_poly(other)

    def __rsub__(self, other):
        return as_poly(other) + -self

    def __mul__(self, other):
        other = as_poly(other)
        d = {}
        for a, av in self.d.items():
            for b, bv in other.d.items():
                key = tuple(x+y for x, y in zip(a, b))
                if key[0] > 2:
                    continue
                d[key] = d.get(key, F(0)) + av*bv
        return Poly(d)

    __rmul__ = __mul__

    def __truediv__(self, other):
        return self * (F(1) / F(other))

    def __pow__(self, order):
        need(isinstance(order, int) and order >= 0, 'nonnegative polynomial power')
        result = Poly(1)
        for _ in range(order):
            result = result*self
        return result

    def coefficient(self, order):
        return Poly({(0,)+k[1:]: v for k, v in self.d.items() if k[0] == order})

    def derivative(self, name):
        index = NAMES.index(name)
        d = {}
        for key, value in self.d.items():
            if key[index]:
                new = list(key)
                new[index] -= 1
                d[tuple(new)] = value*key[index]
        return Poly(d)

    def __eq__(self, other):
        return self.d == as_poly(other).d


def as_poly(value):
    return value if isinstance(value, Poly) else Poly(value)


def symbol(name, power=1):
    exponent = list(ZERO)
    exponent[NAMES.index(name)] = power
    return Poly({tuple(exponent): 1})


def inverse_unit(value):
    delta = value-1
    need(delta.coefficient(0) == 0, 'unit constant term')
    return 1-delta+delta*delta


def gauge_algebra(mutation):
    e, rho, U, Uy, Ui, T, f, fy, fi, q, a = map(symbol, NAMES[:11])
    rhoi = symbol('rho', -1)
    A, B, P = 1-4*e*U, 1+2*e*U+e*e*f/2, 1+e*U
    Ai, Bi, Pi = map(inverse_unit, (A, B, P))
    Ay, By, Py = -4*e*Uy, 2*e*Uy+e*e*fy/2, e*Uy
    Ad, Bd, Pd = -4*e*Ui, 2*e*Ui+e*e*fi/2, e*Ui
    # Exact initial ADM prescription in the actual four-dimensional slice.
    Ndot = rho*P*P*((-2*e*T+e*e*rho*q/2)*Ai
                        +3*(e*T+e*e*rho*q/2)*Bi)
    bydot = rho*rho*P*P*(Ay*Ai*Ai/2-3*By*Ai*Bi/2
                         -4*a*Ai-Py*Ai*Pi+a*Pi*Pi+3*a*Bi)
    bidot = -P*P*Bi*(Ad*Ai/2+Bd*Bi/2+Pd*Pi)
    wanted_N2 = 4*rho*rho*q-24*rho*U*T
    wanted_by2 = rho*rho*(-74*U*Uy-3*fy/2-194*a*U*U-3*a*f)
    wanted_bi2 = 22*U*Ui-fi/2
    if mutation == 'half_lapse_quadratic':
        wanted_N2 = wanted_N2/2
    if mutation == 'omit_radial_fy':
        wanted_by2 = wanted_by2+3*rho*rho*fy/2
    if mutation == 'half_angular_quadratic':
        wanted_bi2 = wanted_bi2/2
    checks = [
        (Ndot.coefficient(1), rho*T, 'first lapse velocity'),
        (bydot.coefficient(1), -6*rho*rho*(Uy+4*a*U), 'first gauge acceleration'),
        (bidot.coefficient(1), 0, 'first angular velocity'),
        (2*Ndot.coefficient(2), wanted_N2, 'second lapse velocity'),
        (2*bydot.coefficient(2), wanted_by2, 'second radial shift velocity'),
        (2*bidot.coefficient(2), wanted_bi2, 'second angular shift velocity'),
    ]
    N = rho*P
    gtt2 = 2*(-N*N).coefficient(2)
    gttdot2 = 2*(-2*N*Ndot).coefficient(2)
    gtydot2 = 2*(A*bydot).coefficient(2)
    gtidot2 = 2*(rho*rho*B*bidot).coefficient(2)
    Kyy = -2*e*T*rhoi+e*e*q/2
    Kangular = e*T*rhoi+e*e*q/2
    if mutation == 'omit_mixed_metric_conversion':
        gtydot2 = 2*bydot.coefficient(2)
    expected_gtt2 = 0 if mutation == 'omit_lapse_square' else -2*rho*rho*U*U
    checks += [
        (gtt2, expected_gtt2, 'second metric lapse square'),
        (gttdot2, -8*rho**3*q+44*rho*rho*U*T, 'second metric tt velocity'),
        (gtydot2, rho*rho*(-26*U*Uy-3*fy/2-2*a*U*U-3*a*f), 'mixed metric conversion'),
        (gtidot2, rho*rho*(22*U*Ui-fi/2), 'angular metric conversion'),
        (2*(2*N*Kyy).coefficient(2), 2*rho*q-8*U*T, 'physical yy velocity'),
        (2*(2*N*Kangular).coefficient(2), 2*rho*q+4*U*T, 'physical angular velocity'),
    ]
    for actual, expected, label in checks:
        need(actual == expected, label)
    return len(checks)


def scalar_cancellation(mutation):
    p, py, Vp, fy, fyy, a = map(symbol, ('p', 'py', 'Vp', 'fy', 'fyy', 'a'))
    pi = symbol('p', -1)
    # Differentiate zeta=-3 f_y/(2p), then add both scalar linear terms.
    zeta = -3*fy*pi/2
    zetay = -3*(fyy*pi-fy*pi*pi*py)/2
    if mutation == 'flip_scalar_compensator':
        zeta, zetay = -zeta, -zetay
    scalar = -2*p*zetay-2*Vp*zeta
    # Enforce the actual background identity p_y=V'-4ap by polynomial substitution.
    substituted = Poly()
    for key, value in scalar.d.items():
        power = key[NAMES.index('py')]
        base_key = list(key)
        base_key[NAMES.index('py')] = 0
        substituted += Poly({tuple(base_key): value})*(Vp-4*a*p)**power
    need(substituted == 3*fyy+12*a*fy, 'scalar compensator cancellation')
    # No admissible sector has the excluded lambda=3 denominator.
    for lam in (0, 8, 24):
        need(2*(lam-3) in (-6, 10, 42), 'three physical sector denominators')
    return 4


def bump_fourth(mutation):
    # Low-to-high coefficients of Q''' / (v exp(4-s)).
    p3 = [0, 0, 0, -12, 24, -10, 1]
    diff = [i*p3[i] for i in range(1, len(p3))]+[0]
    p4 = [0]*9
    for i, value in enumerate(p3):
        p4[i] -= 2*value
        p4[i+1] -= 4*(value-diff[i])
        p4[i+2] += value-diff[i]
    expected = [0, 0, 0, -120, 420, -384, 136, -20, 1]
    if mutation == 'omit_fourth_packet_term':
        p4[-1] = 0
    need(p4 == expected, 'fourth packet derivative identity')
    shifted = []
    for sign in (-1, 1):
        polynomial = [sign*x for x in p4]
        polynomial[8] += 1
        shifted.append([sum(polynomial[i]*comb(i,j)*4**(i-j)
                            for i in range(j,9)) for j in range(9)])
    need(all(x >= 0 for row in shifted for x in row), 'global fourth derivative polynomial envelope')
    # exp(2)>sum_{j=0}^4 2^j/j!=7, hence exp(4)>49.
    bound = F(8**8, 49)
    need(bound < 342393, 'strict global fourth derivative cap')
    return {'polynomial_low_to_high': p4, 'shifted_nonnegative_coefficients': shifted,
            'strict_scaled_fourth_bound': str(bound),
            'domain': 'all real r, including both smooth flat endpoints'}


def canceled_physical_sources(mutation):
    U, v, vy, p, a, h, Vpp, T, vt, E = map(symbol,
        ('U', 'v', 'vy', 'p', 'a', 'h', 'Vpp', 'T', 'vt', 'E'))
    ir, ip = symbol('rho', -1), symbol('p', -1)
    ay = -ir*ir-p*p/3
    uy = -2*a*U-p*v/3
    uyy = (-2*ay+4*a*a)*U+(2*a*p-p*h/3)*v-p*vy/3
    original_A = (-24*uy*uy-288*a*U*uy-24*U*uyy
        -(192*ay+384*a*a+32*p*p)*U*U
        -ir*ir*(336*U*U+12*T*T+2*vt*vt)
        -2*vy*vy-16*p*U*vy-2*Vpp*v*v)
    canceled_A = (-2*vy*vy-8*p*U*vy+8*p*(h+2*a)*U*v
        +(16*p*p-192*ir*ir)*U*U-(2*Vpp+F(8,3)*p*p)*v*v
        -ir*ir*(12*T*T+2*vt*vt))
    need(original_A == canceled_A, 'all-point canceled Hamiltonian source')
    d = 3*ir*ip
    dy = d*(3*a-h)
    ut = E*U
    uty = E*(uy+ir*U)
    vt = E*(v-d*U)
    vty = E*(vy-dy*U-d*uy+ir*(v-d*U))
    Fy = ir*(18*uy*ut+18*U*uty-vy*vt-v*vty
             -a*(18*U*ut-v*vt))
    G = ir*(12*U*uty-48*a*U*ut-6*uy*ut-2*vt*vy)
    original_Pi = (G-Fy)*ip
    if mutation == 'omit_physical_affine_term':
        original_Pi = original_Pi+6*uy*ut*ir*ip
    canceled_Pi = E*ir*ip*(d*U*vy+(d*(h-ir)+10*p)*U*v
                     +(30*a-6*ir)*U*U+(2*ir-a)*v*v)
    need(original_Pi == canceled_Pi, 'all-point canceled physical-affine scalar velocity')
    return 2


CONTROLS = {
    'half_lapse_quadratic': 'second lapse velocity',
    'omit_radial_fy': 'second radial shift velocity',
    'half_angular_quadratic': 'second angular shift velocity',
    'omit_mixed_metric_conversion': 'mixed metric conversion',
    'omit_lapse_square': 'second metric lapse square',
    'flip_scalar_compensator': 'scalar compensator cancellation',
    'omit_fourth_packet_term': 'fourth packet derivative identity',
    'omit_physical_affine_term': 'all-point canceled physical-affine scalar velocity',
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    parser.add_argument('--inject', choices=tuple(CONTROLS))
    parser.add_argument('--run-controls', action='store_true')
    args = parser.parse_args()
    try:
        count = (gauge_algebra(args.inject)+scalar_cancellation(args.inject)
                 +canceled_physical_sources(args.inject))
        bump = bump_fourth(args.inject)
        controls = {}
        if args.run_controls:
            for name, reason in CONTROLS.items():
                command = [sys.executable]
                if not __debug__:
                    command.append('-O')
                command += [str(Path(__file__).resolve()), '--inject', name]
                result = subprocess.run(command, capture_output=True, text=True)
                need(result.returncode == 1 and result.stderr.strip() == reason,
                     'mutation-specific independent control: '+name)
                controls[name] = reason
        report = {'status': 'pass', 'independent_symbolic_identities': count,
                  'fourth_packet_derivative': bump, 'controls': controls,
                  'scope': 'initial formula identities and global packet derivative cap; no evolved response or source approximation error'}
        if args.output:
            need(not args.output.exists(), 'output already exists')
            args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
        else:
            print(json.dumps(report, sort_keys=True))
    except Failure as error:
        print(error, file=sys.stderr)
        raise SystemExit(1)


if __name__ == '__main__':
    main()
