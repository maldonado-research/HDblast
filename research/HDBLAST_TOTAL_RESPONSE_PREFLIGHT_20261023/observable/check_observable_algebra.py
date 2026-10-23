#!/usr/bin/env python3
"""Finite exact algebra checks; standard library only; no physical evaluation.

Commutative polynomials use rational coefficients, truncated only at e^3.
Noncommutative gauge words use a separate minimal free-word algebra.
No frozen producer is imported and no floating-point calculation is used.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import permutations
import json
from math import factorial
from pathlib import Path


class P:
    def __init__(self, value=0):
        if isinstance(value, P):
            self.d = dict(value.d)
        elif isinstance(value, dict):
            self.d = {k: F(v) for k, v in value.items() if v}
        else:
            self.d = {(): F(value)} if value else {}

    def __add__(self, other):
        other = P(other)
        out = dict(self.d)
        for k, v in other.d.items():
            out[k] = out.get(k, F(0)) + v
        return P(out)

    __radd__ = __add__

    def __neg__(self):
        return P({k: -v for k, v in self.d.items()})

    def __sub__(self, other):
        return self + -P(other)

    def __rsub__(self, other):
        return P(other) + -self

    def __mul__(self, other):
        other = P(other)
        out = {}
        for k, v in self.d.items():
            for l, w in other.d.items():
                key = tuple(sorted(k + l))
                if key.count('e') <= 2:
                    out[key] = out.get(key, F(0)) + v*w
        return P(out)

    __rmul__ = __mul__

    def __truediv__(self, other):
        return self * (F(1)/F(other))

    def __pow__(self, degree):
        out = P(1)
        for _ in range(degree):
            out = out*self
        return out

    def coeff(self, name, degree):
        return P({tuple(x for x in k if x != name): v
                  for k, v in self.d.items() if k.count(name) == degree})

    def substitute(self, name, value):
        out = P(0)
        for k, v in self.d.items():
            n = k.count(name)
            out += P({tuple(x for x in k if x != name): v})*P(value)**n
        return out


def var(name):
    return P({(name,): F(1)})


def matrix(n, m, fn):
    return [[P(fn(i, j)) for j in range(m)] for i in range(n)]


def eye(n):
    return matrix(n, n, lambda i, j: i == j)


def add(*ms):
    return matrix(len(ms[0]), len(ms[0][0]),
                  lambda i, j: sum((m[i][j] for m in ms), P(0)))


def scale(c, m):
    return [[c*x for x in row] for row in m]


def mm(a, b):
    return matrix(len(a), len(b[0]),
                  lambda i, j: sum((a[i][k]*b[k][j]
                                    for k in range(len(b))), P(0)))


def transpose(a):
    return [list(row) for row in zip(*a)]


def tr(a):
    return sum((a[i][i] for i in range(len(a))), P(0))


def det(a):
    out = P(0)
    for perm in permutations(range(len(a))):
        inversions = sum(perm[i] > perm[j]
                         for i in range(len(a)) for j in range(i+1, len(a)))
        term = P((-1)**inversions)
        for i, j in enumerate(perm):
            term *= a[i][j]
        out += term
    return out


def generic_symmetric(prefix, n):
    return matrix(n, n, lambda i, j: var(prefix+str(min(i, j))+str(max(i, j))))


CHECKS = []
CONTROLS = []


def zero(label, value):
    if P(value).d:
        raise ValueError('IDENTITY_FAILED:'+label)
    CHECKS.append(label)


def matrix_zero(label, value):
    for row in value:
        for item in row:
            if item.d:
                raise ValueError('MATRIX_IDENTITY_FAILED:'+label)
    CHECKS.append(label)


def nonzero(label, value):
    if not P(value).d:
        raise ValueError('MUTATION_NOT_DETECTED:'+label)
    CONTROLS.append(label)


def average_sphere(poly):
    """Exact uniform S^3 monomial moments in four ambient coordinates."""
    out = P(0)
    for key, coefficient in P(poly).d.items():
        counts = Counter(key)
        powers = [counts.get('n'+str(i), 0) for i in range(4)]
        if any(p % 2 for p in powers):
            continue
        halves = [p//2 for p in powers]
        moment = F(1)
        for h in halves:
            moment *= F(factorial(2*h), 2**h*factorial(h))
        for j in range(sum(halves)):
            moment /= 4+2*j
        keep = tuple(x for x in key if x not in ('n0', 'n1', 'n2', 'n3'))
        out += P({keep: coefficient*moment})
    return out


def run():
    e = var('e')
    # Full generic Lorentzian 4x4 inverse and eikonal, including shift.
    eta = matrix(4, 4, lambda i, j: (1 if i else -1) if i == j else 0)
    m, n = generic_symmetric('m', 4), generic_symmetric('k', 4)
    gamma = add(eta, scale(e, m), scale(e*e/2, n))
    inv = add(eta, scale(-e, mm(mm(eta, m), eta)),
              scale(e*e, add(mm(mm(mm(mm(eta, m), eta), m), eta),
                             scale(F(-1, 2), mm(mm(eta, n), eta)))))
    matrix_zero('generic_metric_inverse_through_second_order',
                add(mm(gamma, inv), scale(-1, eye(4))))
    als, bes = var('alpha_s'), var('beta_s')
    da = [var('alpha_'+str(i)) for i in range(3)]
    db = [var('beta_'+str(i)) for i in range(3)]
    clock = [[1+e*als+e*e*bes/2]] + [[e*da[i]+e*e*db[i]/2] for i in range(3)]
    residual = mm(mm(transpose(clock), inv), clock)[0][0] + 1
    zero('first_clock_equation', residual.coeff('e', 1)+2*als+m[0][0])
    beta_expected = -n[0][0]/2-m[0][0]**2/4+sum(
        ((m[0][i+1]+da[i])**2 for i in range(3)), P(0))
    second = residual.coeff('e', 2).substitute('alpha_s', -m[0][0]/2)
    zero('second_clock_equation_with_shift', second + bes-beta_expected)
    nonzero('omitted_clock_shift_detected',
            beta_expected-(-n[0][0]/2-m[0][0]**2/4+sum((x*x for x in da), P(0))))
    f = var('lapse_f')
    zero('pure_lapse_exact_clock_beta_zero', -(-2*f*f)/2-(-2*f)**2/4)

    # Independent graph pullback from generic metric/time jets and Jacobian.
    a, b, bg0, bg1, bg2 = [var(s) for s in ('a', 'b', 'bg0', 'bg1', 'bg2')]
    ai, bi = [var('a'+str(i)) for i in range(3)], [var('b'+str(i)) for i in range(3)]
    mt = generic_symmetric('mt', 4)
    delta = e*a+e*e*b/2
    graph_gamma = matrix(4, 4, lambda i, j:
        (-1 if i == j == 0 else
         (bg0+delta*bg1+delta*delta*bg2/2 if i == j else 0))
        +e*(m[i][j]+delta*mt[i][j])+e*e*n[i][j]/2)
    jac = [[e*ai[i]+e*e*bi[i]/2 for i in range(3)]] + eye(3)
    pulled = mm(mm(transpose(jac), graph_gamma), jac)
    expected1 = matrix(3, 3, lambda i, j: m[i+1][j+1]+a*bg1*(i == j))
    expected2 = matrix(3, 3, lambda i, j:
        n[i+1][j+1]+2*a*mt[i+1][j+1]+(b*bg1+a*a*bg2)*(i == j)
        +2*(m[0][i+1]*ai[j]+m[0][j+1]*ai[i])-2*ai[i]*ai[j])
    matrix_zero('graph_first_metric', matrix(3, 3,
                lambda i, j: pulled[i][j].coeff('e', 1)-expected1[i][j]))
    matrix_zero('graph_second_metric_all_terms', matrix(3, 3,
                lambda i, j: 2*pulled[i][j].coeff('e', 2)-expected2[i][j]))
    nonzero('omitted_graph_tilt_detected', -2*ai[0]*ai[0])
    # Invert T(S)=tau at the coefficient level.
    al, be = var('alpha'), var('beta')
    zero('clock_inverse_second_coefficient',
         (-be+2*al*als)/2 + (-al)*als+be/2)

    # Generic determinant of a symmetric spatial metric, no isotropic shortcut.
    m3, n3 = generic_symmetric('p', 3), generic_symmetric('r', 3)
    determinant = det(add(eye(3), scale(e, m3), scale(e*e/2, n3)))
    d1, d2 = determinant.coeff('e', 1), determinant.coeff('e', 2)
    v2_from_det = d2-d1*d1/4
    v2_trace = tr(n3)/2+tr(m3)**2/4-tr(mm(m3, m3))/2
    zero('volume_determinant_first', d1/2-tr(m3)/2)
    zero('volume_determinant_second_factorial', v2_from_det-v2_trace)
    nonzero('omitted_anisotropic_determinant_detected', tr(mm(m3, m3))/2)

    # Sphere contractions from ambient projector/Hessian, integrated by moments.
    nv = [[var('n'+str(i))] for i in range(4)]
    projector = add(eye(4), scale(-1, mm(nv, transpose(nv))))
    D = matrix(4, 4, lambda i, j: (1 if i == 0 else -1 if i == 1 else 0) if i == j else 0)
    Y = mm(mm(transpose(nv), D), nv)[0][0]
    grad = scale(2, mm(mm(projector, D), nv))
    hess = add(scale(2, mm(mm(projector, D), projector)), scale(-2*Y, projector))
    zero('sphere_Y_mean', average_sphere(Y))
    zero('sphere_Y_squared', average_sphere(Y*Y)-F(1, 6))
    zero('sphere_gradient_squared', average_sphere(mm(transpose(grad), grad)[0][0])-F(4, 3))
    zero('sphere_hessian_squared', average_sphere(tr(mm(hess, hess)))-8)
    zero('sphere_laplacian_identity_squared', average_sphere((tr(hess)+8*Y)**2))
    E, d = var('E'), var('d')
    B = scale(2, add(scale(E*Y, projector), scale(d, hess)))
    determinant_quad = average_sphere(tr(B)**2/4-tr(mm(B, B))/2)
    zero('quadrupole_determinant_mean', determinant_quad-(E*E/2-F(8, 3)*E*d-F(16, 3)*d*d))

    # Complete graph h2 contraction assembled as tensors, then sphere integrated.
    C, Z, H, Hp, Zp, dp, J, Bclock, V, Ainv = [var(s) for s in
        ('Cclock', 'Z', 'H', 'Hp', 'Zp', 'dp', 'J0', 'Bclock', 'V', 'Ainv')]
    grad2 = mm(grad, transpose(grad))
    h2 = add(scale(J, projector),
             scale(-4*C*Y, add(scale((Zp+2*H*Z)*Y, projector), scale(dp+2*H*d, hess))),
             scale(2*H*Bclock, projector),
             scale(C*C*Y*Y*(2*Hp+4*H*H), projector),
             scale(-4*C*V*Ainv-2*C*C*Ainv*Ainv, grad2))
    contracted = average_sphere(tr(h2))/2+determinant_quad.substitute('E', Z-H*C)
    expected = (F(3, 2)*J+3*H*Bclock
                -C*(3*Zp+6*H*Z-8*dp-16*H*d)/3
                +C*C*(Hp+2*H*H)/2-F(8, 3)*C*V*Ainv-F(4, 3)*C*C*Ainv*Ainv
                +(Z-H*C)**2/2-F(8, 3)*(Z-H*C)*d-F(16, 3)*d*d)
    zero('full_second_volume_mean_formula', contracted-expected)
    nonzero('wrong_J0_trace_normalization_detected', J)
    X, Xp, N0 = [var(s) for s in ('X', 'Xp', 'N0')]
    beta_mean_prime = -N0/2-X*X/6+F(4, 3)*(V+C*Ainv)**2
    Bprime = -beta_mean_prime+(X*X+C*Xp)/3
    zero('clock_mean_second_formula', Bprime-(N0/2+X*X/2+C*Xp/3-F(4, 3)*(V+C*Ainv)**2))

    # Independent mixed pure gauge from the separate literature review:
    # xi=f(s)Y partial_s+c(s)grad Y, f(0)=c(0)=0, n=L_xi^2 gamma0.
    # The reviewer derived N0,J0 directly from Lie derivatives, not from (8).
    # Expand the combinations A^-1 V=cp-f A^-2 and V+C/A=A cp before
    # substitution; this keeps this standard-library check polynomial.
    f, fp, fpp, c, cp, A2 = [var(s) for s in ('gf', 'gfp', 'gfpp', 'gc', 'gcp', 'gA2')]
    pure_N = -(f*fpp+2*fp*fp)/3+F(8, 3)*(-c*fp+A2*cp*cp-f*cp)
    pure_J = (f*f*(Hp+2*H*H)+H*f*fp-F(8, 3)*H*f*c
              +F(32, 3)*c*c-F(8, 3)*f*f*Ainv*Ainv)/3
    pure_beta = f*fp/6+F(4, 3)*f*c
    pure_beta_prime = (fp*fp+f*fpp)/6+F(4, 3)*(fp*c+f*cp)
    zero('independent_mixed_pure_gauge_clock',
         -pure_N/2-fp*fp/6+F(4, 3)*A2*cp*cp-pure_beta_prime)
    pure_terms = {
        'total_spatial_second': F(3, 2)*pure_J,
        'second_clock': 3*H*(-pure_beta+f*fp/3),
        'first_evaluation': -f*(3*(Hp*f+H*fp)+6*H*H*f-8*cp-16*H*c)/3,
        'background_evaluation': f*f*(Hp+2*H*H)/2,
        'shift_graph': -F(8, 3)*f*(cp-f*Ainv*Ainv),
        'clock_gradient': -F(4, 3)*f*f*Ainv*Ainv,
        'quadratic_determinant': -F(16, 3)*c*c,
    }
    pure_total = sum(pure_terms.values(), P(0))
    zero('independent_mixed_pure_gauge_volume_zero', pure_total)
    for name, term in pure_terms.items():
        nonzero('mixed_gauge_deletion_'+name, pure_total-term)

    # Radius and log volume expansions independently reconstructed by cubing
    # and exponentiating their proposed coefficients.
    v1, v2 = var('v1'), var('v2')
    radius = 1+e*v1/3+e*e*(v2/6-v1*v1/9)
    zero('radius_cubed_matches_volume', radius**3-(1+e*v1+e*e*v2/2))
    logvol = e*v1+e*e*(v2-v1*v1)/2
    zero('log_volume_second_coefficient', 1+logvol+logvol*logvol/2-(1+e*v1+e*e*v2/2))

    # Free noncommuting words check the full second-order gauge cancellation.
    def word(symbol): return {(symbol,): F(1)}
    def wa(*terms):
        ans = {}
        for term in terms:
            for key, value in term.items(): ans[key] = ans.get(key, F(0))+value
        return {key: value for key, value in ans.items() if value}
    def ws(coefficient, term): return {key: coefficient*value for key, value in term.items() if coefficient*value}
    def wm(left, right):
        ans = {}
        for k, v in left.items():
            for l, w in right.items(): ans[k+l] = ans.get(k+l, F(0))+v*w
        return {key: value for key, value in ans.items() if value}
    xi, xi2, z, z2, q, k, g = [word(s) for s in ('xi', 'xi2', 'z', 'z2', 'q', 'k', 'g')]
    zp = wa(z, ws(-1, xi))
    z2p = wa(z2, ws(-1, xi2), wm(xi, z), ws(-1, wm(z, xi)))
    qp = wa(q, wm(xi, g))
    kp = wa(k, ws(2, wm(xi, q)), wm(wa(xi2, wm(xi, xi)), g))
    h2p = wa(kp, ws(2, wm(zp, qp)), wm(wa(z2p, wm(zp, zp)), g))
    h2old = wa(k, ws(2, wm(z, q)), wm(wa(z2, wm(z, z)), g))
    if wa(h2p, ws(-1, h2old)):
        raise ValueError('SECOND_GAUGE_CANCELLATION_FAILED')
    CHECKS.append('noncommutative_second_gauge_cancellation')
    wrong_commutator = wm(wa(ws(-2, wm(xi, z)), ws(2, wm(z, xi))), g)
    if not wrong_commutator:
        raise ValueError('GAUGE_COMMUTATOR_MUTATION_NOT_DETECTED')
    CONTROLS.append('wrong_gauge_commutator_sign_detected')

    return {'status': 'EXACT_FINITE_ALGEBRA_CHECKS_PASS',
            'coefficient_convention': 'parameter derivatives, e^2/2',
            'arithmetic': 'Python standard-library Fraction sparse polynomials; no floating point',
            'checks': CHECKS, 'mutation_controls': CONTROLS,
            'not_established': ['actual observable value or sign', 'no-caustic numerical time',
                                'finite-amplitude nonlinear physical family',
                                'numerical input authentication or numerical error propagation']}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    report = run()
    with Path(args.output).open('x', encoding='utf-8') as out:
        json.dump(report, out, sort_keys=True, indent=2)
        out.write('\n')
    print('PASS', len(report['checks']), 'identities;', len(report['mutation_controls']), 'controls')


if __name__ == '__main__':
    main()
