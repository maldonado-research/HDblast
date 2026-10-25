#!/usr/bin/env python3
"""Independent polynomial controls for B2 shell primitives and triangle algebra.

No profile producer or quadrature experiment is run. Polynomial integration
uses exact Lagrange weights on each reflection branch of the frozen kernel.
"""
import argparse
from functools import lru_cache
import hashlib
import importlib.util
import json
import math
import sys
from fractions import Fraction as Q
from pathlib import Path

sys.dont_write_bytecode = True
PARITY = (1, -1, 1, 1, 1)
ZERO = (Q(0),)*5


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def va(a, b):
    return tuple(x+y for x, y in zip(a, b))


def scale(c, a):
    return tuple(c*x for x in a)


def project(a):
    return tuple(x if sign == 1 else Q(0) for x, sign in zip(a, PARITY))


def mv(matrix, vector):
    return tuple(sum((a*b for a, b in zip(row, vector)), Q(0)) for row in matrix)


def primitive(coefficients, x, order):
    return tuple(sum((c*Q(math.factorial(k), math.factorial(k+order))*x**(k+order)
                      for k, c in enumerate(row)), Q(0)) for row in coefficients)


@lru_cache(None)
def integration_weights(degree):
    """Exact integral weights of Lagrange cardinal polynomials on [0,1]."""
    nodes = tuple(Q(j, degree) for j in range(degree+1)) if degree else (Q(0),)
    weights = []
    for j, node in enumerate(nodes):
        coefficients = [Q(1)]
        denominator = Q(1)
        for k, other in enumerate(nodes):
            if k == j:
                continue
            nxt = [Q(0)]*(len(coefficients)+1)
            for power, value in enumerate(coefficients):
                nxt[power] -= other*value
                nxt[power+1] += value
            coefficients = nxt
            denominator *= node-other
        weights.append(sum((value/Q(power+1) for power, value in enumerate(coefficients)), Q(0))/denominator)
    need(sum(weights) == 1, 'CONSTANT_INTEGRATION_CONTROL')
    for power in range(degree+1):
        need(sum((w*x**power for w, x in zip(weights, nodes)), Q(0)) == Q(1, power+1),
             'MONOMIAL_INTEGRATION_CONTROL')
    return nodes, tuple(weights)


def integrate_polynomial(function, left, right, degree):
    if left == right:
        return ZERO
    nodes, weights = integration_weights(degree)
    result = ZERO
    for node, weight in zip(nodes, weights):
        result = va(result, scale((right-left)*weight, function(left+(right-left)*node)))
    return result


def closed(F, G, g, T):
    f1, f2, f3, f4 = (F(T, k) for k in range(1, 5))
    g2, g3, g4 = (G(T, k) for k in range(2, 5))
    gm3, gm4 = G(-T, 3), G(-T, 4)
    boundary2, boundary3 = g(T, 2), g(T, 3)
    H = tuple(T*f2[j]-2*f3[j]+(g3[j]-gm3[j])/2 for j in range(5))
    B = tuple((H[j] if PARITY[j] == 1 else Q(0))-boundary2[j] for j in range(5))
    C = tuple((T*T*f1[j]-T*f2[j]+T*g2[j]-g3[j]-gm3[j]+PARITY[j]*H[j])/4-boundary2[j]/2
              for j in range(5))
    V = tuple((T*T*f2[j]-3*T*f3[j]+3*f4[j]+T*g3[j]-2*g4[j]+gm4[j]
               +PARITY[j]*(T*f3[j]-3*f4[j]+(g4[j]+gm4[j])/2))/4-boundary3[j]/2
              for j in range(5))
    return {'B': B, 'C': C, 'V': V}


def polynomial_controls(kernel):
    z = [[Q(0)] for _ in range(5)]
    cases = []
    cases.append(('F_constant', [[Q((-1)**i*(i+1), 7)] for i in range(5)], z, z))
    cases.append(('G_constant', z, [[Q((-1)**(i+1)*(i+2), 11)] for i in range(5)], z))
    cases.append(('g_constant', z, z, [[Q(i+1, 13) if i != 1 else Q(0)] for i in range(5)]))
    for degree in range(5):
        f = [[Q((-1)**(i+k)*(i+2)*(k+1), 7+k) for k in range(degree+1)] for i in range(5)]
        gg = [[Q((-1)**k*(i+1)*(k+3), 13+k) for k in range(5-degree)] for i in range(5)]
        g = [[Q((-1)**(i+k)*(i+3)*(k+2), 17+k) if i != 1 else Q(0)
              for k in range((degree % 3)+1)] for i in range(5)]
        cases.append(('mixed_degree_'+str(degree), f, gg, g))
    results = []
    negative_controls = {'omit_boundary_lift': False, 'flip_negative_G_argument': False}
    for name, fc, gc, bc in cases:
        F = lambda y, k=0: primitive(fc, y, k)
        G = lambda y, k=0: primitive(gc, y, k)
        g = lambda y, k=0: primitive(bc, y, k)
        degree = max(max(map(len, fc))+1, max(map(len, gc))+1, max(map(len, bc)))
        def Z(x, r):
            result = kernel.reflected_kernel(F, G, x, r)['value']
            return va(result, scale(-1, g(x-r, 1))) if x > r else result
        def direct_C(a):
            return va(integrate_polynomial(lambda s: Z(s, a-s), Q(0), a/2, degree),
                      integrate_polynomial(lambda s: Z(s, a-s), a/2, a, degree))
        for T in (Q(3, 7), Q(9, 8), Q(5, 4)):
            direct = {'B': integrate_polynomial(lambda s: Z(s, Q(0)), Q(0), T, degree),
                      'C': direct_C(T),
                      'V': integrate_polynomial(direct_C, Q(0), T, degree+1)}
            formula = closed(F, G, g, T)
            need(direct == formula, 'CLOSED_PRIMITIVE_MISMATCH:'+name+':'+str(T))
            if any(g(T, 2)):
                broken = closed(F, G, lambda y, k=0: ZERO, T)
                negative_controls['omit_boundary_lift'] |= broken != direct
            if any(G(-T, 3)):
                wrong_G = lambda y, k=0: G(abs(y), k)
                negative_controls['flip_negative_G_argument'] |= closed(F, wrong_G, g, T) != direct
            results.append({'case': name, 'T': str(T), 'branch_polynomial_degree_upper': degree,
                            'B_C_V_exactly_equal': True})
    need(all(negative_controls.values()), 'NEGATIVE_CONTROL_NOT_DETECTED')
    return results, negative_controls


def triangle_controls():
    """Integrate monomial Z and its partial derivatives over the triangle."""
    A = [[Q((-1)**(i+j)*(2*i+j+1), 17) for j in range(5)] for i in range(5)]
    Bmat = [[Q((i+1)*(j+2), 23) for j in range(5)] for i in range(5)]
    D = [[Q((-1)**j*(i+j+3), 29) for j in range(5)] for i in range(5)]
    M = [[Q((i+2)*(j+1), 31) if i != 1 and j != 1 else Q(0) for j in range(5)] for i in range(5)]
    rows = []
    for index in range(6):
        exponents = [(1+(index+i) % 4, (2*index+i) % 4) for i in range(5)]
        coefficients = [Q((-1)**(index+i)*(i+1), index+3) for i in range(5)]
        for T in (Q(2, 3), Q(9, 8), Q(5, 4)):
            def moment(a, b):
                return Q(math.factorial(a)*math.factorial(b), math.factorial(a+b+2))*T**(a+b+2)
            X, R, volume, characteristic, shell = [], [], [], [], []
            for coefficient, (a, b) in zip(coefficients, exponents):
                X.append(coefficient*a*moment(a-1, b))
                R.append(coefficient*b*moment(a, b-1) if b else Q(0))
                volume.append(coefficient*moment(a, b))
                characteristic.append(coefficient*Q(math.factorial(a)*math.factorial(b), math.factorial(a+b+1))*T**(a+b+1))
                shell.append(coefficient*T**(a+1)/Q(a+1) if b == 0 else Q(0))
            need(tuple(X) == tuple(characteristic), 'TRIANGLE_DX')
            need(tuple(R) == tuple(c-b for c, b in zip(characteristic, shell)), 'TRIANGLE_DR')
            direct = va(project(va(va(mv(A, X), mv(Bmat, R)), mv(D, volume))),
                        scale(-1, mv(M, project(shell))))
            formula = va(project(va(va(mv(A, characteristic), mv(Bmat, characteristic)),
                                    va(scale(-1, mv(Bmat, shell)), mv(D, volume)))),
                         scale(-1, mv(M, project(shell))))
            need(direct == formula, 'B2_TRIANGLE_IDENTITY')
            rows.append({'case': index, 'T': str(T), 'exact_triangle_B2_identity': True})
    return rows


def producer_controls(producer, kernel, matrices):
    """Exercise producer algebra against direct kernel polynomial integration."""
    A, Bmat, D, M = (matrices[k] for k in ('A', 'B', 'D', 'M'))
    ApB = [[a+b for a, b in zip(ar, br)] for ar, br in zip(A, Bmat)]
    BmA = [[b-a for a, b in zip(ar, br)] for ar, br in zip(A, Bmat)]
    rows = []
    for index in range(3):
        uc, vc = [], []
        for j in range(5):
            du = (3 if j == 1 else 2)+2*(index % 2)
            dv = (1 if j == 1 else 0)+2*((index+j) % 2)
            uc.append([Q(0)]*du+[Q((-1)**(index+j)*(j+2), 7+index)])
            vc.append([Q(0)]*dv+[Q((-1)**j*(j+1), 13+index)])
        duc = [[k*c for k, c in enumerate(row) if k] for row in uc]
        def reflected(coefficients, y, k, parity):
            positive = primitive(coefficients, abs(y), k)
            return tuple((-1)**k*sign*v for sign, v in zip(parity, positive)) if y < 0 else positive
        def ridge_F(y, k=0):
            du = reflected(duc, y, k, tuple(-s for s in PARITY))
            u = reflected(uc, y, k, PARITY)
            v = reflected(vc, y, k, PARITY)
            pv = reflected(vc, y, k+1, PARITY)
            return scale(Q(1, 2), va(mv(ApB, va(du, v)), mv(D, va(u, pv))))
        def ridge_G(y, k=0):
            du = reflected(duc, y, k, tuple(-s for s in PARITY))
            u = reflected(uc, y, k, PARITY)
            v = reflected(vc, y, k, PARITY)
            pv = reflected(vc, y, k+1, PARITY)
            return scale(Q(1, 2), va(mv(BmA, va(du, scale(-1, v))), mv(D, va(u, scale(-1, pv)))))
        def g(y, k=0):
            return mv(M, project(va(primitive(uc, y, k), primitive(vc, y, k+1))))
        def Z(x, r):
            value = kernel.reflected_kernel(ridge_F, ridge_G, x, r)['value']
            return va(value, scale(-1, g(x-r, 1))) if x > r else value
        degree = max(max(map(len, uc))-1, max(map(len, vc)))+2
        def direct_C(a):
            return va(integrate_polynomial(lambda s: Z(s, a-s), Q(0), a/2, degree),
                      integrate_polynomial(lambda s: Z(s, a-s), a/2, a, degree))
        for T in (Q(3, 7), Q(9, 8), Q(5, 4)):
            U = {k: primitive(uc, T, k) for k in range(5)}
            V = {k: primitive(vc, T, k) for k in range(1, 6)}
            actual = producer.shell_algebra(T, U, V, matrices)
            bh = integrate_polynomial(lambda s: Z(s, Q(0)), Q(0), T, degree)
            ch = direct_C(T)
            vh = integrate_polynomial(direct_C, Q(0), T, degree+1)
            need(tuple(actual['B1star']) == Z(T, Q(0)), 'PRODUCER_B1_KERNEL')
            need(tuple(actual['B1star_shell_history']) == bh, 'PRODUCER_B_HISTORY')
            need(tuple(actual['B1star_characteristic_history']) == ch, 'PRODUCER_C_HISTORY')
            need(tuple(actual['B1star_triangle_history']) == vh, 'PRODUCER_V_HISTORY')
            direct_b2 = va(project(va(va(mv(A, ch), mv(Bmat, va(ch, scale(-1, bh)))), mv(D, vh))),
                           scale(-1, mv(M, project(bh))))
            need(tuple(actual['B2star']) == direct_b2, 'PRODUCER_B2')
            rows.append({'case': index, 'T': str(T), 'producer_B1_B2_and_histories_exactly_equal': True})
    nodes = [Q(0), Q(1, 3), Q(5, 6), Q(3, 2)]
    values = [[Q((-1)**(i+j)*(i+1)*(j+2), 11+j) for j in range(5)] for i in range(len(nodes))]
    data = {'nodes': nodes, 'profiles': {'v': values}}
    velocity_rows = []
    for T in (Q(1, 3), Q(2, 3), Q(9, 8), Q(5, 4)):
        observed = producer.velocity_moments(T, data)
        expected = {k: ZERO for k in range(1, 6)}
        for left, right, vaa, vbb in zip(nodes, nodes[1:], values, values[1:]):
            if T <= left:
                break
            end = min(right, T)
            def v(r):
                return tuple(a+(b-a)*(r-left)/(right-left) for a, b in zip(vaa, vbb))
            for k in range(1, 6):
                expected[k] = va(expected[k], integrate_polynomial(
                    lambda r: scale((T-r)**(k-1)/math.factorial(k-1), v(r)), left, end, k))
        need({k: tuple(row) for k, row in observed.items()} == expected, 'PRODUCER_P1_VELOCITY_MOMENTS')
        velocity_rows.append({'T': str(T), 'moment_orders': [1, 2, 3, 4, 5], 'exactly_equal': True})
    return rows, velocity_rows


class ExactInterval:
    """Unrounded rational endpoints, independent of the producer interval code."""
    def __init__(self, lo=0, hi=None):
        self.lo, self.hi = Q(lo), Q(lo if hi is None else hi)
        need(self.lo <= self.hi, 'REPLAY_INTERVAL_ORDER')
    @staticmethod
    def cv(value):
        return value if isinstance(value, ExactInterval) else ExactInterval(value)
    def __add__(self, other):
        other = self.cv(other)
        return ExactInterval(self.lo+other.lo, self.hi+other.hi)
    __radd__ = __add__
    def __neg__(self):
        return ExactInterval(-self.hi, -self.lo)
    def __sub__(self, other):
        return self+-self.cv(other)
    def __rsub__(self, other):
        return self.cv(other)+-self
    def __mul__(self, other):
        other = self.cv(other)
        products = [a*b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return ExactInterval(min(products), max(products))
    __rmul__ = __mul__
    def __truediv__(self, other):
        other = Q(other)
        need(other != 0, 'REPLAY_DIVIDE_BY_ZERO')
        return self*(1/other)
    def mag(self):
        return max(abs(self.lo), abs(self.hi))


def box(value):
    return ExactInterval(*value) if isinstance(value, list) else ExactInterval(value)


def contained(outer, inner, message):
    outer, inner = box(outer), ExactInterval.cv(inner)
    need(outer.lo <= inner.lo <= inner.hi <= outer.hi, message)


def p1_moments_binomial(T, profile):
    """Independent local-power expansion, rather than producer endpoint powers."""
    nodes = list(map(Q, profile['nodes']))
    values = [list(map(Q, row)) for row in profile['profiles']['v']]
    result = {k: [Q(0)]*5 for k in range(1, 6)}
    for a, b, vaa, vbb in zip(nodes, nodes[1:], values, values[1:]):
        if T <= a:
            break
        width = min(T, b)-a
        for k in range(1, 6):
            for j in range(5):
                slope = (vbb[j]-vaa[j])/(b-a)
                for power in range(k):
                    coefficient = Q(math.comb(k-1, power)*(-1)**power, math.factorial(k-1))*(T-a)**(k-1-power)
                    result[k][j] += coefficient*(vaa[j]*width**(power+1)/Q(power+1)
                                                 +slope*width**(power+2)/Q(power+2))
    return result


def moment_algebra_independent(T, U, V, matrices):
    A, Bmat, D, M = (matrices[k] for k in ('A', 'B', 'D', 'M'))
    ApB = [[a+b for a, b in zip(ar, br)] for ar, br in zip(A, Bmat)]
    BmA = [[b-a for a, b in zip(ar, br)] for ar, br in zip(A, Bmat)]
    def reflected_moment(values, k, parity, y):
        return tuple((-1)**k*sign*x for sign, x in zip(parity, values)) if y < 0 else tuple(values)
    def Ffun(y, k):
        du = reflected_moment(U[k-1], k, tuple(-s for s in PARITY), y)
        u = reflected_moment(U[k], k, PARITY, y)
        v = reflected_moment(V[k], k, PARITY, y)
        pv = reflected_moment(V[k+1], k+1, PARITY, y)
        return scale(Q(1, 2), va(mv(ApB, va(du, v)), mv(D, va(u, pv))))
    def Gfun(y, k):
        du = reflected_moment(U[k-1], k, tuple(-s for s in PARITY), y)
        u = reflected_moment(U[k], k, PARITY, y)
        v = reflected_moment(V[k], k, PARITY, y)
        pv = reflected_moment(V[k+1], k+1, PARITY, y)
        return scale(Q(1, 2), va(mv(BmA, va(du, scale(-1, v))), mv(D, va(u, scale(-1, pv)))))
    def gfun(y, k):
        need(y == T, 'BOUNDARY_POSITIVE_PRIMITIVE')
        return mv(M, project(va(U[k], V[k+1])))
    histories = closed(Ffun, Gfun, gfun, T)
    bh, ch, vh = (histories[name] for name in ('B', 'C', 'V'))
    b0 = project(va(U[0], V[1]))
    bulk = project(va(va(scale(T, Ffun(T, 1)), scale(-1, Ffun(T, 2))),
                      scale(Q(1, 2), va(Gfun(T, 2), Gfun(-T, 2)))))
    boundary = scale(-1, gfun(T, 1))
    b1 = va(bulk, boundary)
    b2 = va(project(va(va(mv(ApB, ch), scale(-1, mv(Bmat, bh))), mv(D, vh))),
            scale(-1, mv(M, project(bh))))
    return {'B0': b0, 'B1star_bulk': bulk, 'B1star_boundary': boundary,
            'B1star': b1, 'B2star': b2, 'two_term': va(b0, b1),
            'three_term': va(va(b0, b1), b2),
            'B1star_shell_history': bh, 'B1star_characteristic_history': ch,
            'B1star_triangle_history': vh}


def replay_report(report_path, producer_path, repo, roots, matrices):
    data = json.loads(report_path.read_text())
    need(data['source_sha256'] == sha(producer_path), 'REPORT_PRODUCER_SOURCE_PIN')
    need(data['schema'] == 'HDBLAST_EXACT_U_B2_SHELL_CANDIDATE_V1', 'REPORT_SCHEMA')
    for rel, pin in data['frozen_source_pins'].items():
        need(sha(repo/rel) == pin, 'REPORT_FROZEN_SOURCE_PIN:'+rel)
    need({k: [list(map(Q, row)) for row in m] for k, m in data['matrices'].items()} == matrices,
         'FROZEN_MATRIX_EQUALITY')
    profile = json.loads((roots['current']/'results/PROFILES.json').read_text())
    frozen = json.loads((roots['I']/'initial/RECENTERED_CANDIDATE_NORMAL.json').read_text())
    coefficients = frozen['representation_background_coefficients_ascending']
    def nominal_polynomial(name, z):
        return sum((Q(c)*z**k for k, c in enumerate(coefficients[name])), Q(0))
    expected_panels = []
    for a, b, count in ((Q(65,64),Q(66,64),16),(Q(66,64),Q(9,8),112),(Q(9,8),Q(5,4),128)):
        step = (b-a)/count
        expected_panels += [(a+k*step, a+(k+2)*step) for k in range(0, count, 2)]
    need(len(data['panel_certificates']) == len(expected_panels) == 128, 'REPORT_PANEL_COUNT')
    for panel, expected in zip(data['panel_certificates'], expected_panels):
        a, b = map(Q, panel['r'])
        need((a,b) == expected and Q(panel['midpoint']) == (a+b)/2, 'REPORT_FIXED_COVERAGE')
        need(set(panel['targets']) == {str(T) for T in (Q(9,8),Q(5,4)) if b <= T}, 'REPORT_PANEL_TARGETS')
        for target in panel['targets'].values():
            need(set(target) == {'1','2','3','4'}, 'REPORT_MOMENT_ORDERS')
            for entries in target.values():
                need(len(entries) == 3, 'REPORT_MOMENT_COMPONENTS')
                for entry in entries:
                    fourth = box(entry['fourth_derivative_interval'])
                    need(Q(entry['remainder_upper']) == (b-a)**5*fourth.mag()/2880,
                         'REPORT_REMAINDER_FORMULA')
    results = []
    for sample in data['samples']:
        T = Q(sample['x'])
        need(T in (Q(9,8),Q(5,4)) and Q(sample['r']) == 0, 'REPORT_SAMPLE_LOCATION')
        U = {int(k): [box(x) for x in row] for k, row in sample['exact_u_moment_enclosures'].items()}
        sums = {k: [ExactInterval(0)]*3 for k in range(1,5)}
        errors = {k: [Q(0)]*3 for k in range(1,5)}
        for panel in data['panel_certificates']:
            if str(T) not in panel['targets']:
                continue
            for k, entries in panel['targets'][str(T)].items():
                for j, entry in enumerate(entries):
                    sums[int(k)][j] = sums[int(k)][j]+box(entry['simpson_interval'])
                    errors[int(k)][j] += Q(entry['remainder_upper'])
        C_exact = Q(3)/(2*Q(1,10**6)*nominal_polynomial('rho',T/Q(10**6))*nominal_polynomial('p',T/Q(10**6)))
        central = {0: U[0]}
        error_input = {0: [ExactInterval(0)]*5}
        for k in range(1,5):
            central[k], error_input[k] = [ExactInterval(0)]*5, [ExactInterval(0)]*5
            for j, component in enumerate((2,3,4)):
                saved_sum = sample['moment_simpson_sum_enclosures_before_scalar_endpoint'][str(k)][j]
                contained(saved_sum, sums[k][j], 'REPORT_OUTWARD_PANEL_SUM')
                need(errors[k][j] == Q(sample['moment_remainder_uppers_N_J_w'][str(k)][j]),
                     'REPORT_REMAINDER_SUM')
                central[k][component] = box(saved_sum)
                if k == 1 and component == 4:
                    central[k][component] += C_exact*U[0][3]
                error_input[k][component] = ExactInterval(-errors[k][j],errors[k][j])
                contained(sample['exact_u_moment_enclosures'][str(k)][component],
                          central[k][component]+error_input[k][component], 'REPORT_MOMENT_ASSEMBLY')
            for component in (0,1):
                need(U[k][component].lo == U[k][component].hi == 0, 'REPORT_ZERO_K_L_MOMENTS')
        V = p1_moments_binomial(T, profile)
        need(V == {int(k): list(map(Q,row)) for k,row in sample['velocity_moments_exact'].items()},
             'REPORT_INDEPENDENT_FROZEN_P1_MOMENTS')
        evaluated = moment_algebra_independent(T,U,V,matrices)
        central_evaluation = moment_algebra_independent(T,central,V,matrices)
        quadrature = moment_algebra_independent(T,error_input,{k:[Q(0)]*5 for k in V},matrices)
        for name, values in evaluated.items():
            for j, value in enumerate(values):
                contained(sample['evaluation'][name][j], value, 'REPORT_FULL_ASSEMBLY:'+name+':'+str(j))
                contained(sample['evaluation_omitting_quadrature_remainders'][name][j],
                          central_evaluation[name][j], 'REPORT_CENTRAL_ASSEMBLY:'+name+':'+str(j))
                need(Q(sample['quadrature_error_component_uppers'][name][j]) >= ExactInterval.cv(quadrature[name][j]).mag(),
                     'REPORT_QUADRATURE_PROPAGATION:'+name+':'+str(j))
        ji = box(sample['evaluation']['three_term'][3])
        results.append({'x':str(T),'all_stored_remainders_sums_moments_and_assemblies_passed':True,
                        'J_three_term_candidate_interval':sample['evaluation']['three_term'][3],
                        'J_candidate_strictly_negative':ji.hi < 0})
    need(len(results)==2, 'REPORT_TWO_SAMPLES')
    return {'schema':'HDBLAST_B2_INDEPENDENT_REPORT_REPLAY_V1',
            'status':'PASS_INDEPENDENT_B2_REPORT_ARITHMETIC_AND_SOURCE_REVIEW',
            'checker_sha256':sha(Path(__file__)), 'producer_source_sha256':sha(producer_path),
            'normal_report_sha256':sha(report_path), 'samples':results,
            'checks':['Every panel midpoint, coverage, target assignment and Simpson remainder formula.',
                      'All twelve moment integrals per target: outward sums, remainder sums and scalar endpoint assembly.',
                      'Independent local-binomial P1 integration verifies velocity moments through order five.',
                      'Unrounded rational interval replay verifies B0, B1, B2, two-term and three-term assembly and all histories.',
                      'Separate replay verifies every propagated pure-quadrature component cap.'],
            'limits':['Candidate arithmetic and source review only; no actual-response or observable sign claim.',
                      'Stored derivative enclosures are replayed algebraically; profile jets and quadrature are not recomputed.']}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo-root', required=True, type=Path)
    p.add_argument('--output', type=Path, default=Path(__file__).with_name('B2_POLYNOMIAL_CONTROLS_FINAL.json'))
    p.add_argument('--report', type=Path, default=Path(__file__).with_name('EXACT_B2_NORMAL.json'))
    p.add_argument('--review-output', type=Path, default=Path(__file__).with_name('B2_INDEPENDENT_REVIEW.json'))
    args = p.parse_args()
    need(not args.output.exists(), 'FRESH_OUTPUT_REQUIRED')
    need(not args.review_output.exists(), 'FRESH_REVIEW_OUTPUT_REQUIRED')
    auth_path = Path(__file__).resolve().parents[1]/'exact_boundary.py'
    auth = load('b2_parent_authentication', auth_path)
    roots = {key: args.repo_root.resolve()/'research'/name for key, name in auth.LABELS.items()}
    need(all(root not in args.output.resolve().parents for root in roots.values()), 'OUTPUT_INSIDE_FROZEN_CAPSULE')
    before = auth.authenticate(roots)
    kernel_path = roots['current']/'reflected_ridge.py'
    kernel = load('b2_authenticated_frozen_kernel', kernel_path)
    primitive_rows, negative_controls = polynomial_controls(kernel)
    triangle_rows = triangle_controls()
    producer_path = Path(__file__).with_name('exact_b2.py')
    producer = load('b2_new_producer_under_review', producer_path)
    profile = json.loads((roots['current']/'results/PROFILES.json').read_text())
    matrices = {k: [list(map(Q, row)) for row in mat] for k, mat in profile['matrices'].items()}
    producer_rows, velocity_rows = producer_controls(producer, kernel, matrices)
    review = replay_report(args.report, producer_path, args.repo_root.resolve(), roots, matrices)
    need(auth.authenticate(roots) == before, 'FROZEN_CLOSURE_CHANGED')
    result = {'schema': 'HDBLAST_B2_INDEPENDENT_POLYNOMIAL_CONTROLS_V1',
              'status': 'PASS_EXACT_B2_PRIMITIVE_AND_TRIANGLE_CONTROLS',
              'checker_sha256': sha(Path(__file__)), 'authentication_source_sha256': sha(auth_path),
              'producer_source_sha256': sha(producer_path),
              'frozen_kernel_sha256': sha(kernel_path), 'authentication': before,
              'complete_closures_authenticated_before_and_after': True,
              'primitive_comparison': 'Exact Lagrange integration of frozen reflected_kernel separately on each reflection branch; outer exact integration for V.',
              'B2_comparison': 'Independent monomial triangle moments for Z_x, Z_r, Z, characteristic history and boundary history; arbitrary coupled rational constant matrices.',
              'arithmetic': 'Exact Fraction only; no floating-point decisions',
              'primitive_controls': primitive_rows, 'negative_controls': negative_controls,
              'triangle_controls': triangle_rows,
              'producer_controls': producer_rows, 'P1_velocity_moment_controls': velocity_rows,
              'manual_producer_source_review': ['Negative G primitive signs correctly reflect u, du, v and Pv before applying the coupled matrices.',
                                                'Scalar W1 has C(T)f(T) endpoint minus integral Cprime f; Wk uses (p{k-1}C-pk Cprime)f for k>=2.',
                                                'Jet orders 2,3,6 supply the values, endpoint f_r, and fourth derivatives needed by the declared formulas.',
                                                'The candidate uses exact registered displacement in both constant-operator feedback terms.',
                                                'Authentication precedes frozen imports, and the complete closures are checked again after evaluation.'],
              'limits': ['Polynomial formula controls only; no actual-profile B2 evaluation, total-response error or observable sign claim.',
                         'No frozen producer, profile stage or boundary quadrature experiment was run.']}
    with args.output.open('x') as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
        handle.write('\n')
    review['polynomial_controls_sha256'] = sha(args.output)
    review['complete_closures_authenticated_before_and_after'] = True
    with args.review_output.open('x') as handle:
        json.dump(review, handle, indent=2, sort_keys=True)
        handle.write('\n')
    print(json.dumps({'status': result['status'], 'primitive_cases': len(primitive_rows),
                      'triangle_cases': len(triangle_rows), 'output': str(args.output),
                      'sha256': sha(args.output), 'review_output': str(args.review_output),
                      'review_sha256': sha(args.review_output)}, sort_keys=True))


if __name__ == '__main__':
    main()
