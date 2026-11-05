#!/usr/bin/env python3
"""Finite exact checks only; no HDBLAST analytic/numerical inputs are supplied."""
from fractions import Fraction as F
from itertools import permutations, product
from pathlib import Path
import argparse
import json


class CheckFailure(Exception):
    pass


CHECKS = []


def require(condition, message):
    if not condition:
        raise CheckFailure(message)


def record(name, condition):
    require(condition, name)
    CHECKS.append(name)


def pa(a, b):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, F(0)) + v
    return {k: v for k, v in out.items() if v}


def scale(a, c):
    return {k: v*c for k, v in a.items() if v*c}


def pm(a, b):
    out = {}
    for k, v in a.items():
        for l, w in b.items():
            t = tuple(x+y for x, y in zip(k, l))
            out[t] = out.get(t, F(0)) + v*w
    return {k: v for k, v in out.items() if v}


def pd(a, i):
    out = {}
    for k, v in a.items():
        if k[i]:
            l = list(k)
            l[i] -= 1
            out[tuple(l)] = v*k[i]
    return out


def pe(a, x):
    return sum((v*prod(t**n for t, n in zip(x, k)) for k, v in a.items()), F(0))


def prod(xs):
    ans = F(1)
    for x in xs:
        ans *= x
    return ans


def psum(xs):
    out = {}
    for x in xs:
        out = pa(out, x)
    return out


def polynomial_certificates():
    one = {(0, 0, 0): F(1)}
    xs = [{tuple(int(i == j) for i in range(3)): F(1)} for j in range(3)]
    s = psum(xs)
    sq = psum([pm(x, x) for x in xs])
    qdet = pa(scale(pm(s, s), F(1, 4)), scale(sq, F(-1, 2)))
    interp = {}
    values = []
    for v in product((-1, 1), repeat=3):
        w = one
        for x, sign in zip(xs, v):
            w = pm(w, scale(pa(one, scale(x, sign)), F(1, 2)))
        fv = pe(qdet, v)
        values.append(fv)
        interp = pa(interp, scale(w, fv))
    record('whole_cube_lower_bound_interpolation_identity',
           qdet == pa(interp, scale(pa(scale(one, 3), scale(sq, -1)), F(1, 4))))
    differences = psum([pm(pa(xs[i], scale(xs[j], -1)), pa(xs[i], scale(xs[j], -1)))
                       for i in range(3) for j in range(i+1, 3)])
    record('whole_cube_upper_bound_sum_of_squares_identity',
           pa(scale(one, F(3, 4)), scale(qdet, -1)) ==
           scale(pa(pa(scale(one, 3), scale(sq, -1)), differences), F(1, 4)))
    record('sharp_absolute_coefficient_5_over_4', min(values) == F(-5, 4)
           and max(values) == F(3, 4) and pe(qdet, (1, 1, -1)) == F(-5, 4))
    fplus = pa(qdet, scale(s, F(-3, 2)))
    for i in range(3):
        other = [j for j in range(3) if j != i]
        cert = scale(psum([pa(one, xs[i]), pa(one, scale(xs[other[0]], -1)),
                           pa(one, scale(xs[other[1]], -1))]), F(-1, 2))
        record('correlated_endpoint_coordinate_monotonicity_' + str(i), pd(fplus, i) == cert)
    record('correlated_endpoint_range_9', pe(fplus, (-1, -1, -1)) == F(21, 4)
           and pe(fplus, (1, 1, 1)) == F(-15, 4))
    # The note proves nonnegative barycentric weights and range convexity.
    for b in (F(-3, 2), F(-1), F(0), F(1, 3), F(3, 2)):
        theta = (b+F(3, 2))/3
        fminus = pa(qdet, scale(s, F(3, 2)))
        record('correlated_affine_interpolation_' + str(b),
               pa(qdet, scale(s, -b)) == pa(scale(fplus, theta), scale(fminus, 1-theta)))


def matmul(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), F(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def transpose(a):
    return list(map(list, zip(*a)))


def madd(*xs):
    return [[sum((x[i][j] for x in xs), F(0)) for j in range(len(xs[0][0]))]
            for i in range(len(xs[0]))]


def mscale(a, c):
    return [[c*x for x in row] for row in a]


def trace(a):
    return sum((a[i][i] for i in range(len(a))), F(0))


def inverse(a):
    n = len(a)
    x = [list(row)+[F(int(i == j)) for j in range(n)] for i, row in enumerate(a)]
    for i in range(n):
        pivot = next((j for j in range(i, n) if x[j][i]), None)
        require(pivot is not None, 'singular rational matrix')
        x[i], x[pivot] = x[pivot], x[i]
        z = x[i][i]
        x[i] = [v/z for v in x[i]]
        for j in range(n):
            if j != i:
                z = x[j][i]
                x[j] = [v-z*w for v, w in zip(x[j], x[i])]
    return [row[n:] for row in x]


def jmul(a, b):
    return [sum((a[k]*b[i-k] for k in range(i+1)), F(0)) for i in range(3)]


def determinant_jet(g, h, k):
    ans = [F(0)]*3
    for perm in permutations(range(3)):
        sign = (-1)**sum(perm[i] > perm[j] for i in range(3) for j in range(i+1, 3))
        term = [F(sign), F(0), F(0)]
        for i, j in enumerate(perm):
            term = jmul(term, [g[i][j], h[i][j], k[i][j]/2])
        ans = [a+b for a, b in zip(ans, term)]
    return ans


def determinant_and_pullback_checks():
    for seed in range(1, 21):
        l = [[F(1), F(0), F(0)], [F(seed, 7), F(1), F(0)],
             [F(-seed, 11), F(seed+1, 13), F(1)]]
        g = matmul(matmul(l, [[F(2), F(0), F(0)], [F(0), F(3), F(0)],
                            [F(0), F(0), F(5)]]), transpose(l))
        h = [[F((i+j+seed)%7-3, 5) for j in range(3)] for i in range(3)]
        k = [[F((i*j+seed)%5-2, 3) for j in range(3)] for i in range(3)]
        a, b = matmul(inverse(g), h), matmul(inverse(g), k)
        if seed == 1:
            record('noncommuting_metric_and_first_derivative_witness', matmul(g, h) != matmul(h, g))
        d = determinant_jet(g, h, k)
        record('determinant_first_jet_' + str(seed), d[1]/d[0] == trace(a))
        record('determinant_second_jet_' + str(seed),
               2*d[2]/d[0] == trace(b)+trace(a)**2-trace(matmul(a, a)))
        record('sqrt_determinant_second_ratio_' + str(seed),
               d[2]/d[0]-F(1, 4)*(d[1]/d[0])**2 ==
               trace(b)/2+trace(a)**2/4-trace(matmul(a, a))/2)
    # Differentiate T(e)^T H(e) T(e) directly with polynomial jets.
    t = [[F(1), F(2), F(-1)], [F(3), F(-1), F(2)]]
    t1, t2 = mscale(t, F(2, 7)), [[F(2), F(1), F(0)], [F(-1), F(3), F(2)]]
    h, h1, h2 = [[F(-2), F(1)], [F(1), F(3)]], [[F(2), F(3)], [F(3), F(5)]], [[F(1), F(-1)], [F(-1), F(4)]]
    jets = [[F(0)]*3 for _ in range(9)]
    for i in range(3):
        for j in range(3):
            for r in range(2):
                for s in range(2):
                    term = jmul(jmul([t[r][i], t1[r][i], t2[r][i]/2],
                                     [h[r][s], h1[r][s], h2[r][s]/2]),
                                [t[s][j], t1[s][j], t2[s][j]/2])
                    jets[3*i+j] = [a+b for a, b in zip(jets[3*i+j], term)]
    tri = lambda left, middle, right: matmul(matmul(transpose(left), middle), right)
    first = madd(tri(t, h1, t), tri(t1, h, t), tri(t, h, t1))
    second = madd(tri(t, h2, t), mscale(tri(t1, h1, t), 2), mscale(tri(t, h1, t1), 2),
                  tri(t2, h, t), tri(t, h, t2), mscale(tri(t1, h, t1), 2))
    record('full_pullback_metric_first_jet', all(jets[3*i+j][1] == first[i][j]
                                               for i in range(3) for j in range(3)))
    record('full_pullback_metric_second_jet', all(2*jets[3*i+j][2] == second[i][j]
                                                for i in range(3) for j in range(3)))
    # g(e,x)=c00+c10 e+c01 x+c20 e^2/2+c11 e x+c02 x^2/2 at x(0)=0.
    c00, c10, c01, c20, c11, c02, z, w = map(F, (2, 3, 5, 7, 11, 13, 17, 19))
    x = [F(0), z, w/2]
    ex = jmul([F(0), F(1), F(0)], x)
    xx = jmul(x, x)
    composition = [c00, c10+c01*z, c20/2+c01*w/2+c11*ex[2]+c02*xx[2]/2]
    record('ambient_composition_clock_and_embedding_second_jet',
           2*composition[2] == c20+2*c11*z+c02*z*z+c01*w)


def normalized_quotient_checks():
    for seed in range(1, 31):
        ys = [F(-1), F(-1, 3), F(0), F(1, 5), F(1)]
        js = [F(seed+i+2, 7) for i in range(5)]
        es = [F((seed+2*i)%9-4, 5) for i in range(5)]
        bs = [F((3*seed+i)%7-3, 11) for i in range(5)]
        a0 = [F(seed-i, 13) for i in range(5)]
        a1 = [F(seed+2*i, 17) for i in range(5)]
        a2 = [F(seed-3*i, 19) for i in range(5)]
        ts = [b+e*e for b, e in zip(bs, es)]
        d = sum(y*y*j for y, j in zip(ys, js))
        nu = [y*y*j/d for y, j in zip(ys, js)]
        mu = sum(n*e for n, e in zip(nu, es))
        nt = sum(n*t for n, t in zip(nu, ts))
        nb = sum(n*b for n, b in zip(nu, bs))
        var = sum(n*(e-mu)**2 for n, e in zip(nu, es))
        n0 = sum(y*j*a for y, j, a in zip(ys, js, a0))
        n1 = sum(y*j*(a+e*b) for y, j, a, e, b in zip(ys, js, a1, es, a0))
        n2 = sum(y*j*(a+2*e*b+t*c) for y, j, a, e, b, t, c in zip(ys, js, a2, es, a1, ts, a0))
        exact = n2/d-2*n1*mu/d-n0*nt/d+2*n0*mu*mu/d
        grouped = sum(y*j/d*(a+2*b*(e-mu)+c*((t-nt)-2*mu*(e-mu)))
                      for y, j, a, b, c, e, t in zip(ys, js, a2, a1, a0, es, ts))
        score = sum(y*j/d*(a+2*b*(e-mu)+c*((e-mu)**2-var+v-nb))
                    for y, j, a, b, c, e, v in zip(ys, js, a2, a1, a0, es, bs))
        record('quotient_grouping_' + str(seed), exact == grouped == score)


def metric_guard(a, adm, p_budget, q_budget, e1, e2, lower, upper, contrast):
    require(0 < a < adm, 'common admitted interval not strictly interior')
    require(min(p_budget, q_budget, e1, e2) >= 0, 'negative norm budget')
    require(0 < lower <= 1 <= upper, 'invalid anchored spectral interval')
    require(contrast > 0 and contrast**2*lower**3 >= upper**3,
            'rational contrast certificate fails')
    p, q = p_budget/lower, q_budget/lower
    h = 2*contrast*(e2+6*p*e1+a*e1*(3*q+9*p*p))
    require(a*h < 106, 'strict endpoint guard fails')
    return h


def rejection_test(name, callback):
    try:
        callback()
    except CheckFailure:
        CHECKS.append(name)
        return
    raise CheckFailure('missing rejection: ' + name)


def guard_checks():
    args = [F(1, 100), F(1, 10), F(2), F(3), F(5), F(7), F(49, 50), F(51, 50), F(11, 10)]
    h = metric_guard(*args)
    record('synthetic_guard_is_not_a_physical_amplitude', args[0]*h < 106)
    for index, replacement, name in [(0, F(1, 10), 'reject_admissibility_equality'),
                                      (2, F(-1), 'reject_negative_budget'),
                                      (6, F(0), 'reject_zero_spectral_gap'),
                                      (8, F(1), 'reject_insufficient_contrast'),
                                      (5, F(100000), 'reject_failed_endpoint_guard')]:
        changed = list(args)
        changed[index] = replacement
        rejection_test(name, lambda changed=changed: metric_guard(*changed))
    # Exact equality at 106 must also fail, with identity metric.
    rejection_test('reject_endpoint_guard_equality',
                   lambda: metric_guard(F(1), F(2), F(0), F(0), F(0), F(53), F(1), F(1), F(1)))
    record('zero_budget_weak_bound_is_zero',
           metric_guard(F(1), F(2), F(0), F(0), F(0), F(0), F(1), F(1), F(1)) == 0)
    x, n = F(3, 5), 4
    exp_envelope = (1-x/n)**(-n)
    record('optional_exponential_envelope_positive_rational', 0 <= x < n and exp_envelope > 1)
    return {'a': str(args[0]), 'a_adm': str(args[1]), 'P': str(args[2]), 'Q': str(args[3]),
            'E1': str(args[4]), 'E2': str(args[5]), 'lambda': str(args[6]),
            'Lambda': str(args[7]), 'C': str(args[8]), 'H': str(h), 'aH': str(args[0]*h),
            'meaning': 'synthetic algebraic test only; no actual HDBLAST budget is evaluated'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    opts = parser.parse_args()
    if opts.output:
        require(opts.output.is_absolute(), '--output requires an absolute path')
    polynomial_certificates()
    determinant_and_pullback_checks()
    normalized_quotient_checks()
    synthetic = guard_checks()
    report = {'status': 'PASS', 'arithmetic': 'Python standard-library fractions.Fraction',
              'checks': len(CHECKS), 'check_names': CHECKS, 'synthetic_guard': synthetic,
              'limitations': ['finite algebra checker; transparent mathematical proof in companion note',
                              'no PDE, observer-domain, metric or scalar analytic budget is certified',
                              'no physical or numerical amplitude claim',
                              'no network, external packages, frozen-file writes or legacy campaign replay']}
    serialized = json.dumps(report, indent=2, sort_keys=True)+'\n'
    if opts.output:
        opts.output.write_text(serialized)
    print(serialized, end='')


if __name__ == '__main__':
    main()
