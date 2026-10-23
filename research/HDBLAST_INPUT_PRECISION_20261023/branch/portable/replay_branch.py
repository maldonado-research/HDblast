#!/usr/bin/env python3
"""Portable exact-rational replay of a retained-branch necessary enclosure.

Python standard library only. No project imports, network, or absolute paths.
Upstream premises: the retained theorem supplies existence and uniqueness in its
specified regular Dirichlet tube, its old beta enclosure, regular-cone asymptotic
data, scalar sign/monotonicity, and the positive-a constraint branch. This script
does not authenticate or reprove those premises. It proves the conditional
polynomial comparison and necessary contraction of that same existing root.
"""
import argparse
from fractions import Fraction as F
import json
import math
from pathlib import Path

# Complete degree4 coefficient listing, x=eta and h=rho^-2. P=x*Z.
Z_ROWS = '''
0 0 14/9
0 1 -42/5
1 0 25/23
0 2 18387/100
1 1 405/253
2 0 -9104/42849
0 3 -1852389/250
1 2 -538731/14168
2 1 -7457194/2531265
3 0 133240/1350537
0 4 14993308737/40000
1 3 37378017/25300
2 2 1042602915537/36382715600
3 1 3530653345/1723285212
4 0 -2911768910/27676554741
'''
A_ROWS = '''
0 0 1/9
0 1 9/2
0 2 -729/8
2 0 -7/27
0 3 59049/16
2 1 7/10
3 0 -25/207
0 4 -23914845/128
3 1 -45/506
4 0 2276/128547
'''

def require(test, stage):
    if not test:
        raise ValueError(stage)

def parse(rows):
    return {(int(i), int(j)): F(c) for i, j, c in
            (line.split() for line in rows.splitlines() if line.strip())}

def clean(p):
    return {m: c for m, c in p.items() if c}

def add(*polys):
    out = {}
    for p in polys:
        for m, c in p.items():
            out[m] = out.get(m, F(0)) + c
    return clean(out)

def scale(p, c):
    return clean({m: c*v for m, v in p.items()})

def mul(p, q):
    out = {}
    for (i, j), c in p.items():
        for (k, l), d in q.items():
            m = (i+k, j+l)
            out[m] = out.get(m, F(0)) + c*d
    return clean(out)

def power(p, n):
    out = {(0, 0): F(1)}
    for _ in range(n):
        out = mul(out, p)
    return out

def derivative(p, axis):
    return clean({(i-(axis == 0), j-(axis == 1)): c*(i if axis == 0 else j)
                  for (i, j), c in p.items() if (i if axis == 0 else j)})

def evaluate(p, x, h):
    return sum((c*x**i*h**j for (i, j), c in p.items()), F(0))

def majorant(p, x, h):
    return sum((abs(c)*x**i*h**j for (i, j), c in p.items()), F(0))

def coarse(p, x, h):
    c = p.get((0, 0), F(0))
    radius = majorant(add(p, {(0, 0): -c}), x, h)
    return c-radius, c+radius

class Interval:
    def __init__(self, lo, hi=None):
        self.lo, self.hi = F(lo), F(lo if hi is None else hi)
        require(self.lo <= self.hi, 'interval ordering')

    def __add__(self, other):
        other = other if isinstance(other, Interval) else Interval(other)
        return Interval(self.lo+other.lo, self.hi+other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)

    def __mul__(self, other):
        other = other if isinstance(other, Interval) else Interval(other)
        products = [self.lo*other.lo, self.lo*other.hi,
                    self.hi*other.lo, self.hi*other.hi]
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def __pow__(self, n):
        require(isinstance(n, int) and n >= 0, 'nonnegative interval power')
        if n == 0:
            return Interval(1)
        if n % 2:
            return Interval(self.lo**n, self.hi**n)
        lower = 0 if self.lo <= 0 <= self.hi else min(self.lo**n, self.hi**n)
        return Interval(lower, max(self.lo**n, self.hi**n))

    def encode(self):
        return [str(self.lo), str(self.hi)]

def interval_value(p, x, h):
    return sum((c*x**i*h**j for (i, j), c in p.items()), Interval(0))

def outward_sqrt(q, places=32):
    require(q >= 0, 'nonnegative square root input')
    base = 10**places
    numerator = math.isqrt((q.numerator*base**2)//q.denominator)
    lo, hi = F(numerator, base), F(numerator+1, base)
    require(lo**2 <= q <= hi**2, 'outward square root')
    return lo, hi

def floor_decimal(q, places):
    base = 10**places
    return F((q.numerator*base)//q.denominator, base)

def ceil_decimal(q, places):
    return -floor_decimal(-q, places)

def replay(control=None):
    Z, A = parse(Z_ROWS), parse(A_ROWS)
    if control == 'coefficient':
        Z[(1, 0)] += F(1, 10**12)
    x, h, one = {(1, 0): F(1)}, {(0, 1): F(1)}, {(0, 0): F(1)}
    W = {(0, 0): F(1, 3), (2, 0): F(1), (3, 0): F(1, 3)}
    U = add(scale(power(derivative(W, 0), 2), F(1, 2)),
            scale(power(W, 2), F(-2, 3)))
    U1 = derivative(U, 0)
    require(all(i > 0 for i, j in U1), 'U prime divisible by x')
    M = {(i-1, j): c for (i, j), c in U1.items()}
    B = add(scale(mul(h, derivative(Z, 1)), 2), scale(Z, -4))
    R = add(M, scale(power(Z, 2), -1),
            scale(mul(mul(x, derivative(Z, 0)), Z), -1), mul(A, B))
    D = add(h, scale(mul(power(x, 2), power(Z, 2)), F(1, 12)),
            scale(U, F(-1, 6)), scale(power(A, 2), -1))
    require(all(i+j >= 5 for i, j in R), 'degree4 residual cancellation')
    require(all(i+j >= 5 for i, j in D), 'degree4 square root cancellation')

    eps, h0 = F(9, 100000), F(1, 100)
    amin, amax, damping = F(1, 9)-F(1, 1000000), F(3, 20), F(37, 20)
    if control == 'damping':
        damping = F(2)
    zlo, zhi = coarse(Z, eps, h0)
    alo, ahi = coarse(A, eps, h0)
    blo, bhi = coarse(B, eps, h0)
    xzx = majorant(mul(x, derivative(Z, 0)), eps, h0)
    require(0 < zlo and zhi < 2, 'polynomial slope range')
    require(alo > 0, 'polynomial square root positivity')
    require(bhi < 0, 'favorable secant correction sign')
    secant_lower = 4*amin+zlo-xzx
    require(secant_lower >= damping, 'secant damping bound')
    require(damping-12*amax > 0, 'memory damping bound')
    Mlo, Mhi = coarse(M, eps, F(0))
    require(0 < Mlo and Mhi < 4, 'regular slope invariant interval')
    Udiff_poly = add(U, {(0, 0): F(2, 27)})
    Udiff_over_x2 = {(i-2, j): c for (i, j), c in Udiff_poly.items()}
    require(coarse(Udiff_over_x2, eps, F(0))[0] > 0, 'U minus U0 nonnegative')
    Udiff = majorant(Udiff_poly, eps, F(0))
    require(-F(2, 27)+Udiff < 0, 'negative potential')
    require(F(1, 81)-Udiff/6 > amin**2, 'actual a lower bound')
    require(h0+F(1, 81)+eps**2/3 < amax**2, 'actual a upper bound')

    # e=z-Z has e'=R+B*D/(a+A)-K*e exactly, with
    # K=4a+z+Z+xZ_x-B*x^2*(z+Z)/(12*(a+A)) >= damping.
    # On the domain, B<0 and z+Z>=0, so the last term is favorable.
    factor = max(abs(blo), abs(bhi))/(amin+alo)
    C = [F(0)]*6
    for p, weight in ((R, F(1)), (D, factor)):
        for (i, j), c in p.items():
            C[min(j, 5)] += weight*abs(c)*eps**i*h0**max(j-5, 0)
    barrier = [C[j]/(damping-2*j*amax) for j in range(6)]

    def E(H):
        return 2*(H/h0)**6 + sum((barrier[j]*H**j for j in range(6)), F(0))

    # Along a regular negative scalar branch, |x| is nondecreasing and h
    # decreases; d(|x|^i*h^j)/dy=(i*z-2*j*a)|x|^i*h^j.
    # Entry |e|<2, K>=damping, j<=5, and damping-12amax>0 prove |e|<=E.
    c, delta = F(96714321425692, 161839258930731), F(1, 1000)
    beta0, oldradius = -F(9, 64)*c*delta, F(799, 250)*delta**2
    old = Interval(beta0-oldradius, beta0+oldradius)
    require(old.hi < 0 and max(abs(old.lo), abs(old.hi)) < eps,
            'retained negative scalar enclosure')
    f = add(scale(one, 1+c), scale(x, c))
    H = add(scale(add(scale(mul(W, f), F(1, 9)),
                      scale(derivative(W, 0), -c/F(12))), delta),
            scale(add(scale(power(f, 2), F(1, 36)),
                      scale(one, -c*c/F(48))), delta**2))
    p = add(scale(derivative(W, 0), -1), scale(one, -delta*c/2))
    a = add(scale(W, F(1, 3)), scale(f, delta/6))
    Hbox = interval_value(H, old, Interval(0))
    require(0 < Hbox.lo < Hbox.hi < h0, 'shell inside comparison domain')
    Hp, Zx, Zh = derivative(H, 0), derivative(Z, 0), derivative(Z, 1)
    W1, W2 = derivative(W, 0), derivative(derivative(W, 0), 0)

    def g(t):
        return t*evaluate(Z, t, evaluate(H, t, F(0))) + evaluate(W1, t, F(0)) + delta*c/2

    def gp(t):
        ht = evaluate(H, t, F(0))
        return (evaluate(Z, t, ht) + t*(evaluate(Zx, t, ht) +
                evaluate(Zh, t, ht)*evaluate(Hp, t, F(0))) + evaluate(W2, t, F(0)))

    seed, seed_steps = beta0, []
    for iteration in range(1, 5):
        require(gp(seed) != 0, 'nonzero seed derivative')
        seed = floor_decimal(seed-g(seed)/gp(seed), 24)
        seed_steps.append({'iteration': iteration, 'beta_seed': str(seed),
                           'polynomial_residual': str(g(seed))})
    gpbox = (interval_value(Z, old, Hbox) +
             old*(interval_value(Zx, old, Hbox) +
                  interval_value(Zh, old, Hbox)*interval_value(Hp, old, Interval(0))) +
             interval_value(W2, old, Interval(0)))
    require(gpbox.lo > 0, 'shell polynomial derivative positivity')
    error = eps*E(Hbox.hi)
    radius = (abs(g(seed))+error)/gpbox.lo
    require(old.lo <= seed-radius and seed+radius <= old.hi,
            'necessary enclosure lies inside retained region')
    # Mean-value bound: the retained existing root beta* satisfies
    # |g(beta*)|<=eps*E(Hbox.hi), hence |beta*-seed|<=radius.
    lo, hi = floor_decimal(seed-radius, 20), ceil_decimal(seed+radius, 20)
    if control == 'rounding':
        lo, hi = ceil_decimal(seed-radius, 20), floor_decimal(seed+radius, 20)
    require(lo <= seed-radius and hi >= seed+radius, 'outward beta rounding')
    beta = Interval(lo, hi)
    newH = interval_value(H, beta, Interval(0))
    newp = interval_value(p, beta, Interval(0))
    newa = interval_value(a, beta, Interval(0))
    rho = [outward_sqrt(1/newH.hi)[0], outward_sqrt(1/newH.lo)[1]]
    checks = {
        'M_bounds': [str(Mlo), str(Mhi)],
        'Z_bounds': [str(zlo), str(zhi)],
        'A_bounds': [str(alo), str(ahi)],
        'B_bounds': [str(blo), str(bhi)],
        'eta_Zeta_abs': str(xzx), 'secant_damping_lower': str(secant_lower),
        'selected_damping': str(damping), 'memory_damping_margin': str(damping-12*amax)}
    return {
        'exact_checks': checks, 'residual_majorant_Cj': [str(v) for v in C],
        'barrier_Ej': [str(v) for v in barrier], 'old_beta': old.encode(),
        'old_H2': Hbox.encode(), 'endpoint_error_z': str(E(Hbox.hi)),
        'endpoint_error_scalar_residual': str(error), 'scalar_polynomial_derivative': gpbox.encode(),
        'seed_steps': seed_steps, 'pre_rounding_radius': str(radius), 'beta': beta.encode(),
        'beta_midpoint': str((lo+hi)/2), 'beta_radius': str((hi-lo)/2),
        'provisional_beta_precision_target_passed': (hi-lo)/2 <= F(1, 10**11),
        'shell_outer_projections': {'H2': newH.encode(), 'p': newp.encode(),
                                    'a': newa.encode(), 'rho': [str(v) for v in rho]}}

def project_branch_report(report):
    keys = ('residual_majorant_Cj', 'barrier_Ej', 'old_beta', 'old_H2',
            'scalar_polynomial_derivative', 'seed_steps', 'beta', 'beta_midpoint',
            'provisional_beta_precision_target_passed', 'shell_outer_projections')
    out = {key: report[key] for key in keys}
    for key in ('endpoint_error_z', 'endpoint_error_scalar_residual',
                'pre_rounding_radius', 'beta_radius'):
        out[key] = report[key]['exact']
    out['exact_checks'] = {}
    for key in ('M_bounds', 'Z_bounds', 'A_bounds', 'B_bounds'):
        out['exact_checks'][key] = [v['exact'] for v in report['exact_checks'][key]]
    for key in ('eta_Zeta_abs', 'secant_damping_lower'):
        out['exact_checks'][key] = report['exact_checks'][key]['exact']
    for key in ('selected_damping', 'memory_damping_margin'):
        out['exact_checks'][key] = report['exact_checks'][key]
    return out

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--negative-control', choices=('coefficient', 'damping', 'rounding'))
    parser.add_argument('--compare-branch-report', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        numerical = replay(args.negative_control)
        reference = json.loads((Path(__file__).parent/'NUMERICAL_REFERENCE.json').read_text())
        require(numerical == reference, 'exact portable numerical reference comparison')
        if args.compare_branch_report:
            branch = json.loads(args.compare_branch_report.read_text())
            require(numerical == project_branch_report(branch), 'exact branch numerical field comparison')
        result = {
            'status': 'PASS_CONDITIONAL_RETAINED_BRANCH_RATIONAL_REPLAY',
            'arithmetic': 'Python fractions.Fraction and integer isqrt; no float arithmetic',
            'upstream_premises': [
                'Existence and uniqueness of the retained regular branch inside its specified Dirichlet tube',
                'The retained beta enclosure beta0 +/- (799/250)*delta^2',
                'Regular cone data and expansions, scalar sign and monotonicity, and positive-a constraint branch'],
            'scope': 'Necessary enclosure of that existing root; no theorem authentication, new endpoint-map inclusion, initial-data error join, or total response claim',
            'numerical_reference_exact_match': True, 'beta': numerical['beta'],
            'beta_radius': numerical['beta_radius'], 'numerical_fields': numerical}
        if args.output:
            args.output.write_text(json.dumps(result, indent=2)+'\n')
        print(json.dumps({key: result[key] for key in ('status', 'numerical_reference_exact_match',
                                                       'beta', 'beta_radius')}, sort_keys=True))
    except (ValueError, KeyError, OSError) as error:
        print(json.dumps({'status': 'FAIL', 'reason': str(error)}, sort_keys=True))
        return 1
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
