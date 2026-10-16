#!/usr/bin/env python3
"""Independent ordinary-derivative and provenance review of the actual oracle.

The checker uses ordinary, binomial-weighted derivatives, rather than the
producer's factorial-normalized Taylor coefficients.  Finite manufactured
witnesses supplement the accompanying analytic all-point review.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from math import comb, factorial
from pathlib import Path

ROOT = Path(__file__).resolve().parent
NUMERIC = ROOT.parent/'numeric'
ORDER = 4
KEYS = tuple((i, n-i) for n in range(ORDER+1) for i in range(n+1))


def need(ok, text):
    if not ok:
        raise ValueError(text)


class D:
    """Exact ordinary derivatives, deliberately not normalized coefficients."""
    def __init__(self, entries=None):
        self.d = {k: F((entries or {}).get(k, 0)) for k in KEYS}

    @staticmethod
    def cv(value):
        return value if isinstance(value, D) else D({(0, 0): value})

    def __add__(self, other):
        other = D.cv(other)
        return D({k: self.d[k]+other.d[k] for k in KEYS})

    __radd__ = __add__

    def __neg__(self):
        return D({k: -v for k, v in self.d.items()})

    def __sub__(self, other):
        return self+-D.cv(other)

    def __rsub__(self, other):
        return D.cv(other)+-self

    def __mul__(self, other):
        other = D.cv(other)
        return D({(i, j): sum((comb(i, k)*comb(j, l)*self.d[k, l]*other.d[i-k, j-l]
                              for k in range(i+1) for l in range(j+1)), F(0))
                  for i, j in KEYS})

    __rmul__ = __mul__

    def inv(self):
        out = {(0, 0): 1/self.d[0, 0]}
        for i, j in KEYS[1:]:
            out[i, j] = -sum((comb(i, k)*comb(j, l)*self.d[k, l]*out[i-k, j-l]
                             for k in range(i+1) for l in range(j+1) if k+l), F(0))/self.d[0, 0]
        return D(out)

    def __truediv__(self, other):
        return self*D.cv(other).inv()

    def __rtruediv__(self, other):
        return D.cv(other)*self.inv()

    def __pow__(self, n):
        if n < 0:
            return self.inv()**(-n)
        value = D.cv(1)
        for _ in range(n):
            value = value*self
        return value

    def deriv(self, x=0, r=0):
        return D({(i, j): self.d.get((i+x, j+r), 0) for i, j in KEYS})


def exp_negative(value, terms=420):
    """Independent direct alternating series, without argument-halving."""
    need(0 <= value <= 70, 'test exponential range')
    n = terms if terms % 2 == 0 else terms+1
    upper = sum(((-value)**k/F(factorial(k)) for k in range(n+1)), F(0))
    lower = upper- value**(n+1)/F(factorial(n+1))
    need(0 < lower <= upper, 'independent alternating exponential')
    return lower, upper


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def contains(box, value):
    return box.lo <= value <= box.hi


def intersects(box, lo, hi):
    return max(box.lo, lo) <= min(box.hi, hi)


def check_pins(a):
    freeze = json.loads((NUMERIC/'ORACLE_FREEZE.json').read_text())
    for path, pin in freeze['files'].items():
        raw = (NUMERIC/path).read_bytes()
        need(len(raw) == pin['bytes'] and hashlib.sha256(raw).hexdigest() == pin['sha256'], 'oracle freeze '+path)
    data, errors, pins = a.inputs()
    for path, pin in pins.items():
        raw = (NUMERIC/path).read_bytes()
        need(len(raw) == pin['bytes'] and hashlib.sha256(raw).hexdigest() == pin['sha256'], 'authenticated local frozen predecessor copy '+path)
    ordinary = json.loads((NUMERIC/'INPUT_PINS.json').read_text())
    need(errors['input_sha256'] == {k: v['sha256'] for k, v in ordinary.items()}, 'error certificate uses same five inputs')
    need(data['coefficients']['input_sha256']['HIGHER_JETS_NORMAL.json'] == ordinary['master']['sha256'], 'master/coefficient cross-pin')
    for field in ('rho', 'eta', 'p', 'a'):
        need(data['coefficients']['first_exit'][field]['shell']['exact'] == data['shell']['boxes'][field+'_b']['exact'], 'shell/coefficient box '+field)
    return data, errors, freeze


def to_t(a, value):
    return a.T({k: v/F(factorial(k[0])*factorial(k[1])) for k, v in value.d.items() if sum(k) <= 3})


def manufactured(a, injection):
    x = D({(0, 0): F(2, 7), (1, 0): 1})
    r = D({(0, 0): F(5, 4), (0, 1): 1})
    ell = a.ELL
    W = 2+x+3*r+5*x*r+7*x*x*r+11*x*r*r+13*r**3
    f = 17+x*r+3*x*x*r+5*x*r*r+7*r**3
    b = 19+2*x+3*r+5*x*x*r+7*x*r*r+11*r**3
    d = 23+3*x+5*r+7*x*x*r+11*x*r*r+13*r**3
    rho = 29+ell*r+ell**2*r*r
    p = -31+ell*r+2*ell**2*r*r
    aa = -rho.deriv(r=1)/(ell*rho*rho)
    R = 37+ell*r+3*ell**2*r*r
    B = -41+2*ell*r+5*ell**2*r*r
    k = 43+3*ell*r+7*ell**2*r*r
    c = 1+ell*x+ell**2*x*x
    H = ell*x+ell**2*x*x
    vtheta = 1-ell*x+3*ell**2*x*x
    u = R*vtheta*W
    v = B*vtheta*(W.deriv(r=1)/ell+k*W)
    # Independently form coordinate components (Eq.4), then divide by
    # the frame factors.  This route differs from the producer Eq.7.
    dy = lambda q: -q.deriv(r=1)/(ell*rho)
    dt = lambda q: q.deriv(x=1)/ell
    raw = {'A': -4*u+2*dy(f),
           'B': (dt(f)-rho*rho*dy(b))/rho,
           'C': (f+rho*rho*c*c*dy(d))/(rho*c),
           'D': -2*rho*rho*(u+aa*f+dt(b))/(rho*rho),
           'E': rho*rho*(c*c*dt(d)-b)/(rho*rho*c),
           'L': 2*rho*rho*c*c*(u+aa*f+H*b)/(rho*rho*c*c),
           'M': 2*rho*rho*c*c*d/(rho*rho*c*c)}
    scalar = v+p*f
    amps, s = a.frame_amplitudes(to_t(a, W), list(map(lambda q: to_t(a, q), (f, b, d))),
        *map(lambda q: to_t(a, q), (rho, p, aa, R, B, k, c, H, vtheta)),
        injection=injection if injection in ('wrong_scalar_shift', 'drop_radial_rho') else 'none')
    count = 0
    suffix = {'': (0, 0), 'y': (0, 1), 't': (1, 0), 'yy': (0, 2), 'yt': (1, 1), 'tt': (2, 0)}
    for name, actual in list(amps.items())+[('s', s)]:
        expected = scalar if name == 's' else raw[name]
        actual_values = a.physical_derivatives(actual, to_t(a, rho), to_t(a, aa))
        for end, (ix, iy) in suffix.items():
            target = expected
            for _ in range(iy):
                target = dy(target)
            for _ in range(ix):
                target = dt(target)
            value = actual_values[end]
            if injection == 'drop_radial_conversion' and end == 'yy':
                value = actual.derivative(0, 2)/(ell*to_t(a, rho).val())**2
            need(contains(value, target.d[0, 0]), 'independent coordinate-to-frame ordinary derivative '+name+end)
            count += 1
    # Pure normalized Taylor quotient derivatives separately reach total3.
    q = (x*x+3*x*r+r**3)/(2+r+x*r)
    tq = (to_t(a, x)**2+3*to_t(a, x)*to_t(a, r)+to_t(a, r)**3)/(2+to_t(a, r)+to_t(a, x)*to_t(a, r))
    for key in a.KEYS:
        target = q.d[key]
        value = tq.derivative(*key)
        if injection == 'wrong_factorial' and key == (2, 1):
            value = tq.c[key]
        need(contains(value, target), 'independent ordinary quotient derivative '+str(key))
        count += 1
    return count


def check_bump_comparisons(a, data):
    count = 0
    for point in (F(65, 64), F(5, 4), F(3, 2), F(7, 4), F(127, 64)):
        s = D({(0, 0): point, (0, 1): 1})
        exponent = 4-1/((s-1)*(2-s))
        q = -exponent.d[0, 0]
        elo, ehi = exp_negative(q)
        # Complete Bell polynomial recursion for exp(g), ordinary jets.
        bell = [F(1)]
        for n in range(3):
            bell.append(sum((comb(n, k)*exponent.d[0, k+1]*bell[n-k] for k in range(n+1)), F(0)))
        for sign in (1, -1):
            boxes = a.qodd(a.I(sign*point), 3, data)
            for n in range(4):
                factor = sign**(n+1)*bell[n]
                lo, hi = sorted((factor*elo, factor*ehi))
                need(intersects(boxes[n], lo, hi), 'analytic compact exponential derivative '+str((point, sign, n)))
                # The independent rational exponential bracket is much
                # narrower than the fixed producer grid.
                need(contains(boxes[n], (lo+hi)/2), 'strict analytic midpoint in packet box')
                count += 1
    for box in ((-8, -2), (-1, 1), (2, 8), (1, 1), (2, 2)):
        need(all(v.lo == v.hi == 0 for v in a.qodd(a.I(*box), 3, data)), 'flat or compact zero '+str(box))
    original = a.qodd
    # Manufacture Qo=sign(s)*s^4 and Pe=s^5 away from zero.  Substitute
    # only the leaf derivative oracle and compare the uncollected Pe-rQo
    # formula with every producer mixed derivative.
    def polynomial_q(box, order, ignored):
        need(box.lo == box.hi and box.lo != 0, 'manufactured point branch')
        point = box.lo
        sg = 1 if point > 0 else -1
        return [a.I(sg*F(factorial(4), factorial(4-n))*point**(4-n)) for n in range(order+1)]
    try:
        a.qodd = polynomial_q
        for xx, rr in ((F(1, 2), F(5, 4)), (F(3, 2), F(1, 4))):
            x = D({(0, 0): xx, (1, 0): 1})
            r = D({(0, 0): rr, (0, 1): 1})
            plus, minus = r+x, r-x
            sg = 1 if rr > xx else -1
            w = (plus**4+sg*minus**4)/2
            p = (6*a.ELL/130)/4*(plus**5+minus**5-r*(plus**4+sg*minus**4))
            wb, pb = a.comparisons(a.I(xx), a.I(rr), data)
            for key in a.KEYS:
                need(contains(wb.derivative(*key), w.d[key]), 'uncollected comparator master '+str(key))
                need(contains(pb.derivative(*key), p.d[key]), 'uncollected Pe-rQo mixed comparator '+str(key))
                count += 2
    finally:
        a.qodd = original
    return count


def check_temporal(a):
    count = 0
    for point in (F(0), F(1, 4), F(3)):
        t = a.ELL*point
        lo, hi = exp_negative(t, 80)
        box = a.core.expminus(t)
        need(contains(box, (lo+hi)/2), 'independent exponential direct Taylor bracket')
        c, h, v = a.temporal(a.I(point))
        # Value identities and differentiated ODEs ensure correct
        # exp/cosh/tanh/sqrt extensions without differentiating a box.
        for expression, constant in ((c*c-h*h*c*c, 1), (v*v*c*c*c, 1)):
            for key in a.KEYS:
                need(contains(expression.derivative(*key), constant if key == (0, 0) else 0), 'temporal differentiated identity')
                count += 1
        need(c.val().lo > 0 and v.val().lo > 0, 'positive square-root branch')
    return count


def check_background(a, data):
    count = 0
    def polynomial(coefficients, value):
        out = D.cv(0)
        for coefficient in reversed(coefficients):
            out = out*value+coefficient
        return out
    for point in (F(0), F(3, 2), F(8)):
        rho, eta, p, aa, hh, R, B, kk, U2, E = a.core.background(a.I(point), data, 'none')
        kg = a.I(*data['coefficients']['alpha_global']['log_derivative']['exact'])
        abox = a.core.expbox(-kg*a.ELL*point)
        state = {name: D({(0, 0): (box.lo+box.hi)/2}) for name, box in
                 (('rho', rho.val()), ('eta', eta.val()), ('p', p.val()), ('a', aa.val()), ('alpha', abox))}
        # Ordinary ODE recursion f^(n+1)=d^n(rhs), separate from the
        # producer's normalized-coefficient integration recurrence.
        for n in range(4):
            r, e, pp, av, alpha = (state[k] for k in ('rho', 'eta', 'p', 'a', 'alpha'))
            h = polynomial(a.core.UP[1], e)/pp
            rhs = {'rho': -r*r*av, 'eta': -r*pp,
                   'p': -r*(polynomial(a.core.UP[1], e)-4*av*pp),
                   'a': 1/r+r*pp*pp/3,
                   'alpha': r*(F(5, 2)*av-h)*alpha}
            for name in state:
                state[name].d[0, n+1] = rhs[name].d[0, n]
        r, e, pp, av, alpha = (state[k] for k in ('rho', 'eta', 'p', 'a', 'alpha'))
        h = polynomial(a.core.UP[1], e)/pp
        direct = (r, e, pp, av, h, alpha/r**2, 3*alpha/(r**3*pp), r*(F(5, 2)*av-h), polynomial(a.core.UP[2], e))
        for actual, expected in zip((rho, eta, p, aa, hh, R, B, kk, U2), direct):
            for n in range(5):
                need(contains(actual.c[n]*factorial(n), expected.d[0, n]), 'independent ordinary actual background derivative '+str(n))
                count += 1
        elo, ehi = exp_negative(a.ELL*point, 80)
        for n in range(5):
            need(contains(E.c[n]*factorial(n), (-1)**n*(elo+ehi)/2), 'independent background exponential derivative')
            count += 1
    return count


def shell_exterior_and_domains(a, data, errors):
    count = 0
    for xx in (a.I(0), a.I(F(5, 4), F(5, 4)+F(1, 4096)), a.I(0, 3)):
        m, s, bg = a.actual_jet(xx, a.I(0))
        need(all(m[key].lo == m[key].hi == 0 for key in ('B', 'Bt', 'Btt', 'C', 'Ct', 'Ctt')), 'homogeneous shell mixed rows through two tangential derivatives')
        count += 6
    # Selectivity test: removing the lift is legal on a whole cell only
    # when r.lo>=x.hi.  Instrument the branch, without changing estimates.
    calls = []
    original = a.exterior_caps
    def instrumented():
        calls.append(True)
        return original()
    try:
        a.exterior_caps = instrumented
        for xx, rr, expected in ((a.I(1, 2), a.I(2, 3), 1), (a.I(1, 2), a.I(F(3, 2), 3), 0), (a.I(1), a.I(1), 1)):
            calls.clear()
            a.actual_jet(xx, rr)
            need(len(calls) == expected, 'whole-cell exterior criterion')
            count += 1
    finally:
        a.exterior_caps = original
    for xx, rr in ((a.I(-1), a.I(0)), (a.I(4), a.I(0)), (a.I(1), a.I(-1)), (a.I(3), a.I(6)), (a.I(2, 3), a.I(5, 6))):
        try:
            a.actual_jet(xx, rr)
        except ValueError:
            count += 1
        else:
            raise ValueError('invalid whole-cell domain accepted')
    # Corner-clipped evaluation is valid on D intersection only.
    a.actual_jet(a.I(2, 3), a.I(5, 6), clip_causal=True)
    return count


def run(injection='none'):
    a = load_module(NUMERIC/'actual_first_jets.py', 'independently_reviewed_oracle')
    data, errors, freeze = check_pins(a)
    count = manufactured(a, injection)
    packet = check_bump_comparisons(a, data)
    temporal = check_temporal(a)
    background = check_background(a, data)
    shell = shell_exterior_and_domains(a, data, errors)
    return {'status': 'PASS_INDEPENDENT_ACTUAL_FIRST_JET_REVIEW',
            'exact_manufactured_ordinary_derivative_checks': count,
            'independent_packet_and_uncollected_comparator_checks': packet,
            'temporal_identity_checks': temporal,
            'independent_actual_background_ordinary_derivative_checks': background,
            'shell_exterior_domain_checks': shell,
            'oracle_freeze_sha256': hashlib.sha256((NUMERIC/'ORACLE_FREEZE.json').read_bytes()).hexdigest(),
            'oracle_source_sha256': freeze['files'],
            'review_source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'scope': 'Analytic composition review plus exact finite support checks; inherited smooth actual PDE solutions and derivative certificates remain authenticated conditional premises. No C0 error differentiation, second solution, or force claim.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--inject-error', default='none', choices=('none', 'wrong_factorial', 'wrong_scalar_shift', 'drop_radial_rho', 'drop_radial_conversion'))
    args = parser.parse_args()
    need(not args.output.exists(), 'refuse overwrite')
    report = run(args.inject_error)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(report['status'])


if __name__ == '__main__':
    main()
