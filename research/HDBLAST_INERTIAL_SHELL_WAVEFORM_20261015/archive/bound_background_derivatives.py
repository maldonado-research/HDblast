#!/usr/bin/env python3
"""Fresh rational interval bounds for the wider candidate's actual background.

No archived producer is imported or executed. This standard-library program
reads authenticated theorem/certificate data and the retained shell box.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
COMMIT = 'c490f17139f1266426995595f5ba0eaefb14c636'
PINS = {
    'OCT11_JACOBI_ENDPOINT_MONOTONICITY.md': 'f1e2571a551ec122dd7755a13145bae73bd92bd7f7426bdc3002e2f52f8dbfc8',
    'OCT11_RESCALED_BOUND_CERTIFICATE.json': '22c6cd4e505ee2c7e549420355b94a5b3c079ebb6945314f5b285c13bb67835d',
    'PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json': '8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8',
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


class I:
    def __init__(self, lo, hi=None):
        self.lo, self.hi = F(lo), F(lo if hi is None else hi)
        require(self.lo <= self.hi, 'Reversed interval')

    @staticmethod
    def box(x):
        return x if isinstance(x, I) else I(x)

    def __add__(self, other):
        other = self.box(other)
        return I(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -self.box(other)

    def __rsub__(self, other):
        return self.box(other) + -self

    def __mul__(self, other):
        other = self.box(other)
        products = [a*b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return I(min(products), max(products))

    __rmul__ = __mul__

    def recip(self):
        require(not self.lo <= 0 <= self.hi, 'Reciprocal of a zero-containing interval')
        return I(1/self.hi, 1/self.lo)

    def __truediv__(self, other):
        return self*self.box(other).recip()

    def __rtruediv__(self, other):
        return self.box(other)*self.recip()

    def square(self):
        return I(0 if self.lo <= 0 <= self.hi else min(self.lo**2, self.hi**2),
                 max(self.lo**2, self.hi**2))

    def hull(self, other):
        return I(min(self.lo, other.lo), max(self.hi, other.hi))

    def mag(self):
        return max(abs(self.lo), abs(self.hi))

    def record(self):
        den = 10**12
        outward_lo = F(self.lo.numerator*den//self.lo.denominator, den)
        outward_hi = F(-((-self.hi.numerator*den)//self.hi.denominator), den)
        require(outward_lo <= self.lo <= self.hi <= outward_hi, 'Outward display conversion failed')
        return {'exact': [str(self.lo), str(self.hi)],
                'outward_rational_12_places': [str(outward_lo), str(outward_hi)],
                'display_only': [float(self.lo), float(self.hi)]}


def add(a, b):
    return [(a[k] if k < len(a) else F(0))+(b[k] if k < len(b) else F(0))
            for k in range(max(len(a), len(b)))]


def mul(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def scale(a, x):
    return [v*x for v in a]


def derivative(a):
    return [k*a[k] for k in range(1, len(a))]


def evaluate(a, x):
    total = 0
    for v in reversed(a):
        total = total*x+v
    return total


def trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def sqrt_box(q, bits=128):
    require(q.lo > 0, 'Positive square-root interval required')
    den = 1 << bits
    nlo = isqrt((q.lo.numerator << (2*bits))//q.lo.denominator)
    nhi = isqrt((q.hi.numerator << (2*bits))//q.hi.denominator)+1
    lo, hi = F(nlo, den), F(nhi, den)
    require(lo*lo <= q.lo <= q.hi <= hi*hi, 'Exact square-root bracket failed')
    return I(lo, hi)


def main():
    for name, pin in PINS.items():
        require(sha256((ROOT/'inputs'/name).read_bytes()).hexdigest() == pin, 'Input pin differs: '+name)
    old = json.loads((ROOT/'inputs/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json').read_text())
    cert = json.loads((ROOT/'inputs/OCT11_RESCALED_BOUND_CERTIFICATE.json').read_text())
    c, delta = F(96714321425692, 161839258930731), F(1, 1000)
    ell, Z, T = F(1, 10**6), F(8, 10**6), F(3, 10**6)
    require(old['parameters_exact'] == {'c': str(c), 'delta': str(delta)}, 'Actual branch parameter mismatch')
    W = [F(1, 3), 0, 1, F(1, 3)]
    W1, W2 = derivative(W), derivative(derivative(W))
    U = add(scale(mul(W1, W1), F(1, 2)), scale(mul(W, W), F(-2, 3)))
    U1, U2 = derivative(U), derivative(derivative(U))
    U3 = derivative(U2)
    f = [1+c, c]
    sigma = add(scale(W, 2), scale(f, delta))
    ppoly = scale(derivative(sigma), F(-1, 2))
    apoly = scale(sigma, F(1, 6))
    A = add(scale(mul(W, f), F(1, 9)), scale(W1, -c/12))
    Bfun = add(scale(mul(f, f), F(1, 36)), [-c*c/48])
    H2 = add(scale(A, delta), scale(Bfun, delta*delta))
    direct = add(add(scale(mul(sigma, sigma), F(1, 36)),
                     scale(mul(derivative(sigma), derivative(sigma)), F(-1, 48))), scale(U, F(1, 6)))
    require(trim(direct) == trim(H2), 'Shell Hamiltonian polynomial identity failed')
    constants = cert['primary_uniform_certificate']['exact_constants']
    require(trim(A) == trim([F(x) for x in constants['A']])
            and trim(Bfun) == trim([F(x) for x in constants['B_fun']]), 'Retained certificate shell polynomial differs')
    outer = {'rho': I(129, 131), 'eta': I(F(-9, 10**5), F(-8, 10**5)),
             'p': I(F(-14, 10**5), F(-12, 10**5)), 'a': I(F(111, 1000), F(112, 1000))}
    initial = {key: I(*old['boxes'][name]['exact']) for key, name in
               [('rho', 'rho_b'), ('eta', 'eta_b'), ('p', 'p_b'), ('a', 'a_b')]}
    rho, eta, p, a = (outer[k] for k in ('rho', 'eta', 'p', 'a'))
    u1, u2, u3 = (evaluate(poly, eta) for poly in (U1, U2, U3))
    rhs = {'rho': -rho.square()*a, 'eta': -rho*p,
           'p': -rho*(u1-4*a*p), 'a': 1/rho+rho*p.square()/3}
    margins = {}
    for name, box in outer.items():
        tube = initial[name]+I(0, Z)*rhs[name]
        require(box.lo < tube.lo <= tube.hi < box.hi, 'Wider collar first-exit inclusion failed: '+name)
        margins[name] = {'outer': box.record(), 'rhs': rhs[name].record(), 'flow_integral': tube.record(),
                         'strict_lower_margin': str(tube.lo-box.lo), 'strict_upper_margin': str(box.hi-tube.hi)}
    require(131*Z < 9, 'The wider inward collar could reach the regular endpoint')
    h = u1/p
    hz = -rho*(u2-h.square()+4*a*h)
    Q = 2*h.square()-10*a*h+F(35, 4)*a.square()-u2-p.square()/6
    Qz = (4*h-10*a)*hz+(-10*h+F(35, 2)*a)*rhs['a']+rho*p*u3-p*rhs['p']/3
    Wz = -2*rho.square()*rho*a*Q+rho.square()*Qz
    Mz, Mt = F(3*10**7), F(35, 2)*T
    require(Wz.mag() < Mz, 'Requested strict W_z bound did not pass')
    temporal_variation = F(35, 4)*T*T
    total_variation = Mz*Z+temporal_variation
    eta_b = initial['eta']
    # Rational scalar subdivisions preserve shell correlation. They enclose
    # every eta in the retained analytic interval, not sampled trajectories.
    shell = {}
    panels = 128
    for i in range(panels):
        x = I(eta_b.lo+(eta_b.hi-eta_b.lo)*i/panels,
              eta_b.lo+(eta_b.hi-eta_b.lo)*(i+1)/panels)
        hp, pp, aa, ww2 = evaluate(H2, x), evaluate(ppoly, x), evaluate(apoly, x), evaluate(W2, x)
        rr = sqrt_box(hp).recip()
        hh = evaluate(U1, x)/pp
        bb = hh-4*aa+ww2
        gamma = rr*(hh-F(5, 2)*aa)
        mu = 1/(rr*bb)
        beta0 = gamma+F(5, 2)*mu
        gap = rr*(ww2-F(3, 2)*aa)
        gap0 = gap-F(5, 2)*mu
        shell_Q = 2*hh.square()-10*aa*hh+F(35, 4)*aa.square()-evaluate(U2, x)-pp.square()/6
        values = {'rho_b': rr, 'p_b': pp, 'a_b': aa, 'h_b': hh, 'B_b': bb,
                  'mu': mu, 'gamma': gamma, 'beta0': beta0,
                  'inverse_mu_minus_gamma': gap, 'inverse_mu_minus_beta0': gap0,
                  'W_zero_zero': rr.square()*shell_Q}
        require(pp.hi < 0 and hp.lo > 0 and bb.lo > 0 and mu.lo > 0 and gamma.lo > 0
                and gap.lo > 0 and gap0.lo > 0, 'Correlated shell signs failed')
        for name, value in values.items():
            shell[name] = value if name not in shell else shell[name].hull(value)
    require(shell['mu'].lo > F(1, 1000) and shell['mu'].hi < F(1, 400), 'Retained simple inertia bounds failed')
    require(shell['inverse_mu_minus_gamma'].lo > 238
            and shell['inverse_mu_minus_gamma'].hi < 239
            and shell['inverse_mu_minus_beta0'].lo > 238
            and shell['inverse_mu_minus_beta0'].hi < 239, 'Positive-tail simple bounds failed')
    compact_shell = {
        'mu': I(F(2059814, 10**9), F(2271616, 10**9)),
        'gamma': I(F(202092285, 10**6), F(247360282, 10**6)),
        'beta0': I(F(202097960, 10**6), F(247365436, 10**6)),
        'inverse_mu_minus_gamma': I(F(2381203, 10**4), F(2381231, 10**4)),
        'inverse_mu_minus_beta0': I(F(2381151, 10**4), F(2381174, 10**4)),
        'W_zero_zero': I(F(283488885, 10**4), F(691450116, 10**4)),
    }
    require(all(box.lo < shell[name].lo <= shell[name].hi < box.hi for name, box in compact_shell.items()),
            'Compact strict shell enclosure failed')
    output = {
        'status': 'PASS_EXACT_RATIONAL_WIDER_COLLAR_AND_DERIVATIVE_BOUNDS',
        'immutable_input_main': COMMIT, 'input_pins': PINS,
        'candidate_state_exact': {'delta': str(delta), 'c': str(c), 'ell': str(ell), 'Z': str(Z), 'T': str(T)},
        'scope': 'Actual retained near-phi=+1 background, wider newly specified candidate. Conditional analytic bounds; no archived producer execution, numerical trajectory or wave response claim.',
        'first_exit': {'status': 'PASS', 'component_margins': margins, 'proper_radial_distance_upper': str(131*Z),
                       'retained_proper_shell_radius_lower': '9'},
        'outer_derived': {'h': h.record(), 'h_z': hz.record(), 'Q': Q.record(), 'Q_z': Qz.record(), 'W_z': Wz.record()},
        'strict_simple_bounds': {'absolute_W_z': str(Mz), 'conservative_absolute_W_z': '100000000',
                                 'absolute_W_t_on_window': str(Mt),
                                 'absolute_beta_t_on_window': str(Mt/400)},
        'variation_bounds': {
            'pointwise_formula': f'|W(t,z)-W(0,0)| <= {Mz}*z+(35/4)*t^2',
            'temporal_only_upper': str(temporal_variation),
            'whole_rectangle_upper': str(total_variation),
            'beta_temporal_only_upper': str(temporal_variation/400),
            'beta_identity': 'beta(t)=beta0-(35/4)*mu*tanh(t)^2'},
        'shell_correlation': {'scalar_panels': panels, 'same_retained_eta_interval': eta_b.record(),
                              'parameters': {name: value.record() for name, value in shell.items()},
                              'compact_strict_enclosures': {name: value.record() for name, value in compact_shell.items()}},
        'exact_identities': {
            'B_b': 'h_b-4*a_b+W_second(eta_b)',
            'gamma': 'rho_b*(h_b-5*a_b/2)',
            'inverse_mu': 'rho_b*B_b',
            'beta0': 'gamma+5*mu/2',
            'positive_tail_gap': '1/mu-gamma=rho_b*(W_second(eta_b)-3*a_b/2)',
            'beta0_tail_gap': '1/mu-beta0=rho_b*(W_second(eta_b)-3*a_b/2)-5*mu/2',
            'h_z': '-rho*(U_second-h^2+4*a*h)',
            'Q_z': '(4*h-10*a)*h_z+(-10*h+35*a/2)*a_z+rho*p*U_third-p*p_z/3',
            'W_z': '-2*rho^3*a*Q+rho^2*Q_z',
            'W_t': '-(35/2)*tanh(t)*sech(t)^2'},
    }
    (ROOT/'BACKGROUND_DERIVATIVE_BOUNDS.json').write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps({'status': output['status'], 'strict_simple_bounds': output['strict_simple_bounds'],
                      'W_z_display_only': Wz.record()['display_only'],
                      'shell_display_only': {name: value.record()['display_only'] for name, value in shell.items()},
                      'variation_bounds': output['variation_bounds']}, indent=2))


if __name__ == '__main__':
    main()
