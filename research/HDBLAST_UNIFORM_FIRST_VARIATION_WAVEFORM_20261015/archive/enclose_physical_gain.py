#!/usr/bin/env python3
"""New exact-rational shell-gain audit; archived programs are input data only."""
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parent
PINS = {
    'OCT11_JACOBI_ENDPOINT_MONOTONICITY.md': 'f1e2571a551ec122dd7755a13145bae73bd92bd7f7426bdc3002e2f52f8dbfc8',
    'OCT11_RESCALED_BOUND_CERTIFICATE.json': '22c6cd4e505ee2c7e549420355b94a5b3c079ebb6945314f5b285c13bb67835d',
    'PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json': '8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8',
    'REFERENCE_ONLY_prior_derive_shell_box.py': '3a4d2a71847114b8dce9f910679a1f4ce14898212001b0355279ec633ebed110',
    'PRIOR_ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md': '0832a16c4bdf18cc188861a2c2860e7e6ff736dec39426d6ca7c589dc0434d5c',
    'PRIOR_SHELL_RESPONSE_ERROR_CERTIFICATE.md': '5f01be536c77123e9f667f4cd5514b52067641b2b389ba4e6eb655cb24654022',
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

    def __add__(self, o):
        o = self.box(o)
        return I(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, o):
        return self + -self.box(o)

    def __rsub__(self, o):
        return self.box(o) + -self

    def __mul__(self, o):
        o = self.box(o)
        values = [a * b for a in (self.lo, self.hi) for b in (o.lo, o.hi)]
        return I(min(values), max(values))

    __rmul__ = __mul__

    def recip(self):
        require(not self.lo <= 0 <= self.hi, 'Zero-containing reciprocal')
        return I(1 / self.hi, 1 / self.lo)

    def __truediv__(self, o):
        return self * self.box(o).recip()

    def __rtruediv__(self, o):
        return self.box(o) * self.recip()

    def square(self):
        return I(0 if self.lo <= 0 <= self.hi else min(self.lo**2, self.hi**2),
                 max(self.lo**2, self.hi**2))

    def record(self):
        den = 10**12
        low = F(self.lo.numerator * den // self.lo.denominator, den)
        high = F(-((-self.hi.numerator * den) // self.hi.denominator), den)
        require(low <= self.lo <= self.hi <= high, 'Outward decimal bracket failed')
        return {'exact': [str(self.lo), str(self.hi)],
                'outward_rational_at_12_decimal_places': [str(low), str(high)],
                'display_only': [float(self.lo), float(self.hi)]}


def padd(a, b):
    return [(a[k] if k < len(a) else F(0)) + (b[k] if k < len(b) else F(0))
            for k in range(max(len(a), len(b)))]


def pmul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def scale(a, x):
    return [v * x for v in a]


def derivative(a):
    return [k * a[k] for k in range(1, len(a))]


def peval(a, x):
    total = 0
    for v in reversed(a):
        total = total * x + v
    return total


def normalized(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def sqrt_interval(q, bits=128):
    require(q > 0, 'Positive square-root input required')
    den = 1 << bits
    n = isqrt((q.numerator << (2 * bits)) // q.denominator)
    lo, hi = F(n, den), F(n + 1, den)
    require(lo**2 <= q <= hi**2, 'Integer square-root enclosure failed')
    return I(lo, hi)


def atan_interval(x, terms=40):
    require(0 < x < 1 and terms % 2 == 0, 'Even alternating-series length required')
    total = sum((-1)**k * x**(2*k + 1) / (2*k + 1) for k in range(terms))
    return I(total, total + x**(2*terms + 1) / (2*terms + 1))


def pi_interval():
    a, b = atan_interval(F(1, 5)), atan_interval(F(1, 239))
    t2 = 2 * F(1, 5) / (1 - F(1, 5)**2)
    t4 = 2 * t2 / (1 - t2**2)
    tan_difference = (t4 - F(1, 239)) / (1 + t4 * F(1, 239))
    require(tan_difference == 1 and (4*a-b).lo > 0 and (4*a-b).hi < 1,
            'Machin identity/range check failed')
    # 0 < 4 atan(1/5)-atan(1/239) < 1 < pi/2 selects the principal
    # angle with tangent one. The alternating remainders are rigorous.
    out = 16*a - 4*b
    require(out.lo > 3 and out.hi < F(22, 7), 'Pi sanity bracket failed')
    return out


def main():
    for name, pin in PINS.items():
        require(sha256((ROOT/'inputs'/name).read_bytes()).hexdigest() == pin,
                'Reference input changed: ' + name)
    prior = json.loads((ROOT/'inputs/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json').read_text())
    cert = json.loads((ROOT/'inputs/OCT11_RESCALED_BOUND_CERTIFICATE.json').read_text())
    delta, c = F(1, 1000), F(96714321425692, 161839258930731)
    require(prior['parameters_exact'] == {'delta': str(delta), 'c': str(c)}, 'State parameter mismatch')
    require(prior['collar']['status'] == 'PASS_EXACT_FIRST_EXIT_TUBE_AND_COEFFICIENT_BOUNDS'
            and F(prior['collar']['inward_depth_exact']) == F(8, 10**10), 'Retained collar mismatch')
    W, f = [F(1, 3), 0, 1, F(1, 3)], [1+c, c]
    Wp, Wpp = derivative(W), derivative(derivative(W))
    sigma = padd(scale(W, 2), scale(f, delta))
    sigma1 = derivative(sigma)
    U = padd(scale(pmul(Wp, Wp), F(1, 2)), scale(pmul(W, W), F(-2, 3)))
    A = padd(scale(pmul(W, f), F(1, 9)), scale(Wp, -c/12))
    B = padd(scale(pmul(f, f), F(1, 36)), [-c*c/48])
    H2 = padd(scale(A, delta), scale(B, delta*delta))
    H2_direct = padd(padd(scale(pmul(sigma, sigma), F(1, 36)),
                          scale(pmul(sigma1, sigma1), F(-1, 48))), scale(U, F(1, 6)))
    ppoly = scale(sigma1, F(-1, 2))
    require(normalized(H2) == normalized(H2_direct), 'Exact shell-constraint cancellation failed')
    require(normalized(ppoly) == normalized(padd(scale(Wp, -1), [-delta*c/2])),
            'Scalar-junction polynomial failed')
    exact_constants = cert['primary_uniform_certificate']['exact_constants']
    require(normalized(A) == normalized([F(x) for x in exact_constants['A']])
            and normalized(B) == normalized([F(x) for x in exact_constants['B_fun']]),
            'Independently derived shell polynomials differ from retained certificate')
    pi = pi_interval()
    pifactor = 3 * pi.square()
    hp, pp = derivative(H2), derivative(ppoly)
    qg = padd(scale(pmul(hp, ppoly), F(3, 2)), scale(pmul(H2, pp), -1))
    qk = padd(pmul(hp, ppoly), scale(pmul(H2, pp), -2))
    eta_prior = I(*prior['boxes']['eta_b']['exact'])
    eta_center = -9*c*delta/64
    radius = F(exact_constants['K_eta']) * delta*delta
    eta_sharper = I(eta_center-radius, eta_center+radius)
    require(eta_prior.lo <= eta_sharper.lo < eta_sharper.hi <= eta_prior.hi,
            'Exact-certificate refinement must be inside retained shell box')

    def endpoints(eta):
        pbox, hbox, dh = peval(ppoly, eta), peval(H2, eta), peval(hp, eta)
        sign_g, sign_k = peval(qg, eta), peval(qk, eta)
        require(pbox.hi < 0 and hbox.lo > 0 and dh.hi < 0
                and sign_g.lo > 0 and sign_k.lo > 0,
                'Monotonicity/sign verification failed')
        def gain(x):
            h, p = peval(H2, x), peval(ppoly, x)
            return 3*h*sqrt_interval(h)/p
        gl, gu = gain(eta.lo), gain(eta.hi)
        # G is increasing, K decreasing: denominator p^3 changes the K sign.
        G = I(gl.lo, gu.hi)
        kl = pifactor * peval(H2, eta.hi) / peval(ppoly, eta.hi)**2
        ku = pifactor * peval(H2, eta.lo) / peval(ppoly, eta.lo)**2
        K = I(kl.lo, ku.hi)
        H2box = I(peval(H2, eta.hi), peval(H2, eta.lo))
        pexact = I(peval(ppoly, eta.hi), peval(ppoly, eta.lo))
        require(G.hi < 0 and K.lo > 0, 'Physical gain signs failed')
        return {'eta_interval': eta.record(), 'p_endpoint_interval': pexact.record(),
                'H_squared_endpoint_interval': H2box.record(), 'G': G.record(), 'K': K.record(),
                'derivative_sign_witnesses': {'H_squared_eta_derivative': dh.record(),
                    'G_positive_numerator': sign_g.record(), 'K_negative_derivative_numerator': sign_k.record()},
                'monotonicity': {'G': 'strictly increasing in eta', 'K': 'strictly decreasing in eta'}}

    rho_old, p_old = (I(*prior['boxes'][name]['exact']) for name in ('rho_b', 'p_b'))
    independent_G = 3 / (rho_old.square()*rho_old*p_old)
    independent_K = pifactor/(rho_old.square()*p_old.square())
    output = {
        'status': 'PASS_EXACT_RATIONAL_PHYSICAL_GAIN_ENCLOSURE',
        'state_exact': {'delta': str(delta), 'c': str(c), 'ell': '1/10000000000',
                        'observation_T1': '3/10000000000', 'inward_collar': '1/1250000000'},
        'reference_input_pins': PINS,
        'definitions': {'G': '3/(rho_b^3 p_b)', 'K': '3 pi^2/(rho_b^2 p_b^2)',
                        'S(t)': 'G cosh(t)^(-3/2) [w_z(t,0)-(C_b/2)w(t,0)]',
                        'leading_C': 'K integral_0^T1 [w_z(t,0)-(C_b/2)w(t,0)]^2 dt'},
        'independent_old_shell_box': {'G': independent_G.record(), 'K': independent_K.record()},
        'correlated_retained_eta_box': endpoints(eta_prior),
        'correlated_exact_certificate_eta_box': endpoints(eta_sharper),
        'pi_interval': pi.record(),
        'pi_method': 'Machin identity with exact tangent/range check and 40-term alternating arctangent enclosures',
        'sqrt_method': '128-bit dyadic bounds verified by squaring exact rational endpoints',
        'polynomial_identities_verified': ['p=-sigma_eta/2=-W_eta-delta*c/2',
            'rho_b^-2=sigma^2/36-sigma_eta^2/48+U/6=delta*A+delta^2*B',
            'independently derived A and B equal retained exact certificate coefficients'],
        'scope': 'Conditional consequences of retained analytic branch bounds. No trajectory, wave PDE, effective epsilon or cosmological scale is validated here.',
        'archived_programs_executed': False,
    }
    dest = ROOT/'PHYSICAL_GAIN_ENCLOSURE.json'
    dest.write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps({'status': output['status'],
        'correlated_retained_G': output['correlated_retained_eta_box']['G']['outward_rational_at_12_decimal_places'],
        'correlated_retained_K': output['correlated_retained_eta_box']['K']['outward_rational_at_12_decimal_places'],
        'sharper_G': output['correlated_exact_certificate_eta_box']['G']['outward_rational_at_12_decimal_places'],
        'sharper_K': output['correlated_exact_certificate_eta_box']['K']['outward_rational_at_12_decimal_places']}, indent=2))


if __name__ == '__main__':
    main()
