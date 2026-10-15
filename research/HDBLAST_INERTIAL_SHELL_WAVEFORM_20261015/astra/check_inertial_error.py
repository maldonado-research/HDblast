#!/usr/bin/env python3
"""Independent exact-rational inertial-comparison error and signed-tail audit.

Analytic identities and the retained PDE/ODE theorem are mathematical inputs.
No archived producer is imported or executed. Decisions survive Python -O.
"""
import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

PINS = {
    'inputs/previous/archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json':
    '8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8',
    'inputs/previous/archive/ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md':
    '0832a16c4bdf18cc188861a2c2860e7e6ff736dec39426d6ca7c589dc0434d5c',
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    def __init__(self, lo, hi=None):
        object.__setattr__(self, 'lo', F(lo))
        object.__setattr__(self, 'hi', F(lo if hi is None else hi))
        require(self.lo <= self.hi, 'reversed interval')

    def __add__(self, other):
        other = interval(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)
    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -interval(other)

    def __rsub__(self, other):
        return interval(other) + -self

    def __mul__(self, other):
        other = interval(other)
        endpoints = [a*b for a in (self.lo, self.hi)
                     for b in (other.lo, other.hi)]
        return Interval(min(endpoints), max(endpoints))
    __rmul__ = __mul__

    def __truediv__(self, other):
        other = interval(other)
        require(not other.lo <= 0 <= other.hi, 'denominator contains zero')
        return self * Interval(1/other.hi, 1/other.lo)

    def __rtruediv__(self, other):
        return interval(other) / self

    def __pow__(self, exponent):
        require(isinstance(exponent, int) and exponent >= 0,
                'invalid nonnegative integer power')
        value = Interval(1)
        for unused in range(exponent):
            value = value * self
        return value

    def magnitude(self):
        return max(abs(self.lo), abs(self.hi))


def interval(value):
    return value if isinstance(value, Interval) else Interval(value)


def derivative_polynomials(eta):
    # Independent Horner evaluation of the exactly expanded potential.
    up = eta*(F(28,9)+eta*(F(50,9)+eta*(F(-2,3)
         +eta*(F(-20,9)-F(4,9)*eta))))
    upp = F(28,9)+eta*(F(100,9)+eta*(-2+eta*(F(-80,9)-F(20,9)*eta)))
    uppp = F(100,9)+eta*(-4+eta*(F(-80,3)-F(80,9)*eta))
    return up, upp, uppp


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, Interval):
        return {'exact': [str(value.lo), str(value.hi)],
                'display_only': [float(value.lo), float(value.hi)]}
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True,
                        help='Frozen prior capsule containing inputs/previous/archive')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    checks, rejected = [], []

    def check(name, condition):
        require(condition, name)
        checks.append(name)

    def reject(name, function):
        try:
            function()
        except ValueError:
            rejected.append(name)
        else:
            raise ValueError('negative control accepted: ' + name)

    for name, digest in PINS.items():
        check('input_pin:' + name,
              sha256((args.source_root/name).read_bytes()).hexdigest() == digest)
    source = json.loads((args.source_root/next(iter(PINS))).read_text())
    check('same_exact_background_parameters', source['parameters_exact'] ==
          {'delta': '1/1000', 'c': '96714321425692/161839258930731'})

    boxes = source['collar']['first_exit_inclusions']
    outer = {key: Interval(*item['outer_box']['exact'])
             for key, item in boxes.items()}
    initial = {key: Interval(*item['initial_box']['exact'])
               for key, item in boxes.items()}
    rho, eta, p, a = (outer[key] for key in ('rho', 'eta', 'p', 'a'))
    up, upp, uppp = derivative_polynomials(eta)
    rhs = {'rho': -rho**2*a, 'eta': -rho*p,
           'p': -rho*(up-4*a*p), 'a': 1/rho+rho*p**2/3}
    ell = F(1, 1000000)
    T, Z = 3*ell, 8*ell
    margins = {}
    for key in outer:
        margins[key] = [initial[key].lo + min(F(0), Z*rhs[key].lo)-outer[key].lo,
                        outer[key].hi-initial[key].hi-max(F(0), Z*rhs[key].hi)]
        check('strict_wider_first_exit:' + key, min(margins[key]) > 0)

    h = up/p
    hz = rho*(h**2-4*a*h-upp)
    potential_factor = 2*h**2-10*a*h+F(35,4)*a**2-upp-p**2/6
    factor_z = (4*h-10*a)*hz+(-10*h+F(35,2)*a)*rhs['a'] \
               -uppp*rhs['eta']-p*rhs['p']/3
    spatial_potential = rho**2*potential_factor
    W = spatial_potential + Interval(F(-35,4), 0)
    Wz = 2*rho*rhs['rho']*potential_factor + rho**2*factor_z

    sb = source['boxes']
    shell_eta, shell_rho, shell_p, shell_a = (
        Interval(*sb[key]['exact']) for key in ('eta_b', 'rho_b', 'p_b', 'a_b'))
    shell_up, unused, unused2 = derivative_polynomials(shell_eta)
    shell_h = shell_up/shell_p
    mu = 1/(shell_rho*(shell_h-4*shell_a+2+2*shell_eta))
    gamma = shell_rho*(2*shell_h-5*shell_a)/2
    beta0 = gamma+F(5,2)*mu
    beta = gamma+mu*Interval(F(-25,4), F(5,2))

    m, mu_max, M, B, Gamma, Kz = F(1,500), F(1,400), F(102000), F(250), F(250), F(30000000)
    check('positive_boundary_inertia', m < mu.lo <= mu.hi < mu_max)
    check('boundary_gamma_positive_and_bounded', 200 < gamma.lo <= gamma.hi < Gamma)
    check('positive_constant_beta0', 0 < beta0.lo <= beta0.hi < B)
    check('variable_beta_bound', beta.magnitude() < B)
    check('potential_bound', W.magnitude() < M)
    check('spatial_potential_derivative_bound', Wz.magnitude() < Kz)

    U = 1+T
    r, s = M*T*T/2, B*T*T/(2*m)
    determinant = (1-r)*(1-s)-T*r/m
    check('majorant_positive', r < 1 and s < 1 and determinant > 0)
    A = U*(1-s+T/m)/determinant
    V = U/(m*determinant)
    check('majorant_first_row', (1-r)*A-T*V == U)
    check('majorant_second_row', (1-s)*V-r*A/m == U/m)

    # Energy supplies |G'|<=1/mu and |G|<=t/mu. The exact factor mu
    # cancels in the beta-variation contribution before using m.
    memory_error = (1+Gamma*T)*M*A*T*T/(2*m)
    beta_variation_error = F(35,16)*(1+Gamma*T)*V*T**4
    global_error = ell*(M*T*A+memory_error+beta_variation_error)

    # Safe paired-q bound; the unused Q reflection symmetry can tighten it.
    paired_q = (Kz+F(35,2)*T)*ell**2/2
    velocity_part = M*ell**2
    boundary_lift = M*V*(T-ell)**2/4
    bulk_remainder = M*M*A*T**3/6
    I_tail = paired_q+velocity_part+boundary_lift+bulk_remainder
    tail_error = ell*(I_tail+memory_error+beta_variation_error)

    # Independent uniform sign proof for the inertial comparator tail.
    # L(0,s)=1/mu^2-(gamma+exp(-s))/mu >=59600 by monotonicity in mu.
    leading_lower = 1/mu_max**2-(Gamma+1)/mu_max
    check('leading_kernel_mu_monotonicity', (Gamma+1)*mu_max < 2)
    Gpp_bound = (1+B*T)/m**2
    Gppp_bound = (Gpp_bound+B/m)/m
    L_derivative_bound = Gppp_bound+(Gamma+1)*Gpp_bound+Gamma/m
    tail_kernel_lower = leading_lower-T*L_derivative_bound
    check('positive_tail_kernel', tail_kernel_lower > 0)
    # On [4/3,5/3], Q>=exp(-1/2)>1/2, so integral Q>1/6.
    bump_area_lower = F(1,6)
    comparator_tail_lower = ell**2*bump_area_lower*tail_kernel_lower
    actual_tail_lower = comparator_tail_lower-tail_error

    check('global_error_rounded', global_error < F('0.0000003067'))
    check('tail_error_rounded', tail_error < F('0.000000000297'))
    check('comparator_tail_lower_rounded', comparator_tail_lower > F('0.0000000098'))
    check('actual_tail_uniform_positive', actual_tail_lower > F('0.0000000095'))
    check('new_width_no_narrow_response_transfer', ell != F(1,10**10))

    reject('zero_denominator_is_rejected', lambda: 1/Interval(-1,1))
    reject('missing_inertial_feedback_fails_majorant', lambda: require(
        (1-r)*(U*(1-s+T/m)/((1-r)*(1-s)))
        -T*(U/(m*((1-r)*(1-s)))) >= U,
        'dropping determinant feedback fails first majorant inequality'))
    reject('tail_sign_not_from_global_error', lambda: require(
        comparator_tail_lower-global_error > 0,
        'global potential bound is too large to certify this tail sign'))
    reject('false_zero_tail_rejected', lambda: require(actual_tail_lower <= 0,
                                                       'positive lower bound excludes zero'))
    reject('insufficient_Kz_rejected', lambda: require(Wz.magnitude() < 1,
                                                    'derivative bound is not one'))

    receipt = {
        'status': 'PASS_CONDITIONAL_INERTIAL_COMPARATOR_AND_POSITIVE_TAIL',
        'optimized_mode': bool(sys.flags.optimize),
        'input_pins': PINS, 'checks': checks, 'negative_controls_rejected': rejected,
        'state': {'ell': ell, 'T': T, 'Z': Z},
        'background_first_exit_margins': margins,
        'independent_coefficient_intervals': {'mu': mu, 'gamma': gamma,
          'beta0': beta0, 'beta_all_real_t': beta, 'W_all_real_t': W, 'W_z': Wz,
          'U_etaetaeta': uppp, 'h_z': hz},
        'rounded_constants': {'m': m, 'mu_upper': mu_max, 'M': M,
                              'B': B, 'Gamma': Gamma, 'K_z': Kz},
        'majorant': {'r': r, 's': s, 'determinant': determinant, 'A': A, 'V': V},
        'error': {'global': global_error, 'tail': tail_error,
          'incoming_tail': I_tail, 'paired_initial_q': paired_q,
          'initial_velocity': velocity_part, 'boundary_lift': boundary_lift,
          'bulk_Duhamel': bulk_remainder, 'history_convolution': memory_error,
          'beta_variation': beta_variation_error},
        'independent_tail_sign': {'leading_kernel_lower': leading_lower,
          'kernel_time_derivative_bound': L_derivative_bound,
          'kernel_lower': tail_kernel_lower, 'bump_area_lower': bump_area_lower,
          'comparator_H_lower': comparator_tail_lower, 'actual_H_lower': actual_tail_lower},
        'display_only': {'global_error': float(global_error),
          'tail_error': float(tail_error), 'comparator_tail_lower': float(comparator_tail_lower),
          'actual_tail_lower': float(actual_tail_lower)},
        'limitations': ['Conditional first variation only; no effective nonlinear remainder',
          'Retained analytic branch and coupled PDE identities are mathematical inputs',
          'Exact-family inverse, gauge, support, reflection and evolution hypotheses remain unproved inputs',
          'Finite rational checks do not formally verify the analytic proof',
          'No response coefficient transferred from the narrow state'],
    }
    output = json.dumps(encode(receipt), indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(output)
    else:
        print(output, end='')


if __name__ == '__main__':
    main()
