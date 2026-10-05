"""Exact algebra and conditional analytic bounds; no physical source evaluation."""
import argparse
from decimal import Decimal, localcontext, ROUND_CEILING
from fractions import Fraction as F
import hashlib
import json
from math import factorial
from pathlib import Path

NAMES = ('k', 'L', 'X', 'Z', 'T', 'h', 'h1', 'h2', 'eps', 'df',
         'rX', 'rZ', 'rT', 'K', 'q', 's', 'c', 'A')
ZERO = (0,) * len(NAMES)


class Poly:
    def __init__(self, terms=None):
        self.terms = {m: F(v) for m, v in (terms or {}).items() if v}

    @staticmethod
    def value(value):
        return value if isinstance(value, Poly) else Poly({ZERO: F(value)})

    def __add__(self, other):
        other = self.value(other)
        result = dict(self.terms)
        for m, value in other.terms.items():
            result[m] = result.get(m, F(0)) + value
        return Poly(result)

    __radd__ = __add__

    def __neg__(self):
        return Poly({m: -v for m, v in self.terms.items()})

    def __sub__(self, other):
        return self + -self.value(other)

    def __rsub__(self, other):
        return self.value(other) - self

    def __mul__(self, other):
        other = self.value(other)
        result = {}
        for a, av in self.terms.items():
            for b, bv in other.terms.items():
                m = tuple(x + y for x, y in zip(a, b))
                result[m] = result.get(m, F(0)) + av * bv
        return Poly(result)

    __rmul__ = __mul__

    def __pow__(self, n):
        result = self.value(1)
        for _ in range(n):
            result = result * self
        return result

    def partial(self, name):
        index = NAMES.index(name)
        result = {}
        for m, value in self.terms.items():
            if m[index]:
                new = list(m)
                new[index] -= 1
                result[tuple(new)] = value * m[index]
        return Poly(result)

    def derivative(self, directions):
        return sum((self.partial(name) * value for name, value in directions.items()), Poly())

    def replace(self, replacements):
        result = Poly()
        for m, value in self.terms.items():
            term = self.value(value)
            for name, exponent in zip(NAMES, m):
                term = term * self.value(replacements.get(name, variables[name])) ** exponent
            result = result + term
        return result


variables = {}
for index, name in enumerate(NAMES):
    m = list(ZERO)
    m[index] = 1
    variables[name] = Poly({tuple(m): F(1)})
globals().update(variables)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def upward(value):
    with localcontext() as context:
        context.prec = 8
        context.rounding = ROUND_CEILING
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def verify():
    checks = []
    def check(label, value):
        require(value, label)
        checks.append({'name': label, 'pass': True})

    # Operators below are separately given by the action, never solved from a Ward equation.
    # Denominators 2*k*eps are cleared. All variables are formal real symbols.
    g = 4 * L**2 * h - 2 * L * h1 - h2
    R_mode = (2 * k**2 + 3 * L**2) * X - k * T - L * Z
    P_mode = (F(2, 3) * k**2 - L**2) * X - k * T - L * Z
    R_contact = eps * (L * h1 - 4 * k**2 * h - 4 * L**2 * h)
    P_contact = eps * (L * h1 - F(4, 3) * k**2 * h)
    baseline_work = eps * h1 * (4 * k**2 + 3 * L**2)
    directions = {'L': L**2, 'h': h1, 'h1': h2,
                  'X': Z + rX, 'Z': -2 * k * T - eps * (g + df) + rZ,
                  'T': 2 * k * Z + rT}
    full_R, full_P = R_mode + R_contact, P_mode + P_contact
    actual = full_R.derivative(directions) - L * (full_R - 3 * full_P) + baseline_work
    residual = (2 * k**2 + 3 * L**2) * rX - k * rT - L * rZ + eps * L * df
    check('full_action_residual_identity', not (actual - residual).terms)
    check('on_shell_matched_source_ward_cancellation', not actual.replace({'rX': 0, 'rZ': 0, 'rT': 0, 'df': 0}).terms)
    check('forcing_mismatch_is_explicit', not (actual.replace({'rX': 0, 'rZ': 0, 'rT': 0}) - eps * L * df).terms)
    G = 3 * L**2 * X - L * Z
    canonical_derivative = (R_mode - G).derivative(directions)
    check('canonical_drift_from_two_residuals', not (canonical_derivative - 2 * k**2 * rX + k * rT).terms)
    check('negative_control_G_only_drops_nonzero_canonical_defect', bool(canonical_derivative.terms))
    wrong_pressure = full_R.derivative(directions) - L * (full_R - 3 * P_mode) + baseline_work
    check('negative_control_missing_pressure_contact', bool((wrong_pressure - residual).terms))
    check('negative_control_missing_baseline_work', bool((actual - baseline_work - residual).terms))
    check('negative_control_wrong_forcing_sign', bool((actual - residual + 2 * eps * L * df).terms))

    # For a real unit forcing impulse, 2*k times the raw (X,Z,T) response is
    # (-eps*sin, -2*k*eps*cos, -2*k*eps*sin). Verify the direct operator kernels.
    impulse = {'X': -eps*s, 'Z': -2*k*eps*c, 'T': -2*k*eps*s}
    check('density_forcing_kernel_from_direct_operator', not
          (R_mode.replace(impulse) - eps*(2*k*L*c - 3*L**2*s)).terms)
    check('pressure_forcing_kernel_from_direct_operator', not
          (P_mode.replace(impulse) - eps*(F(4, 3)*k**2*s + L**2*s + 2*k*L*c)).terms)

    # A smooth compactly supported spectral amplitude between probes is a
    # first-order Bogoliubov state perturbation, not a physical source sample.
    state = {'X': A*c, 'Z': -2*k*A*s, 'T': 2*k*A*c}
    phase_derivative = {'L': L**2, 's': 2*k*c, 'c': -2*k*s}
    check('hidden_state_first_mode_equation', not (state['X'].derivative(phase_derivative)-state['Z']).terms)
    check('hidden_state_real_second_mode_equation', not (state['Z'].derivative(phase_derivative)+2*k*state['T']).terms)
    check('hidden_state_imaginary_second_mode_equation', not (state['T'].derivative(phase_derivative)-2*k*state['Z']).terms)
    check('hidden_state_preserves_linear_wronskian', not (2*k*state['X']-state['T']).terms)
    hidden_R, hidden_P = R_mode.replace(state), P_mode.replace(state)
    check('hidden_state_obeys_direct_ward_identity', not
          (hidden_R.derivative(phase_derivative)-L*(hidden_R-3*hidden_P)).terms)
    check('hidden_state_has_nonzero_positive_initial_density', not
          (hidden_R.replace({'s': 0, 'c': 1})-3*L**2*A).terms and bool((3*L**2*A).terms))

    # q = 2*tau, s=sin(q*K), c=cos(q*K). Prove the closed primitives by differentiating.
    dK = {'K': Poly.value(1), 's': q * c, 'c': -q * s}
    C_num = q * K * s + c - 1
    S0_num = 1 - c
    S2_num = -q**2 * K**2 * c + 2 * q * K * s + 2 * c - 2
    for name, numerator, denominator, integrand in (
            ('k_cos', C_num, q**2, K * c),
            ('sin', S0_num, q, s),
            ('k_squared_sin', S2_num, q**3, K**2 * s)):
        check(name + '_exact_derivative', not (numerator.derivative(dK) - denominator * integrand).terms)
        check(name + '_zero_lower_limit', not numerator.replace({'K': 0, 's': 0, 'c': 1}).terms)
    check('negative_control_sine_squared_kernel_sign', bool(((-S2_num).derivative(dK) - q**3 * K**2 * s).terms))

    # Exact Taylor coefficients inspect removable zero-lag poles without evaluating any sine/source.
    series = []
    for power in range(16):
        cosine = F((-1)**(power//2), factorial(power)) if power % 2 == 0 else F(0)
        sine = F((-1)**((power-1)//2), factorial(power)) if power % 2 else F(0)
        series.append((sine, cosine))
    def coefficient(n, index):
        return series[n][index] if 0 <= n < len(series) else F(0)
    for n in range(10):
        c1 = coefficient(n+1, 0) + coefficient(n+2, 1)
        s0 = -coefficient(n+1, 1)
        s2 = -coefficient(n+1, 1) + 2*coefficient(n+2, 0) + 2*coefficient(n+3, 1)
        expected_c1 = coefficient(n, 1) / (n+2)
        expected_s0 = coefficient(n, 0) / (n+1)
        expected_s2 = coefficient(n, 0) / (n+3)
        check('entire_kernel_series_coefficient_'+str(n), (c1, s0, s2) == (expected_c1, expected_s0, expected_s2))
    check('zero_lag_extensions', (F(1, 2), F(0), F(0)) ==
          (coefficient(1, 0)+coefficient(2, 1), -coefficient(1, 1),
           -coefficient(1, 1)+2*coefficient(2, 0)+2*coefficient(3, 1)))

    M, radius, halfwidth, degree, panels = F(64), F(1, 8), F(1, 128), 24, 64
    ratio = halfwidth / radius
    tail = M * ratio**(degree+1) / (1-ratio)
    length = panels * 2 * halfwidth
    alpha = length * tail / (degree+2)
    centers = [F(-9, 2) + F(2*j+1, 128) for j in range(panels)]
    check('registered_cauchy_uniform_tail', tail == F(1, 15*2**90))
    check('integrated_absolute_tail', alpha == F(1, 390*2**90))
    check('unit_interval_and_reflection_paired_centers', length == 1 and all(centers[j]+centers[-j-1] == -8 for j in range(panels)))
    check('represented_pi_satisfies_lower_bound', F(14488038916154245685, 4611686018427387904) >= 3)
    ell, pi_lower = F(2, 7), F(3)
    table = []
    for cutoff in (64, 128, 256):
        # Only exact analytic truncation; no coefficient, contact or trajectory error is covered.
        mass = F(cutoff**2) / (8*pi_lower**2)
        density_kernel = alpha * (ell*cutoff**2 + 3*ell**2*cutoff)/(8*pi_lower**2)
        density_moment = alpha * mass*(ell+F(3, 2)*ell**2)
        density = min(density_kernel, density_moment)
        pressure = alpha*(F(cutoff**3, 18)/pi_lower**2 + ell*cutoff**2/(8*pi_lower**2) + ell**2*cutoff/(8*pi_lower**2))
        source_work = alpha*mass*ell
        row = {'K': cutoff, 'density_truncation_bound': str(density),
               'density_kernel_bound': str(density_kernel), 'density_first_moment_bound': str(density_moment),
               'pressure_truncation_bound': str(pressure), 'integrated_source_mismatch_bound': str(source_work),
               'upward_decimal_bounds': {'density': upward(density), 'pressure': upward(pressure), 'source_mismatch': upward(source_work)}}
        table.append(row)
    check('strict_K256_density_bound', F(table[-1]['density_truncation_bound']) < F(55, 10**29))
    check('strict_K256_pressure_bound', F(table[-1]['pressure_truncation_bound']) < F(22, 10**26))
    check('strict_K256_source_mismatch_bound', F(table[-1]['integrated_source_mismatch_bound']) < F(55, 10**29))
    return {'schema_version': 1, 'status': 'PASS_EXACT_WARD_AND_FINITE_BAND_TRUNCATION_COMPONENT',
            'checks': checks, 'checks_count': len(checks), 'tail_L1': str(alpha),
            'tail_first_time_moment': str(alpha/2), 'tail_uniform': str(tail),
            'bounds': table, 'conditions': {'source_degree': degree, 'panels': panels,
            'interval': ['-9/2', '-7/2'], 'Pi_lower': '3', 'L_upper': '2/7',
            'contact': 'Exact separately defined action/subtraction contacts; otherwise add their residual',
            'source_polynomial': 'Exact degree24 Taylor coefficients; rounding/coefficients require separate bounds'},
            'scope': 'Conditional analytic source-truncation component, exact formal Ward identities and exact finite-band kernel algebra only.',
            'full_twelve_case_pressure_contact_certificate': 'UNRESOLVED',
            'external_mathematical_novelty': 'NOT_ASSESSED', 'real_source_evaluations': 0,
            'inherited_quantum_arrays_decoded': 0, 'numerical_trajectory_cases_run': 0}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = verify()
    raw = (json.dumps(result, indent=2, sort_keys=True) + '\n').encode()
    if args.output:
        with args.output.open('xb') as handle:
            handle.write(raw)
    print(json.dumps({'status': result['status'], 'checks': result['checks_count'],
                      'result_sha256': hashlib.sha256(raw).hexdigest(), 'real_source_evaluations': 0}))
