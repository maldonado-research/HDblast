"""Exact rational audit of an independently derived inertial residual bound.

Analytic identities are proved in INDEPENDENT_INERTIAL_BOUND.md; arithmetic
checks do not formally verify PDE hypotheses. No archived program is run.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys

PINS = {
    'archive/BACKGROUND_DERIVATIVE_BOUNDS.json':
    'fd6fd277176033703a4e75a2f66ab94894318aac49c83264ef979e9a27f0a02c',
    'archive/inputs/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json':
    '8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8',
}


def require(ok, label):
    if not ok:
        raise RuntimeError(label)


def encode(x):
    if isinstance(x, F):
        return str(x)
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    return x


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, default=Path(__file__).resolve().parent.parent)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    checks = []

    def check(label, ok):
        require(ok, label)
        checks.append(label)

    data = {}
    for name, pin in PINS.items():
        raw = (args.source_root/name).read_bytes()
        check('authenticated_input:' + name, sha256(raw).hexdigest() == pin)
        data[name] = json.loads(raw)
    background = data['archive/BACKGROUND_DERIVATIVE_BOUNDS.json']
    shell = data['archive/inputs/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json']
    ell, T, Z = F(1, 10**6), F(3, 10**6), F(8, 10**6)
    check('same_actual_retained_background_and_new_state', background['candidate_state_exact'] == {
        'delta': '1/1000', 'c': '96714321425692/161839258930731',
        'ell': str(ell), 'Z': str(Z), 'T': str(T)})
    check('four_strict_wider_collar_first_exit_margins', len(background['first_exit']['component_margins']) == 4 and all(
        F(r['strict_lower_margin']) > 0 and F(r['strict_upper_margin']) > 0
        for r in background['first_exit']['component_margins'].values()))
    m, M, B, Gamma, K, a0 = F(1, 500), F(102000), F(250), F(250), F(3*10**7), F(35, 4)
    qlo, qhi = map(F, background['outer_derived']['Q']['exact'])
    check('potential_bound_from_actual_outer_rectangle', qlo > 0 and 129**2*qlo-F(35, 4) > 0 and 131**2*qhi < M)
    wzlo, wzhi = map(F, background['outer_derived']['W_z']['exact'])
    check('strict_spatial_derivative_bound', -K < wzlo <= wzhi < K)
    mulo, muhi = map(F, shell['boxes']['mu']['exact'])
    clo, chi = map(F, shell['boxes']['C_b']['exact'])
    glo, ghi = clo/2, chi/2
    check('positive_actual_boundary_inertia', m < mulo <= muhi)
    check('actual_gamma_bound', 0 < glo <= ghi < Gamma)
    check('constant_and_variable_beta_bounds', 0 < glo-F(25, 4)*muhi and ghi+F(5, 2)*muhi < B)

    U = 1+T
    r, s = M*T*T/2, B*T*T/(2*m)
    det = (1-r)*(1-s)-T*r/m
    check('positive_coupled_actual_field_majorant', r < 1 and s < 1 and det > 0)
    A = U*(1-s+T/m)/det
    V = U/(m*det)
    check('first_actual_field_majorant_equation', (1-r)*A-T*V == U)
    check('second_actual_field_majorant_equation', (1-s)*V-r*A/m == U/m)
    check('support_moment_integral', ((2*ell)**2-ell**2)/2 == F(3, 2)*ell**2)
    R0 = F(3, 8)*(K+a0*T)*ell**2+F(3, 4)*M*ell**2
    Itail = R0+M*V*T*T/4+M*M*A*T**3/6
    J1 = F(3, 4)*M*ell**2+T*R0+M*V*T**3/12+M*M*A*T**4/24
    check('history_integral_stricter_than_supremum_product', 0 < J1 < T*(M*ell+Itail))
    Ev = (J1/m+a0*V*T**4/4)/(1-s)
    check('constant_boundary_error_feedback_equation', (1-s)*Ev == J1/m+a0*V*T**4/4)
    global_error = ell*(M*ell+Itail+(1+Gamma*T)*Ev)
    tail_error = ell*(Itail+(1+Gamma*T)*Ev)
    check('global_error_below_1_022e_minus7', global_error < F('1.022e-7'))
    check('tail_error_below_1_645e_minus10', tail_error < F('1.645e-10'))

    # An independent coarse free-tail lower bound avoids relying on samples.
    # Positive-energy bounds give |k'|<=1/mu and |k|<=t/mu.
    # At r=0, P>=a^2-(Gamma+1)a; use mu<=1/400 and its monotonicity.
    mu_upper = F(1, 400)
    check('actual_mu_below_coarse_upper', muhi < mu_upper)
    check('free_kernel_initial_lower_monotonicity', (Gamma+1)*mu_upper < 2)
    P0 = 1/mu_upper**2-(Gamma+1)/mu_upper
    kpp = (1+B*T)/m**2
    kppp = (kpp+B/m)/m
    Pr = kppp+(Gamma+1)*kpp+Gamma/m
    Pmin = P0-T*Pr
    check('continuous_free_tail_kernel_lower_positive', Pmin > 0)
    # Q>=exp(-1/2)>1/2 on [4/3,5/3]; area>1/6.
    # exp(1/2)<sum (1/2)^n=2 proves the strict elementary inequality.
    free_lower = ell**2*Pmin/6
    actual_lower = free_lower-tail_error
    check('continuous_actual_normalized_tail_positive', actual_lower > F('9.64e-9'))

    out = {
        'status': 'PASS_INDEPENDENT_CONDITIONAL_INERTIAL_RESIDUAL_BOUND',
        'python_optimization': sys.flags.optimize,
        'check_count': len(checks), 'checks': checks, 'input_pins': PINS,
        'state': {'ell': ell, 'T': T, 'Z': Z},
        'constants': {'m': m, 'M': M, 'B': B, 'Gamma': Gamma, 'Kz': K},
        'majorant': {'r': r, 's': s, 'determinant': det, 'A': A, 'V': V},
        'residual': {'R0': R0, 'incoming_tail': Itail, 'history_integral': J1,
                     'boundary_error_velocity': Ev, 'global_H_error': global_error,
                     'tail_H_error': tail_error},
        'independent_continuous_tail_sign': {'free_kernel_lower': Pmin,
                    'bump_area_strict_lower': F(1, 6), 'free_H_tail_lower': free_lower,
                    'actual_H_tail_lower': actual_lower},
        'display_only': {'global_H_error': float(global_error), 'tail_H_error': float(tail_error),
                         'actual_H_tail_lower': float(actual_lower)},
        'own_source_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
        'scope': 'Independent exact rational arithmetic for the proved conditional first-variation estimator. No nonlinear existence or remainder is established.',
    }
    text = json.dumps(encode(out), indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(text)
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
