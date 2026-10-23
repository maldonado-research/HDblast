#!/usr/bin/env python3
"""Standalone rational implications of explicitly assumed scientific premises.

This authenticates only its two small local input files. It does not
authenticate or replay any upstream PDE, source, initial-data or gauge proof.
"""
import argparse
import hashlib
import json
from pathlib import Path

PINS = {
    'LOCAL_TOTAL_REVIEW.md': '69a7a4741d599aad07d952a9b48423bf9b06f18c2b2431f3de8e821134db610c',
    'EXPLICIT_PREMISES.json': '436305e19f375ddead2958dd44ea48dd080302852ef056631efd3ae2e6b433c1',
}


def need(ok, label):
    if not ok:
        raise ValueError(label)


def authenticate(root):
    raw = {}
    for name, pin in PINS.items():
        path = root / name
        need(not path.is_symlink(), 'LOCAL_INPUT_SYMLINK_' + name)
        raw[name] = path.read_bytes()
        need(hashlib.sha256(raw[name]).hexdigest() == pin, 'LOCAL_INPUT_PIN_' + name)
    return raw


def run(root):
    before = authenticate(root)
    from fractions import Fraction as Q
    p = json.loads(before['EXPLICIT_PREMISES.json'])
    A, B, D, M = [Q(p['coefficient_norm_upper'][key]) for key in ('A', 'B', 'D', 'M')]
    K, L = [Q(p['kernel_norm_upper'][key]) for key in ('K', 'L')]
    F, G = [Q(p['source_norm_upper'][key]) for key in ('F', 'G')]
    cap, e0, e1, target = [Q(p[key]) for key in ('exact_total_X_upper', 'design_e0', 'design_e1', 'design_target')]
    q = K*(A+B+Q(3,4)*D)+L*Q(1,4)*M
    q_old = K*(A+B+3*D)+L*3*M
    forcing = K*F+L*G
    tail = q*q*cap
    budget = tail+(1+q)*e0+e1
    need(q == Q(3409,8000000), 'WEDGE_CONTRACTION')
    need(q_old == Q(3609,2000000) and 0<q<q_old<1, 'CONTRACTION_ORDER')
    need(forcing == Q(52162547172519,2), 'WHOLE_FORCING_ALLOWANCE')
    need(forcing < Q(2609)*10**10, 'FORCING_ROUNDING')
    need(tail < Q(1870300)*10**20, 'TAIL_ROUNDING')
    need(budget < target, 'DESIGN_BUDGET_ARITHMETIC')

    distinctions = {}
    t, free, source = Q(1,10), Q(1), Q(1,5)
    total, hom, W = (free+source)/(1-t), free/(1-t), source/(1-t)
    need(total == hom+W and total != free+W, 'TOTAL_VS_FREE')
    distinctions['homogeneous_is_not_free'] = [str(total), str(free+W)]
    b, B0, B1 = free+source, Q(5,4), Q(3,25)
    exact_error = total-B0-B1
    rhs = t*t*total+(1+t)*(b-B0)+(t*B0-B1)
    wrong = t*t*total+(b-B0)+(t*B0-B1)
    need(exact_error == rhs and exact_error != wrong, 'E0_FEEDBACK_FACTOR')
    distinctions['two_term_e0_feedback_factor'] = [str(exact_error), str(wrong)]
    residual = B0-b-t*B0
    need(abs(total-B0) == abs(residual)/(1-t), 'FULL_RESIDUAL_ALREADY_COMPLETE')
    distinctions['full_residual_already_complete'] = str(abs(residual)/(1-t))
    rho, a, cov, cov_y = Q(2), Q(3), Q(5), Q(15)
    coordinate_residual = a*rho*cov+rho*cov_y-2*a*rho*cov
    frame_residual = rho*(cov_y-a*cov)
    wrong_frame = rho*(cov_y-2*a*cov)
    need(coordinate_residual == frame_residual == 0 and wrong_frame != 0, 'COVECTOR_FRAME')
    distinctions['covector_frame_one_warp_factor'] = [str(frame_residual), str(wrong_frame)]
    need(Q(7,1000000) != Q(7), 'SCALED_INITIAL_VELOCITY')
    distinctions['ordinary_t_velocity_needs_ell'] = ['7/1000000', '7']
    need(p['mean_gram'][1]*3**2 == 18 and 18 != 9, 'MIXED_GRAM')
    distinctions['mixed_metric_gram_two'] = ['18', '9']
    # For the constant V=1, X=1/3 and the false wedge bound would be 1/4.
    need(Q(1) > Q(3,4)*Q(1,3), 'ZERO_WEDGE_HYPOTHESIS_NECESSARY')
    distinctions['constant_function_without_zero_wedge'] = ['1', '1/4']
    need(-2*Q(3)**2 == -18 and -18 != 0, 'LAPSE_SQUARE')
    distinctions['zero_second_lapse_not_zero_gtt2'] = ['-18', '0']
    need(authenticate(root) == before, 'LOCAL_INPUTS_CHANGED')
    return {
        'status': 'PASS_PORTABLE_CONDITIONAL_ARITHMETIC_ONLY',
        'input_sha256': PINS,
        'local_inputs_authenticated_before_math_and_after': True,
        'upstream_scientific_capsules_authenticated_by_this_run': False,
        'q_zero_wedge': str(q), 'q_previous_local': str(q_old),
        'resolvent_zero_wedge': str(1/(1-q)),
        'forcing_allowance': str(forcing), 'two_term_tail': str(tail),
        'design_error_bound': str(budget), 'design_target': str(target),
        'design_error_is_not_an_evaluated_error': True,
        'synthetic_formula_distinctions': distinctions,
        'distinction_count': len(distinctions),
        'limits': p['limits'],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-root', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = run(args.source_root.resolve())
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(report['status'])


if __name__ == '__main__':
    main()
