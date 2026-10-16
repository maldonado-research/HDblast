#!/usr/bin/env python3
"""Exact prospective budget arithmetic; this computes no response candidate."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def need(ok, text):
    if not ok:
        raise ValueError(text)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    need(not args.output.exists(), 'refuse overwrite')
    inputs = ROOT/'followon-response/inputs'
    paths = {'source': inputs/'ACTUAL_SOURCES_NORMAL.json',
             'retained_response': inputs/'REFINED_RESPONSE_NORMAL.json',
             'evolved_signs': inputs/'EVOLVED_SCALAR_SIGNS_NORMAL.json'}
    source = json.loads(paths['source'].read_text())
    retained = json.loads(paths['retained_response'].read_text())
    q = Q(retained['retained_join']['q'])
    B = Q(retained['response']['X_simple_upper'])
    F = Q(source['global_cover']['summary']['positive_bundle_bulk_forcing_H_bound']['F_H_upper'])
    G = Q(retained['retained_join']['G_H'])
    need(q == Q(2109, 10**6) and B == Q('1.03e33'), 'same retained X contraction and solution cap')
    need(F == 17387513156573 and G == 3851400 and G < 4*10**6, 'complete current bulk and inherited full boundary sources')
    hypothetical = Q(10**29)
    tail = q*q*B
    source_lift = 3*F+G
    omitted = (1+q)*source_lift
    mean_full = 3*(tail+hypothetical)
    mean_omitted = 3*(tail+omitted+hypothetical)
    need(tail == Q('4.58131743e27'), 'independent two-term exact Neumann tail')
    need(mean_full == Q('3.1374395229e29'), 'probability mean norm one and X value weight three')
    need(mean_omitted < Q('3.138e29'), 'explicitly retained true source omission budget')
    need(omitted < 10**14, 'source lift scale in current scaled X units')
    signs = json.loads(paths['evolved_signs'].read_text())
    margin = Q('1.85e30')
    need(all(Q(row['signed_magnitude_exceeds_common']) == margin for row in signs['persistence']['cells']), 'same evolved mean sign margin')
    increment = source_lift/(1-q)
    increment_mean = 3*increment
    need(increment_mean < margin and increment_mean/margin < Q('8.477e-17'), 'auxiliary zero-source mean signs persist')
    report = {'status': 'PASS_INDEPENDENT_PROSPECTIVE_RESPONSE_BUDGET',
              'q': str(q), 'retained_U_X_bound': str(B), 'current_bulk_F_H': str(F),
              'inherited_full_boundary_G_H': str(G), 'K_C0_to_X_bound': '3', 'L_C0_to_X_bound': '1',
              'exact_N1_tail_X_upper': str(tail), 'hypothetical_total_numerical_X_error': str(hypothetical),
              'scalar_probability_mean_error_full_N1_upper': str(mean_full),
              'source_lift_X_upper': str(source_lift), 'propagated_true_source_omission_X_upper': str(omitted),
              'scalar_probability_mean_error_omitted_source_V_upper': str(mean_omitted),
              'auxiliary_forced_increment_X_upper': str(increment),
              'auxiliary_forced_increment_mean_upper': str(increment_mean),
              'increment_mean_cap_to_existing_sign_margin_ratio': str(increment_mean/margin),
              'forced_increment_two_term_tail_X_upper': str(q*q*increment),
              'interpretation': 'Uhom=(I-T)^-1 Ufree retains the full coupled reduced bulk and Robin operator. W=U-Uhom=(I-T)^-1(KF+LG) is an auxiliary source increment. The huge current mean signs survive this auxiliary source removal and do not diagnose the forced increment. Neither split term is asserted separately to satisfy the actual quadratic physical constraints or gauge completion. A total-response error certificate allowing 1e29 in X does not certify resolution of W at its much smaller upper-bound scale; this does not prove the actual approximation error is large or resolution impossible.',
              'candidate_evaluated': False, 'continuous_candidate_X_error_certified': False,
              'limits': 'Arithmetic for conditional prospective bounds only. A global continuous C1 reconstructed candidate, valid H-operator approximation and complete numerical/input/quadrature error budget are still required. No angular pointwise factor enters the scalar mean functional.',
              'input_sha256': {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in paths.items()},
              'review_program_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True)+'\n')
    print(report['status'])


if __name__ == '__main__':
    main()
