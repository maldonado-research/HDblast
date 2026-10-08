"""Exact sign corollary from a checked uniform branch certificate.

This does not certify an ODE trajectory. It propagates the explicit analytic
remainder bound and records its source certificate by SHA256.
"""
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path
import sys


def main():
    if len(sys.argv) != 3:
        raise SystemExit("usage: derive_curvature_suppression.py CONSTANTS.json OUT.json")
    inp, out = map(Path, sys.argv[1:])
    raw = inp.read_bytes()
    record = json.loads(raw)
    if record['status'] != 'PASS_EXACT_SUFFICIENT_INEQUALITIES_CONDITIONAL_ON_ANALYTIC_PROOF':
        raise ValueError('Expected passing exact sufficient-inequality certificate')
    rows = []
    for model in record['models']:
        if len(model['checks']) != 12 or not all(c['pass_check'] for c in model['checks']):
            raise ValueError('Missing sufficient inequalities')
        params = {key: Q(value) for key, value in model['parameters'].items()}
        k, w, f1 = params['k'], params['w'], params['f1']
        if not (k > 0 and w > 4*k and f1 != 0):
            raise ValueError('Strict suppression requires positive model gap and nonzero scalar coupling')
        coefficient = k*f1*f1/(24*(w-2*k))
        remainder = Q(model['exact_constants']['K_H'])
        delta_star = Q(model['delta_star'])
        if not (remainder > 0 and delta_star > 0):
            raise ValueError('Invalid certified interval')
        delta_sign = min(delta_star, coefficient/(2*remainder))
        margin = coefficient-remainder*delta_sign
        if not margin >= coefficient/2 > 0:
            raise ArithmeticError('Strict sign implication failed')
        rows.append(dict(
            name=model['name'], delta_sign=str(delta_sign),
            covers_entire_certified_interval=delta_sign == delta_star,
            positive_scalar_coefficient=str(coefficient),
            uniform_remainder_coefficient=str(remainder),
            minimum_suppression_coefficient=str(margin),
            maximum_suppression_coefficient=str(coefficient+remainder*delta_sign),
            display_only=dict(delta_sign=float(delta_sign),
                              minimum_suppression_coefficient=float(margin),
                              maximum_suppression_coefficient=float(coefficient+remainder*delta_sign))
        ))
    result = dict(
        status='PASS_EXACT_CURVATURE_SUPPRESSION_COROLLARY',
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256=hashlib.sha256(raw).hexdigest(),
        theorem_sha256=record['theorem_sha256'],
        arithmetic='Exact fractions.Fraction; display_only fields are approximate.',
        implication='For every 0 < delta <= delta_sign, the regular local branch satisfies '
                    '-maximum_suppression_coefficient*delta^2 <= H^2-H_metric^2 '
                    '<= -minimum_suppression_coefficient*delta^2 < 0, '
                    'where H_metric^2=k*f0*delta/3+f0^2*delta^2/36.',
        scope='Depends on the reviewed analytic theorem and exact sufficient-inequality certificate. '
              'The metric-only reference fails the scalar junction when f1 is nonzero. '
              'This is a static comparison, without a stability, novelty or cosmological-origin claim.',
        models=rows)
    out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(status=result['status'], models=[dict(name=x['name'],
                      covers_entire_certified_interval=x['covers_entire_certified_interval'],
                      **x['display_only']) for x in rows]), indent=2))


if __name__ == '__main__':
    main()
