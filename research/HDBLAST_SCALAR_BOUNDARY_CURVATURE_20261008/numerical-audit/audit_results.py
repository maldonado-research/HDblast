"""Independently recompute diagnostics from the preserved and replayed outputs.

This checks recorded arithmetic and cross-environment reproduction. It neither
reruns the solver nor turns sampled residuals into a continuum error bound.
"""
import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--reference', type=Path, required=True)
    ap.add_argument('--replay', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    models = []
    maxima = dict(cross_environment_H2=0., cross_environment_eta_b=0.,
                  setting_H2=0., junction_residual=0., sampled_constraint=0.,
                  identity_difference=0., observed_remainder_over_delta3=0.,
                  observed_scalar_shift_error_over_delta2=0.)
    count = fitted = 0
    all_pass = True
    for filename in ('family-benchmark-results.json', 'matched-pair-results.json'):
        original = json.loads((args.reference / filename).read_text())
        replay = json.loads((args.replay / filename).read_text())
        if not (original['complete'] and replay['complete']):
            raise ValueError('Incomplete inputs')
        if len(original['cases']) != len(replay['cases']):
            raise ValueError('Case count changed')
        for oc, rc in zip(original['cases'], replay['cases'], strict=True):
            p = rc['parameters']
            if p != oc['parameters']:
                raise ValueError('Model parameters changed')
            k,w,f0,f1 = (Fraction.from_float(float(p[n])) for n in ('k','w','f0','f1'))
            q = -k*f1*f1/(24*(w-2*k))
            a = -f1/(4*(w-2*k))
            rows = []
            coefficient_errors = []
            for old, new in zip(oc['rows'], rc['rows'], strict=True):
                delta = new['delta']
                if delta != old['delta']:
                    raise ValueError('Detuning changed')
                d = Fraction.from_float(delta)
                pred = k*f0*d/3 + (f0*f0/36+q)*d*d
                metric = k*f0*d/3 + f0*f0*d*d/36
                for setting in ('standard','refined'):
                    o,n = old[setting],new[setting]
                    count += 1
                    fitted += p['f1'] != 0
                    e = Fraction.from_float(n['eta_b'])
                    v,f2 = (Fraction.from_float(float(p[x])) for x in ('v','f2'))
                    W, Wp = 3*k+w*e*e/2+v*e**3/6, w*e+v*e*e/2
                    f, fp = f0+f1*e+f2*e*e/2, f1+f2*e
                    boundary = d*(W*f/9-Wp*fp/12)+d*d*(f*f/36-fp*fp/48)
                    independently_recalculated_pass = (
                        max(abs(x) for x in n['junction_residual']) < 2e-11
                        and n['max_sampled_constraint_residual'] < 2e-10
                        and abs(Fraction.from_float(n['H2'])-boundary) < Fraction(2,10**11)
                        and math.isfinite(n['H2']) and n['H2'] > 0)
                    all_pass &= independently_recalculated_pass and n['numerical_checks_pass'] and n['solver_reported_success']
                    maxima['cross_environment_H2'] = max(maxima['cross_environment_H2'],abs(n['H2']-o['H2']))
                    maxima['cross_environment_eta_b'] = max(maxima['cross_environment_eta_b'],abs(n['eta_b']-o['eta_b']))
                    maxima['junction_residual'] = max(maxima['junction_residual'],*(abs(x) for x in n['junction_residual']))
                    maxima['sampled_constraint'] = max(maxima['sampled_constraint'],n['max_sampled_constraint_residual'])
                    maxima['identity_difference'] = max(maxima['identity_difference'],abs(n['H2_identity_difference']))
                    maxima['observed_remainder_over_delta3'] = max(maxima['observed_remainder_over_delta3'],abs(float((Fraction.from_float(n['H2'])-pred)/d**3)))
                    maxima['observed_scalar_shift_error_over_delta2'] = max(maxima['observed_scalar_shift_error_over_delta2'],abs(float((Fraction.from_float(n['eta_b'])-a*d)/d**2)))
                diff = abs(new['standard']['H2']-new['refined']['H2'])
                maxima['setting_H2'] = max(maxima['setting_H2'],diff)
                h = Fraction.from_float(new['refined']['H2'])
                estimate = (h-metric)/d**2
                error = abs(estimate-q)
                if q:
                    coefficient_errors.append(error)
                rows.append(dict(delta=delta,coefficient_estimate=float(estimate),
                                 coefficient_error=float(error),relative_coefficient_error=float(error/abs(q)) if q else None,
                                 remainder_over_delta3=float((h-pred)/d**3),setting_difference=diff))
            ratios = [float(a/b) for a,b in zip(coefficient_errors,coefficient_errors[1:])] if q else []
            models.append(dict(name=p['name'],parameters=p,predicted_coefficient=float(q),
                               predicted_scalar_shift=float(a),leading_scalar_junction_inverse=float(1/(2*(w-2*k))),
                               coefficient_error_halving_ratios=ratios,rows=rows))
    algebra_original=(args.reference/'portable-algebra-results.json').read_bytes()
    algebra_replay=(args.replay/'portable-algebra-results.json').read_bytes()
    algebra=json.loads(algebra_replay)
    result = dict(status='PASS_REPLAY_AND_RECALCULATED_DIAGNOSTICS' if all_pass and algebra['pass_all'] else 'FAIL',
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  symbolic_checks=len(algebra['checks']),symbolic_replay_byte_identical=algebra_original==algebra_replay,
                  configurations_at_two_settings=count//2,final_integrations=count,fitted_final_integrations=fitted,
                  exact_zero_scalar_final_integrations=count-fitted,
                  original_runtime=dict(python='3.12.14',numpy='2.3.5',scipy='1.18.1',sympy='1.14.0'),
                  replay_runtime=json.loads((args.replay/'family-benchmark-results.json').read_text()) | {},
                  maxima=maxima,models=models,
                  limitations=['Cross-environment rerun uses the same solver source, not independent numerical software.',
                               'The solver performs additional trial integrations; final_integrations does not count them.',
                               'All reported errors, ratios, residuals and condition estimates are floating-point diagnostics, not certified solution errors.',
                               'No interval enclosure of the entire boundary-value problem or stability calculation follows.'])
    result['replay_runtime']={k:result['replay_runtime'][k] for k in ('python','numpy','scipy')}
    args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('models','limitations')},indent=2))
    if result['status'].startswith('FAIL'):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
