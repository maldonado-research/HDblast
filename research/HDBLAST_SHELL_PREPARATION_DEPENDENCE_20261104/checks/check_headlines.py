#!/usr/bin/env python3
"""Join the preparation certificates and audit published outward displays exactly."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
CAP = HERE.parent
PINS = {
    'numeric/references/SCIENTIFIC_REFERENCE.json': 'fb31050c93ddeeb7c58f576e23a9d982ff813671f434f480dce832bcd9e7a7f1',
    'theory/references/NORMAL.json': '7a43ac153d6a6af3208657ad2b362f00c054617568844313107f491eb7eec3eb',
}
OLD = 'research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/response/references/FINAL_NORMAL.json'
OLD_PIN = '1a9ebffb6bee9332bc0fd0df0f5d1d9f7f2a796825c2a563e4c644974ba4d8a3'


class Rejected(Exception):
    pass


def need(ok, witness):
    if not ok:
        raise Rejected(witness)


def regular(path):
    return path.is_file() and all(not p.is_symlink() for p in (path, *path.absolute().parents))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path, pin=None):
    need(regular(path), 'REGULAR_INPUT')
    need(pin is None or digest(path) == pin, 'PINNED_INPUT_BYTES')
    return json.loads(path.read_bytes())


def scalar(value):
    return F(value['exact'])


def interval(value):
    lo, hi = map(F, value['exact'])
    need(lo <= hi, 'ORDERED_RATIONAL_INTERVAL')
    return lo, hi


def product(left, right):
    values = [x*y for x in left for y in right]
    return min(values), max(values)


def scale(value, factor):
    return product(value, (factor, factor))


def add(left, right):
    return left[0]+right[0], left[1]+right[1]


def widened(value, error):
    return value[0]-error, value[1]+error


def display(value, bounds, witness, records):
    lo, hi = interval(value)
    left, right = map(F, bounds)
    need(left <= lo <= hi <= right, witness)
    records[witness] = {'outward_display': bounds, 'exact': value['exact']}


def finite_rule(base_v, base_t, slope_v, slope_t, lam, threshold):
    need(max(*base_v, *base_t, *slope_v, *slope_t) < 0, 'STRICT_NEGATIVE_BASE_AND_SLOPES')
    expected = max((-base_v[0])/(-slope_v[1]), (-base_t[0])/(-slope_t[1]))
    need(threshold == expected, 'EXACT_LAMBDA_THRESHOLD')
    power = F(1)
    while power <= expected:
        power *= 10
    need(lam == -power == -F(10**21), 'LEAST_POWER_OF_TEN_LAMBDA')
    return add(base_v, scale(slope_v, lam)), add(base_t, scale(slope_t, lam))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', required=True, type=Path)
    parser.add_argument('--numeric-report', required=True, type=Path)
    parser.add_argument('--raw-witness', required=True, type=Path)
    parser.add_argument('--theory-report', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    repo = args.repo_root.resolve(strict=True)
    output = args.output.absolute()
    need(all(not p.is_symlink() for p in (output, *output.parents)), 'OUTPUT_SYMLINK')
    output = output.resolve()
    need(not output.exists(), 'FRESH_OUTPUT')
    need(output.parent.is_dir(), 'OUTPUT_PARENT')
    need(all(base != output and base not in output.parents for base in (repo, CAP)), 'OUTPUT_OUTSIDE_CODE_AND_SOURCE')
    before = {name: digest(CAP/name) for name in PINS}
    reference = read(CAP/'numeric/references/SCIENTIFIC_REFERENCE.json', PINS['numeric/references/SCIENTIFIC_REFERENCE.json'])
    coarse_ref = read(CAP/'theory/references/NORMAL.json', PINS['theory/references/NORMAL.json'])
    base = read(repo/OLD, OLD_PIN)
    report, raw, coarse = read(args.numeric_report), read(args.raw_witness), read(args.theory_report)
    science = {k: v for k, v in report.items() if k not in ('source_authentication', 'portable_provenance')}
    need(science == reference and coarse == coarse_ref, 'FRESH_CERTIFICATES_MATCH_FROZEN_SCIENCE')
    need(report['campaign'] == 'ONE_DIRECTION_P1_N1024_T5_4' and scalar(report['target']) == F(5, 4), 'ONE_REGISTERED_DIRECTION_AND_TARGET')
    q, n, error = map(scalar, (report['q'], report['N_h'], report['E_h']))
    need(q == F(1201, 3200000) and error == q*n/(1-q), 'FULL_HOMOGENEOUS_RESOLVENT')
    gram = report['initial_Gram_norm_caps']
    need(n == max((scalar(gram['D0'])+F(5, 4)*scalar(gram['Dv']))/3, scalar(gram['Dr'])+scalar(gram['Dv'])), 'COMPLETE_INITIAL_NORM')
    caps = report['initial_component_caps']
    need(set(caps) == {'U0', 'U0_r', 'V0'} and all(len(row) == 5 for row in caps.values()), 'FIVE_COMPONENT_ADAPTER_ROSTER')
    for name in ('U0', 'U0_r'):
        row = list(map(scalar, caps[name]))
        need(row[0] == row[3] and row[1] == row[2] == 0 and row[4] > 0, 'METRIC_AND_SCALAR_ADAPTER_RETAINED')
    velocity = list(map(scalar, caps['V0']))
    need(velocity[1] > 0 and all(velocity[j] == 0 for j in (0, 2, 3, 4)), 'MIXED_METRIC_VELOCITY_RETAINED')
    delta, eta, ell = F(1, 4), F(1, 10**11), F(1, 10**6)
    need(scalar(report['ell']) == ell, 'FIXED_SCALE')
    need(interval(report['M2_slope']) == widened(interval(report['free_scalar_target']), delta*error), 'SCALAR_TARGET_FULL_FEEDBACK')
    for kind, factor in [('volume', delta**3/6), ('expansion', delta**2/2)]:
        expected_error = factor*error + eta*scalar(report['A_'+kind])
        need(scalar(report[kind+'_history_error']) == expected_error, 'COMPLETE_KERNEL_HISTORY_ERROR:'+kind)
        need(interval(report['kernel_'+kind+'_history']) == widened(interval(report['J_'+kind]), expected_error), 'COMPLETE_HISTORY:'+kind)
    rho = interval(raw['shared_inputs']['rho_b'])
    clock = scale(rho, ell)
    fp = scale(interval(report['Fprime']), F(1, 2))
    need(interval(report['volume_slope']) == product(product(clock, clock), product(fp, interval(report['kernel_volume_history']))), 'EXACT_VOLUME_PROPER_TIME_FACTORS')
    need(interval(report['expansion_slope']) == product(clock, product(fp, interval(report['kernel_expansion_history']))), 'EXACT_EXPANSION_PROPER_TIME_FACTORS')
    witness = report['finite_lambda_witness']
    bv, bt = interval(base['proper_v2']), interval(base['proper_Theta2'])
    need(bv == interval(witness['base_volume']) and bt == interval(witness['base_expansion']), 'PRESERVED_BASE_REPORT_JOIN')
    nv, nt = finite_rule(bv, bt, interval(report['volume_slope']), interval(report['expansion_slope']), scalar(witness['lambda']), scalar(witness['threshold']))
    need(nv == interval(witness['new_volume_second_coefficient']) and nt == interval(witness['new_expansion_second_coefficient']) and min(*nv, *nt) > 0, 'EXACT_AFFINE_SIMULTANEOUS_REVERSAL')
    coarse_bv = (F('-5.4520e16'), F('-2.6521e16'))
    coarse_bt = (F('-9.0208e21'), F('-4.4458e21'))
    need(coarse_bv[0] <= bv[0] <= bv[1] <= coarse_bv[1] and coarse_bt[0] <= bt[0] <= bt[1] <= coarse_bt[1], 'COARSE_OUTWARD_BASE_HEADLINES')
    cnv, cnt = finite_rule(coarse_bv, coarse_bt, interval(coarse['slopes']['v2']), interval(coarse['slopes']['Theta2']), F(coarse['deterministic_choice']['lambda']), scalar(coarse['deterministic_choice']['threshold']))
    need(cnv == interval(coarse['new_coefficients']['v2']) and cnt == interval(coarse['new_coefficients']['Theta2']), 'COARSE_AFFINE_REVERSAL')
    need(cnv[0] > F('7e16') and cnt[0] > F('5e21'), 'COARSE_SHORTENED_POSITIVE_LOWER_BOUNDS')
    records = {}
    for key, bounds in [('M2_slope', ['1.5393e9', '1.7068e9']), ('volume_slope', ['-2.1094e-4', '-1.5795e-4']), ('expansion_slope', ['-22.9410', '-17.4038']), ('kernel_volume_history', ['2.1313e7', '2.5493e7']), ('kernel_expansion_history', ['3.0511e8', '3.6022e8'])]:
        display(report[key], bounds, key, records)
    for key, bounds in [('base_volume', ['-5.4520e16', '-2.6521e16']), ('base_expansion', ['-9.0208e21', '-4.4458e21']), ('new_volume_second_coefficient', ['1.0343e17', '1.8442e17']), ('new_expansion_second_coefficient', ['8.3831e21', '1.8496e22'])]:
        display(witness[key], bounds, key, records)
    for key, bounds in [('M2', ['1.269e9', '2.198e10']), ('v2', ['-0.005883', '-0.0001282']), ('Theta2', ['-359.5', '-14.21'])]:
        display(coarse['slopes'][key], bounds, 'coarse_'+key, records)
    need(n < F('4.5253e10') and error < F('1.6991e7'), 'OUTWARD_NORM_AND_ERROR_UPPER_DISPLAYS')
    need(report['oracle_calls']['background'] == report['oracle_calls']['evolved_actual_first_field'] == 0, 'NO_NEW_BACKGROUND_OR_PDE_ORACLE_CALLS')
    need(before == {name: digest(CAP/name) for name in PINS} and digest(repo/OLD) == OLD_PIN, 'SOURCE_POSTCHECK')
    result = {'schema': 'HDBLAST_PREPARATION_DEPENDENCE_HEADLINES_V1', 'status': 'PASS_EXACT_PREPARATION_DEPENDENCE_HEADLINES',
              'outward_displays_checked': records, 'norm_upper_display': '4.5253e10', 'resolvent_error_upper_display': '1.6991e7',
              'lambda': str(-10**21), 'both_methods_same_finite_lambda': True,
              'full_five_component_initial_adapter_checked': True, 'full_coupled_homogeneous_resolvent_checked': True,
              'proper_time_factors_checked_exactly': True, 'same_preserved_base_report': True,
              'prepared_directions': 1, 'independent_two_observable_control_established': False, 'rank_two_established': False,
              'conditional_identity_trace': {
                  'initial_constraints': 'DERIVATION.md sections 2-3: the exact zero moment and compact metric correction preserve the complete second constraint equation; first tangent and second TT pair persist.',
                  'five_fields': ['K', 'L', 'N', 'J', 'w'],
                  'dynamics': 'DERIVATION.md sections 4-5: the same operator and physical/normal-gauge shell rows retain every coupled channel; unchanged quadratic forcing cancels in the difference.',
                  'anchored_observables': 'DERIVATION.md section 4: identical first proper-clock, measure and first-only Gauss terms cancel only in the difference; they remain in each absolute response.',
                  'family_scope': 'DERIVATION.md section 8: each fixed finite lambda has a qualitative family-dependent small-parameter interval; no effective nonlinear radius is computed.'},
              'identity_trace_is_conditional_analytic_reasoning_not_a_proof_from_json': True,
              'new_background_or_evolved_first_field_oracle_calls': 0,
              'source_sha256': {**PINS, OLD: OLD_PIN}}
    with output.open('xb') as handle:
        handle.write((json.dumps(result, sort_keys=True, indent=2)+'\n').encode())
    print(json.dumps({'status': 'PASS', 'witness': 'EXACT_PREPARATION_DEPENDENCE_HEADLINES'}, sort_keys=True))


if __name__ == '__main__':
    try:
        main()
    except (Rejected, OSError, ValueError, KeyError, TypeError) as error:
        print(json.dumps({'status': 'FAIL', 'witness': str(error) if isinstance(error, Rejected) else 'HEADLINE_INPUT_FAILURE'}, sort_keys=True))
        raise SystemExit(2)
