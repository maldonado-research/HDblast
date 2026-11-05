"""Independent exact guards for odd-part graph control and past sign transfer.

Seven selected input copies are authenticated. No project executable, numerical
PDE solver or field/source/bump oracle is imported. The retained analytic and
future-certificate arguments remain premises; finite checks do not prove them.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent


def need(condition, label):
    if not condition:
        raise ValueError(label)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def pair(raw):
    need(isinstance(raw, list) and len(raw) == 2, 'interval needs two endpoints')
    lo, hi = map(F, raw)
    need(lo <= hi, 'interval endpoints reversed')
    return lo, hi


def multiply(a, b):
    products = [x * y for x in a for y in b]
    return min(products), max(products)


def strings(interval):
    return [str(x) for x in interval]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    out = args.output
    need(out.is_absolute(), 'absolute output required')
    out = out.resolve()
    need(not out.exists() and not out.is_relative_to(ROOT), 'fresh output outside component required')
    need(out.parent.is_dir() and all(not p.is_symlink() for p in [out, *out.parents]), 'regular output ancestry required')
    pin_bytes = (ROOT / 'SOURCE_PINS.json').read_bytes()
    pins = json.loads(pin_bytes)
    need(pins['schema'] == 'HDBLAST_PAST_FUTURE_INDEPENDENT_SOURCE_PINS_V1', 'source schema')
    input_bytes = {}
    for name, row in pins['files'].items():
        rel = Path(name)
        need(not rel.is_absolute() and '..' not in rel.parts and name == rel.as_posix(), 'relative input path')
        path = ROOT / rel
        need(path.is_file() and all(not p.is_symlink() for p in [path, *path.parents]), 'regular source input')
        raw = path.read_bytes()
        need(len(raw) == row['bytes'] and digest(raw) == row['sha256'], 'source input mismatch: ' + name)
        need(hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest() == row['git_blob_sha1'], 'input blob mismatch')
        input_bytes[name] = raw
    checks = []

    def gate(condition, label):
        need(label not in checks, 'duplicate guard name')
        need(condition, label)
        checks.append(label)

    gate(len(input_bytes) == 7, 'seven_selected_primary_input_copies_authenticated')
    old = json.loads(input_bytes['inputs/OLD_FUTURE_REPORT.json'])
    contact = json.loads(input_bytes['inputs/CONTACT_NORMAL.json'])
    shell = json.loads(input_bytes['inputs/SHELL_BOX.json'])
    bg = json.loads(input_bytes['inputs/WIDE_BACKGROUND.json'])
    ell = F(1, 1000000)
    gate(old['schema'] == 'hdblast-inertial-shell-waveform-v1' and old['status'] == 'PASS', 'retained_actual_future_certificate_schema_status')
    gate(old['state']['ell'] == str(ell) and old['state']['initial_velocity'] == 'exp(-z) Q(z/ell)', 'same_width_and_nonzero_registered_initial_velocity')
    gate(F(old['state']['T']) == 3 * ell and old['state']['domain'] == [0, 3], 'old_actual_future_certificate_covers_full_comparison_interval')
    gate(bg['candidate_state_exact']['ell'] == str(ell) and F(bg['candidate_state_exact']['T']) == 3 * ell and F(bg['candidate_state_exact']['Z']) == 8 * ell, 'same_raw_wide_background_interval')
    gate(bg['candidate_state_exact']['delta'] == shell['parameters_exact']['delta'] == old['state']['delta'] == '1/1000', 'same_branch_detuning')
    gate(bg['candidate_state_exact']['c'] == shell['parameters_exact']['c'] == '96714321425692/161839258930731', 'same_registered_model_parameter')
    gate(old['mathematical_inputs']['background']['sha256'] == pins['files']['inputs/WIDE_BACKGROUND.json']['sha256'], 'old_actual_report_pins_same_wide_background')
    gate(old['mathematical_inputs']['shell']['sha256'] == pins['files']['inputs/SHELL_BOX.json']['sha256'], 'old_actual_report_pins_same_raw_shell')
    gate(contact['status'] == 'PASS_CONDITIONAL_PAST_SHELL_CONTACT_SCALARS' and contact['physical_scalar_trace']['normalization_coefficients'] == {'scalar_factor': '3', 'q_factor': '5/2', 'graph_ell_power': -2}, 'physical_scalar_endpoint_normalization_retained')
    mu = pair(contact['raw_derived_boundary_intervals']['mu_b'])
    gamma = pair(contact['raw_derived_boundary_intervals']['C_b_over_two'])
    rho = pair(shell['boxes']['rho_b']['exact'])
    pb = pair(shell['boxes']['p_b']['exact'])
    gate(F(1, 1000) < mu[0] <= mu[1] < F(1, 400), 'actual_positive_endpoint_mass_ceiling')
    gate(rho[0] >= 129 and pb[1] <= -F(3, 25000), 'physical_scalar_denominator_lower_guards')
    gate(gamma[0] / mu[1] + F(5, 2) > 5500, 'static_boundary_coercivity_including_endpoint_mass')
    gate(F(129)**2 / 3 > 5500 and F(6000) > 5500, 'actual_and_auxiliary_static_bulk_coercivity')
    epsilon = F(35, 4) * ell**2 * (3 * ell)**2 / (5500 * ell**2)
    gate(epsilon == F(63, 4400) * ell**2 and 0 < epsilon < 1, 'even_temporal_perturbation_graph_coercivity_loss')
    gate(F(2835, 8) * ell**4 < 1, 'retained_common_operator_wave_Volterra_contraction')
    # O = (Z_plus - Z_minus)/2 has zero initial displacement and velocity g.
    initial_O_energy = ell**2 / 2
    O_trace_square = 4 * initial_O_energy / 500
    delta_trace_square = 4 * O_trace_square
    gate(delta_trace_square == 2 * ell**2 / 125 and delta_trace_square < ell**2, 'odd_part_energy_and_single_doubling_give_delta_trace_below_ell')
    A = (16 + ell)**2 + 103000 * ell**2
    O_graph_square = ell**2 * A / (1 - epsilon)
    delta_graph_square = 4 * O_graph_square
    gate(delta_graph_square < (33 * ell)**2, 'independent_odd_part_graph_estimate_delta_norm_below_33ell')
    gate(mu[1] * ell < F(1, 20000)**2, 'endpoint_H_mass_to_boundary_graph_component')
    B = F(25000, 129**3) * (F(33, 20000) / ell + ell / 160)
    gate(B == F(264000000001, 13738809600) and B < 20, 'physical_first_scalar_difference_exact_absolute_ceiling')
    G = pair(old['physical_gain_interval'])
    gate(G[1] < 0, 'old_physical_gain_has_expected_negative_sign')
    points = {}
    for name, row in old['sign_witnesses'].items():
        x = F(row['x'])
        H0, H = pair(row['H0']), pair(row['H'])
        E = F(row['error_upper'])
        gate(E >= 0 and H[0] <= H0[0] - E and H[1] >= H0[1] + E, name + '_actual_H_contains_comparator_plus_retained_error')
        ellS = pair(row['ell_S'])
        exact_product = multiply(multiply(G, (1 - (ell*x)**2, F(1))), H)
        gate(ellS[0] <= exact_product[0] and ellS[1] >= exact_product[1], name + '_ell_S_contains_physical_gain_cosh_actual_H_product')
        future = tuple(y / ell for y in ellS)
        past = future[0] - B, future[1] + B
        points[name] = {'positive_scaled_time': str(x), 'future_actual_S_interval': strings(future),
                        'past_actual_S_enclosure': strings(past)}
    rising = tuple(map(F, points['rising']['past_actual_S_enclosure']))
    falling = tuple(map(F, points['falling']['past_actual_S_enclosure']))
    gate(rising[1] < -37376 and falling[0] > 37379, 'strict_past_lobe_signs_survive_absolute_transport')
    gate(-F(old['sign_witnesses']['falling']['x']) < -F(old['sign_witnesses']['rising']['x']), 'coarse_past_zero_bracket_is_open_and_correctly_ordered')
    for name in ['central_left', 'central_right', 'late_left', 'late_right']:
        iv = tuple(map(F, points[name]['past_actual_S_enclosure']))
        gate(iv[0] < 0 < iv[1], name + '_old_narrow_sign_witness_does_not_transfer')
    tail_H = pair(old['tail_theorem']['actual_H_interval'])
    tail_ellS = multiply(multiply(G, (1 - (3 * ell)**2, F(1))), tail_H)
    tail_future = tuple(v / ell for v in tail_ellS)
    tail_past = tail_future[0] - B, tail_future[1] + B
    gate(tail_future[1] < 0 and tail_past[0] < 0 < tail_past[1], 'negative_future_tail_sign_not_implied_in_past_by_this_ceiling')
    gate(all((ROOT / name).read_bytes() == raw for name, raw in input_bytes.items()), 'selected_input_bytes_unchanged_after_independent_arithmetic')
    report = {
        'schema': 'HDBLAST_PAST_FUTURE_INDEPENDENT_ARITHMETIC_V1',
        'status': 'PASS_CONDITIONAL_ODD_PART_GRAPH_BOUND_AND_COARSE_PAST_ZERO',
        'source_pins_sha256': digest(pin_bytes),
        'verifier_sha256': digest(Path(__file__).read_bytes()),
        'proof_sha256': digest((ROOT / 'PROOF.md').read_bytes()),
        'finite_guards': len(checks), 'check_names': checks,
        'ell': str(ell), 'graph_coercivity_loss': str(epsilon),
        'delta_master_trace_square_upper': str(delta_trace_square),
        'delta_graph_norm_square_upper': str(delta_graph_square),
        'physical_scalar_difference_strict_upper': str(B),
        'safe_difference_strict_upper': '20', 'selected_point_transfers': points,
        'tail_future_actual_S_enclosure': strings(tail_future),
        'tail_past_actual_S_enclosure': strings(tail_past),
        'coarse_past_zero_scaled_bracket': ['-7/4', '-5/4'],
        'conditional_zero_argument': 'Actual first scalar continuity plus strict opposite lobe signs; at least one zero, no uniqueness or count assertion.',
        'premises': ['Retained regular branch, shell box and full coupled physical first reconstruction.',
                     'Common positive endpoint-mass operator and smooth unique Wentzell developments.',
                     'Older actual future first-scalar enclosure theorem; no legacy waveform campaign rerun.'],
        'limits': ['Absolute first-linear model-unit difference only; no relative error, individual amplitude below20, finite nonlinear window, radiation or cosmological-origin claim.',
                   'No transfer of old narrow zero brackets, tiny tail sign or quadratic response integral.',
                   'Existing general time reflection, energy and interval transport principles are not new mathematics.'],
        'activity_this_process': {'selected_inputs_authenticated': 7, 'legacy_or_project_producer_imports': 0,
                                  'numerical_PDE_solves': 0, 'source_background_bump_or_waveform_oracle_queries': 0}
    }
    with out.open('x') as stream:
        stream.write(json.dumps(report, sort_keys=True, indent=2) + '\n')
    print(json.dumps({'status': report['status'], 'finite_guards': len(checks), 'sha256': digest(out.read_bytes())}, sort_keys=True))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError) as error:
        print(json.dumps({'status': 'FAIL', 'witness': str(error)}, sort_keys=True))
        sys.exit(2)
