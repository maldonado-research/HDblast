"""Exact scalar guards and retained-witness transport for one registered packet.

Only standard-library rational arithmetic and selected-source authentication.
The continuous arguments and inherited enclosure theorem remain explicit.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PIN = '6b19c5ed2bbf4ac895d2fef24a1e926f9c656cd8965dfd5b99feb2c608283816'


def need(ok, reason):
    if not ok:
        raise ValueError(reason)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def regular(path):
    return path.is_file() and all(not p.is_symlink() for p in [path, *path.parents])


def authenticate():
    path = HERE / 'SOURCE_PINS.json'
    need(regular(path), 'REGULAR_SOURCE_PINS')
    raw = path.read_bytes()
    need(sha(raw) == PIN, 'SOURCE_PIN_FILE')
    pins = json.loads(raw)
    need(pins['schema'] == 'HDBLAST_PAST_FUTURE_PRIMARY_SOURCE_PINS_V1', 'SOURCE_PIN_SCHEMA')
    need(len(pins['files']) == 7, 'SEVEN_SELECTED_INPUTS')
    for name, row in pins['files'].items():
        relative = Path(name)
        need(not relative.is_absolute() and '..' not in relative.parts and
             relative.as_posix() == name, 'SAFE_RELATIVE_INPUT')
        source = HERE / relative
        need(regular(source), 'REGULAR_INPUT:' + name)
        content = source.read_bytes()
        need(len(content) == row['bytes'] and sha(content) == row['sha256'], 'INPUT_BYTES:' + name)
        blob = hashlib.sha1(b'blob ' + str(len(content)).encode() + b'\0' + content).hexdigest()
        need(blob == row['git_blob_sha1'], 'INPUT_GIT_BLOB:' + name)
        need(len(row['source_commit']) == 40 and
             set(row['source_commit']) <= set('0123456789abcdef'), 'SOURCE_COMMIT')
        origin = row['source_repository_path']
        need(not origin.startswith('/') and '..' not in Path(origin).parts and
             '/workspace/' not in origin, 'PUBLIC_SOURCE_LOCATION')
    return pins


def run():
    pins = authenticate()
    load = lambda name: json.loads((HERE / 'inputs' / name).read_bytes())
    prior = load('CONTACT_NORMAL.json')
    future = load('OLD_FUTURE_REPORT.json')
    shell = load('SHELL_BOX.json')
    background = load('WIDE_BACKGROUND.json')
    checks = []

    def check(ok, name):
        need(ok, name)
        checks.append(name)

    check(prior['status'] == 'PASS_CONDITIONAL_PAST_SHELL_CONTACT_SCALARS' and
          prior['finite_checks'] == 27, 'retained_contact_certificate_status')
    check(future['schema'] == 'hdblast-inertial-shell-waveform-v1' and
          future['status'] == 'PASS', 'retained_actual_future_certificate_status')
    ell = F(1, 10**6)
    state = future['state']
    model = shell['parameters_exact']
    candidate = background['candidate_state_exact']
    check(F(state['ell']) == ell == F(candidate['ell']) and
          F(state['T']) == 3*ell == F(candidate['T']) and
          F(state['Z']) == 8*ell == F(candidate['Z']), 'same_width_time_and_collar_state')
    check(F(state['delta']) == F(model['delta']) == F(candidate['delta']) == F(1, 1000) and
          F(model['c']) == F(candidate['c']) == F(96714321425692, 161839258930731),
          'same_registered_background_parameters')
    check(state['initial_velocity'] == 'exp(-z) Q(z/ell)' and
          state['time_coordinate'] == 'x=t/ell=tau/(rho_b ell)',
          'nonzero_registered_velocity_and_proper_clock')
    for role, name in [('shell', 'SHELL_BOX.json'), ('background', 'WIDE_BACKGROUND.json')]:
        row = future['mathematical_inputs'][role]
        data = (HERE / 'inputs' / name).read_bytes()
        check(row['bytes'] == len(data) and row['sha256'] == sha(data),
              'retained_future_uses_same_raw_' + role)
    loss = F(prior['energy']['relative_energy_change_upper'])
    check(loss == F(63, 4400)*ell**2 and 0 < loss < 1 and
          F(prior['energy']['scaled_potential_lower']) == 5500*ell**2,
          'static_coercivity_and_relative_temporal_loss')
    mass = list(map(F, prior['raw_derived_boundary_intervals']['mu_b']))
    old_mass = list(map(F, future['shell_parameter_intervals']['mu']))
    check(F(1, 1000) < mass[0] <= mass[1] < F(1, 400) and
          old_mass[0] <= mass[0] <= mass[1] <= old_mass[1],
          'common_endpoint_mass_enclosures_and_safe_caps')
    check(prior['energy']['backward_elapsed_slab'] == ['0', '3'] and
          F(prior['energy']['support_outer_depth']) == 5 <
          F(prior['energy']['actual_static_collar_depth']) == 8,
          'actual_packet_support_stays_inside_collar')
    caps = prior['graph_energy']['unit_bump_L2_derivative_caps']
    check(F(caps['Q']) == 1 and F(caps['Q_r']) == 16 and
          F(prior['energy']['all_real_W_upper']) < 103000,
          'retained_unit_bump_and_static_potential_majorants')
    physical = prior['physical_scalar_trace']
    check(physical['normalization_coefficients'] ==
          {'graph_ell_power': -2, 'q_factor': '5/2', 'scalar_factor': '3'} and
          physical['background_proper_clock'] == 'tau=rho_b*t' and
          physical['actual_selected_first_linear_coefficient_bound'] is True,
          'same_embedding_included_physical_scalar_normalization')
    # Exact temporal-reflection cancellation: omega is even in scaled time.
    a0 = F(35, 4)
    check(a0 * 3**2 / 5500 == F(63, 4400),
          'three_time_coercivity_loss_factor')
    E0 = 2*ell**2
    check(F(1, 2)*(-2)**2*ell**2 == E0,
          'zero_displacement_minus_twice_velocity_initial_energy')
    q_square = E0 / 125
    check(q_square < ell**2 and F(2, 1000)*2**2 == F(1, 125),
          'complete_endpoint_mass_and_contact_trace_integral')
    twice_half_initial = 4*ell**2*((16+ell)**2 + 103000*ell**2)
    graph_square = twice_half_initial/(1-loss)
    check(graph_square < (33*ell)**2,
          'complete_half_graph_energy_norm_below_33ell')
    sqrt_factor = F(1, 20000)
    check(mass[1]*ell < sqrt_factor**2,
          'weighted_endpoint_graph_factor')
    rho_guard = F(physical['rho_b_lower_guard'])
    p_guard = F(physical['minus_p_b_lower_guard'])
    check(rho_guard == 129 and p_guard == F(3, 25000),
          'physical_reconstruction_denominator_guards')
    prefactor = F(3)/(rho_guard**3*p_guard)
    q_term = F(5, 2)*F(1, 400)*ell
    check(prefactor == F(25000, 129**3) and q_term == ell/160,
          'scalar_reconstruction_q_term_and_prefactor')
    B = prefactor*(33/(20000*ell)+q_term)
    check(B == F(264000000001, 13738809600) and B < 20,
          'uniform_absolute_past_future_difference_below_twenty')

    transported = {}
    for label, row in future['sign_witnesses'].items():
        x = F(row['x'])
        H0 = list(map(F, row['H0']))
        H = list(map(F, row['H']))
        error = F(row['error_upper'])
        check(error > 0 and H[0] <= H0[0]-error <= H0[1]+error <= H[1],
              'actual_error_enlargement_present_' + label)
        exported = list(map(F, row['ell_S']))
        check(exported[0] <= exported[1] and 0 <= x <= 3,
              'ordered_exported_actual_witness_' + label)
        future_S = [v/ell for v in exported]
        past_S = [future_S[0]-B, future_S[1]+B]
        sign = 'negative' if past_S[1] < 0 else 'positive' if past_S[0] > 0 else 'not_determined'
        transported[label] = {
            'matched_positive_scaled_time': str(x),
            'future_exported_ell_S_interval': list(map(str, exported)),
            'future_actual_S_interval': list(map(str, future_S)),
            'conditional_past_actual_S_interval': list(map(str, past_S)),
            'past_sign_from_this_bound': sign}
    rising, falling = transported['rising'], transported['falling']
    check(F(rising['matched_positive_scaled_time']) == F(5, 4) and
          F(rising['conditional_past_actual_S_interval'][1]) < -37376,
          'strict_negative_past_lobe_margin_at_minus_five_quarters')
    check(F(falling['matched_positive_scaled_time']) == F(7, 4) and
          F(falling['conditional_past_actual_S_interval'][0]) > 37379,
          'strict_positive_past_lobe_margin_at_minus_seven_quarters')
    check(F(-7, 4) < F(-5, 4) < -1,
          'distinct_past_contact_times_and_open_zero_bracket')
    check(F(rising['conditional_past_actual_S_interval'][0]) > -41261 and
          F(falling['conditional_past_actual_S_interval'][1]) < 41264,
          'readable_outward_integer_lobe_enclosures')
    for label in ['central_left', 'central_right', 'late_left', 'late_right']:
        low, high = map(F, transported[label]['conditional_past_actual_S_interval'])
        check(low < 0 < high, 'fine_old_witness_does_not_keep_sign_' + label)
    gain = list(map(F, future['physical_gain_interval']))
    tail_H = list(map(F, future['tail_theorem']['actual_H_interval']))
    tail_absolute_upper = (-gain[0])*tail_H[1]/ell
    check(gain[0] <= gain[1] < 0 < tail_H[0] <= tail_H[1] and
          0 < tail_absolute_upper < B,
          'old_small_tail_margin_does_not_survive_difference_ceiling')
    check(authenticate() == pins, 'all_selected_inputs_unchanged_after_arithmetic')
    return {
        'schema': 'HDBLAST_PAST_FUTURE_PRIMARY_V1',
        'status': 'PASS_CONDITIONAL_PAST_FUTURE_PHYSICAL_SCALAR_TRANSFER',
        'client_date': '2026-11-04', 'client_timezone': 'America/Los_Angeles',
        'source_pins_sha256': PIN, 'proof_sha256': sha((HERE/'PROOF.md').read_bytes()),
        'producer_sha256': sha(Path(__file__).read_bytes()),
        'selected_public_inputs_authenticated': len(pins['files']),
        'finite_checks': len(checks), 'check_names': checks,
        'state': {'ell': str(ell), 'scaled_elapsed_interval': ['0', '3'],
                  'same_nonzero_registered_velocity': True,
                  'coordinate_times': '(-ell*s,+ell*s)',
                  'proper_times': '(-rho_b*ell*s,+rho_b*ell*s)',
                  'same_embedding_included_quadrupole_scalar': True},
        'difference_energy': {
            'initial_displacement': '0', 'initial_velocity': '-2ell*exp(-ell*r)*Q(r)',
            'initial_basic_energy_upper': str(E0),
            'boundary_master_difference_square_upper': str(q_square),
            'boundary_master_difference_strict_upper': str(ell),
            'relative_temporal_loss_upper': str(loss),
            'twice_initial_half_graph_energy_upper': str(twice_half_initial),
            'A0_delta_state_squared_upper': str(graph_square),
            'A0_delta_state_norm_strict_upper': str(33*ell),
            'endpoint_terms_retained': True,
            'nonlocal_half_power_compact_support_claimed': False,
            'physical_Noether_or_injected_energy_claimed': False},
        'physical_scalar_difference': {
            'extra_time_reflection_minus_sign': False,
            'prefactor_upper': str(prefactor), 'weighted_endpoint_sqrt_factor': str(sqrt_factor),
            'q_term_upper': str(q_term), 'derived_uniform_upper': str(B),
            'safe_strict_uniform_upper': '20', 'units': 'fixed model units of first scalar coefficient',
            'individual_amplitude_below_twenty_claimed': False,
            'finite_nonlinear_parameter_window_supplied': False},
        'selected_point_transfers': transported,
        'coarse_zero': {'negative_coefficient_at_scaled_past_time': '-5/4',
                        'strict_negative_upper': '-37376',
                        'positive_coefficient_at_scaled_past_time': '-7/4',
                        'strict_positive_lower': '37379',
                        'open_scaled_past_interval': ['-7/4', '-5/4'],
                        'at_least_one_by_retained_continuity': True,
                        'uniqueness_claimed': False},
        'tail_nontransfer': {'retained_actual_tail_H_interval': list(map(str, tail_H)),
                             'future_physical_absolute_upper': str(tail_absolute_upper),
                             'below_difference_ceiling': True,
                             'past_sign_not_established': True},
        'retained_premises': [
            'Actual branch and shell enclosures, positive endpoint-mass operator and smooth unique linked-domain linear development.',
            'Coupled embedding-included scalar reconstruction, proper-clock interpretation and gauge equivalence.',
            'The older same-state actual future pointwise enclosure theorem; its producer, gain proof and continuous-cell campaign are not rerun.',
            'Continuity needed for the intermediate-value conclusion; exact arithmetic alone does not prove the analytic premises.'],
        'prior_second_order_scope': 'Earlier conditional formal complete second Taylor coefficients through future contact exist. No quantitative particular past second-order observable/baseline or exact nonlinear defect bound is supplied here.',
        'limits': [
            'Difference and coarse zero for one registered first scalar coefficient; not a scalar monopole or globally signed angular field.',
            'No fine old zero bracket, late-tail sign, relative-error, unique zero or quadratic response-integral transport.',
            'No physical energy, remote source selection, radiation, SI calibration or cosmological-origin inference.',
            'Time reflection, graph energy and interval widening are standard methods, not claimed as new general mathematics.'],
        'activity': {'producer_processes_this_mode': 1, 'project_scientific_imports': 0,
                     'legacy_producer_or_checker_invocations': 0, 'new_background_ODE_images': 0,
                     'background_field_source_or_bump_queries': 0, 'numerical_PDE_solves': 0}}


def main():
    class Parser(argparse.ArgumentParser):
        def error(self, message):
            raise ValueError('CLI_ARGUMENTS')
    parser = Parser(add_help=False)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    need(args.output.is_absolute(), 'ABSOLUTE_OUTPUT_REQUIRED')
    need(all(not path.is_symlink() for path in [args.output, *args.output.parents]), 'OUTPUT_SYMLINK')
    output = args.output.resolve()
    need(not output.exists() and output.parent.is_dir() and not output.is_relative_to(HERE),
         'FRESH_EXTERNAL_OUTPUT')
    result = run()
    raw = (json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n').encode()
    with output.open('xb') as stream:
        stream.write(raw)
    print(json.dumps({'status': result['status'], 'finite_checks': result['finite_checks'],
                      'sha256': sha(raw)}, sort_keys=True))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, ArithmeticError, KeyError, TypeError, IndexError, AttributeError) as error:
        print(json.dumps({'status': 'REJECTED', 'witness': str(error)}, sort_keys=True))
        raise SystemExit(2)
