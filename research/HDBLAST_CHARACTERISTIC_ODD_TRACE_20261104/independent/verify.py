"""Exact selected-source and scalar guards for the direct odd-trace proof."""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
import math
from pathlib import Path, PurePosixPath
import sys

ROOT = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise ValueError(message)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def product_interval(*intervals):
    values = [math.prod(values) for values in itertools.product(*intervals)]
    return min(values), max(values)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    need(args.output.is_absolute(), 'absolute output required')
    need(all(not p.is_symlink() for p in [args.output, *args.output.parents]), 'regular output ancestry required')
    out = args.output.resolve()
    need(not out.exists() and out.parent.is_dir() and not out.is_relative_to(ROOT), 'fresh external output required')
    pin_bytes = (ROOT/'SOURCE_PINS.json').read_bytes()
    pins = json.loads(pin_bytes)
    names = {'inputs/CONTACT_PROOF.md', 'inputs/CONTACT_NORMAL.json',
             'inputs/REGISTERED_COMPARISON_PROOF.md', 'inputs/OLD_FUTURE_PROOF.md',
             'inputs/OLD_FUTURE_REPORT.json', 'inputs/PHYSICAL_SHELL_RESPONSE.md'}
    need(pins['schema'] == 'HDBLAST_ODD_TRACE_INDEPENDENT_SOURCE_PINS_V1' and set(pins['files']) == names,
         'exact six selected input records required')
    raw_inputs = {}
    for name, row in pins['files'].items():
        path = PurePosixPath(name)
        need(not path.is_absolute() and path.as_posix() == name and '..' not in path.parts and '\\' not in name,
             'safe selected input path required')
        p = ROOT/name
        need(p.is_file() and all(not a.is_symlink() for a in [p, *p.parents]), 'regular selected input required')
        raw = p.read_bytes()
        blob = hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
        need(len(raw) == row['bytes'] and digest(raw) == row['sha256'] and blob == row['git_blob_sha1'],
             'selected input authentication failed: '+name)
        need(row['source_commit'] == '74d0455f1446dbb10dd2f63721d966349d092162' and
             row['source_repository_path'].startswith('research/HDBLAST_PAST_FUTURE_SCALAR_TRANSFER_20261104/') and
             '/workspace/' not in row['source_repository_path'], 'public selected source identity required')
        raw_inputs[name] = raw
    contact = json.loads(raw_inputs['inputs/CONTACT_NORMAL.json'])
    old = json.loads(raw_inputs['inputs/OLD_FUTURE_REPORT.json'])
    checks = []

    def check(ok, label):
        need(ok, label)
        checks.append(label)

    check(len(raw_inputs) == 6, 'six_selected_public_inputs_authenticated')
    check(contact['status'] == 'PASS_CONDITIONAL_PAST_SHELL_CONTACT_SCALARS', 'retained_contact_certificate_status')
    check(old['schema'] == 'hdblast-inertial-shell-waveform-v1' and old['status'] == 'PASS', 'retained_actual_future_certificate_status')
    ell = F(1, 10**6)
    T = 3*ell
    check(F(old['state']['ell']) == ell and F(old['state']['T']) == T and
          old['state']['initial_velocity'] == 'exp(-z) Q(z/ell)', 'same_registered_width_time_and_nonzero_velocity')
    m, M, B, Gamma, Gstar = F(1,1000), F(103000), F(300), F(300), F(3,250)
    mu_lo, mu_hi = map(F, contact['raw_derived_boundary_intervals']['mu_b'])
    beta_lo, beta_hi = map(F, contact['raw_derived_boundary_intervals']['beta_all_real_t'])
    gamma_lo, gamma_hi = map(F, contact['raw_derived_boundary_intervals']['C_b_over_two'])
    check(m < mu_lo < mu_hi, 'positive_actual_boundary_mass_lower_bound')
    check(0 < beta_lo < beta_hi < B, 'complete_dynamic_spring_coarse_bound')
    check(0 < gamma_lo < gamma_hi < Gamma, 'actual_spatial_conjugation_constant_bound')
    check(0 < F(contact['energy']['all_real_W_lower']) < F(contact['energy']['all_real_W_upper']) < M,
          'actual_bulk_potential_coarse_bound')
    rho = F(contact['physical_scalar_trace']['rho_b_lower_guard'])
    abs_p = F(contact['physical_scalar_trace']['minus_p_b_lower_guard'])
    check(rho >= 129 and abs_p >= F(3,25000) and 3/(rho**3*abs_p) < Gstar,
          'physical_gain_denominator_and_safe_absolute_upper')
    check(F(contact['energy']['support_outer_depth']) == 5 < F(contact['energy']['actual_static_collar_depth']) == 8,
          'actual_first_support_inside_retained_collar')
    check(contact['physical_scalar_trace']['normalization_coefficients']['scalar_factor'] == '3' and
          contact['physical_scalar_trace']['background_proper_clock'] == 'tau=rho_b*t',
          'retained_physical_scalar_and_clock_normalization')

    # The two continuous analytic inequalities are proved in PROOF.md.
    free_amplitude = ell
    input_history = 2*ell
    r = M*T*T/2
    s = B*T*T/(2*m)
    d = (1-r)*(1-s)-T*r/m
    check(0 < r < 1 and 0 < s < 1 and d > 0, 'positive_two_by_two_majorant_inverse')
    A = (free_amplitude*(1-s)+T*input_history/m)/d
    V = ((1-r)*input_history/m+(r/m)*free_amplitude)/d
    check((1-r)*A-T*V == ell and -(r/m)*A+(1-s)*V == 2*ell/m,
          'both_complete_majorant_rows_close_exactly')
    check(A == ell*(1-s+2*T/m)/d and V == ell*(2-r)/(m*d),
          'independent_closed_form_majorants_agree')
    check(A == F(40239946000,39999927404405029) and V == F(79999981460000,39999927404405029),
          'printed_actual_field_and_boundary_velocity_majorants')
    C = M*T*A+(1+Gamma*T)*V
    check(C == F(80084415586628,39999927404405029), 'incoming_potential_and_boundary_history_terms_retained')
    global_bound = Gstar*(2+C)
    check(global_bound == F(120119908836595029,4999990925550628625) < F(1,40),
          'uniform_actual_scalar_difference_below_one_fortieth')
    x_late = F(49,25)
    exponent = 1/((x_late-1)*(2-x_late))-4
    check(exponent == F(529,24) > 22, 'late_bump_exponent_exact_lower_bound')
    exp_two_partial = sum(F(2)**k/math.factorial(k) for k in range(6))
    check(exp_two_partial == F(109,15) > 7 and 7**11 > 10**9,
          'analytic_exponential_lower_series_certifies_late_bump_upper')
    late_bound = Gstar*(F(2,10**9)+C)
    check(late_bound == F(120126743379724213215087,4999990925550628625000000000) < F(1,40000),
          'late_actual_scalar_difference_below_one_forty_thousandth')
    check(x_late < 2 < 3 and 1-T*T > 0, 'complete_late_and_tail_domain_and_positive_cosh_lower')

    gain = tuple(map(F, old['physical_gain_interval']))
    check(gain[0] < gain[1] < 0 and max(map(abs,gain)) < Gstar,
          'retained_future_negative_physical_gain_matches_safe_bound')
    transfers = {}
    order = ['central_left','central_right','late_left','late_right']
    caps = [F(-3,2), F(1), F(1,1000), F(-29,100000)]
    for j, key in enumerate(order):
        witness = old['sign_witnesses'][key]
        x = F(witness['x'])
        H0 = tuple(map(F, witness['H0']))
        H = tuple(map(F, witness['H']))
        error = F(witness['error_upper'])
        check(H[0] <= H0[0]-error and H[1] >= H0[1]+error,
              key+'_uses_actual_field_error_enlargement')
        ell_S = tuple(map(F, witness['ell_S']))
        physical_product = product_interval(gain, (1-(ell*x)**2,F(1)), H)
        check(ell_S[0] <= physical_product[0] and ell_S[1] >= physical_product[1],
              key+'_ell_S_encloses_negative_gain_cosh_actual_H_product')
        future = (ell_S[0]/ell, ell_S[1]/ell)
        epsilon = F(1,40) if j < 2 else F(1,40000)
        past = (future[0]-epsilon, future[1]+epsilon)
        correct = past[1] < caps[j] if j in [0,3] else past[0] > caps[j]
        check(correct, key+'_past_strict_sign_margin_after_ell_conversion')
        transfers[key] = {'matched_positive_x':str(x), 'difference_safe_upper':str(epsilon),
                          'future_S':[str(v) for v in future], 'past_S':[str(v) for v in past]}
    first = (-F(old['sign_witnesses']['central_right']['x']),-F(old['sign_witnesses']['central_left']['x']))
    second = (-F(old['sign_witnesses']['late_right']['x']),-F(old['sign_witnesses']['late_left']['x']))
    check(first == (F(-149999,100000),F(-74999,50000)) and second == (F(-197,100),F(-49,25)) and
          second[1] < first[0] and first[0] < first[1], 'two_distinct_ordered_open_past_zero_brackets')
    check(first[1]-first[0] == F(1,100000) and F(1,2048)/(first[1]-first[0]) == F(3125,64),
          'central_width_improvement_compares_only_prior_past_cell_bracket')

    Htail = tuple(map(F,old['tail_theorem']['actual_H_interval']))
    H0tail = tuple(map(F,old['tail_theorem']['free_H0_interval']))
    Etail = F(old['error_certificate']['tail_error'])
    check(Htail[0] <= H0tail[0]-Etail and Htail[1] >= H0tail[1]+Etail,
          'retained_whole_tail_encloses_actual_comparator_error')
    future_tail = product_interval(gain,(1-T*T,F(1)),Htail)
    future_tail = tuple(v/ell for v in future_tail)
    past_tail = (future_tail[0]-F(1,40000),future_tail[1]+F(1,40000))
    check(F(-607,10**6) < past_tail[0] < past_tail[1] < F(-296,10**6) < 0,
          'whole_past_postpulse_physical_coefficient_strictly_negative')

    delta_H_norm = ell*(2+2*C)
    check(delta_H_norm == F(40080011819991657,19999963702202514500000) < 3*ell,
          'normalized_derivative_difference_L2_with_local_and_history_terms')
    I_plus = tuple(map(F,old['quadrature']['actual_H_integral_interval']))
    delta_I = 18*ell+9*ell*ell
    check(0 < I_plus[0] < I_plus[1] < 9 and delta_I == F(18000009,10**12),
          'squared_integral_difference_retains_the_Cauchy_cross_term')
    K = tuple(map(F,old['physical_norm_prefactor_interval']))
    delta_C = K[1]/ell*delta_I
    check(0 < K[0] < K[1] and delta_C == F(1021088587031956677959217,500000000000000000),
          'same_positive_background_prefactor_scales_integral_difference_by_inverse_ell')
    C_plus = tuple(map(F,old['quadrature']['C_interval']))
    C_minus = (C_plus[0]-delta_C,C_plus[1]+delta_C)
    check(F(606,100)*10**11 < C_minus[0] < C_minus[1] < F(744,100)*10**11,
          'past_quadratic_first_signal_response_coefficient_safe_interval')
    check(all((ROOT/name).read_bytes() == raw for name,raw in raw_inputs.items()),
          'selected_input_bytes_unchanged_after_arithmetic')
    report = {
        'schema':'HDBLAST_ODD_TRACE_INDEPENDENT_CHARACTERISTIC_V1',
        'status':'PASS_CONDITIONAL_DIRECT_ODD_TRACE_FINE_ZEROS_AND_NEGATIVE_TAIL',
        'ell':str(ell),'T':str(T),'finite_guards':len(checks),'check_names':checks,
        'source_pins_sha256':digest(pin_bytes),'proof_sha256':digest((ROOT/'PROOF.md').read_bytes()),
        'verifier_sha256':digest(Path(__file__).read_bytes()),
        'constants':{k:str(v) for k,v in [('m',m),('M',M),('B',B),('Gamma',Gamma),('Gstar',Gstar)]},
        'majorants':{k:str(v) for k,v in [('r',r),('s',s),('determinant',d),('field_A',A),('boundary_velocity_V',V),('trace_history_C',C)]},
        'difference_bounds':{'derived_global_upper':str(global_bound),'global_safe_strict_upper':'1/40',
                             'derived_late_upper':str(late_bound),'late_safe_strict_upper':'1/40000',
                             'late_scaled_time_interval':['49/25','3']},
        'selected_actual_witness_transfers':transfers,
        'past_zero_brackets_in_units_of_ell':[[str(v) for v in first],[str(v) for v in second]],
        'central_bracket_width_in_units_of_ell':'1/100000',
        'central_width_reduction_against_previous_past_cell_bracket':'3125/64',
        'past_whole_tail':{'matched_positive_x_interval':['2','3'],'actual_past_S_enclosure':[str(v) for v in past_tail],
                          'safe_strict_S_interval':['-607/1000000','-296/1000000']},
        'quadratic_first_signal_response':{'delta_H_L2_upper':str(delta_H_norm),'safe_delta_H_L2_upper':str(3*ell),
            'squared_integral_difference_upper':str(delta_I),'coefficient_difference_upper':str(delta_C),
            'past_C_interval':[str(v) for v in C_minus],'safe_strict_past_C_interval':['606000000000','744000000000'],
            'cross_term_retained':True,'physical_injected_energy_interpretation':False,
            'common_exact_nonlinear_family_on_past_interval_proved':False},
        'continuous_premises':['Admitted branch, full reconstruction, smooth first development and finite propagation retained.',
                               'Zero-displacement free Dirichlet cancellation and complete damped boundary majorant proved in PROOF.md.',
                               'Same-state physical clock/normal identity and old actual future waveform/tail theorems retained.'],
        'limits':['First linear scalar coefficient in fixed model units; not an individual amplitude upper, relative error or nonlinear signal.',
                  'Two distinct past zeros established conditionally, without uniqueness or a total zero count.',
                  'Negative coefficient tail is not a globally signed angular field, radiation or physical energy.',
                  'No effective nonlinear interval, incoming source law, SI calibration or cosmological-origin evidence.'],
        'activity_this_process':{'selected_public_input_copies':6,'fresh_exact_fraction_producer_processes':1,
                                 'legacy_producer_or_checker_executions':0,'scientific_project_imports':0,
                                 'numerical_PDE_solves':0,'background_field_source_or_bump_oracle_queries':0}
    }
    with out.open('x') as stream:
        json.dump(report,stream,sort_keys=True,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps({'status':report['status'],'finite_guards':len(checks)},sort_keys=True))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError) as error:
        print(json.dumps({'status':'FAIL','witness':str(error)},sort_keys=True))
        sys.exit(2)
