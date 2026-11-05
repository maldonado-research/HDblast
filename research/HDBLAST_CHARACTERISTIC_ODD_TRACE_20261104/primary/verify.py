"""Exact selected-source and scalar guards for the velocity-only trace proof.

No scientific project import, old campaign, PDE solver or numerical oracle.
The continuous proof and retained old theorems remain stated premises.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent


def need(ok, witness):
    if not ok:
        raise ValueError(witness)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    out = args.output
    need(out.is_absolute(), 'ABSOLUTE_OUTPUT_REQUIRED')
    need(all(not p.is_symlink() for p in [out, *out.parents]), 'REGULAR_OUTPUT_ANCESTRY')
    out = out.resolve()
    need(not out.exists() and out.parent.is_dir() and not out.is_relative_to(ROOT), 'FRESH_EXTERNAL_OUTPUT')
    checks = []

    def gate(ok, name):
        need(ok, name)
        need(name not in checks, 'UNIQUE_GUARD_NAME')
        checks.append(name)

    pins_path = ROOT / 'SOURCE_PINS.json'
    pins_raw = pins_path.read_bytes()
    pins = json.loads(pins_raw)
    gate(pins['schema'] == 'HDBLAST_ODD_TRACE_SELECTED_SOURCE_PINS_V1' and len(pins['files']) == 5,
         'five_selected_public_source_bindings')
    sources = {}
    for name, row in pins['files'].items():
        rel = Path(name)
        need(not rel.is_absolute() and rel.as_posix() == name and '..' not in rel.parts,
             'SAFE_SELECTED_SOURCE_PATH')
        path = ROOT / rel
        need(path.is_file() and all(not q.is_symlink() for q in [path, *path.parents]),
             'REGULAR_SELECTED_SOURCE')
        raw = path.read_bytes()
        blob = hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest()
        gate(len(raw) == row['bytes'] and digest(raw) == row['sha256'] and blob == row['git_blob_sha1'],
             'selected_source_authenticated:' + name)
        gate(row['source_commit'] == '4bfc15d0dee151aab288f33f34cf1ffd32e4c870' and
             row['source_repository_path'].startswith('research/HDBLAST_PAST_FUTURE_SCALAR_TRANSFER_20261104/') and
             '/workspace/' not in row['source_repository_path'],
             'selected_public_origin:' + name)
        sources[name] = raw
    contact = json.loads(sources['inputs/CONTACT_NORMAL.json'])
    old = json.loads(sources['inputs/OLD_FUTURE_REPORT.json'])
    gate(contact['status'] == 'PASS_CONDITIONAL_PAST_SHELL_CONTACT_SCALARS' and
         old['schema'] == 'hdblast-inertial-shell-waveform-v1' and old['status'] == 'PASS',
         'retained_actual_contact_and_future_theorem_reports')
    ell = F(1, 1000000)
    T = 3*ell
    m, M, B, Gamma, a0 = F(1,1000), F(103000), F(300), F(300), F(35,4)
    gate(F(old['state']['ell']) == ell and F(old['state']['T']) == T and
         old['state']['initial_velocity'] == 'exp(-z) Q(z/ell)',
         'same_registered_width_velocity_and_horizon')
    gate(old['mathematical_inputs']['shell']['sha256'] ==
         '8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8' and
         old['mathematical_inputs']['background']['sha256'] ==
         'fd6fd277176033703a4e75a2f66ab94894318aac49c83264ef979e9a27f0a02c',
         'retained_same_shell_and_wide_background_bindings')
    mu_lo, mu_hi = map(F, contact['raw_derived_boundary_intervals']['mu_b'])
    gamma_lo, gamma_hi = map(F, contact['raw_derived_boundary_intervals']['C_b_over_two'])
    gate(m < mu_lo <= mu_hi < F(1,400), 'coarse_positive_endpoint_mass')
    gate(0 < gamma_lo <= gamma_hi < Gamma, 'coarse_actual_spatial_coefficient')
    gate(gamma_hi + F(5,2)*mu_hi < B, 'positive_constant_spring_upper')
    gate(gamma_lo - F(25,4)*mu_hi > 0, 'actual_spring_stays_positive')
    gate(0 < F(contact['energy']['all_real_W_lower']) <=
         F(contact['energy']['all_real_W_upper']) < M, 'coarse_actual_bulk_potential')
    Gstar = F(25000,129**3)
    gate(F(contact['physical_scalar_trace']['rho_b_lower_guard']) == 129 and
         F(contact['physical_scalar_trace']['minus_p_b_lower_guard']) == F(3,25000) and
         Gstar < F(3,250), 'physical_scalar_gain_guard')
    gate(contact['physical_scalar_trace']['background_proper_clock'] == 'tau=rho_b*t' and
         contact['physical_scalar_trace']['normalization_coefficients'] ==
         {'scalar_factor':'3','q_factor':'5/2','graph_ell_power':-2},
         'retained_physical_clock_and_reconstruction')

    r = M*T*T/2
    s = B*T*T/(2*m)
    d = (1-r)*(1-s)-T*r/m
    gate(0 < r < 1 and 0 < s < 1 and d > 0, 'positive_velocity_mass_majorant_inverse')
    A = ell*(1-s+T/m)/d
    V = ell/(m*d)
    gate((1-r)*A-T*V == ell, 'field_majorant_row_exact')
    gate((1-s)*V-r*A/m == ell/m, 'boundary_speed_majorant_row_exact')
    I = F(3,4)*M*ell**2 + M*V*T*T/4 + M*M*A*T**3/6
    Ev = (T*I/m+a0*V*T**4/4)/(1-s)
    R = I+(1+Gamma*T)*Ev
    memory = (1+Gamma*T)*ell/m
    gate(0 < I < F(7749,10**11), 'global_velocity_only_incoming_residual')
    gate(0 < Ev < F(2325,10**13), 'boundary_error_keeps_complete_prior_history')
    gate(0 < R < F(7772,10**11), 'actual_normal_trace_comparator_error')
    gate(memory == F(10009,10**7), 'constant_inertia_memory_mass_factor')
    global_derived = 2*Gstar*(1+memory+R)
    global_safe = F(1,40)
    gate(global_derived < global_safe, 'uniform_actual_past_future_difference_below_one_fortieth')
    late_start = F(49,25)
    exponent = 4-1/((late_start-1)*(2-late_start))
    gate(exponent == -F(529,24) and exponent < -22, 'late_bump_exponent_analytic')
    gate(F(8,3)**22 > 10**9, 'positive_exponential_series_integer_guard')
    late_derived = 2*Gstar*(F(1,10**9)+memory+R)
    late_safe = F(1,40000)
    gate(late_derived < late_safe, 'late_actual_past_future_difference_below_one_forty_thousandth')

    point_transfers = {}
    for name, row in old['sign_witnesses'].items():
        x = F(row['x'])
        low, high = map(F, row['ell_S'])
        bound = late_safe if x >= late_start else global_safe
        future = [low/ell, high/ell]
        past = [future[0]-bound, future[1]+bound]
        gate(future[0] <= future[1] and past[0] < past[1], 'outward_actual_ell_S_conversion:' + name)
        point_transfers[name] = {'x':str(x), 'difference_upper':str(bound),
                                 'future_S':list(map(str,future)), 'past_S':list(map(str,past))}
    gate(F(point_transfers['central_left']['past_S'][1]) < 0 <
         F(point_transfers['central_right']['past_S'][0]), 'two_central_actual_past_signs_survive')
    gate(F(point_transfers['late_right']['past_S'][1]) < 0 <
         F(point_transfers['late_left']['past_S'][0]), 'two_late_actual_past_signs_survive')
    central = [-F(old['sign_witnesses']['central_right']['x']),
               -F(old['sign_witnesses']['central_left']['x'])]
    late = [-F(old['sign_witnesses']['late_right']['x']),
            -F(old['sign_witnesses']['late_left']['x'])]
    gate(central == [-F(149999,100000),-F(74999,50000)] and
         late == [-F(197,100),-F(49,25)] and late[0]<late[1]<central[0]<central[1],
         'disjoint_open_past_zero_brackets')
    central_width = central[1]-central[0]
    improvement = F(1,2048)/central_width
    gate(central_width == F(1,100000) and improvement == F(3125,64),
         'central_past_bracket_width_comparison_only')
    gain_lo, gain_hi = map(F, old['physical_gain_interval'])
    tail_h_lo, tail_h_hi = map(F, old['tail_theorem']['actual_H_interval'])
    gate(gain_lo <= gain_hi < 0 < tail_h_lo <= tail_h_hi and 1-T*T>0,
         'retained_actual_continuous_tail_and_negative_gain')
    tail_future = [gain_lo*tail_h_hi/ell, gain_hi*tail_h_lo*(1-T*T)/ell]
    tail_past = [tail_future[0]-late_safe,tail_future[1]+late_safe]
    gate(-F(1,1500) < tail_past[0] <= tail_past[1] < -F(1,4000),
         'uniform_actual_past_negative_tail_with_readable_strict_margins')

    deltaH_sup = 2*ell*(1+memory+R)
    deltaH_norm_square = 2*deltaH_sup**2
    epsilon_H = F(3,10**6)
    gate(deltaH_norm_square < epsilon_H**2, 'incoming_trace_to_H_difference_norm_with_front_support')
    I_lo, I_hi = map(F, old['quadrature']['actual_H_integral_interval'])
    gate(0 < I_lo <= I_hi < 9, 'retained_actual_future_H_norm_below_three')
    integral_error = 6*epsilon_H+epsilon_H**2
    gate(integral_error == F(18000009,10**12), 'quadratic_integral_cross_term_retained')
    past_I = [I_lo-integral_error,I_hi+integral_error]
    gate(0 < past_I[0] <= past_I[1], 'past_H_squared_integral_positive_interval')
    K_lo,K_hi = map(F,old['physical_norm_prefactor_interval'])
    gate(0 < K_lo <= K_hi, 'retained_positive_proper_time_norm_prefactor')
    C_past = [K_lo/ell*past_I[0],K_hi/ell*past_I[1]]
    gate(F(606,100)*10**11 < C_past[0] <= C_past[1] < F(744,100)*10**11,
         'past_first_field_quadratic_norm_strict_interval')
    gate(all((ROOT/name).read_bytes()==raw for name,raw in sources.items()) and
         pins_path.read_bytes()==pins_raw, 'all_selected_sources_unchanged_after_arithmetic')

    result = {
        'schema':'HDBLAST_ODD_TRACE_CHARACTERISTIC_COMPARISON_V1',
        'status':'PASS_CONDITIONAL_SHARP_ODD_TRACE_AND_PAST_WAVEFORM',
        'client_date':'2026-11-04','client_timezone':'America/Los_Angeles',
        'ell':str(ell),'T':str(T),'selected_public_inputs_authenticated':len(sources),
        'source_pins_sha256':digest(pins_raw),'proof_sha256':digest((ROOT/'PROOF.md').read_bytes()),
        'verifier_sha256':digest(Path(__file__).read_bytes()),
        'finite_guards':len(checks),'check_names':checks,
        'majorants':{k:str(v) for k,v in {'m':m,'M':M,'B':B,'Gamma':Gamma,'r':r,'s':s,'d':d,
                    'velocity_mass':ell,'A':A,'V':V,'Istar':I,'error_velocity':Ev,'trace_error':R,
                    'kernel_memory':memory,'physical_gain_upper':Gstar}.items()},
        'difference':{'derived_global_upper':str(global_derived),'safe_strict_global_upper':str(global_safe),
                      'derived_late_upper':str(late_derived),'safe_strict_late_upper':str(late_safe),
                      'late_scaled_interval':[str(late_start),'3'],'prefront_difference_zero':True,
                      'quantity':'S_current(-t)-S_current(+t), paired proper clocks -rho_b*t,+rho_b*t'},
        'selected_point_transfers':point_transfers,
        'past_zero_brackets_in_units_of_ell':list(map(lambda row:list(map(str,row)),[central,late])),
        'at_least_two_distinct_past_amplitude_zeros':True,
        'central_width_in_units_of_ell':str(central_width),
        'width_reduction_factor_against_previous_ell_over_2048_past_bracket':str(improvement),
        'tail':{'scaled_positive_elapsed_interval':['2','3'],'future_actual_S':list(map(str,tail_future)),
                'past_actual_S':list(map(str,tail_past)),
                'readable_strict_past_interval':['-1/1500','-1/4000']},
        'quadratic_norm':{'deltaH_uniform_upper':str(deltaH_sup),
                          'deltaH_L2_square_upper':str(deltaH_norm_square),
                          'safe_deltaH_L2_strict_upper':str(epsilon_H),
                          'H_squared_integral_difference_upper':str(integral_error),
                          'past_H_squared_integral_interval':list(map(str,past_I)),
                          'past_C_interval':list(map(str,C_past)),
                          'readable_strict_past_C_interval':['606000000000','744000000000'],
                          'physical_energy_claimed':False,'finite_amplitude_family_law_claimed':False},
        'retained_premises':['Analytic branch and complete smooth linear development/reconstruction/uniqueness',
                             'Actual old same-state point, tail, gain and integral enclosure theorems'],
        'limits':['Absolute first-coefficient differences in model units; no individual amplitude below the ceilings or relative-error law.',
                  'No uniqueness/total zero count, physical energy, finite nonlinear amplitude interval or cosmological-origin evidence.',
                  'Selected-source hashes and new Fraction guards do not rerun the old analytic/numerical waveform campaign.',
                  'Prior conditional formal complete future second Taylor constructions are preserved.'],
        'new_activity':{'this_process_exact_fraction_postprocessors':1,'legacy_producers_or_checkers':0,
                        'background_ODE_or_trajectory_queries':0,'numerical_PDE_solves':0,
                        'field_source_bump_oracles':0,'scientific_project_imports':0}
    }
    raw = (json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n').encode()
    with out.open('xb') as stream:
        stream.write(raw)
    print(json.dumps({'status':result['status'],'finite_guards':len(checks)},sort_keys=True))


if __name__ == '__main__':
    try:
        main()
    except (ValueError,OSError,KeyError,TypeError) as error:
        print(json.dumps({'status':'FAIL','witness':str(error)},sort_keys=True))
        sys.exit(2)
