"""Authenticated fresh prescribed quadratic-forcing calculations; the full U2 response remains unresolved."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.metadata
import json
from pathlib import Path, PurePosixPath
import re
import subprocess
import sys
import time

# BEGIN EMBEDDED APPROVED CONFIGURATION
CAPSULE = 'HDBLAST_QUADRATIC_SHELL_FORCING_20261015'
PYTHON_VERSION = '3.12.14'
DEPENDENCY_VERSIONS = {'sympy': '1.14.0'}
CONFIG_SHA256 = '0ecf96a7373bfdce0b0728379951518d010ed35b75430b030164f952d7a25679'
CHECKERS = [{'name': 'generator',
  'script': 'compute_quadratic_forcing.py',
  'dependencies': ['inputs/prior/INERTIAL_WAVEFORM_CELLS.csv',
                   'inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json',
                   'inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json',
                   'inputs/prior/PRIOR_PHYSICAL_GAIN_ENCLOSURE.json'],
  'arguments': ['--output-dir', 'generated'],
  'output': 'generated/QUADRATIC_FORCING_ENCLOSURE.json',
  'golden': 'computed/QUADRATIC_FORCING_ENCLOSURE.json',
  'status': 'PASS',
  'byte_equal': True,
  'additional_outputs': [['generated/QUADRATIC_FORCING_CELLS.csv',
                          'computed/QUADRATIC_FORCING_CELLS.csv']],
  'required_fields': [{'path': ['cell_data', 'rows'], 'equals': 4353},
                      {'path': ['scope', 'full_U2_response_computed'], 'equals': False}],
  'status_path': ['status'],
  'mode_fields': [],
  'checks': None,
  'embedded_checks': None,
  'check_count_path': None,
  'embedded_count_path': None,
  'check_names': [],
  'embedded_control_names': [],
  'groups': 0,
  'embedded': 0},
 {'name': 'root',
  'script': 'check_quadratic_forcing.py',
  'dependencies': ['inputs/prior/INERTIAL_WAVEFORM_CELLS.csv',
                   'inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json',
                   'inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json',
                   'inputs/prior/PRIOR_PHYSICAL_GAIN_ENCLOSURE.json',
                   'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                   'computed/QUADRATIC_FORCING_CELLS.csv'],
  'arguments': ['--output', 'checked.json'],
  'output': 'checked.json',
  'golden': 'computed/CHECK_QUADRATIC_FORCING.json',
  'status': 'PASS_QUADRATIC_MATERIAL_TENSION_ENCLOSURE',
  'checks': {'path': ['checks'], 'kind': 'names'},
  'check_count_path': ['check_count'],
  'additional_outputs': [],
  'byte_equal': False,
  'status_path': ['status'],
  'mode_fields': [],
  'embedded_checks': None,
  'required_fields': [],
  'embedded_count_path': None,
  'check_names': ['authenticated complete prior scientific inputs',
                  'quadratic material tension report identity',
                  'same registered state and full proper shell window',
                  'second derivative factorial convention',
                  'correct scalar SUM and per-copy tension sources',
                  'actual registered positive tension Hessian and cubic derivative',
                  'physical gain and proper measure normalization',
                  'exact unit three-sphere moments',
                  'orthogonal squared angular power shares',
                  'exact projected scalar and constitutive coefficients',
                  'outward arithmetic and retained positive clock inequality',
                  'complete quadratic CSV authenticated',
                  'all continuous coefficient cells and semantic columns present',
                  'every quadratic continuous cell independently reconstructed',
                  'entire inherited second moment reconstructed',
                  'entire weighted fourth moment independently reconstructed',
                  'first-order proper-time norm identity',
                  'signed forcing integral and L1 identity',
                  'fourth-moment proper-time L2 Jacobian and normalization',
                  'constitutive tension integrals retain actual Hessian',
                  'strict monopole bound on full continuous tail',
                  'background proper shell measure explicit',
                  'explicit contribution scope and unresolved second response retained'],
  'embedded_control_names': [],
  'groups': 23,
  'embedded': 0},
 {'name': 'astra',
  'script': 'astra/check_astra_material_forcing.py',
  'dependencies': [],
  'arguments': [],
  'output': None,
  'golden': 'astra/ASTRA_NORMAL.json',
  'optimized_golden': 'astra/ASTRA_OPTIMIZED.json',
  'status': 'PASS',
  'checks': {'path': ['check_names'], 'kind': 'names'},
  'check_count_path': ['checks'],
  'mode_fields': [{'path': ['python_optimized'], 'normal': False, 'optimized': True}],
  'additional_outputs': [],
  'byte_equal': False,
  'status_path': ['status'],
  'embedded_checks': None,
  'required_fields': [],
  'embedded_count_path': None,
  'check_names': ['registered sigma second derivative',
                  'registered sigma third derivative',
                  'outward sum background shell law',
                  'first scalar sum row',
                  'second scalar sum derivative',
                  'factorial coefficient is half the derivative',
                  'identical copy division',
                  'tension second derivative',
                  'lower-index stress product',
                  'inverse-metric cancellation in mixed stress',
                  'four-dimensional mixed stress trace',
                  'material Israel product',
                  'clock covariance using differentiated first junction',
                  'tension response clock transport',
                  'degree two homogeneous harmonic',
                  'degree four homogeneous harmonic',
                  'homogeneous square decomposition',
                  'sphere mean Y squared',
                  'sphere mean Y fourth',
                  'degree two norm mean',
                  'degree four norm mean',
                  'degree two zero mean',
                  'degree four zero mean',
                  'cross degree orthogonality',
                  'Y orthogonal to its square',
                  'degree two projection coefficient',
                  'degree four projection coefficient',
                  'projection squared norms sum',
                  'proper-volume cosh and width cancellation',
                  'integrated tension retains unknown measure and response',
                  'squared response removes second jet and measure at order two',
                  'tension lower rational product',
                  'tension upper rational product'],
  'embedded_control_names': [],
  'groups': 33,
  'embedded': 0},
 {'name': 'archive',
  'script': 'archive/check_registered_shell_coefficients.py',
  'dependencies': ['inputs/INPUT_PINS.json',
                   'inputs/prior/INERTIAL_WAVEFORM_CELLS.csv',
                   'inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json',
                   'inputs/prior/INERTIAL_SHELL_WAVEFORM.md',
                   'inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json',
                   'inputs/prior/PRIOR_PHYSICAL_GAIN_ENCLOSURE.json',
                   'inputs/prior/FINAL_INERTIAL_WAVEFORM_REVIEW.json',
                   'inputs/prior/FINAL_INERTIAL_WAVEFORM_REVIEW.md',
                   'inputs/reference/SCALAR_SHELL_BOUNDARY_SYSTEM.md',
                   'inputs/reference/MATERIAL_SECOND_ORDER_JUNCTIONS.md'],
  'arguments': ['--output', 'checked.json'],
  'output': 'checked.json',
  'golden': 'archive/REGISTERED_SHELL_COEFFICIENT_AUDIT_NORMAL.json',
  'optimized_golden': 'archive/REGISTERED_SHELL_COEFFICIENT_AUDIT_OPTIMIZED.json',
  'status': 'PASS',
  'checks': {'path': ['groups'], 'kind': 'names'},
  'check_count_path': ['positive_check_groups'],
  'embedded_checks': {'path': ['embedded_control_names'], 'kind': 'names'},
  'embedded_count_path': ['embedded_wrong_formula_controls_rejected'],
  'mode_fields': [{'path': ['python_optimized'], 'normal': False, 'optimized': True}],
  'additional_outputs': [],
  'byte_equal': False,
  'status_path': ['status'],
  'required_fields': [],
  'check_names': ['all_nine_selected_public_inputs_match_pins',
                  'same_registered_exact_delta_and_model_c',
                  'actual_registered_cubic_tension_and_one_copy_convention_retained',
                  'registered_polynomial_three_derivatives',
                  'same_near_positive_phi_branch_and_strict_sigma_prime_monotonicity',
                  'correlated_monotone_sigma_prime_enclosure_refines_retained_box',
                  'sigma_second_exact_affine_enclosure_matches_retained_box',
                  'retained_shell_identity_sigma_prime_equals_minus_two_p',
                  'unit_normal_second_metric_derivative',
                  'unit_normal_second_scalar_derivative',
                  'complete_Eulerian_metric_source_after_first_shell_rows',
                  'complete_Eulerian_scalar_source_after_first_shell_rows',
                  'factorial_second_scalar_tension_chain_rule',
                  'one_copy_material_scalar_source_has_correct_factor_and_sign',
                  'eight_embedded_factor_sign_and_wrong_model_controls_rejected'],
  'embedded_control_names': ['quadratic_tension_wrong_third_derivative',
                             'wrong_sign_third_derivative',
                             'doubled_third_derivative',
                             'wrong_material_scalar_sign',
                             'lost_Taylor_factor_in_material_source',
                             'scalar_source_omits_normal_terms',
                             'wrong_Eulerian_metric_v_squared_factor',
                             'wrong_Eulerian_scalar_u_squared_sign'],
  'groups': 15,
  'embedded': 8},
 {'name': 'sol',
  'script': 'sol/check_quadratic_material_forcing.py',
  'dependencies': ['sol/SOURCE_PINS.json',
                   'sol/QUADRATIC_MATERIAL_FORCING.md',
                   'inputs/INPUT_PINS.json',
                   'inputs/prior/INERTIAL_SHELL_WAVEFORM.md',
                   'inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json',
                   'inputs/prior/INERTIAL_WAVEFORM_CELLS.csv',
                   'inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json',
                   'inputs/reference/SCALAR_SHELL_BOUNDARY_SYSTEM.md',
                   'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                   'computed/QUADRATIC_FORCING_CELLS.csv'],
  'arguments': ['--check-computed',
                'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                'computed/QUADRATIC_FORCING_CELLS.csv',
                '--output',
                'checked.json'],
  'output': 'checked.json',
  'golden': 'sol/QUADRATIC_MATERIAL_NORMAL.json',
  'optimized_golden': 'sol/QUADRATIC_MATERIAL_OPTIMIZED.json',
  'status': 'PASS',
  'checks': {'path': ['check_names'], 'kind': 'names'},
  'check_count_path': ['check_group_count'],
  'mode_fields': [{'path': ['runtime', 'optimization'], 'normal': 0, 'optimized': 1}],
  'additional_outputs': [],
  'byte_equal': False,
  'status_path': ['status'],
  'embedded_checks': None,
  'required_fields': [],
  'embedded_count_path': None,
  'check_names': ['portable_input_byte_authentication_and_external_manifest_identity',
                  'ambient_harmonics_unit_S3_spectrum_exact_moments_norms_and_power_shares',
                  'material_sum_copy_factorial_stress_product_cubic_algebra_and_source_sign_counterexample',
                  'proper_time_background_measure_quadratic_quartic_norm_identities_and_Jacobians',
                  'retained_4353_continuous_cells_H2_H4_clock_norm_and_tail_projection_consequences',
                  'new_csv_all_4353_cells_reconstructed_with_exact_byte_equality',
                  'new_report_twelve_numeric_moment_parameter_integral_norm_and_tail_paths',
                  'new_report_state_angular_conventions_scope_metadata_and_strict_tail',
                  'independent_sharper_clock_integral_and_norm_inclusion'],
  'embedded_control_names': [],
  'groups': 9,
  'embedded': 0},
 {'name': 'referee',
  'script': 'referee/check_material_source.py',
  'dependencies': ['inputs/prior/INERTIAL_WAVEFORM_CELLS.csv',
                   'inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json',
                   'inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json',
                   'inputs/prior/PRIOR_PHYSICAL_GAIN_ENCLOSURE.json',
                   'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                   'computed/QUADRATIC_FORCING_CELLS.csv'],
  'arguments': ['--output', 'checked.json'],
  'output': 'checked.json',
  'golden': 'referee/MATERIAL_SOURCE_AUDIT_NORMAL.json',
  'optimized_golden': 'referee/MATERIAL_SOURCE_AUDIT_OPTIMIZED.json',
  'status': 'PASS_INDEPENDENT_QUADRATIC_MATERIAL_SOURCE_AUDIT',
  'checks': {'path': ['checks'], 'kind': 'names'},
  'check_count_path': ['check_count'],
  'mode_fields': [{'path': ['python_optimization'], 'normal': 0, 'optimized': 1}],
  'additional_outputs': [],
  'byte_equal': False,
  'status_path': ['status'],
  'embedded_checks': None,
  'required_fields': [],
  'embedded_count_path': None,
  'check_names': ['ambient_harmonic_B_Z2_H4',
                  'exact_ambient_degree_0_2_4_decomposition',
                  'normalized_sphere_volume_convention',
                  'quadrupole_second_moment',
                  'quadrupole_fourth_moment',
                  'degree_two_norm',
                  'degree_four_norm',
                  'orthogonality_of_three_degrees',
                  'full_scalar_source_squared_angular_mean',
                  'all_three_orthogonal_squared_power_shares',
                  'actual_registered_tension_second_and_third_derivatives',
                  'second_derivative_chain_rule_without_spurious_half',
                  'sum_and_per_copy_source_factors',
                  'authenticated_prior_input:prior/INERTIAL_WAVEFORM_CELLS.csv',
                  'authenticated_prior_input:prior/INERTIAL_WAVEFORM_ENCLOSURE.json',
                  'authenticated_prior_input:prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json',
                  'authenticated_prior_input:prior/PRIOR_PHYSICAL_GAIN_ENCLOSURE.json',
                  'new_report_and_complete_csv_authenticated',
                  'same_wider_state_and_derivative_normal_conventions',
                  'all_4353_source_and_target_cells_present',
                  'actual_positive_tension_hessian_interval',
                  'physical_gain_and_norm_prefactor_use_authenticated_background',
                  'all_six_harmonic_coefficient_ranges_rebuilt_at_every_cell',
                  'fourth_moment_retains_positive_cosh_minus3_weight',
                  'same_state_C_identity_not_cross_width_transfer',
                  'signed_integral_and_L1_are_exact_multiples_of_C',
                  'full_scalar_L2_squared_with_correct_clock_and_inverse_cubic_width',
                  'explicit_constitutive_term_integral_and_L2_squared',
                  'continuous_tail_monopole_negative_and_strict_bounds',
                  'no_full_second_order_or_total_stress_claim',
                  'inputs_and_outputs_unchanged_during_audit'],
  'embedded_control_names': [],
  'groups': 31,
  'embedded': 0}]
CONTROLS = [{'checker': 'root',
  'name': 'wrong_sum_source',
  'reason': 'correct scalar SUM and per-copy tension sources',
  'mechanism': 'copied input/certificate corruption',
  'failure_format': 'fail_line',
  'target': 'computed/QUADRATIC_FORCING_ENCLOSURE.json',
  'json_path': ['definitions', 'scalar_material_source'],
  'value': 'Qphi=-2 psi1^2'},
 {'checker': 'root',
  'name': 'wrong_tension_hessian',
  'reason': 'actual registered positive tension Hessian and cubic derivative',
  'mechanism': 'copied input/certificate corruption',
  'failure_format': 'fail_line',
  'target': 'computed/QUADRATIC_FORCING_ENCLOSURE.json',
  'json_path': ['parameters', 'sigma_second_interval'],
  'value': ['4', '4']},
 {'checker': 'root',
  'name': 'wrong_angular_power_shares',
  'reason': 'orthogonal squared angular power shares',
  'mechanism': 'copied input/certificate corruption',
  'failure_format': 'fail_line',
  'target': 'computed/QUADRATIC_FORCING_ENCLOSURE.json',
  'json_path': ['angular', 'orthogonal_squared_power_shares_l0_l2_l4'],
  'value': ['1/3', '1/3', '1/3']},
 {'checker': 'root',
  'name': 'wrong_scalar_l2_jacobian',
  'reason': 'fourth-moment proper-time L2 Jacobian and normalization',
  'mechanism': 'copied input/certificate corruption',
  'failure_format': 'fail_line',
  'target': 'computed/QUADRATIC_FORCING_ENCLOSURE.json',
  'json_path': ['quadrature', 'scalar_L2_squared_interval'],
  'factor': '1/1000000000000'},
 {'checker': 'root',
  'name': 'wrong_weighted_fourth_moment',
  'reason': 'entire weighted fourth moment independently reconstructed',
  'mechanism': 'copied input/certificate corruption',
  'failure_format': 'fail_line',
  'target': 'computed/QUADRATIC_FORCING_ENCLOSURE.json',
  'json_path': ['quadrature', 'cosh_minus3_H_fourth_integral_interval'],
  'factor': '2'},
 {'checker': 'root',
  'name': 'shortened_tail_domain',
  'reason': 'strict monopole bound on full continuous tail',
  'mechanism': 'copied input/certificate corruption',
  'failure_format': 'fail_line',
  'target': 'computed/QUADRATIC_FORCING_ENCLOSURE.json',
  'json_path': ['tail_theorem', 'domain_x'],
  'value': [2, 2.9]},
 {'checker': 'root',
  'name': 'claim_full_bulk_source',
  'reason': 'explicit contribution scope and unresolved second response retained',
  'mechanism': 'copied input/certificate corruption',
  'failure_format': 'fail_line',
  'target': 'computed/QUADRATIC_FORCING_ENCLOSURE.json',
  'json_path': ['scope', 'complete_bulk_coordinate_quadratic_source'],
  'value': True},
 {'checker': 'root',
  'name': 'alter_complete_csv',
  'reason': 'complete quadratic CSV authenticated',
  'mechanism': 'copied input/certificate corruption',
  'failure_format': 'fail_line',
  'target': 'computed/QUADRATIC_FORCING_CELLS.csv',
  'append': '\n'},
 {'checker': 'astra',
  'name': 'reverse_one_normal',
  'reason': 'second scalar sum derivative: nonzero symbolic residue 2*p_minus2',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject', 'reverse_one_normal'],
  'failure_format': 'json',
  'failure_path': ['failure']},
 {'checker': 'astra',
  'name': 'double_tension',
  'reason': 'second scalar sum derivative: nonzero symbolic residue -4*(eta_b*psi2 + psi1**2 + '
            'psi2)',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject', 'double_tension'],
  'failure_format': 'json',
  'failure_path': ['failure']},
 {'checker': 'astra',
  'name': 'half_hessian',
  'reason': 'second scalar sum derivative: nonzero symbolic residue 2*psi1**2',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject', 'half_hessian'],
  'failure_format': 'json',
  'failure_path': ['failure']},
 {'checker': 'astra',
  'name': 'wrong_projection',
  'reason': 'homogeneous square decomposition: nonzero symbolic residue (X1**2 + X2**2 - X3**2 - '
            'X4**2)*(X1**2 + X2**2 + X3**2 + X4**2)/12',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject', 'wrong_projection'],
  'failure_format': 'json',
  'failure_path': ['failure']},
 {'checker': 'astra',
  'name': 'omit_unknown_response',
  'reason': 'second scalar sum derivative: nonzero symbolic residue 4*psi2*(eta_b + 1)',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject', 'omit_unknown_response'],
  'failure_format': 'json',
  'failure_path': ['failure']},
 {'checker': 'astra',
  'name': 'omit_clock_transport',
  'reason': 'clock covariance using differentiated first junction: nonzero symbolic residue '
            '-8*clock_velocity*dpsi1*(eta_b + 1)',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject', 'omit_clock_transport'],
  'failure_format': 'json',
  'failure_path': ['failure']},
 {'checker': 'astra',
  'name': 'drop_metric_stress',
  'reason': 'lower-index stress product: nonzero symbolic residue -2*h1*psi1*(c*delta + 2*eta_b**2 '
            '+ 4*eta_b)',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject', 'drop_metric_stress'],
  'failure_format': 'json',
  'failure_path': ['failure']},
 {'checker': 'archive',
  'name': 'third-sign',
  'reason': 'registered_polynomial_three_derivatives',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--output', 'checked.json', '--inject-error', 'third-sign'],
  'failure_format': 'traceback',
  'exception': 'RuntimeError',
  'diagnostic_prefix': 'REJECTED: registered_polynomial_three_derivatives\n'},
 {'checker': 'archive',
  'name': 'material-factor',
  'reason': 'one_copy_material_scalar_source_has_correct_factor_and_sign',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--output', 'checked.json', '--inject-error', 'material-factor'],
  'failure_format': 'traceback',
  'exception': 'RuntimeError',
  'diagnostic_prefix': 'REJECTED: one_copy_material_scalar_source_has_correct_factor_and_sign\n'},
 {'checker': 'archive',
  'name': 'omit-normal',
  'reason': 'complete_Eulerian_scalar_source_after_first_shell_rows',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--output', 'checked.json', '--inject-error', 'omit-normal'],
  'failure_format': 'traceback',
  'exception': 'RuntimeError',
  'diagnostic_prefix': 'REJECTED: complete_Eulerian_scalar_source_after_first_shell_rows\n'},
 {'checker': 'archive',
  'name': 'metric-factor',
  'reason': 'complete_Eulerian_metric_source_after_first_shell_rows',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--output', 'checked.json', '--inject-error', 'metric-factor'],
  'failure_format': 'traceback',
  'exception': 'RuntimeError',
  'diagnostic_prefix': 'REJECTED: complete_Eulerian_metric_source_after_first_shell_rows\n'},
 {'checker': 'sol',
  'name': 'sum_vs_copy',
  'reason': 'ValueError: outward SUM source and derivative convention',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--check-computed',
                'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                'computed/QUADRATIC_FORCING_CELLS.csv',
                '--output',
                'checked.json',
                '--inject-error',
                'sum_vs_copy'],
  'failure_format': 'json',
  'failure_path': ['error']},
 {'checker': 'sol',
  'name': 'factorial_half',
  'reason': 'ValueError: outward SUM source and derivative convention',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--check-computed',
                'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                'computed/QUADRATIC_FORCING_CELLS.csv',
                '--output',
                'checked.json',
                '--inject-error',
                'factorial_half'],
  'failure_format': 'json',
  'failure_path': ['error']},
 {'checker': 'sol',
  'name': 'wrong_sign',
  'reason': 'ValueError: outward SUM source and derivative convention',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--check-computed',
                'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                'computed/QUADRATIC_FORCING_CELLS.csv',
                '--output',
                'checked.json',
                '--inject-error',
                'wrong_sign'],
  'failure_format': 'json',
  'failure_path': ['error']},
 {'checker': 'sol',
  'name': 'sphere_dimension',
  'reason': 'ValueError: unit-S3 eigenvalue',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--check-computed',
                'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                'computed/QUADRATIC_FORCING_CELLS.csv',
                '--output',
                'checked.json',
                '--inject-error',
                'sphere_dimension'],
  'failure_format': 'json',
  'failure_path': ['error']},
 {'checker': 'sol',
  'name': 'unit_norm',
  'reason': 'ValueError: Y is not unit norm',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--check-computed',
                'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                'computed/QUADRATIC_FORCING_CELLS.csv',
                '--output',
                'checked.json',
                '--inject-error',
                'unit_norm'],
  'failure_format': 'json',
  'failure_path': ['error']},
 {'checker': 'sol',
  'name': 'wrong_fourth_moment',
  'reason': 'ValueError: fourth moment',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--check-computed',
                'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                'computed/QUADRATIC_FORCING_CELLS.csv',
                '--output',
                'checked.json',
                '--inject-error',
                'wrong_fourth_moment'],
  'failure_format': 'json',
  'failure_path': ['error']},
 {'checker': 'sol',
  'name': 'omit_degree_two',
  'reason': 'ValueError: complete harmonic decomposition',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--check-computed',
                'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                'computed/QUADRATIC_FORCING_CELLS.csv',
                '--output',
                'checked.json',
                '--inject-error',
                'omit_degree_two'],
  'failure_format': 'json',
  'failure_path': ['error']},
 {'checker': 'sol',
  'name': 'wrong_harmonic_constant',
  'reason': 'ValueError: ambient harmonicity',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--check-computed',
                'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                'computed/QUADRATIC_FORCING_CELLS.csv',
                '--output',
                'checked.json',
                '--inject-error',
                'wrong_harmonic_constant'],
  'failure_format': 'json',
  'failure_path': ['error']},
 {'checker': 'sol',
  'name': 'wrong_norm_jacobian',
  'reason': 'ValueError: quadratic observable Jacobian',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--check-computed',
                'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                'computed/QUADRATIC_FORCING_CELLS.csv',
                '--output',
                'checked.json',
                '--inject-error',
                'wrong_norm_jacobian'],
  'failure_format': 'json',
  'failure_path': ['error']},
 {'checker': 'sol',
  'name': 'cancel_cosh_in_quartic',
  'reason': 'ValueError: quartic cosh weight survives',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--check-computed',
                'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                'computed/QUADRATIC_FORCING_CELLS.csv',
                '--output',
                'checked.json',
                '--inject-error',
                'cancel_cosh_in_quartic'],
  'failure_format': 'json',
  'failure_path': ['error']},
 {'checker': 'sol',
  'name': 'tension_remainder_factor',
  'reason': 'ValueError: exact constitutive cubic remainder',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--check-computed',
                'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                'computed/QUADRATIC_FORCING_CELLS.csv',
                '--output',
                'checked.json',
                '--inject-error',
                'tension_remainder_factor'],
  'failure_format': 'json',
  'failure_path': ['error']},
 {'checker': 'sol',
  'name': 'response_sign_from_source',
  'reason': 'ValueError: negative source need not give a negative response',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--check-computed',
                'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                'computed/QUADRATIC_FORCING_CELLS.csv',
                '--output',
                'checked.json',
                '--inject-error',
                'response_sign_from_source'],
  'failure_format': 'json',
  'failure_path': ['error']},
 {'checker': 'sol',
  'name': 'computed_cell',
  'reason': 'ValueError: new coefficient CSV differs from independent cell reconstruction',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--check-computed',
                'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                'computed/QUADRATIC_FORCING_CELLS.csv',
                '--output',
                'checked.json',
                '--inject-error',
                'computed_cell'],
  'failure_format': 'json',
  'failure_path': ['error']},
 {'checker': 'sol',
  'name': 'computed_J4',
  'reason': 'ValueError: new report numeric discrepancy: '
            'quadrature.cosh_minus3_H_fourth_integral_interval',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--check-computed',
                'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                'computed/QUADRATIC_FORCING_CELLS.csv',
                '--output',
                'checked.json',
                '--inject-error',
                'computed_J4'],
  'failure_format': 'json',
  'failure_path': ['error']},
 {'checker': 'sol',
  'name': 'computed_norm',
  'reason': 'ValueError: new report numeric discrepancy: quadrature.scalar_L2_squared_interval',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--check-computed',
                'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                'computed/QUADRATIC_FORCING_CELLS.csv',
                '--output',
                'checked.json',
                '--inject-error',
                'computed_norm'],
  'failure_format': 'json',
  'failure_path': ['error']},
 {'checker': 'sol',
  'name': 'computed_tail',
  'reason': 'ValueError: new report numeric discrepancy: '
            'tail_theorem.ell2_scalar_monopole_interval',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--check-computed',
                'computed/QUADRATIC_FORCING_ENCLOSURE.json',
                'computed/QUADRATIC_FORCING_CELLS.csv',
                '--output',
                'checked.json',
                '--inject-error',
                'computed_tail'],
  'failure_format': 'json',
  'failure_path': ['error']}]
# END EMBEDDED APPROVED CONFIGURATION

def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def authenticate(root, expected):
    """Stdlib bootstrap, before any capsule module import or execution."""
    root = Path(root).absolute()
    require(not any(p.is_symlink() for p in [root, *root.parents]), 'Real capsule path required')
    require(re.fullmatch('[0-9a-f]{64}', expected) is not None, 'Complete external manifest SHA256 required')
    manifest = root / 'MANIFEST.json'
    require(manifest.is_file() and not manifest.is_symlink(), 'Real regular manifest required')
    raw = manifest.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == expected, 'External manifest pin differs')
    data = json.loads(raw)
    require(data['capsule'] == root.name == CAPSULE, 'Capsule directory identity differs')
    records = {}
    for row in data['files']:
        name = row['path']
        relative = PurePosixPath(name)
        require(name and '\\' not in name and not relative.is_absolute() and '..' not in relative.parts
                and relative.as_posix() == name and name not in records and name != 'MANIFEST.json',
                'Noncanonical, reserved or duplicate source path')
        require(type(row['bytes']) is int and row['bytes'] >= 0
                and re.fullmatch('[0-9a-f]{64}', row['sha256']) is not None, 'Invalid source record')
        records[name] = row
        read_source(root, records, name)
    actual = set()
    for path in root.rglob('*'):
        require(not path.is_symlink(), 'Symlink in frozen tree')
        if path.is_file():
            actual.add(path.relative_to(root).as_posix())
        else:
            require(path.is_dir(), 'Special file in frozen tree')
    require(records and actual == set(records) | {'MANIFEST.json'}, 'Frozen tree contains missing or extra files')
    return {'status': 'PASS', 'files_verified': len(records), 'manifest_sha256': expected}, records


def read_source(root, inventory, name):
    require(name in inventory, 'Missing declared source dependency: ' + name)
    path = root / name
    require(not any(p.is_symlink() for p in [path, *path.parents]) and path.is_file()
            and path.resolve().is_relative_to(root.resolve()), 'Real source inside capsule required: ' + name)
    raw = path.read_bytes()
    item = inventory[name]
    require(len(raw) == item['bytes'] and hashlib.sha256(raw).hexdigest() == item['sha256'],
            'Frozen bytes differ: ' + name)
    return raw


def read_generated(path, work):
    require(path.is_file() and not any(p.is_symlink() for p in [path, *path.parents])
            and path.resolve().is_relative_to(work.resolve()),
            'Fresh generated output must be a real regular file inside its workspace')
    return path.read_bytes()


def record(path, base):
    raw = read_generated(path, base)
    return {'path': path.relative_to(base).as_posix(), 'bytes': len(raw),
            'sha256': hashlib.sha256(raw).hexdigest()}


def json_equal(left, right):
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return left.keys() == right.keys() and all(json_equal(left[k], right[k]) for k in left)
    if isinstance(left, list):
        return len(left) == len(right) and all(json_equal(a, b) for a, b in zip(left, right))
    return left == right


def at_path(data, path):
    for key in path:
        data = data[key]
    return data


def inventory_names(data, specification):
    if specification is None:
        return []
    value = at_path(data, specification['path'])
    if specification['kind'] == 'names':
        require(isinstance(value, list) and all(isinstance(x, str) and x for x in value), 'Named list required')
        names = value
    else:
        require(specification['kind'] == 'true_map' and isinstance(value, dict)
                and all(isinstance(k, str) and k and v is True for k, v in value.items()), 'Passing named boolean map required')
        names = list(value)
    require(len(names) == len(set(names)), 'Duplicate named group/control')
    return names


def positive_inventory(checker, data, mode):
    require(at_path(data, checker['status_path']) == checker['status'], 'Positive status differs: '+checker['name'])
    for field in checker['mode_fields']:
        value = at_path(data, field['path'])
        expected = field[mode]
        require(type(value) is type(expected) and value == expected, 'Execution-mode metadata differs')
    checks = inventory_names(data, checker['checks'])
    embedded = inventory_names(data, checker['embedded_checks'])
    require(checks == checker['check_names'] and embedded == checker['embedded_control_names'],
            'Exact approved group/control identities differ')
    for item in checker['required_fields']:
        value = at_path(data, item['path'])
        require(json_equal(value, item['equals']), 'Required positive output field differs')
    if checker.get('check_count_path') is not None:
        require(type(at_path(data, checker['check_count_path'])) is int
                and at_path(data, checker['check_count_path']) == len(checks), 'Reported group count differs')
    if checker.get('embedded_count_path') is not None:
        require(type(at_path(data, checker['embedded_count_path'])) is int
                and at_path(data, checker['embedded_count_path']) == len(embedded), 'Reported embedded count differs')
    return len(checks), len(embedded)


def validate_failure_diagnostic(control, diagnostic):
    if control['failure_format'] == 'fail_line':
        require(diagnostic == 'FAIL: '+control['reason']+'\n', 'Unrelated or ambiguous failure diagnostic')
    elif control['failure_format'] == 'json':
        failure = json.loads(diagnostic)
        require(at_path(failure, control.get('status_path', ['status'])) == control.get('status', 'FAIL')
                and at_path(failure, control['failure_path']) == control['reason'],
                'JSON failure does not contain the exact expected named reason')
    else:
        require(control['failure_format'] == 'traceback'
                and diagnostic.startswith(control.get('diagnostic_prefix', '')+'Traceback (most recent call last):\n')
                and diagnostic.endswith(control['exception']+': '+control['reason']+'\n')
                and diagnostic.splitlines()[-1] == control['exception']+': '+control['reason'],
                'Traceback did not end with the exact named expected failure')


def execute(root, inventory, output, checker, mode, control=None):
    work = output / mode / checker['name'] / (control['name'] if control else 'positive')
    work.mkdir(parents=True)
    copied = []
    for name in [checker['script'], *checker['dependencies']]:
        raw = read_source(root, inventory, name)
        target = work / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(raw)
        require(target.read_bytes() == raw, 'Fresh copied bytes differ: '+name)
        copied.append({'source': inventory[name], 'copy': record(target, output)})
    mutation = None
    if control and 'target' in control:
        target = work / control['target']
        require(control['target'] != checker['script']
                and control['target'] in checker['dependencies'], 'Corruption must target a declared copied dependency, preserving the checker')
        original = target.read_bytes()
        if 'append' in control:
            mutated = original+control['append'].encode()
        elif 'json_path' in control:
            document = json.loads(original)
            parent = document
            for key in control['json_path'][:-1]:
                parent = parent[key]
            key = control['json_path'][-1]
            if 'factor' in control:
                value = parent[key]
                parent[key] = [str(F(v)*F(control['factor'])) for v in value] if isinstance(value, list) else str(F(value)*F(control['factor']))
            else:
                parent[key] = control['value']
            mutated = (json.dumps(document, indent=2)+'\n').encode()
        else:
            require(original.decode().count(control['before']) == 1, 'Formula replacement must be unique')
            mutated = original.decode().replace(control['before'], control['after'], 1).encode()
        require(mutated != original, 'Corruption must alter copied bytes')
        target.write_bytes(mutated)
        mutation = {'source': inventory[control['target']], 'mutated': record(target, output),
                    'recipe': control, 'independent_checker_unchanged': True}
    arguments = control.get('arguments', checker.get('arguments', [])) if control else checker.get('arguments', [])
    flags = ['-I', '-B']+(['-O'] if mode == 'optimized' else [])
    result_path = work / checker['output'] if checker['output'] else work/'execution.log'
    if checker['output']:
        require(not result_path.exists(), 'Stale output present before execution')
    for path, _ in checker.get('additional_outputs', []):
        require(not (work/path).exists(), 'Stale additional generated output present')
    started = time.monotonic()
    with (work/'execution.log').open('w', encoding='utf-8') as log:
        result = subprocess.run([sys.executable, *flags, str(work/checker['script']), *arguments], cwd=work,
                                stdout=log, stderr=subprocess.STDOUT, timeout=300)
    if control:
        diagnostic = read_generated(work/'execution.log', work).decode()
        require(result.returncode == 1, 'Corruption did not fail with exit one: '+control['name'])
        require(not checker['output'] or not (result_path.exists() or result_path.is_symlink()), 'A corrupted process emitted a passing output file')
        validate_failure_diagnostic(control, diagnostic)
        data = {'status': 'FAIL', 'failure': control['reason'], 'exit_code': result.returncode,
                'passing_output_absent': True, 'complete_diagnostic': diagnostic,
                'result_origin': 'Runner-parsed exact named failure; original diagnostic retained in execution.log'}
        raw = (json.dumps(data, indent=2)+'\n').encode()
    else:
        require(result.returncode == 0 and result_path.is_file(), 'Successful fresh result missing: '+checker['name'])
        raw = read_generated(result_path, work)
        data = json.loads(raw)
    fresh = work/'fresh.json'
    fresh.write_bytes(raw)
    row = {'kind': 'external_corruption' if control else 'positive', 'checker': checker['name'], 'mode': mode,
           'script': checker['script'], 'flags': flags, 'arguments': arguments, 'exit_code': result.returncode,
           'wall_seconds': round(time.monotonic()-started, 3), 'copied_inputs': copied,
           'copied_inputs_snapshot': 'Authenticated initial copies before any separately recorded corruption',
           'log': record(work/'execution.log', output), 'output': record(fresh, output), 'additional_outputs': [],
           'output_origin': 'Exact named failure diagnostic' if control else checker['output'] or 'complete stdout JSON'}
    require((work/checker['script']).read_bytes() == read_source(root, inventory, checker['script']),
            'Independent executed checker script changed')
    if control:
        row.update({'mutation': mutation, 'control': control['name'],
            'control_mechanism': control['mechanism'],
            'expected_reason': control['reason'], 'reported_failure': data['failure'],
            'expected_reason_detected': True, 'expected_failure_exit_code': 1})
    else:
        golden_name = checker.get('optimized_golden', checker['golden']) if mode == 'optimized' else checker['golden']
        golden_raw = read_source(root, inventory, golden_name)
        require(json_equal(data, json.loads(golden_raw)), 'Complete scientific JSON differs from frozen golden: '+checker['name'])
        if checker.get('byte_equal'):
            require(raw == golden_raw, 'Complete generated JSON bytes differ: '+checker['name'])
        groups, embedded = positive_inventory(checker, data, mode)
        require(groups == checker['groups'] and embedded == checker['embedded'], 'Positive finite inventory differs')
        for path, golden in checker.get('additional_outputs', []):
            fresh_path = work/path
            require(read_generated(fresh_path, work) == read_source(root, inventory, golden),
                    'Complete generated data bytes differ: '+path)
            row['additional_outputs'].append({'output': record(fresh_path, output), 'golden': inventory[golden],
                                              'complete_bytes_equal': True})
        row.update({'golden': inventory[golden_name], 'scientific_output_matches': True,
                    'complete_generated_json_bytes_equal': checker.get('byte_equal', False),
                    'finite_positive_check_groups': groups, 'embedded_wrong_formula_controls': embedded,
                    'finite_positive_check_names': checker['check_names'],
                    'embedded_control_names': checker['embedded_control_names'],
                    'ignored_golden_fields': []})
    return row


def roster_counts():
    return {'successful_checker_runs': 2*len(CHECKERS),
            'external_corruption_processes_rejected': 2*len(CONTROLS),
            'embedded_wrong_formula_controls_rejected': 2*sum(c['embedded'] for c in CHECKERS),
            'finite_positive_check_groups_each_mode': sum(c['groups'] for c in CHECKERS)}


def validate_declared_roster(declared):
    counts = roster_counts()
    expected_positive = [{'name': c['name'], 'script': c['script'],
                          'finite_positive_groups_each_mode': c['groups'],
                          'embedded_controls_each_mode': c['embedded'],
                          'finite_positive_check_names': c['check_names'],
                          'embedded_control_names': c['embedded_control_names']} for c in CHECKERS]
    expected_controls = [{'checker': c['checker'], 'name': c['name'], 'expected_reason': c['reason'],
                          'mechanism': c['mechanism']} for c in CONTROLS]
    require(declared['schema'] == 'hdblast-replay-roster-v1' and declared['status'] == 'APPROVED_FOR_FROZEN_REPLAY'
            and declared['capsule'] == CAPSULE and declared['modes'] == ['normal', 'optimized']
            and declared['counts'] == counts and declared['positive_processes'] == expected_positive
            and declared['external_controls'] == expected_controls,
            'Frozen reviewed roster differs from actual runner')
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest-sha256', required=True)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    require(CONFIG_SHA256 is not None and CHECKERS, 'An approved embedded runner configuration is required')
    root = Path(__file__).absolute().parent
    authentication, inventory = authenticate(root, args.manifest_sha256)
    verified = subprocess.run([sys.executable, '-I', '-S', '-B', str(root/'verify_manifest.py'),
                               '--expected', args.manifest_sha256], capture_output=True, text=True, timeout=30)
    require(verified.returncode == 0 and json.loads(verified.stdout) == authentication,
            'Independent source verifier disagrees with bootstrap')
    declared = json.loads(read_source(root, inventory, 'REPLAY_ROSTER.json'))
    counts = validate_declared_roster(declared)
    requested = args.output_dir.absolute()
    require(not any(p.is_symlink() for p in [requested, *requested.parents]), 'Real output path required')
    output = requested.resolve()
    require(not output.is_relative_to(root.resolve()), 'Output must be outside frozen source')
    require(not output.exists() or (output.is_dir() and not any(output.iterdir())), 'New or empty output required')
    require(sys.version.split()[0] == PYTHON_VERSION, 'Tested Python version differs')
    dependency_versions = {name: importlib.metadata.version(name) for name in DEPENDENCY_VERSIONS}
    require(dependency_versions == DEPENDENCY_VERSIONS, 'Pinned dependency version differs')
    output.mkdir(parents=True, exist_ok=True)
    rows = []
    checkers_by_name = {c['name']: c for c in CHECKERS}
    require(len(checkers_by_name) == len(CHECKERS), 'Duplicate checker identity')
    for mode in declared['modes']:
        for checker in CHECKERS:
            rows.append(execute(root, inventory, output, checker, mode))
        for control in CONTROLS:
            rows.append(execute(root, inventory, output, checkers_by_name[control['checker']], mode, control))
    agreements = []
    for checker in CHECKERS:
        pair = [r for r in rows if r['kind'] == 'positive' and r['checker'] == checker['name']]
        require(len(pair) == 2 and [r['mode'] for r in pair] == declared['modes'], 'Complete mode pair required')
        data = [json.loads((output/r['output']['path']).read_bytes()) for r in pair]
        normalized = [field['path'] for field in checker['mode_fields']]
        for item in data:
            for path in normalized:
                parent = item
                for key in path[:-1]:
                    parent = parent[key]
                del parent[path[-1]]
        require(json_equal(data[0], data[1]), 'Normal and optimized scientific outputs differ')
        agreements.append({'checker': checker['name'], 'scientific_json_equal': True, 'only_normalized_fields': normalized})
    actual = {'successful_checker_runs': sum(r['kind'] == 'positive' for r in rows),
        'external_corruption_processes_rejected': sum(r['kind'] == 'external_corruption' for r in rows),
        'embedded_wrong_formula_controls_rejected': sum(r.get('embedded_wrong_formula_controls', 0) for r in rows)}
    require(all(actual[k] == counts[k] for k in actual), 'Incomplete process inventory')
    groups = {m: sum(r.get('finite_positive_check_groups', 0) for r in rows if r['mode'] == m)
              for m in declared['modes']}
    require(all(v == counts['finite_positive_check_groups_each_mode'] for v in groups.values()), 'Per-mode group count differs')
    require(authenticate(root, args.manifest_sha256) == (authentication, inventory), 'Frozen source changed during replay')
    receipt = {'schema': 'hdblast-quadratic-shell-forcing-fresh-replay-v1', 'status': 'PASS', 'capsule': CAPSULE,
        'authentication': authentication, 'runner': inventory['reproduce.py'],
        'integrity_verifier': inventory['verify_manifest.py'], 'requirements': inventory['requirements.txt'],
        'roster': inventory['REPLAY_ROSTER.json'], 'python': {'version': sys.version, 'executable': sys.executable},
        'dependency_versions': dependency_versions, 'embedded_configuration_sha256': CONFIG_SHA256,
        **counts, 'finite_positive_check_groups_by_mode': groups,
        'mode_agreements': agreements, 'steps': rows,
        'scope': 'Fresh finite arithmetic for prescribed quadratic material shell forcing and known constitutive stress terms only. The full second-order U2 response remains unresolved quantitatively. Complete generated JSON and CSV reproduce byte for byte. All golden JSON fields are compared; only individually verified execution-mode metadata is normalized across modes. Executions do not prove analytic PDE/exact-family hypotheses, compute a finite-amplitude nonlinear waveform, set an effective nonlinear amplitude or observational scale, or replace external peer review.'}
    (output/'REPRODUCTION_RECEIPT.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', **counts}, indent=2))


if __name__ == '__main__':
    main()
