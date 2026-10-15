"""Authenticated fresh coordinate shell-source calculations; the full U2 response remains unresolved."""
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
CAPSULE = 'HDBLAST_COORDINATE_SHELL_SOURCES_20261015'
PYTHON_VERSION = '3.12.14'
DEPENDENCY_VERSIONS = {}
CONFIG_SHA256 = 'a163290f91ce4e95f6da78242d7d9f71a0d779c9275ed8f3c0cdd3c74163bb36'
CHECKERS = [{'name': 'generator',
  'script': 'compute_coordinate_sources.py',
  'dependencies': ['inputs/INPUT_PINS.json',
                   'inputs/prior/INERTIAL_WAVEFORM_CELLS.csv',
                   'inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json',
                   'inputs/prior/INERTIAL_SHELL_WAVEFORM.md',
                   'inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json',
                   'inputs/prior/PRIOR_PHYSICAL_GAIN_ENCLOSURE.json',
                   'inputs/prior/FINAL_INERTIAL_WAVEFORM_REVIEW.json',
                   'inputs/prior/FINAL_INERTIAL_WAVEFORM_REVIEW.md',
                   'inputs/reference/SCALAR_SHELL_BOUNDARY_SYSTEM.md',
                   'inputs/reference/MATERIAL_SECOND_ORDER_JUNCTIONS.md',
                   'inputs/reference/BACKGROUND_DERIVATIVE_BOUNDS.json',
                   'inputs/reference/INERTIAL_WAVEFORM_MANIFEST.json',
                   'inputs/reference/PHYSICAL_SHELL_RESPONSE.md',
                   'inputs/reference/PHYSICAL_SHELL_RESPONSE_MANIFEST.json',
                   'inputs/reference/FULL_CONTACT_MANIFEST.json',
                   'inputs/reference/GH_BOUNDARY_SOURCE_AND_CONSTRAINT_COMPLETION.md',
                   'inputs/reference/FULL_FINITE_ANGULAR_SECOND_ORDER_CONTACT.md',
                   'inputs/reference/COMPLETE_L4_SCALAR_BULK.md',
                   'inputs/reference/FULL_L4_SCALAR_SOURCE_AND_OPERATOR.md'],
  'arguments': ['--output-dir', 'generated'],
  'output': 'generated/COORDINATE_SHELL_SOURCE_ENCLOSURE.json',
  'golden': 'computed/COORDINATE_SHELL_SOURCE_ENCLOSURE.json',
  'optimized_golden': 'computed/COORDINATE_SHELL_SOURCE_ENCLOSURE_OPTIMIZED.json',
  'status': 'PASS',
  'byte_equal': True,
  'additional_outputs': [['generated/COORDINATE_SHELL_SOURCE_CELLS.csv',
                          'computed/COORDINATE_SHELL_SOURCE_CELLS.csv']],
  'mode_fields': [{'path': ['runtime', 'optimization'], 'normal': 0, 'optimized': 1}],
  'required_fields': [{'path': ['cell_data', 'rows'], 'equals': 4353},
                      {'path': ['scope', 'full_U2_evaluated'], 'equals': False},
                      {'path': ['scope', 'EN_scope'],
                       'equals': 'Retained degree-four reduced Robin forcing only.'}],
  'status_path': ['status'],
  'checks': None,
  'embedded_checks': None,
  'check_count_path': None,
  'embedded_count_path': None,
  'check_names': [],
  'embedded_control_names': [],
  'groups': 0,
  'embedded': 0},
 {'name': 'astra',
  'script': 'astra/check_coordinate_sources.py',
  'dependencies': ['inputs/INPUT_PINS.json',
                   'inputs/prior/INERTIAL_WAVEFORM_CELLS.csv',
                   'inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json',
                   'inputs/prior/INERTIAL_SHELL_WAVEFORM.md',
                   'inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json',
                   'inputs/reference/PHYSICAL_SHELL_RESPONSE.md',
                   'inputs/reference/SCALAR_SHELL_BOUNDARY_SYSTEM.md',
                   'inputs/reference/PHYSICAL_SHELL_RESPONSE_MANIFEST.json',
                   'inputs/reference/FULL_CONTACT_MANIFEST.json',
                   'inputs/reference/INERTIAL_WAVEFORM_MANIFEST.json',
                   'inputs/reference/COMPLETE_L4_SCALAR_BULK.md',
                   'computed/COORDINATE_SHELL_SOURCE_ENCLOSURE.json',
                   'computed/COORDINATE_SHELL_SOURCE_CELLS.csv'],
  'arguments': ['--output', 'checked.json'],
  'output': 'checked.json',
  'golden': 'astra/ASTRA_COORDINATE_NORMAL.json',
  'optimized_golden': 'astra/ASTRA_COORDINATE_OPTIMIZED.json',
  'status': 'PASS',
  'byte_equal': True,
  'checks': {'path': ['positive_groups'], 'kind': 'names'},
  'check_count_path': ['positive_group_count'],
  'additional_outputs': [],
  'status_path': ['status'],
  'mode_fields': [],
  'embedded_checks': None,
  'required_fields': [],
  'embedded_count_path': None,
  'check_names': ['pin:inputs/prior/INERTIAL_WAVEFORM_CELLS.csv',
                  'registry:inputs/prior/INERTIAL_WAVEFORM_CELLS.csv',
                  'pin:inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json',
                  'registry:inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json',
                  'pin:inputs/prior/INERTIAL_SHELL_WAVEFORM.md',
                  'registry:inputs/prior/INERTIAL_SHELL_WAVEFORM.md',
                  'pin:inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json',
                  'registry:inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json',
                  'pin:inputs/reference/PHYSICAL_SHELL_RESPONSE.md',
                  'registry:inputs/reference/PHYSICAL_SHELL_RESPONSE.md',
                  'pin:inputs/reference/SCALAR_SHELL_BOUNDARY_SYSTEM.md',
                  'registry:inputs/reference/SCALAR_SHELL_BOUNDARY_SYSTEM.md',
                  'pin:inputs/reference/PHYSICAL_SHELL_RESPONSE_MANIFEST.json',
                  'registry:inputs/reference/PHYSICAL_SHELL_RESPONSE_MANIFEST.json',
                  'pin:inputs/reference/FULL_CONTACT_MANIFEST.json',
                  'registry:inputs/reference/FULL_CONTACT_MANIFEST.json',
                  'pin:inputs/reference/INERTIAL_WAVEFORM_MANIFEST.json',
                  'registry:inputs/reference/INERTIAL_WAVEFORM_MANIFEST.json',
                  'pin:inputs/reference/COMPLETE_L4_SCALAR_BULK.md',
                  'registry:inputs/reference/COMPLETE_L4_SCALAR_BULK.md',
                  'manifest_membership:PHYSICAL_SHELL_RESPONSE.md',
                  'manifest_membership:boundary/SCALAR_SHELL_BOUNDARY_SYSTEM.md',
                  'manifest_membership:astra-bulk/COMPLETE_L4_SCALAR_BULK.md',
                  'manifest_membership:INERTIAL_SHELL_WAVEFORM.md',
                  'manifest_membership:computed/INERTIAL_WAVEFORM_CELLS.csv',
                  'manifest_membership:computed/INERTIAL_WAVEFORM_ENCLOSURE.json',
                  'manifest_membership:archive/inputs/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json',
                  'unit_normal_series',
                  'derivative_taylor_coefficient',
                  'extrinsic_derivative_coefficient',
                  'one_copy_scalar_source',
                  'metric_coordinate_source',
                  'robin_join',
                  'degree_two_algebraic_row',
                  'degree_two_robin_join',
                  'harmonic_decomposition',
                  'harmonic_extensions',
                  'harmonic_orthogonality',
                  'harmonic_piece_norms',
                  'angular_total_norm',
                  'proper_volume_power',
                  'dynamic_boundary_gamma_cancellation',
                  'registered_state',
                  'shell_signs',
                  'rounded_mu_encloses_exact_box',
                  'whole_cell_partition',
                  'zero_prefront',
                  'actual_trace_envelope',
                  'scaled_boundary_force',
                  'window_operator_bound',
                  'volterra_resolvent',
                  'metric_reconstruction_radius',
                  'proper_clock_chain_rule',
                  'central_positive_convolution',
                  'positive_central_inverse_mu',
                  'central_metric_positive',
                  'central_actual_scalar_zero',
                  'both_central_source_signs',
                  'central_endpoint_nonzero:central_left',
                  'central_endpoint_source_signs:central_left',
                  'central_endpoint_nonzero:central_right',
                  'central_endpoint_source_signs:central_right',
                  'tail_cells_in_retained_theorem',
                  'strict_tail_scalar_sign',
                  'strict_tail_metric_sign',
                  'registered_robin_symmetrizer_positive',
                  'new_source_table_identity',
                  'norm_pi_interval',
                  'norm_background_identity',
                  'norm_full_Y_squared_prefactor',
                  'norm_EN_degree_four_prefactor',
                  'norm_factor_join:Sp',
                  'norm_factor_join:SK',
                  'norm_factor_join:EN',
                  'norm_factor_join:dSp',
                  'norm_factor_join:dSK',
                  'norm_factor_join:dEN',
                  'norm_factor_join:Sp_material',
                  'norm_factor_join:SK_material',
                  'norm_factor_join:EN_material',
                  'norm_ratio_interval_upper',
                  'norm_ratio_report_identity',
                  'norm_ratio_strict_nine_e_minus_thirteen'],
  'embedded_control_names': [],
  'groups': 83,
  'embedded': 0},
 {'name': 'sol',
  'script': 'sol/check_coordinate_shell_sources.py',
  'dependencies': ['compute_coordinate_sources.py',
                   'inputs/INPUT_PINS.json',
                   'inputs/prior/INERTIAL_WAVEFORM_CELLS.csv',
                   'inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json',
                   'inputs/prior/INERTIAL_SHELL_WAVEFORM.md',
                   'inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json',
                   'inputs/prior/PRIOR_PHYSICAL_GAIN_ENCLOSURE.json',
                   'inputs/prior/FINAL_INERTIAL_WAVEFORM_REVIEW.json',
                   'inputs/prior/FINAL_INERTIAL_WAVEFORM_REVIEW.md',
                   'inputs/reference/SCALAR_SHELL_BOUNDARY_SYSTEM.md',
                   'inputs/reference/MATERIAL_SECOND_ORDER_JUNCTIONS.md',
                   'inputs/reference/BACKGROUND_DERIVATIVE_BOUNDS.json',
                   'inputs/reference/INERTIAL_WAVEFORM_MANIFEST.json',
                   'inputs/reference/PHYSICAL_SHELL_RESPONSE.md',
                   'inputs/reference/PHYSICAL_SHELL_RESPONSE_MANIFEST.json',
                   'inputs/reference/FULL_CONTACT_MANIFEST.json',
                   'inputs/reference/GH_BOUNDARY_SOURCE_AND_CONSTRAINT_COMPLETION.md',
                   'inputs/reference/FULL_FINITE_ANGULAR_SECOND_ORDER_CONTACT.md',
                   'inputs/reference/COMPLETE_L4_SCALAR_BULK.md',
                   'inputs/reference/FULL_L4_SCALAR_SOURCE_AND_OPERATOR.md',
                   'computed/COORDINATE_SHELL_SOURCE_ENCLOSURE.json',
                   'computed/COORDINATE_SHELL_SOURCE_CELLS.csv'],
  'arguments': ['--output', 'checked.json'],
  'output': 'checked.json',
  'golden': 'sol/SOL_NORMAL.json',
  'optimized_golden': 'sol/SOL_OPTIMIZED.json',
  'status': 'PASS',
  'byte_equal': True,
  'checks': {'path': ['finite_groups'], 'kind': 'names'},
  'mode_fields': [{'path': ['optimized'], 'normal': False, 'optimized': True}],
  'required_fields': [{'path': ['cell_count'], 'equals': 4353},
                      {'path': ['primary_authentication', 'status'], 'equals': 'PASS'},
                      {'path': ['primary_authentication', 'scientific_endpoints_authenticated'],
                       'equals': 113178}],
  'additional_outputs': [],
  'status_path': ['status'],
  'embedded_checks': None,
  'check_count_path': None,
  'embedded_count_path': None,
  'check_names': ['input identity and coverage',
                  'signed Volterra primitive and resolvent',
                  'clock/radius derivative reconstruction',
                  'complete coefficient and harmonic rows',
                  'whole-cell source norms'],
  'embedded_control_names': [],
  'groups': 5,
  'embedded': 0},
 {'name': 'referee',
  'script': 'referee/check_full_source_formulation.py',
  'dependencies': ['inputs/INPUT_PINS.json',
                   'inputs/prior/INERTIAL_WAVEFORM_CELLS.csv',
                   'inputs/prior/INERTIAL_WAVEFORM_ENCLOSURE.json',
                   'inputs/prior/INERTIAL_SHELL_WAVEFORM.md',
                   'inputs/prior/PRIOR_ANALYTIC_SHELL_BOX_DELTA_1E3.json',
                   'inputs/prior/PRIOR_PHYSICAL_GAIN_ENCLOSURE.json',
                   'inputs/prior/FINAL_INERTIAL_WAVEFORM_REVIEW.json',
                   'inputs/prior/FINAL_INERTIAL_WAVEFORM_REVIEW.md',
                   'inputs/reference/SCALAR_SHELL_BOUNDARY_SYSTEM.md',
                   'inputs/reference/MATERIAL_SECOND_ORDER_JUNCTIONS.md',
                   'inputs/reference/BACKGROUND_DERIVATIVE_BOUNDS.json',
                   'inputs/reference/INERTIAL_WAVEFORM_MANIFEST.json',
                   'inputs/reference/PHYSICAL_SHELL_RESPONSE.md',
                   'inputs/reference/PHYSICAL_SHELL_RESPONSE_MANIFEST.json',
                   'inputs/reference/FULL_CONTACT_MANIFEST.json',
                   'inputs/reference/GH_BOUNDARY_SOURCE_AND_CONSTRAINT_COMPLETION.md',
                   'inputs/reference/FULL_FINITE_ANGULAR_SECOND_ORDER_CONTACT.md',
                   'inputs/reference/COMPLETE_L4_SCALAR_BULK.md',
                   'inputs/reference/FULL_L4_SCALAR_SOURCE_AND_OPERATOR.md',
                   'computed/COORDINATE_SHELL_SOURCE_ENCLOSURE.json'],
  'arguments': ['--output', 'checked.json'],
  'output': 'checked.json',
  'golden': 'referee/FULL_SOURCE_NORMAL.json',
  'optimized_golden': 'referee/FULL_SOURCE_OPTIMIZED.json',
  'status': 'PASS',
  'byte_equal': True,
  'checks': {'path': ['check_names'], 'kind': 'names'},
  'check_count_path': ['check_count'],
  'required_fields': [{'path': ['primitive_check_count'], 'equals': 557}],
  'additional_outputs': [],
  'status_path': ['status'],
  'mode_fields': [],
  'embedded_checks': None,
  'embedded_count_path': None,
  'check_names': ['authenticate_portable_input_inventory',
                  'coordinate_source_scope_and_degree_four_only_Robin_norms',
                  'actual_shell_Robin_symmetrizer_positive_definite',
                  'generic_5D_first_connection',
                  'generic_5D_second_connection_factorial',
                  'all_fifteen_Einstein_Hessian_components',
                  'scalar_Hessian_all_terms',
                  'lower_wave_map_Hessian_cross_term',
                  'coordinate_scalar_and_Israel_normalization',
                  'degree_four_Robin_join',
                  'distinct_degree_two_Robin_join',
                  'ambient_clock_proportional_initial_profile',
                  'background_proper_time_metric_derivative_chain_rule',
                  'actual_wider_collar_first_exit_and_potential_bound',
                  'interior_Dirichlet_Volterra_operator_bound',
                  'independent_actual_H_boundary_sup_and_interior_value_error',
                  'backward_closed_causal_domain_contains_complete_first_support'],
  'embedded_control_names': [],
  'groups': 17,
  'embedded': 0}]
CONTROLS = [{'checker': 'astra',
  'name': 'normal_quadratic',
  'reason': 'unit_normal_series',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject-error', 'normal_quadratic'],
  'failure_format': 'json',
  'failure_path': ['failed_gate']},
 {'checker': 'astra',
  'name': 'one_copy_doubled',
  'reason': 'one_copy_scalar_source',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject-error', 'one_copy_doubled'],
  'failure_format': 'json',
  'failure_path': ['failed_gate']},
 {'checker': 'astra',
  'name': 'taylor_half',
  'reason': 'derivative_taylor_coefficient',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject-error', 'taylor_half'],
  'failure_format': 'json',
  'failure_path': ['failed_gate']},
 {'checker': 'astra',
  'name': 'sk_metric_sign',
  'reason': 'metric_coordinate_source',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject-error', 'sk_metric_sign'],
  'failure_format': 'json',
  'failure_path': ['failed_gate']},
 {'checker': 'astra',
  'name': 'angular_norm',
  'reason': 'angular_total_norm',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject-error', 'angular_norm'],
  'failure_format': 'json',
  'failure_path': ['failed_gate']},
 {'checker': 'astra',
  'name': 'volume_power',
  'reason': 'proper_volume_power',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject-error', 'volume_power'],
  'failure_format': 'json',
  'failure_path': ['failed_gate']},
 {'checker': 'astra',
  'name': 'scaled_force',
  'reason': 'scaled_boundary_force',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject-error', 'scaled_force'],
  'failure_format': 'json',
  'failure_path': ['failed_gate']},
 {'checker': 'astra',
  'name': 'drop_resolvent',
  'reason': 'volterra_resolvent',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject-error', 'drop_resolvent'],
  'failure_format': 'json',
  'failure_path': ['failed_gate']},
 {'checker': 'astra',
  'name': 'central_mu',
  'reason': 'positive_central_inverse_mu',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject-error', 'central_mu'],
  'failure_format': 'json',
  'failure_path': ['failed_gate']},
 {'checker': 'astra',
  'name': 'tail_drop_square',
  'reason': 'strict_tail_scalar_sign',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject-error', 'tail_drop_square'],
  'failure_format': 'json',
  'failure_path': ['failed_gate']},
 {'checker': 'astra',
  'name': 'proper_derivative',
  'reason': 'proper_clock_chain_rule',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject-error', 'proper_derivative'],
  'failure_format': 'json',
  'failure_path': ['failed_gate']},
 {'checker': 'astra',
  'name': 'use_comparator',
  'reason': 'actual_trace_envelope',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject-error', 'use_comparator'],
  'failure_format': 'json',
  'failure_path': ['failed_gate']},
 {'checker': 'astra',
  'name': 'uncancelled_gamma',
  'reason': 'dynamic_boundary_gamma_cancellation',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject-error', 'uncancelled_gamma'],
  'failure_format': 'json',
  'failure_path': ['failed_gate']},
 {'checker': 'astra',
  'name': 'metric_radius',
  'reason': 'metric_reconstruction_radius',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject-error', 'metric_radius'],
  'failure_format': 'json',
  'failure_path': ['failed_gate']},
 {'checker': 'astra',
  'name': 'endpoint_drop_square',
  'reason': 'central_endpoint_source_signs:central_left',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject-error', 'endpoint_drop_square'],
  'failure_format': 'json',
  'failure_path': ['failed_gate']},
 {'checker': 'astra',
  'name': 'erase_central_metric',
  'reason': 'both_central_source_signs',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject-error', 'erase_central_metric'],
  'failure_format': 'json',
  'failure_path': ['failed_gate']},
 {'checker': 'astra',
  'name': 'ratio_denominator',
  'reason': 'norm_ratio_interval_upper',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject-error', 'ratio_denominator'],
  'failure_format': 'json',
  'failure_path': ['failed_gate']},
 {'checker': 'sol',
  'name': 'wrong_scaled_force',
  'reason': 'wrong_scaled_force: ell/mu factor missing',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject', 'wrong_scaled_force'],
  'failure_format': 'fail_line'},
 {'checker': 'sol',
  'name': 'omit_triangular_half',
  'reason': 'omit_triangular_half: integral of 2*(x-s) is x^2',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject', 'omit_triangular_half'],
  'failure_format': 'fail_line'},
 {'checker': 'sol',
  'name': 'wrong_negative_mu',
  'reason': 'wrong_negative_mu: signed endpoint witness excluded',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject', 'wrong_negative_mu'],
  'failure_format': 'fail_line'},
 {'checker': 'sol',
  'name': 'omit_resolvent',
  'reason': 'omit_resolvent: contraction residual not closed',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject', 'omit_resolvent'],
  'failure_format': 'fail_line'},
 {'checker': 'sol',
  'name': 'wrong_metric_radius',
  'reason': 'wrong_metric_radius: u=rho^-2*C*y',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject', 'wrong_metric_radius'],
  'failure_format': 'fail_line'},
 {'checker': 'sol',
  'name': 'wrong_angular_normalization',
  'reason': 'wrong_angular_normalization: sphere Y^4 mean is 3/40',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject', 'wrong_angular_normalization'],
  'failure_format': 'fail_line'},
 {'checker': 'sol',
  'name': 'omit_measure_clock',
  'reason': 'omit_measure_clock: positive Taylor term of cosh^3 excluded',
  'mechanism': 'formula comparison CLI injection',
  'arguments': ['--inject', 'omit_measure_clock'],
  'failure_format': 'fail_line'},
 {'checker': 'sol',
  'name': 'last_cell_upper_endpoint',
  'reason': 'primary all-cell endpoint mismatch: row 4352 Sp_hi',
  'mechanism': 'copied full-cell CSV endpoint corruption',
  'failure_format': 'fail_line',
  'target': 'computed/COORDINATE_SHELL_SOURCE_CELLS.csv',
  'before': '12272,12288,3711244296534467708204190678860515061206,4212260265796515538231137610302858946969,-41664289715150165891938049036107609301,37827112094574399906771060098370538863,219853727149239712188758681431306228,249536001027035446766653683513546380,-18997261592405819043187990258542619120,17247660966548305037541048742329820263,-648408065873721032668802644424987825,-345797802844931071768657089463851287,-12969827882902510795117372876259328000000,-7161688047994915977414340390158336000000,14322126682635687754542175021014718277844,25937558852033753674781559029761800824960,-15087332600219872476156434706120130073,-4600390954289341794090099922998804761,1533141634663730324794926399645294988,5028317363946325972328184239593708086,6132566442471534565036661008116932100,20113269416689971510920402698861120901,-1160921687886823311059421525915170,-564783442578431166661232332372934,-3731859256503434691647755971,-2896853253988405958901839396,-111110823760156783357055271736,-50682745394345958095121069024\n',
  'after': '12272,12288,3711244296534467708204190678860515061206,4212260265796515538231137610302858946969,-41664289715150165891938049036107609301,37827112094574399906771060098370538863,219853727149239712188758681431306228,249536001027035446766653683513546380,-18997261592405819043187990258542619120,17247660966548305037541048742329820263,-648408065873721032668802644424987825,-345797802844931071768657089463851287,-12969827882902510795117372876259328000000,-7161688047994915977414340390158336000000,14322126682635687754542175021014718277844,25937558852033753674781559029761800824960,-15087332600219872476156434706120130073,-4600390954289341794090099922998804760,1533141634663730324794926399645294988,5028317363946325972328184239593708086,6132566442471534565036661008116932100,20113269416689971510920402698861120901,-1160921687886823311059421525915170,-564783442578431166661232332372934,-3731859256503434691647755971,-2896853253988405958901839396,-111110823760156783357055271736,-50682745394345958095121069024\n'},
 {'checker': 'sol',
  'name': 'last_norm_upper_endpoint',
  'reason': 'JSON mathematical value mismatch: /exact/norms/EN_material/squared_L2/1',
  'mechanism': 'copied full-report norm endpoint corruption',
  'failure_format': 'fail_line',
  'target': 'computed/COORDINATE_SHELL_SOURCE_ENCLOSURE.json',
  'json_path': ['exact', 'norms', 'EN_material', 'squared_L2', 1],
  'value': '431207408323303005560211703625189650878930376603830807573270036146804993403053830906606988109521193352361914856419554991361392235066030155659763122351666228385085187753424921374908812427939799668753728162389495275520974768610910837019866355261134268551837861995169771665502825395308987054032434324001921240056117469751/319079212577344121102871000369765878705533547757495644843624752374731681691256507918990011941561869509415352680962962089874536209342578682187020188091315069962662161923101874959091608250348082347682780314473222133031123599097750682673785459116409367253980778433973675141665652736000000000000000000'},
 {'checker': 'referee',
  'name': 'half_ricci',
  'reason': 'Einstein Hessian (0, 0): nonzero rational residual -8683970233/311040000',
  'failure_format': 'json',
  'failure_path': ['reason'],
  'arguments': ['--output', 'rejected.json', '--inject-error', 'half_ricci'],
  'mechanism': 'formula comparison CLI injection'},
 {'checker': 'referee',
  'name': 'drop_potential_cross',
  'reason': 'Einstein Hessian (0, 2): nonzero rational residual -16/3',
  'failure_format': 'json',
  'failure_path': ['reason'],
  'arguments': ['--output', 'rejected.json', '--inject-error', 'drop_potential_cross'],
  'mechanism': 'formula comparison CLI injection'},
 {'checker': 'referee',
  'name': 'drop_lower_gauge_cross',
  'reason': 'lower wave-map Hessian: nonzero rational residual -767167/86400',
  'failure_format': 'json',
  'failure_path': ['reason'],
  'arguments': ['--output', 'rejected.json', '--inject-error', 'drop_lower_gauge_cross'],
  'mechanism': 'formula comparison CLI injection'},
 {'checker': 'referee',
  'name': 'wrong_robin_sign',
  'reason': 'actual symmetrized Robin matrix positivity failed',
  'failure_format': 'json',
  'failure_path': ['reason'],
  'arguments': ['--output', 'rejected.json', '--inject-error', 'wrong_robin_sign'],
  'mechanism': 'formula comparison CLI injection'},
 {'checker': 'referee',
  'name': 'wrong_metric_clock',
  'reason': 'background proper-time metric derivative: nonzero rational residual '
            '-39312004752/1716715',
  'failure_format': 'json',
  'failure_path': ['reason'],
  'arguments': ['--output', 'rejected.json', '--inject-error', 'wrong_metric_clock'],
  'mechanism': 'formula comparison CLI injection'},
 {'checker': 'referee',
  'name': 'drop_bulk_half_factor',
  'reason': 'interior Volterra contraction: nonzero rational residual 459/1000000000',
  'failure_format': 'json',
  'failure_path': ['reason'],
  'arguments': ['--output', 'rejected.json', '--inject-error', 'drop_bulk_half_factor'],
  'mechanism': 'formula comparison CLI injection'},
 {'checker': 'referee',
  'name': 'wrong_EN_scope_input',
  'reason': 'coordinate report scientific scope failed',
  'failure_format': 'json',
  'failure_path': ['reason'],
  'arguments': ['--output', 'rejected.json'],
  'mechanism': 'copied coordinate-report scope corruption',
  'target': 'computed/COORDINATE_SHELL_SOURCE_ENCLOSURE.json',
  'json_path': ['scope', 'EN_scope'],
  'value': 'All scalar degrees including monopole.'}]
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
    receipt = {'schema': 'hdblast-coordinate-shell-sources-fresh-replay-v1', 'status': 'PASS', 'capsule': CAPSULE,
        'authentication': authentication, 'runner': inventory['reproduce.py'],
        'integrity_verifier': inventory['verify_manifest.py'], 'requirements': inventory['requirements.txt'],
        'roster': inventory['REPLAY_ROSTER.json'], 'python': {'version': sys.version, 'executable': sys.executable},
        'dependency_versions': dependency_versions, 'embedded_configuration_sha256': CONFIG_SHA256,
        **counts, 'finite_positive_check_groups_by_mode': groups,
        'mode_agreements': agreements, 'steps': rows,
        'scope': 'Fresh finite arithmetic for the approved coordinate shell boundary trace, scalar and isotropic Israel source coefficients only. The full second-order U2 response remains unresolved quantitatively. Complete generated JSON and CSV reproduce byte for byte. All golden JSON fields are compared; only individually verified execution-mode metadata is normalized across modes. Executions do not prove analytic PDE/exact-family hypotheses, compute a finite-amplitude nonlinear waveform, set an effective nonlinear amplitude or observational scale, or replace external peer review.'}
    (output/'REPRODUCTION_RECEIPT.json').write_text(json.dumps(receipt, indent=2)+'\n')
    print(json.dumps({'status': 'PASS', **counts}, indent=2))


if __name__ == '__main__':
    main()
