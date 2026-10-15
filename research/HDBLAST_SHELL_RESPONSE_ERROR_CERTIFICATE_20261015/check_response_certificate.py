"""Finite arithmetic checks, separate from the conditional analytic proof."""
from fractions import Fraction as F
import importlib.util
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('response_certificate', ROOT / 'certify_response.py')
cert = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cert)
groups = []
failures = []


def check(name, condition):
    groups.append(name)
    if not condition:
        failures.append(name)


for value in [F(0), F(1), F(2), F(1, 3), F(123456789, 7654321)]:
    low, high = cert.sqrt_bounds(value)
    check('sqrt_outward_' + str(value), 0 <= low <= high and low * low <= value <= high * high)
check('sqrt_exact_zero', cert.sqrt_bounds(0) == (0, 0))
check('sqrt_exact_square', cert.sqrt_bounds(81) == (9, 9))
check('exp_positive_bound', cert.exp_upper(F(1, 2)) >= F(13, 8))
check('exp_zero', cert.exp_upper(0) == 1)
pil, piu = cert.pi_bounds()
check('machin_pi_bracket', F(333, 106) < pil <= piu < F(355, 113))
check('machin_precision', piu - pil < F(1, 10**40))
data = cert.comparator_inputs()
result = cert.response_certificate(data)
check('positive_inertia_forcing', result['forcing_squared_upper'] ==
      data['F_l2_sup'] ** 2 + 1000 * data['G_abs_sup'] ** 2)
expected_trace = (data['T'] * (data['F_l2_sup'] + data['W_abs'] * result['energy_norm_upper'])
                  + cert.sqrt_bounds(data['T'] * (1000 + 250 ** 2))[1]
                  * result['energy_norm_upper'])
check('physical_shift_trace', result['trace_error_upper'] == expected_trace)
check('norm_triangle_lower', result['true_D_norm_lower'] ==
      cert.sqrt_bounds(data['trial_D_sq_lower'])[0] - result['trace_error_upper'])
check('proper_time_prefactor', result['C_lower'] ==
      3 * pil ** 2 / (17161 * F(49, 2500000000)) * result['true_D_norm_lower'] ** 2)
check('selected_state_positive', result['positive_lower_bound'] is True)
check('published_outward_lower', result['C_lower'] > 10**15)
check('published_outward_upper', result['C_upper'] < 4 * 10**17)
check('error_below_eight', result['trace_error_upper'] < 8)
for name, updates in [('negative_mass', {'mu_lower': '-1'}),
                      ('zero_scalar_denominator', {'p_sq_lower': '0'}),
                      ('missing_domain_margin', {'L': str(data['T'])}),
                      ('reversed_norm_interval', {'trial_D_sq_upper': '0'}),
                      ('floating_input', {'T': 0.01})]:
    trial = dict(data)
    trial.update(updates)
    try:
        cert.response_certificate(trial)
    except ValueError:
        rejected = True
    else:
        rejected = False
    check('reject_' + name, rejected)
try:
    cert.comparator_inputs(F(1, 10**8))
except ValueError:
    rejected = True
else:
    rejected = False
check('reject_unvalidated_collar', rejected)
receipt = {'schema': 'hdblast-response-certificate-arithmetic-v1',
           'status': 'FAIL' if failures else 'PASS',
           'finite_check_groups': len(groups), 'groups': groups, 'failures': failures,
           'inputs_exact': cert.rational_json(data),
           'enclosure_exact': cert.rational_json(result),
           'published_coarse_interval_exact': ['1000000000000000', '400000000000000000'],
           'scope': 'Exact-rational arithmetic consequences of the stated conditional energy/trace and analytic-background hypotheses. No waveform solver, formal PDE certification or effective nonlinear amplitude interval.'}
print(json.dumps(receipt, indent=2))
raise SystemExit(1 if failures else 0)
