"""Replay the localized second-order checkpoint and validate fresh receipts.

Python 3.9+, sympy 1.14.0, mpmath 1.3.0. Sources are always regenerated
from text code; uploaded or pre-existing pickle caches are never trusted.
--core omits only the independent radial quadrature check for CI use.
The optional --full mode also runs the original slow exact verifier.
"""
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import pickle
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
PY = sys.executable
STEPS = [
    ('new', 'kids_tensor_only.py', 'kids_tensor_only.json'),
    ('new', 'kids_symbolic_l01.py', 'kids_symbolic_l01.json'),
    ('new', 'l0_full_sector.py', 'l0_full_sector.json'),
    ('new', 'wgg_canonical.py', 'wgg_canonical.json'),
    ('new', 'exact_random_point.py', 'exact_random_point.json'),
    ('new', 'fast_verify.py', 'fast_verify.json'),
    ('new', 'angular_exact.py', 'angular_exact.json'),
    ('independent', 'indep_check.py', 'indep_check.json'),
    ('independent', 'indep_l0.py', 'indep_l0.json'),
    ('independent', 'compare_f0.py', 'compare_f0.json'),
]

def finite_below(value, limit):
    number = float(value)
    return math.isfinite(number) and abs(number) < limit

def expression_equal(actual, expected):
    import sympy as sp
    return sp.simplify(sp.sympify(actual) - sp.sympify(expected)) == 0

def check_tensor_kernel(data):
    for label, points in data.items():
        if set(points) != {'point0', 'point1'}:
            return False
        for result in points.values():
            rank = 2 if label in ('l0', 'l1') else 3
            if result['rank'] != rank or len(result['nullspace']) != 3 - rank:
                return False
            if label == 'l0' and result.get('nullspace_is_Omega_jet') is not False:
                return False
            if label == 'l1' and result.get('nullspace_is_Omega_jet') is not True:
                return False
    return set(data) == {'l0', 'l1', 'l2', 'l2b', 'l3', 'l4'}

def check_symbolic_kernel(data):
    return (data['l1']['rank'] == 2 and data['l1']['nullspace_dim'] == 1
            and data['l1']['nullspace_is_Omega_jet'] == [True]
            and expression_equal(data['l0']['compatibility_residual'],
                                 '-F1**2*M0*O0**2/(3*(O0+2*O1*s_)**2)'))

def check_l0(data):
    reduced = data['weighted_reduced_density']
    return (data['linear'] is True and data['gauge_residual'] == '0'
            and expression_equal(reduced['alpha'], '-9*O0*sqrt(s_)/2')
            and expression_equal(reduced['a'], '-3*O0*s_**(3/2)/2')
            and reduced['Psi'] == '0')

def check_exact_samples(data):
    return (data['all_identities_zero'] is True and data['all_controls_detected'] is True
            and len(data['trials']) == 3 and all(
                trial['H_residual_exact'] == '0' and trial['M_residual_exact'] == ['0'] * 4
                and trial['control_sign_flipped_tf_metric_detected'] is True
                and trial['control_doubled_scalar_momentum_detected'] is True
                and trial['control_omit_l0_subtraction_detected'] is True
                for trial in data['trials']))

def check_independent(data):
    if not all(finite_below(value, 1e-30) for name in
               ('background_beta_ODE_residual_samples', 'background_dU_consistency_samples')
               for value in data[name]):
        return False
    cases = data['cases']
    if set(cases) != {'commuting_A_eq_B', 'pure_velocity_G0'}:
        return False
    count = 0
    for label, case in cases.items():
        if len(case['spheres']) != 2:
            return False
        for sphere in case['spheres']:
            if len(sphere['points']) != (3 if label == 'commuting_A_eq_B' else 2):
                return False
            if not finite_below(sphere['remainder_spread_rel'], 1e-20):
                return False
            for point in sphere['points']:
                count += 1
                if not all(finite_below(point[name], 1e-30) for name in
                           ('background_H', 'first_order_residual_max', 'momentum_residual_rel')):
                    return False
            if label == 'pure_velocity_G0':
                measured = float(sphere['f0_from_remainder'])
                expected = float(sphere['f0_velocity_hand_formula'])
                if not (math.isfinite(measured) and math.isfinite(expected) and expected > 0
                        and abs(measured - expected) / expected < 1e-18):
                    return False
    return count == 10

def check_full(data):
    controls = data['wrong_formula_controls_detected']
    return (data['momentum_exact'] is True
            and data['hamiltonian_exact_given_l0_ode'] is True
            and data['algebraic_coefficients_vanish_off_profile_support'] is True
            and set(controls) == {'omit_scalar_momentum_correction', 'wrong_sign_l0_source',
                                 'drop_tangential_tf_metric'}
            and all(value is True for value in controls.values()))

CHECKS = {
    'kids_tensor_only.json': check_tensor_kernel,
    'kids_symbolic_l01.json': check_symbolic_kernel,
    'l0_full_sector.json': check_l0,
    'wgg_canonical.json': lambda data: data['KK_matches_hand_derivation'] is True
        and expression_equal(data['c2'], '6*s_**(13/2)/(F1**2*O0)'),
    'exact_random_point.json': check_exact_samples,
    'fast_verify.json': lambda data: data['no_radial_derivatives_of_coefficients'] is True
        and data['momentum_exact'] is True
        and data['wrong_formula_controls_detected'] == {'doubled_scalar_momentum_coefficient': True},
    'angular_exact.json': lambda data: data['normalization_basis_checks'] == [True, True]
        and data['coefficient_reconstruction_checks'] == [True, True, True]
        and data['source_first_order_ok'] is True
        and data['hamiltonian_identity_exact'] is True
        and data['wrong_sign_C1_control_detected'] is True,
    'indep_check.json': check_independent,
    'indep_l0.json': lambda data: finite_below(data['weighted_identity_rel_diff'], 1e-15)
        and len(data['gauge_residual_rel']) == 3
        and all(finite_below(value, 1e-25) for value in data['gauge_residual_rel'])
        and finite_below(data['psi_only_weighted_integral_rel'], 1e-20),
    'compare_f0.json': lambda data: len(data) == 4
        and all(finite_below(value['rel_diff'], 1e-18) for value in data.values()),
    'verify_cc.json': check_full,
}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def prepare_full_caches():
    """Recreate old verifier input structures from text JSON, never load old pickle."""
    import sympy as sp
    x = sp.symbols('x1:5')
    mons = [x[0]**(2*a) * x[1]**(2*b) for a in range(4) for b in range(4) if a+b <= 3]
    for part, prefixes in [('M', ('b', 'n')), ('H', ('c',))]:
        data = json.loads((HERE / 'basis/results' / ('construct_cc_' + part + '.json')).read_text())
        values = {sp.Symbol(prefix + str(i)): sp.sympify(data['nonzero_coefficients'].get(prefix + str(i), '0'))
                  for prefix in prefixes for i in range(len(mons))}
        with (HERE / 'basis/results' / ('construct_cc_' + part + '.pkl')).open('wb') as output:
            pickle.dump({'mons': mons, 'vals': values, 'f0': sp.sympify(data['f0'])}, output)

def main():
    import sympy, mpmath
    full = '--full' in sys.argv
    core = '--core' in sys.argv
    if full and core:
        print('--core and --full are mutually exclusive', flush=True)
        return 2
    selected_steps = [step for step in STEPS if not (core and step[1] == 'indep_l0.py')]
    receipt_path = HERE / 'REPRODUCTION_RECEIPT.json'
    receipt = {'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'python': sys.version.split()[0], 'sympy': sympy.__version__, 'mpmath': mpmath.__version__,
               'pickle_policy': 'Sources regenerated from text; full verifier caches recreated from text JSON.',
               'mode': 'core_subset' if core else ('default_plus_original_exact' if full else 'default'),
               'selected_checks': [str(Path(directory) / script) for directory, script, _ in selected_steps]
                   + (['basis/code/verify_cc.py'] if full else []),
               'excluded_checks': [str(Path(directory) / script) for directory, script, _ in STEPS
                                   if (directory, script, _) not in selected_steps],
               'all_checks_passed_scope': 'Only the checks selected by the recorded mode.',
               'steps': [], 'all_checks_passed': False}
    # Do not use a caller-supplied source/cache directory during a clean replay.
    environment = os.environ.copy()
    environment.pop('SO_RESULTS', None)
    def save():
        receipt_path.write_text(json.dumps(receipt, indent=2))
    def run(directory, script, output_name, args=(), validate=True):
        output_path = HERE / directory / output_name
        output_path.unlink(missing_ok=True)
        command = [PY, script, *args]
        t0 = time.monotonic()
        record = {'script': str(Path(directory) / script), 'command': command,
                  'source_sha256': sha(HERE / directory / script)}
        try:
            result = subprocess.run(command, cwd=HERE / directory, env=environment,
                                    capture_output=True, text=True)
            record['returncode'] = result.returncode
            ok = result.returncode == 0 and output_path.is_file()
            if ok and validate:
                ok = bool(CHECKS[output_path.name](json.loads(output_path.read_text())))
            record['passed'] = ok
            record['output_sha256'] = sha(output_path) if output_path.is_file() else None
            if not ok:
                record['stdout_tail'] = result.stdout[-2000:]
                record['stderr_tail'] = result.stderr[-2000:]
        except Exception as error:
            ok = False
            record['passed'] = False
            record['error'] = str(error)
        record['seconds'] = round(time.monotonic() - t0, 1)
        receipt['steps'].append(record)
        save()
        print(script, 'ok' if ok else 'FAILED', flush=True)
        if not ok:
            print(json.dumps(record, indent=2), flush=True)
        return ok
    save()
    if not run('basis/code', 'second_order_sources.py', '../results/sources_cc.pkl', args=('cc',), validate=False):
        return 1
    # This file was removed immediately before the reviewed generator ran;
    # only the newly generated local pickle is loaded here.
    try:
        with (HERE / 'basis/results/sources_cc.pkl').open('rb') as generated:
            first_order_ok = pickle.load(generated).get('first_order_ok') is True
    except Exception as error:
        first_order_ok = False
        receipt['steps'][-1]['error'] = str(error)
    receipt['steps'][-1]['first_order_constraints_vanish'] = first_order_ok
    if not first_order_ok:
        receipt['steps'][-1]['passed'] = False
        save()
        print('regenerated source failed first-order constraints', flush=True)
        return 1
    receipt['regenerated_sources_cc_pkl'] = True
    for directory, script, output in selected_steps:
        if not run(directory, script, output):
            return 1
    if full:
        try:
            prepare_full_caches()
        except Exception as error:
            receipt['steps'].append({'script': 'prepare_full_caches', 'passed': False, 'error': str(error)})
            save()
            print('full verifier cache generation FAILED', str(error), flush=True)
            return 1
        if not run('basis/code', 'verify_cc.py', '../results/verify_cc.json'):
            return 1
    receipt['all_checks_passed'] = True
    receipt['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    save()
    print('done', flush=True)
    return 0

if __name__ == '__main__':
    sys.exit(main())
