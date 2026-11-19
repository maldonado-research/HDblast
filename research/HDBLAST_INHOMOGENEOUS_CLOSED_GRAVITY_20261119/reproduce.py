"""Authenticate the closed-gravity study and optionally recompute its numerical results."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile


HERE = Path(__file__).resolve().parent
HISTORICAL_REVISION = '7cfe063ff0238e81107284d57169b83c025d6e30'
SOURCES_DIGEST = '6896f915803628e6c02e4b198d2e7884fb29bea5c967a59d64fa34df22144a00'
PRODUCTS = {'INITIAL': 'initial/INITIAL.json', 'TIDAL': 'tidal/TIDAL.json',
            'VOLUME': 'perturbative/VOLUME.json', 'EVOLUTION': 'evolution/results.json',
            'BOUNCE': 'checks/BOUNCE.json', 'COMPARISON': 'perturbative/COMPARISON.json'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def authenticate():
    inventory = HERE/'SHA256SUMS'
    listed = set()
    for line in inventory.read_text().splitlines():
        checksum, name = line.split('  ', 1)
        p = HERE/name
        require(p.is_relative_to(HERE) and Path(name).as_posix() == name
                and not Path(name).is_absolute() and p == p.resolve() and p != inventory and p not in listed,
                f'Invalid inventory entry {name}')
        require(digest(p) == checksum, f'Study checksum mismatch for {name}')
        listed.add(p)
    present = {p for p in HERE.rglob('*') if p.is_file() and p != inventory
               and not (p.parent.name == '__pycache__' and p.suffix == '.pyc')}
    require(listed == present, f'Study inventory differs: {present ^ listed}')
    require(digest(HERE/'SOURCES.json') == SOURCES_DIGEST, 'Historical source manifest mismatch')
    pins = read(HERE/'SOURCES.json')
    require(pins['base'] == HISTORICAL_REVISION and pins['previous_doi'] == '10.5281/zenodo.24378057',
            'Historical source revision or DOI differs')
    inputs = set()
    for entry in pins['sources']:
        path = HERE/entry['bundled']
        require(path in listed and path not in inputs and path.stat().st_size == entry['bytes']
                and digest(path) == entry['sha256'], f'Historical input differs: {entry["path"]}')
        inputs.add(path)
    require(inputs == {p for p in listed if p.is_relative_to(HERE/'previous')}, 'Historical input coverage differs')
    print(f'Verified {len(listed)} study files, including {len(inputs)} historical inputs.', flush=True)


def tolerance(path):
    # Physical numbers have relative comparisons; diagnostic cancellation, discretization
    # differences and initial proper means near zero need dimension-specific absolute floors.
    name, field = path[0], path[-1]
    relative, absolute = 1e-7, 0.
    if name == 'INITIAL':
        if 'psi_even_U_coefficients' in path:
            absolute = 5e-17
        elif 'coefficient_errors' in path:
            absolute = 5e-9
        elif field in ('proper_mean_q', 'fine_solution_sup_difference_bound'):
            absolute = 2e-16
        elif 'residual' in field or 'Hamiltonian' in field or field in ('top', 'removed_u_slope'):
            absolute = 8e-14
    elif name == 'TIDAL' and field == 'max_two_formula_difference':
        absolute = 5e-14
    elif name == 'VOLUME' and ('absolute_ODE_minus_quadrature_error' in path or 'absolute_error' in field):
        absolute = 2e-11
    elif name == 'EVOLUTION':
        if 'constraints' in path or 'max_sampled_constraints' in path or 'max_abs_regular_momentum' in field:
            absolute = 1e-10 if 'regular_momentum' in field else 8e-12
        elif ('stop_differences_from_primary' in path or 'maximum_spatial_refinement_differences' in path
              or 'maximum_temporal_refinement_differences' in path):
            absolute = 8e-12
        elif 'homogeneous_FRW_check' in path or 'probe_differences_from_primary' in path or field in ('gram_error', 'maximum_local_expansion_difference'):
            absolute = 2e-14
        elif field in ('mean', 'mean_change', 'round_mean', 'odd_mean', 'even_initial_mean', 'maximum_abs_centered_initial_mean'):
            absolute = 3e-16
    elif name == 'BOUNCE':
        if field == 'mean_laplacian' or field == 'mean' or field.endswith('tic_error'):
            absolute = 4e-17
        elif field == 'max_initial_hamiltonian':
            absolute = 8e-14
        elif field == 'quartic_error_relative_to_change':
            absolute = 6e-9
        elif field.endswith('quartic_error_ratio'):
            absolute = .04
    elif name == 'COMPARISON':
        if field == 'coupled_half_difference_per_A2' or field == 'relative_difference':
            absolute = 3e-8
        elif 'ratio_of_relative_errors_upon_halving' in path:
            absolute = .004
        elif field.startswith('max_sampled_'):
            absolute = 2e-14
    return relative, absolute


def compare(actual, archived, path):
    label = '/'.join(path)
    if isinstance(archived, dict):
        require(isinstance(actual, dict) and actual.keys() == archived.keys(), f'Object fields differ at {label}')
        for key in archived:
            compare(actual[key], archived[key], (*path, key))
    elif isinstance(archived, list):
        require(isinstance(actual, list) and len(actual) == len(archived), f'Array length differs at {label}')
        for index, (new, old) in enumerate(zip(actual, archived)):
            compare(new, old, (*path, str(index)))
    elif isinstance(archived, float):
        relative, absolute = tolerance(path)
        require(type(actual) in (int, float) and math.isfinite(actual) and math.isfinite(archived)
                and math.isclose(actual, archived, rel_tol=relative, abs_tol=absolute),
                f'{label}: computed {actual!r}, archived {archived!r}')
    elif path[-1] == 'number_rhs_evaluations':
        # An adaptive ODE solve can legitimately choose different internal timesteps.
        require(type(actual) is int and actual > 0 and abs(actual-archived) <= max(30, .15*archived),
                f'Adaptive evaluation count differs at {label}: {actual} versus {archived}')
    else:
        require(type(actual) is type(archived) and actual == archived, f'Scalar value differs at {label}')


def near(a, b, label, rel=1e-8, absolute=0.):
    require(math.isclose(a, b, rel_tol=rel, abs_tol=absolute), f'{label}: {a} versus {b}')


def conclusions(results):
    historical = read(HERE/'previous/nonlinear/NONLINEAR.json')['branches']
    for key in ('phi_m=0.5', 'phi_m=-0.5'):
        sign = 1 if key == 'phi_m=0.5' else -1
        initial = results['INITIAL']['branches'][key]
        evolved = results['EVOLUTION']['branches'][key]
        volume = results['VOLUME']['branches'][key]['times']
        near(initial['analytic_proper_volume_shift_over_A2'], volume[0]['volume_fractional_coefficient'], 'Initial volume')
        near(initial['analytic_volume_acceleration_shift_over_A2'], volume[0]['acceleration_over_radius_difference_coefficient'], 'Initial acceleration')
        near(initial['mu2_over_u_top'], results['TIDAL']['branches'][key]['beta'], 'Tidal curvature input')
        for group in ('runs','proper_mean_centered_runs'):
            for name, run in evolved[group].items():
                amplitude = run['amplitude']
                static = initial[group][str(float(amplitude))]
                start, probe = run['samples']['0'], run['at_probe']
                near(static['proper_mean_q'], start['mean'], 'Initial proper mean', absolute=2e-16)
                near(static['proper_volume_over_round'], start['volume_over_background'], 'Initial proper volume', rel=1e-13)
                near(static['proper_mean_q_acceleration'], results['BOUNCE'][key][group][name]['mean_xx'],
                     'Initial independent proper mean acceleration', absolute=2e-12)
                require(max(run['max_sampled_constraints'].values()) < 6e-10, 'Coupled constraints exceed tolerance')
                if amplitude == 0:
                    require(probe['mean'] == 0 and probe['maximum_field'] == 0 and probe['volume_fractional_difference'] == 0,
                            'Stationary control departed from the round solution')
                    continue
                require(sign*probe['mean_change'] > 0 and 0 < probe['H_volume'] and probe['H_volume_minus_background'] < 0,
                        'Mean or volume expansion has the wrong sign')
                if group == 'proper_mean_centered_runs':
                    require(abs(start['mean']) < 2e-16, 'Initial physical centering failed')
                else:
                    near(run['fixed_geometry_samples']['0.12']['mean'], historical[key]['runs'][name]['at_probe']['mean'],
                         'Pinned prescribed-background mean', rel=2e-7)
        extension = evolved['controlled_extensions']
        for run in extension['runs'].values():
            require(run['stop'] == 'peak_field', 'Extension did not reach its prescribed field endpoint')
            near(run['at_stop']['maximum_field'], extension['peak_field_cap'], 'Peak endpoint field', rel=2e-10)
            require(max(run['max_sampled_constraints'].values()) < 6e-10, 'Extension constraint exceeds tolerance')
        pairs = evolved['amplitude_pairs_at_probe']
        ratio = pairs['0.001']['even_mean_change']/pairs['0.0005']['even_mean_change']
        require(abs(ratio-4) < .015, 'Amplitude-halving test of the proper mean failed')
        for amplitude in (.001, .0005):
            pos, neg = (evolved['runs'][str(a)]['at_probe'] for a in (amplitude,-amplitude))
            coefficient = (pos['volume_fractional_difference'] + neg['volume_fractional_difference'])/(2*amplitude**2)
            near(coefficient, volume[-1]['volume_fractional_coefficient'], 'Volume at x=0.12', rel=.002)


def replay(jobs):
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    with tempfile.TemporaryDirectory(prefix='hdblast-gravity-') as tmp, ThreadPoolExecutor(max_workers=jobs) as pool:
        fresh = {name: Path(tmp)/(name+'.json') for name in PRODUCTS}
        def run(name, script, *args):
            print(f'Recomputing {script}', flush=True)
            subprocess.run([sys.executable, '-I', '-B', str(HERE/script), *map(str, args), '--output', str(fresh[name])],
                           env=env, cwd=HERE, check=True)
        previous = ['--previous', HERE/'previous']
        sources = {'INITIAL':'initial/constraint.py', 'TIDAL':'tidal/tidal.py',
                   'VOLUME':'perturbative/volume.py', 'EVOLUTION':'evolution/evolve.py'}
        tasks = [pool.submit(run, name, script, *previous) for name, script in sources.items()]
        for task in tasks:
            task.result()
        dependents = [pool.submit(run, 'BOUNCE', 'checks/bounce.py', *previous, '--initial', fresh['INITIAL'], '--evolution', fresh['EVOLUTION']),
                      pool.submit(run, 'COMPARISON', 'perturbative/compare.py', '--volume', fresh['VOLUME'], '--coupled', fresh['EVOLUTION'])]
        for task in dependents:
            task.result()
        results = {name: read(product) for name, product in fresh.items()}
        for name, result in results.items():
            compare(result, read(HERE/PRODUCTS[name]), (name,))
        conclusions(results)
        print(f'Verified {len(results)} numerical JSON products and independent physical comparisons.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true', help='Regenerate numerical results and compare all archived fields')
    parser.add_argument('--jobs', type=int, default=2, help='Maximum simultaneous numerical processes')
    args = parser.parse_args()
    require(args.jobs >= 1, '--jobs must be positive')
    authenticate()
    if args.verify:
        replay(args.jobs)
