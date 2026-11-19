"""Authenticate this study and optionally regenerate its numerical products."""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory


HERE = Path(__file__).resolve().parent
REVISION = '148d545467542c5d14244ff6700d219f7e0fe533'
SOURCES_DIGEST = '7bdbcbc4c4f6dd7491ef01b64a40cb2b02c12e5526c977b727efbf6564b7ee72'
PRODUCTS = {'rollout': 'rollout/ROLLOUT.json', 'patches': 'patches/PATCHES.json',
            'spectrum': 'spectrum/SPECTRUM.json', 'checks': 'checks/SPECTRAL.json',
            'intrinsic': 'intrinsic/FREQUENCIES.json', 'dynamic': 'dynamic/DYNAMIC.json',
            'causality': 'causality/NULL.json', 'interpretation': 'interpretation/UNFOLDING.json',
            'passage': 'passage/PASSAGE.json'}
SCRIPTS = {'rollout': 'rollout/rollout.py', 'patches': 'patches/patches.py',
           'spectrum': 'spectrum/spectrum.py', 'checks': 'checks/spectral.py',
           'intrinsic': 'intrinsic/frequencies.py', 'dynamic': 'dynamic/spectator.py',
           'causality': 'causality/null_cones.py', 'interpretation': 'interpretation/unfolding.py',
           'passage': 'passage/passage.py'}
CSV_PRODUCT = 'spectrum/SPECTRUM.csv'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read(path):
    return json.loads(path.read_bytes())


def authenticate():
    inventory = HERE / 'SHA256SUMS'
    listed = set()
    for line in inventory.read_text().splitlines():
        checksum, name = line.split('  ', 1)
        path = HERE / name
        require(path.is_relative_to(HERE) and Path(name).as_posix() == name
                and not Path(name).is_absolute() and path == path.resolve()
                and path != inventory and path not in listed, f'Invalid inventory entry {name}')
        require(digest(path) == checksum, f'Checksum mismatch for {name}')
        listed.add(path)
    present = {p for p in HERE.rglob('*') if p.is_file() and p != inventory
               and not (p.parent.name == '__pycache__' and p.suffix == '.pyc')}
    require(listed == present, f'Study inventory differs: {sorted(listed ^ present)}')
    source_manifest = HERE / 'SOURCES.json'
    require(digest(source_manifest) == SOURCES_DIGEST, 'Historical source manifest mismatch')
    pins = read(source_manifest)
    require(pins['base'] == REVISION and pins['previous_doi'] == '10.5281/zenodo.24378925',
            'Historical revision or parent DOI differs')
    historical = set()
    for entry in pins['sources']:
        path = HERE / entry['bundled']
        require(path in listed and path not in historical and path.stat().st_size == entry['bytes']
                and digest(path) == entry['sha256'], f'Historical input mismatch: {entry["path"]}')
        historical.add(path)
    supplied = {p for p in listed if p.is_relative_to(HERE/'gravity') or p.is_relative_to(HERE/'modulus')}
    require(supplied == historical, 'Historical input coverage differs')
    for name in [*PRODUCTS.values(), *SCRIPTS.values(), CSV_PRODUCT]:
        require(HERE / name in listed, f'Missing numerical product or replay program: {name}')
    print(f'Verified {len(listed)} study files, including {len(historical)} historical inputs.', flush=True)


def tolerance(path, guarded):
    product, field = path[0], path[-1]
    if product in ('spectrum', 'csv', 'checks'):
        if ('error' in field or 'residual' in field or 'fractional_change' in field
                or 'fractional_difference' in field):
            return 1e-6, 2e-10
        if field in ('fractional_pair_detuning', 'detuning_times_duration'):
            return 2e-8, 3e-13
        if 'occupation' in field:
            return 2e-5, 1e-24
        if product == 'spectrum' and 'density' in field:
            return 5e-7, 1e-19
        return 2e-7, 1e-27
    if product == 'rollout':
        if field in ('ham_max', 'momentum_max', 'max_step_sampled_ham', 'max_step_sampled_momentum'):
            return .03, 2e-8
        if field == 'constraint_guard' or (guarded and field == 'x'):
            return 0., 1e-3
        if guarded:
            # Coarse constraint-triggered stops can move in time with BLAS summation.
            return .02, .005
        if 'old_cap_differences' in path or 'verification' in path:
            return 1e-5, 5e-9
        if field == 'x' or 'event_times' in path:
            return 0., 3e-7
        return 3e-6, 2e-8
    if product == 'patches':
        if 'refinement' in path or 'maximum_sampled_constraints' in path:
            return 2e-6, 4e-9
        if field in ('fraction', 'integrated_outward_energy_flux'):
            return 2e-6, 3e-13
        return 2e-7, 2e-8
    if product == 'intrinsic':
        if 'error' in field:
            return 2e-5, 1.5e-8 if field == 'max_finite_difference_error' else 5e-11
        return 2e-7, 3e-11
    if product == 'causality':
        return 2e-6, 1e-8
    if product == 'interpretation':
        if 'error' in field or 'difference' in field or 'refinement' in field:
            return 1e-5, 4e-11
        return 2e-8, 1e-16
    if product in ('dynamic', 'passage') and (any('constraints' in part for part in path)
                                               or field == 'maximum_sampled_normalized_physical_constraint'):
        # Spherical differentiation at the poles amplifies roundoff in this diagnostic.
        return .05, 4e-8 if 'regular_momentum' in field else 2e-9
    # The finite-band solvers also enforce method-specific absolute tolerances.
    return 2e-6, 2e-11


def compare(actual, archived, path, guarded=False):
    label = '/'.join(path)
    if isinstance(archived, dict):
        require(isinstance(actual, dict) and actual.keys() == archived.keys(), f'Object fields differ at {label}')
        for key, value in archived.items():
            guard_ending = key == 'end' and archived.get('stop') == 'constraint_guard'
            compare(actual[key], value, (*path, key), guarded or guard_ending)
    elif isinstance(archived, list):
        require(isinstance(actual, list) and len(actual) == len(archived), f'Array length differs at {label}')
        for index, (new, old) in enumerate(zip(actual, archived)):
            compare(new, old, (*path, str(index)), guarded)
    elif path[0] == 'csv' and path[-1] == 'exact_occupation':
        # The ODE controls absolute error in beta even for a near-zero UV mode.
        require(type(actual) in (int, float) and math.isfinite(actual) and actual >= 0
                and archived >= 0 and math.isclose(math.sqrt(actual), math.sqrt(archived),
                                                    rel_tol=1e-5, abs_tol=3e-11),
                f'{label}: computed {actual!r}, archived {archived!r}')
    elif isinstance(archived, float):
        rel, absolute = tolerance(path, guarded)
        require(type(actual) in (int, float) and math.isfinite(actual) and math.isfinite(archived)
                and math.isclose(actual, archived, rel_tol=rel, abs_tol=absolute),
                f'{label}: computed {actual!r}, archived {archived!r}')
    elif path[-1] in ('time_steps', 'hard_domain_retries', 'ode_evaluations'):
        require(type(actual) is int and actual >= 0 and abs(actual-archived) <= max(8, .06*archived),
                f'Adaptive integrator count differs at {label}: {actual} versus {archived}')
    else:
        require(type(actual) is type(archived) and actual == archived, f'Scalar value differs at {label}')


def read_rows(path):
    integers = {'cycles', 'n', 'degeneracy'}
    with path.open(newline='') as stream:
        return [{k: (v if k in ('branch', 'radius_role') else int(v) if k in integers else float(v))
                 for k, v in row.items()} for row in csv.DictReader(stream)]


def near(a, b, label, relative, absolute=0.):
    require(math.isclose(a, b, rel_tol=relative, abs_tol=absolute), f'{label}: {a} versus {b}')


def conclusions(data, rows):
    spectrum, reference = data['spectrum'], data['checks']
    for cycles, row in spectrum['continua'].items():
        independent = reference['cases'][cycles]
        # Oscillator quadrature cuts off at k/m=4. The reference integrates to infinity.
        tail = independent['continuum_energy_tail_upper_above_k4']
        near(row['born']['energy_density_over_m4'], independent['energy'], 'Continuum energy',
             8e-12, tail)
        near(row['born']['number_density_over_m3'], independent['number'], 'Parseval number', 5e-10)
        near(row['flat_width_energy_over_m4'], independent['harmonic_width_energy'], 'Width exposure', 2e-13)
    groups = defaultdict(list)
    for row in rows:
        groups[row['branch'], row['radius_role'], row['cycles']].append(row)
    require(len(groups) == len(spectrum['cases']), 'CSV case count differs')
    for case in spectrum['cases']:
        group = groups[case['branch'], case['radius_role'], case['cycles']]
        require([r['n'] for r in group] == list(range(1, case['n_max']+1)), 'Incomplete shell list')
        for row in group:
            n = row['n']
            require(row['degeneracy'] == (n+1)**2, 'Wrong real-harmonic multiplicity')
            near(row['omega_over_m'], math.sqrt(n*(n+2))/case['radius_m_inverse'], 'Shell frequency', 1e-13)
        for method in ('born', 'exact'):
            key = method + '_occupation'
            volume = case['volume_m_inverse_cubed']
            number = math.fsum(r['degeneracy']*r[key] for r in group) / volume
            energy = math.fsum(r['degeneracy']*r[key]*r['omega_over_m'] for r in group) / volume
            near(case[method]['number_density_over_m3'], number, 'CSV number sum', 2e-12)
            near(case[method]['energy_density_over_m4'], energy, 'CSV energy sum', 2e-12)
    rollout = data['rollout']
    for branch in rollout['branches'].values():
        for levels in branch['runs'].values():
            for run in levels.values():
                residual = max(run['max_step_sampled_ham'], run['max_step_sampled_momentum'])
                require(residual <= rollout['constraint_limit'] * (1+2e-8),
                        'A sampled trajectory exceeded its terminal constraint limit')
    for label, traces in data['causality'].items():
        endpoint = rollout['branches']['phi_m=0.5']['runs'][label]['near_edge']['end']['x']
        for trace in traces.values():
            near(trace['guard_x'], endpoint, 'Causal computation uses the refined rollout endpoint', 0., 2e-12)
            reach = trace['coordinate_null_reach_upper']
            near(trace['initial_domain_gap_lower'], math.pi/2-2*reach, 'Null footprint gap', 2e-12)
            require(trace['initial_domain_gap_lower'] > 0, 'Opposite causal footprints can overlap')
    for case in data['interpretation']['cases'].values():
        near(case['early_log_slope'], 1.5, 'Early isolated quadratic contact', 0., 1e-5)
        near(case['late_log_slope'], .5, 'Late angularly averaged contact', 0., .002)
    for patch in data['patches']['runs']:
        case = data['rollout']['branches'][data['patches']['branch']]['runs'][f'{patch["amplitude"]:g}']['fine']
        crossing = case['events']['first_minimum']
        near(patch['crossing_x'], crossing['x'], 'Independent first-passage time', 0., 4e-7)
        near(patch['normal_velocity'], crossing['first_worldline_Pi'], 'First-passage velocity', 3e-6)
    for branch in data['intrinsic'].values():
        for frozen in branch['cases']['0.0']:
            for sector in ('geometric', 'matter_weighted'):
                for n, group in frozen[sector].items():
                    omega = math.sqrt(int(n)*(int(n)+2)) / frozen['volume_radius']
                    for measured in group['by_ell']:
                        near(measured, omega, 'Frozen homogeneous radial mode', 3e-11)

    # Closed de Sitter transfer evaluated in the original-field instantaneous
    # Hamiltonian convention admits a compact direct occupation formula.
    for name, branch in data['dynamic']['branches'].items():
        h = 1 / data['intrinsic'][name]['cases']['0.0'][0]['volume_radius']
        for run in branch['runs']['0'].values():
            for sample in run['samples'].values():
                s = math.sinh(h * sample['x'])
                eta = math.atan(s)
                for label, values in sample['low_modes'].items():
                    n = int(label)
                    ell = n+1
                    sin_term = math.sin(ell*eta)/ell
                    cos_term = math.cos(ell*eta)
                    number = ((sin_term-s*cos_term)**2/(n*(n+2)) + (s*sin_term)**2)/4
                    near(values['occupation'], number, 'Closed de Sitter instantaneous occupation',
                         5e-9, 2e-14)
    require(data['dynamic']['checks']['passed'], 'Finite-band dynamical controls failed')

    passage = data['passage']
    require(passage['checks']['passed'], 'Passage field or geometry controls failed')
    branch = rollout['branches'][passage['branch']]
    for label, record in passage['runs'].items():
        independent = branch['runs'][label]['fine']['events']['first_minimum']['x']
        times = record['sample_times']
        require(times['before'] < independent < times['after'], 'Test-field samples missed the first passage')
        near(record['crossing_x'], independent, 'Independent test-field crossing time', 0., 3e-7)
        near(times['crossing'], record['crossing_x'], 'Central sample time', 0., 2e-11)
        for side in ('before', 'after'):
            near(abs(times[side]-times['crossing']), passage['sample_half_width'], 'Passage bracket width',
                 0., 2e-11)
        for distance in ('signed_worldline_distance_to_minimum', 'global_signed_distance_to_minimum'):
            geometry = record['geometry_samples']
            require(geometry['before'][distance] > 0 > geometry['after'][distance],
                    'Passage bracket failed the signed field test')
            near(geometry['crossing'][distance], 0., 'Passage contact field distance', 0., 1e-10)
        for scheme in ('matter', 'geometric'):
            result = record['kinetic_controls'][scheme]
            for phase, sample in result['samples'].items():
                near(sample['x'], times[phase], 'Field evolution sample time', 0., 2e-11)
                for n in (2, 4, 8):
                    index = result['labels'].index(n)
                    beta_row = sample['beta_abs_squared'][index]
                    mode = sample['low_modes'][str(n)]
                    require(all(x >= 0 for x in beta_row), 'Negative pair-mixing weight')
                    near(mode['occupation'], math.fsum(beta_row), 'Passage occupation from mixing matrix',
                         3e-12, 2e-20)
                    near(mode['diagonal'], beta_row[index], 'Passage diagonal mode contribution',
                         3e-12, 2e-20)
                    near(mode['from_other_initial_modes'], math.fsum(beta_row[:index]+beta_row[index+1:]),
                         'Passage off-diagonal mode contribution', 3e-10, 3e-17)


def replay(jobs):
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    with TemporaryDirectory(prefix='hdblast-crossing-') as tmp:
        root = Path(tmp)
        def generate(name):
            relative = Path(PRODUCTS[name])
            directory = root / relative.parent
            directory.mkdir(exist_ok=True)
            output = directory if name == 'spectrum' else root / relative
            command = [sys.executable, '-I', '-B', str(HERE/SCRIPTS[name])]
            if name not in ('checks', 'interpretation'):
                command.extend(['--gravity', str(HERE/'gravity')])
            if name == 'causality':
                command.extend(['--rollout', str(root / PRODUCTS['rollout'])])
            command.extend(['--output', str(output)])
            print(f'Recomputing {SCRIPTS[name]}', flush=True)
            completed = subprocess.run(command, cwd=HERE, env=env, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, text=True)
            if completed.returncode:
                raise RuntimeError(f'{name} failed ({completed.returncode})\n{completed.stdout[-4500:]}')
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            tasks = {name: pool.submit(generate, name) for name in PRODUCTS if name != 'causality'}
            tasks['rollout'].result()
            tasks['causality'] = pool.submit(generate, 'causality')
            for task in tasks.values():
                task.result()
        regenerated = {name: read(root/path) for name, path in PRODUCTS.items()}
        for name, result in regenerated.items():
            compare(result, read(HERE/PRODUCTS[name]), (name,))
        rows = read_rows(root/CSV_PRODUCT)
        compare(rows, read_rows(HERE/CSV_PRODUCT), ('csv',))
        conclusions(regenerated, rows)
    print('Verified nine numerical JSON products, the per-mode CSV and independent physical comparisons.', flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true', help='Regenerate all numerical products')
    parser.add_argument('--jobs', type=int, default=2, help='Maximum simultaneous numerical processes')
    args = parser.parse_args()
    require(args.jobs >= 1, '--jobs must be positive')
    authenticate()
    if args.verify:
        replay(args.jobs)
