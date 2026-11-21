"""Check the study files, or recreate the field, source and particle calculations."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(__file__).resolve().parent
PARENT = ROOT / "parent"
PARENT_DIGEST = "582302b0a445bf697a21c89b9aa4ac1f2afc9e380830f7731761e30d8ca7c1ec"
PRODUCTS = (
    "fields/FIELDS.json",
    "fields/STORAGE_SCAN.json",
    "fields/IMPULSE.json",
    "source/SOURCE.json",
    "sensitivity/CHECK.json",
    "orbits/PAIRS.json",
    "bound/WEAK.json",
    "bound/MIXTURES.json",
    "fields/SPIN.json",
    "orbits/LIOUVILLE.json",
    "orbits/COARSE.json",
    "feedback/FEEDBACK.json",
)


def digest(path):
    with path.open("rb") as file:
        return hashlib.file_digest(file, "sha256").hexdigest()


def inventory(directory):
    records = {}
    for line in (directory / 'SHA256SUMS').read_text().splitlines():
        value, name = line.split('  ', 1)
        if name in records or Path(name).is_absolute() or '..' in Path(name).parts:
            raise ValueError(f'Invalid manifest path {name}')
        records[name] = value
    return records


def authenticate():
    registered = inventory(ROOT)
    files = {file.relative_to(ROOT).as_posix(): file for file in ROOT.rglob('*')
             if file.is_file() and '__pycache__' not in file.parts and file != ROOT/'SHA256SUMS'}
    if files.keys() != registered.keys():
        raise ValueError(f'Unexpected study files {files.keys() ^ registered.keys()}')
    for name, file in files.items():
        if file.is_symlink() or digest(file) != registered[name]:
            raise ValueError(f'File differs from study manifest: {name}')
    if digest(PARENT/'SHA256SUMS') != PARENT_DIGEST:
        raise ValueError('Prebirth study manifest differs')
    prior = inventory(PARENT)
    inputs = [name for name in registered if name.startswith('parent/') and name != 'parent/SHA256SUMS']
    for name in inputs:
        if registered[name] != prior.get(name.removeprefix('parent/')):
            raise ValueError(f'Prebirth input differs: {name}')
    print(f'Verified {len(registered)} study files and {len(inputs)} predecessor inputs.', flush=True)


def run(script, *args):
    print('Replaying', script, flush=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    subprocess.run([sys.executable, '-I', '-B', str(ROOT/script), *map(str,args)], check=True, env=env)


def regenerate_model(tmp, jobs, predecessor):
    # Copying the small, pinned reference inputs permits replacing the compact histories
    # with freshly regenerated samples for all downstream calculations.
    model = tmp/'model'
    shutil.copytree(PARENT, model)
    native = tmp/'native'
    run('parent/parent/parent/history/export.py', '--prior', model/'parent/parent/prior',
        '--outdir', native, '--jobs', jobs)
    run('parent/parent/history/extract.py', '--native', native, '--output-dir', model/'parent/history')
    predecessor['compare_history'](model/'parent/history')
    shutil.rmtree(native)
    return model


def compare_product(product, actual, reference, compare):
    if product == 'source/SOURCE.json':
        # The absolute sheet-balance residual subtracts accumulated energies.
        # CPU contraction changes its roundoff floor. Keep both original JSON files.
        # Each run must satisfy the limit regardless of cross-platform agreement.
        for key, record in reference.items():
            expected = record['checks'].pop('sheet_balance_error')
            value = actual[key]['checks'].pop('sheet_balance_error')
            if not (0 <= value < 1e-7 and math.isclose(value, expected, rel_tol=2e-6, abs_tol=1e-8)):
                raise ValueError(f'{product}/{key}: sheet-balance residual {value}, reference {expected}')
    compare(actual, reference, product)


def replay(dest, jobs):
    dest = dest.resolve()
    if dest.is_relative_to(ROOT) or (dest.exists() and (not dest.is_dir() or any(dest.iterdir()))):
        raise ValueError('Choose a new or empty output directory outside the study')
    dest.mkdir(parents=True, exist_ok=True)
    prior = runpy.run_path(str(PARENT/'study.py'))
    with tempfile.TemporaryDirectory(prefix='hdblast-polarization-history-') as work:
        model = regenerate_model(Path(work), jobs, prior)
        history = model/'parent/history'
        fields = dest/'fields'
        run('fields/launch.py', '--parent', model, '--history', history,
            '--output-dir', fields, '--jobs', jobs)
        scan = fields/'scan'
        run('fields/launch.py', '--parent', model, '--history', history,
            '--output-dir', scan, '--jobs', jobs, '--scan')
        run('fields/impulse.py', '--fields', fields, '--scan', scan, '--history', history,
            '--output-dir', fields)
        run('source/backreaction.py', '--parent', model, '--profiles', fields/'FIELDS.json',
            '--output-dir', dest/'source', '--jobs', jobs)
        for component in ('orbits', 'feedback', 'bound', 'sensitivity'):
            (dest/component).mkdir()
        run('sensitivity/energy.py', '--parent', model, '--source', dest/'source/SOURCE.json',
            '--output', dest/'sensitivity/CHECK.json', '--jobs', jobs)
        run('orbits/pairs.py', '--history', history, '--fields', fields/'traces',
            '--scan', scan/'traces', '--source', dest/'source/source_traces',
            '--output', dest/'orbits/PAIRS.json', '--jobs', jobs)
        run('bound/weak_pressure.py', '--history', history, '--fields', fields/'traces',
            '--output', dest/'bound/WEAK.json')
        run('bound/mixtures.py', '--history', history, '--fields', fields/'traces',
            '--gram', fields/'FIELDS.json', '--output', dest/'bound/MIXTURES.json')
        run('feedback/currents.py', '--history', history, '--fields', fields/'traces',
            '--output', dest/'feedback/FEEDBACK.json', '--jobs', jobs)
        run('fields/spin.py', '--history', history, '--scan', scan, '--output-dir', fields)
        run('orbits/liouville.py', '--history', history, '--fields', fields/'traces', '--scan', scan/'traces',
            '--source', dest/'source/source_traces', '--pairs', dest/'orbits/PAIRS.json',
            '--output', dest/'orbits/LIOUVILLE.json')
        run('orbits/coarse.py', '--history', history, '--fields', fields/'traces', '--scan', scan/'traces',
            '--pairs', dest/'orbits/PAIRS.json', '--output', dest/'orbits/COARSE.json', '--jobs', jobs)
        (scan/'STORAGE_SCAN.json').replace(fields/'STORAGE_SCAN.json')
    for product in PRODUCTS:
        compare_product(product, json.loads((dest/product).read_bytes()),
                        json.loads((ROOT/product).read_bytes()), prior['compare'])
        print('Reproduced', product, flush=True)
    print_benchmarks(dest)
    print(f'Rebuilt all four reference histories and all {len(PRODUCTS)} numerical products.', flush=True)


def print_benchmarks(output):
    for name, values in reported_values(output).items():
        row = json.dumps(values, sort_keys=True, allow_nan=False)
        print('Benchmark', name, row, flush=True)
        if os.environ.get('GITHUB_ACTIONS') == 'true':
            print(f'::notice file=study.py,line=1,title=HDblast {name}::{row}', flush=True)


def reported_values(output):
    mixtures = json.loads((output/'bound/MIXTURES.json').read_bytes())
    impulse = json.loads((output/'fields/IMPULSE.json').read_bytes())
    orbits = json.loads((output/'orbits/PAIRS.json').read_bytes())
    coarse = json.loads((output/'orbits/COARSE.json').read_bytes())['runs']
    feedback = json.loads((output/'feedback/FEEDBACK.json').read_bytes())['runs']
    benchmarks = {}
    for interaction in ('dedicated', 'shared'):
        label = interaction + '_half'
        mixture = mixtures[label]['tests']['2,5,7,40,258']['0.5']['0.04']
        fraction = mixtures[label]['particle_fraction_at_sites_below_lateral_target']['0.04']['0.9']
        high_modes = impulse['preparation_scan'][label]['localized_high_mode_cdf']['launch_amplitudes']['1.0']
        uniform = next(run['metrics'] for run in orbits['runs']
                       if run['case'] == 'uniform_local_pressure_match' and run['history'] == label)
        fitted_runs = (run for run in coarse if run['case'] == 'uniform_local_pressure_match' and run['history'] == label)
        counts = max(fitted_runs, key=lambda run: run['pairs'])
        pooled = next(run for run in orbits['conditioned']
                      if run['case'] == 'pooled_pressure_match' and run['history'] == label)
        current = next(run for run in feedback if run['history'] == label and run['case'] == 'two_pressure_fit'
                       and run['feedback'] and run['steps'] == 512)
        benchmarks[interaction] = {
            'Total variation in 20 bins after feedback from both transverse currents, percentage points':
                100 * current['coarse_cell_total_variation'],
            'Total variation distance in 20 energy-direction bins for the uniform pressure-matched control, percentage points':
                100 * counts['coarse_total_variation'],
            'Local smaller pressure when the two global pressures match, percent of thermal pressure':
                100 * pooled['local_lateral_pressure_over_thermal'][0],
            'Transverse count discrepancy for the uniform pressure-matched control, percentage points':
                100 * uniform['ks_transverse'],
            'Emitted energy of the uniform pressure-matched control, multiples of the reference':
                orbits['source_control_emitted_energy_over_reference'][label],
            'Minimum transverse count discrepancy for 160 high-mode source coefficients at unit RMS, percentage points':
                100 * high_modes['guaranteed_transverse_cdf_excess'],
            'Minimum particle fraction at sites with weak pressure below 90% of thermal at 4% source RMS, percent':
                100 * fraction,
            'Weak-direction pressure upper bound for 70 source coefficients on half the sphere at 4% RMS, percent of thermal pressure':
                100 * mixture['pressure_upper_over_thermal'],
        }
    source = json.loads((output/'source/SOURCE.json').read_bytes())
    early_birth = source['simultaneous_time_control']['events']['birth']
    delayed = source['delayed_matched_half']
    end = delayed['events']['end']
    if 'birth' in delayed['events']:
        raise ValueError('Equal-energy delayed source produced daughters in the recorded interval')
    doses = json.loads((output/'sensitivity/CHECK.json').read_bytes())['runs']
    extended = next(trial['one_late'] for trial in doses if trial['mesh']=='half' and trial['input_dose_scale']==1)
    nearby = [trial for trial in doses if trial['input_dose_scale']!=1]
    delayed_count = sum(trial['one_late']['first_birth'] is not None for trial in doses)
    nearby_births = sum(trial[timing]['first_birth'] is not None for trial in nearby for timing in ('both_early','one_late'))
    benchmarks['source'] = {
        'End time of the extended delayed benchmark with no births': extended['source_end'],
        'Delayed trials producing any daughters by the extended end time, three doses across three grids': delayed_count,
        'Nearby-dose trials producing daughters through the extended end time, both schedules across three grids': nearby_births,
        'First birth time, early equal-energy source': early_birth['time'],
        'Time reached with no births, delayed equal-energy source': end['time'],
        'Fraction of emitted energy reserved for later delivery, percent': 100 * end['emitted2'] / end['emitted'],
    }
    return benchmarks


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true', help='Rebuild the histories and calculations')
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--jobs', type=int, default=2)
    args = parser.parse_args()
    if args.jobs < 1 or (args.output_dir and not args.verify):
        parser.error('--jobs must be positive and --output-dir requires --verify')
    authenticate()
    if args.verify and args.output_dir:
        replay(args.output_dir, args.jobs)
    elif args.verify:
        with tempfile.TemporaryDirectory(prefix='hdblast-polarization-results-') as tmp:
            replay(Path(tmp), args.jobs)
