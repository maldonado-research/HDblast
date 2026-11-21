"""Reconstruct the charged encounters, photon emissivities and formation tests."""
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

ROOT = Path(__file__).resolve().parent
RESEARCH = ROOT.parent
POLARIZATION = RESEARCH/'HDBLAST_POLARIZATION_20261121'
SCATTERING = RESEARCH/'HDBLAST_SCATTERING_20261121'
AXIAL = RESEARCH/'HDBLAST_AXIAL_CURRENT_20261121'
PRODUCTS = (
    'events/EVENTS.json',
    'soft/SOFT.json',
    'formation/FORMATION.json',
    'annihilation/ANNIHILATION.json',
    'hard/HARD.json',
    'audit/RADIATORS.json',
)


def check_manifest(base, manifest):
    names = set()
    for line in manifest.read_text().splitlines():
        digest, name = line.split('  ', 1)
        expected = base/name
        if name in names or not expected.resolve().is_relative_to(base) or expected.is_symlink():
            raise ValueError('Invalid checksum path '+name)
        with expected.open('rb') as stream:
            if hashlib.file_digest(stream, 'sha256').hexdigest() != digest:
                raise ValueError('Input checksum failure '+name)
        names.add(name)
    return names


def authenticate():
    files = check_manifest(ROOT, ROOT/'SHA256SUMS')
    available = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*')
                 if p.is_file() and p != ROOT/'SHA256SUMS'}
    if files != available:
        raise ValueError('Scientific file inventory differs from the manifest')
    inputs = check_manifest(RESEARCH.parent, ROOT/'parent-inputs.sha256')
    print(f'Authenticated {len(files)} study files and {len(inputs)} pinned predecessor inputs.', flush=True)


def run(script, *args):
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    print('Replaying '+script.parent.name+'/'+script.name, flush=True)
    subprocess.run([sys.executable, '-I', '-B', str(script), *map(str, args)], env=env, check=True)


def compare(actual, expected, path=()):
    if isinstance(expected, dict) and isinstance(actual, dict) and actual.keys() == expected.keys():
        for key in expected:
            yield from compare(actual[key], expected[key], (*path, key))
    elif isinstance(expected, list) and isinstance(actual, list) and len(actual) == len(expected):
        for index, (observed, reference) in enumerate(zip(actual, expected)):
            yield from compare(observed, reference, (*path, index))
    elif type(expected) is float and type(actual) in (int, float):
        metric = next((piece for piece in reversed(path) if isinstance(piece, str)), '')
        # Only normalized numerical accuracy checks can have an absolute floor.
        # Predicted photon amounts, even below 1e-28, must match relatively.
        accuracy = metric.endswith(('_error', '_residual', '_relative_to_spinor')) or metric in {
            'coarse_grid_relative_change_hard_count', 'coarse_grid_relative_change_total_energy',
            'relative_component_difference', 'relative_interference', 'charge_even_change',
            'reduced_vs_raw_relative', 'ward_raw_relative', 'relative_error',
            'maximum_relative_fsr_difference', 'maximum_relative_isr_difference',
            'half_nodes_energy_relative_change',
        }
        atol = 1e-10 if accuracy else 0.0
        if not math.isfinite(actual) or not math.isclose(actual, expected, rel_tol=2e-5, abs_tol=atol):
            yield (path, actual, expected)
    elif type(actual) != type(expected) or actual != expected:
        yield (path, actual, expected)


def validate(output):
    bad = []
    for name in PRODUCTS:
        observed = json.loads((output/name).read_bytes())
        reference = json.loads((ROOT/name).read_bytes())
        errors = list(compare(observed, reference))
        if errors:
            for path, actual, expected in errors[:6]:
                print(f"{name}/{'/'.join(map(str,path))}: rebuilt {actual}, reference {expected}", flush=True)
            bad.append(f'{name} ({len(errors)} discrepancies)')
        else:
            print('Reproduced '+name, flush=True)
    if bad:
        raise ValueError('Reproduction failed for '+', '.join(bad))


def replay(output, jobs):
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError('Choose an empty or new output directory')
    if any(output.is_relative_to(tree) for tree in (ROOT, AXIAL, POLARIZATION, SCATTERING)):
        raise ValueError('Write outputs outside the research inputs')
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='hdblast-photon-histories-') as work:
        temporary = Path(work)
        history = temporary/'dedicated_half.npz'
        tapes = temporary/'events'
        # Both independent collision instruments use the same reconstructed wall.
        run(AXIAL/'fields/extract.py', '--study', POLARIZATION, '--end', 7.0505,
            '--cases', 'dedicated_half', '--jobs', jobs, '--output-dir', temporary,
            '--output', temporary/'fields.json')
        args = ('--history', history, '--polarization', POLARIZATION, '--scattering', SCATTERING)
        instruments = [
            (ROOT/'events/replay.py', '--axial', AXIAL, *args, '--artifacts', tapes,
             '--jobs', max(1, jobs-jobs//2), '--output', output/'events/EVENTS.json'),
            (ROOT/'annihilation/integrate.py', '--study', AXIAL, *args,
             '--jobs', max(1, jobs//2), '--output', output/'annihilation/ANNIHILATION.json'),
        ]
        with ThreadPoolExecutor(max_workers=min(jobs, len(instruments))) as pool:
            tasks = [pool.submit(run, *instrument) for instrument in instruments]
            for task in tasks:
                task.result()
        run(ROOT/'soft/soft_photons.py', '--artifacts', tapes, '--output', output/'soft/SOFT.json')
        hard = output/'hard/HARD.json'
        run(ROOT/'hard/emission.py', '--artifacts', tapes, '--jobs', jobs, '--output', hard)
        run(ROOT/'formation/measure.py', '--events', tapes, '--hard', hard,
            '--output', output/'formation/FORMATION.json')
        run(ROOT/'audit/radiators.py', '--hard', ROOT/'hard', '--artifacts', tapes,
            '--reference', hard, '--output', output/'audit/RADIATORS.json')
    validate(output)
    for label, values in benchmarks(output).items():
        payload = json.dumps(values, sort_keys=True, allow_nan=False)
        print(label, payload, flush=True)
        if os.environ.get('GITHUB_ACTIONS') == 'true':
            print(f'::notice file=study.py,line=1,title=HDblast {label}::{payload}', flush=True)


def benchmarks(output):
    def fine(product, key='runs'):
        rows = json.loads((output/product).read_bytes())[key]
        return next(row for row in rows
                    if 't7.05000-p14-' in row.get('label', row.get('tape',''))
                    or (row.get('time') == 7.05 and row.get('seed') == 91283))

    events = fine('events/EVENTS.json')
    soft = fine('soft/SOFT.json')
    hard = fine('hard/HARD.json', 'collision_histories')
    formed = fine('formation/FORMATION.json')
    pairs = fine('annihilation/ANNIHILATION.json')
    audit = fine('audit/RADIATORS.json')
    band = next(row for row in formed['band_energy_fractions'] if row['xi_lo'] == .001)
    return {
        'Collisions on the shell at time 7.05': {
            'Recorded collisions': events['events'],
            'Encounters per particle': events['mean_encounters'],
            'Fraction inside the running resonance width': events['running_width_bands']['within']['event_fraction'],
            'Processed pair energy over final gas energy': formed['processed_pair_energy_over_gas']},
        'Photons per charged particle and radiated energy fractions': {
            'Electromagnetic fine-structure parameter': soft['alpha'],
            'Soft photons per charged particle per decade': soft['expected_photons_per_particle_per_decade'],
            'Single-photon energy over final gas energy': hard['total_energy_over_gas'],
            'Fraction of single-photon energy above one quarter endpoint': hard['bands'][-1]['fraction_of_endpoint_gas']/hard['total_energy_over_gas'],
            'One-dimensional radiator audit energy over gas': audit['total_energy_over_gas'],
            'Annihilation energy over final gas energy': pairs['photon_energy_over_final_gas'],
            'Annihilation energy ceiling over final gas energy': pairs['energy_bound_over_final_gas']},
        'Formation and energy-transfer diagnostics': {
            'Soft 0.1 to 1 percent band energy passing time separation': band['ten_radian_time'],
            'Soft 0.1 to 1 percent band energy also passing field criterion': band['integrated_force_envelope_and_ten_radians'],
            'Mean processed pair fraction required for ten percent of gas energy': formed['average_shell_fraction_for_ten_percent'],
            'Energy above one quarter endpoint passing both criteria': formed['above_quarter_hard_band']['joint_both_sides_fraction']}}



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true', help='Regenerate every numerical product')
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--jobs', type=int, default=2)
    args = parser.parse_args()
    if args.jobs < 1 or (args.output_dir and not args.verify):
        parser.error('Choose a positive job count and supply --verify with an output directory.')
    authenticate()
    if args.verify:
        if args.output_dir:
            replay(args.output_dir.resolve(), args.jobs)
        else:
            with tempfile.TemporaryDirectory(prefix='hdblast-photon-results-') as folder:
                replay(Path(folder), args.jobs)


if __name__ == '__main__':
    main()
