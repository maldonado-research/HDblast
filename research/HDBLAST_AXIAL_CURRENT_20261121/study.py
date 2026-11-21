"""Reconstruct the prescribed source, kinetic current and angular diagnostics."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import runpy
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
POLARIZATION = ROOT.parent / 'HDBLAST_POLARIZATION_20261121'
SCATTERING = ROOT.parent / 'HDBLAST_SCATTERING_20261121'
DIGESTS = {
    POLARIZATION: '38c06ac62f996fbfbbd2e106e43c67de2896607500787630c617ef610e42916b',
    SCATTERING: 'cff83423b0b81e2bd54ffa353131545b9a2b6acc81246e50776f18853a89a160',
}
PRODUCTS = (
    'fields/FIELDS.json',
    'live/LIVE.json',
    'live/CROSSINGS.json',
    'live/COHORTS.json',
    'flavor/FLAVOR.json',
    'flavor/EXTENDED.json',
    'independent/INDEPENDENT.json',
    'independent/LATE.json',
    'freeze/FREEZE.json',
    'triad/TRIAD.json',
    'diagnostic/CURRENT.json',
    'diagnostic/MEMORY.json',
    'return/RETURN.json',
    'return/EXTENDED.json',
)


def authenticate():
    # Authenticate the checker itself before reusing it on the current tree.
    for path, expected in DIGESTS.items():
        manifest = (path / 'SHA256SUMS').read_bytes()
        if hashlib.sha256(manifest).hexdigest() != expected:
            raise ValueError('Historical manifest differs in ' + path.name)
        for line in manifest.decode().splitlines():
            digest, name = line.split('  ', 1)
            if Path(name).is_absolute() or '..' in Path(name).parts:
                raise ValueError('Invalid historical path')
            with (path / name).open('rb') as stream:
                if hashlib.file_digest(stream, 'sha256').hexdigest() != digest:
                    raise ValueError('Historical input differs: ' + name)
    check = runpy.run_path(str(SCATTERING / 'study.py'))['authenticate']
    for path in (*DIGESTS, ROOT):
        print(f'Authenticated {check(path)} files in {path.name}.', flush=True)


def run(script, *args):
    print('Replaying ' + script, flush=True)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    subprocess.run([sys.executable, '-I', '-B', str(ROOT / script), *map(str, args)],
                   check=True, env=env)


def replay(output, jobs):
    output = output.resolve()
    if any(output.is_relative_to(path) for path in (ROOT, *DIGESTS)):
        raise ValueError('Replay output must be outside the input studies')
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError('Choose a new or empty replay output directory')
    output.mkdir(parents=True, exist_ok=True)
    for folder in {Path(product).parent for product in PRODUCTS}:
        (output / folder).mkdir()

    with tempfile.TemporaryDirectory(prefix='hdblast-current-inputs-') as temp:
        inputs = Path(temp)
        run('fields/extract.py', '--study', POLARIZATION, '--output-dir', inputs,
            '--output', output/'fields/FIELDS.json', '--jobs', jobs)
        run('diagnostic/current_bounds.py', '--inputs', inputs,
            '--output', output/'diagnostic/CURRENT.json')
        run('live/evolve.py', '--polarization', POLARIZATION, '--scattering', SCATTERING,
            '--history', inputs/'dedicated_half.npz', '--output', output/'live/LIVE.json', '--jobs', jobs)
        run('flavor/shared.py', '--polarization', POLARIZATION, '--scattering', SCATTERING,
            '--history', inputs/'shared_half.npz', '--output', output/'flavor/FLAVOR.json', '--jobs', jobs)
        run('independent/run.py', '--polarization', POLARIZATION, '--scattering', SCATTERING,
            '--history', inputs/'shared_half.npz', '--output', output/'independent/INDEPENDENT.json',
            '--jobs', jobs)
        extended = inputs/'extended'
        run('fields/extract.py', '--study', POLARIZATION, '--end', 7.0505,
            '--output-dir', extended, '--output', extended/'fields.json', '--jobs', jobs)
        run('live/follow.py', '--polarization', POLARIZATION, '--scattering', SCATTERING,
            '--history', extended/'dedicated_half.npz', '--output', output/'live/CROSSINGS.json',
            '--jobs', jobs)
        run('live/cohorts.py', '--polarization', POLARIZATION, '--scattering', SCATTERING,
            '--history', extended/'dedicated_half.npz', '--output', output/'live/COHORTS.json',
            '--jobs', jobs)
        run('flavor/shared.py', '--endpoint', 7.05, '--polarization', POLARIZATION,
            '--scattering', SCATTERING, '--history', extended/'shared_half.npz',
            '--output', output/'flavor/EXTENDED.json', '--jobs', jobs)
        run('independent/run.py', '--extended', '--polarization', POLARIZATION, '--scattering', SCATTERING,
            '--history', extended/'shared_half.npz', '--output', output/'independent/LATE.json', '--jobs', jobs)
        snapshots = output/'live/snapshots'
        run('freeze/compare.py', '--polarization', POLARIZATION, '--scattering', SCATTERING,
            '--snapshots', snapshots, '--extended', output/'live/cross-snapshots',
            '--output', output/'freeze/FREEZE.json', '--jobs', jobs)
        run('diagnostic/memory.py', '--inputs', inputs, '--snapshots', snapshots,
            '--output', output/'diagnostic/MEMORY.json')
        geometry = inputs/'geometry'
        run('return/geometry.py', '--polarization', POLARIZATION, '--output-dir', geometry)
        returned = output/'return/snapshots'
        run('return/coupled.py', '--polarization', POLARIZATION, '--scattering', SCATTERING,
            '--history', inputs/'dedicated_half.npz', '--live-output', output/'live',
            '--geometry', geometry/'dedicated.npz', '--snapshots', returned,
            '--output', output/'return/RETURN.json', '--jobs', jobs)
        late_geometry = extended/'geometry'
        run('return/geometry.py', '--polarization', POLARIZATION, '--end', 7.05,
            '--output-dir', late_geometry)
        run('return/coupled.py', '--extended', '--polarization', POLARIZATION, '--scattering', SCATTERING,
            '--history', extended/'dedicated_half.npz', '--live-output', output/'live',
            '--geometry', late_geometry/'dedicated.npz', '--snapshots', returned,
            '--output', output/'return/EXTENDED.json', '--jobs', jobs)
        run('triad/diagnose.py', '--snapshots', snapshots, '--shared-snapshots', output/'flavor/snapshots',
            '--cross-snapshots', output/'live/cross-snapshots', '--return-snapshots', returned,
            '--output', output/'triad/TRIAD.json')
    compare = runpy.run_path(str(POLARIZATION / 'parent/study.py'))['compare']
    for product in PRODUCTS:
        regenerated = json.loads((output/product).read_bytes())
        archived = json.loads((ROOT/product).read_bytes())
        compare_product(product, regenerated, archived, compare)
        print('Reproduced', product, flush=True)

    for label, values in benchmarks(output).items():
        payload = json.dumps(values, sort_keys=True, allow_nan=False)
        print(label, payload, flush=True)
        if os.environ.get('GITHUB_ACTIONS') == 'true':
            print(f'::notice file=study.py,line=1,title=HDblast {label}::{payload}', flush=True)


def compare_product(product, regenerated, archived, compare):
    # Solve_ivp's internal step selection varies slightly with the host. The
    # independent solve enforces its own 5e-6 relative-error bound in both runs.
    if product == 'diagnostic/CURRENT.json':
        for label in ('dedicated_half', 'shared_half'):
            key = 'independent_solver_scaled_error'
            if abs(regenerated[label][key]-archived[label][key]) > 1.e-8:
                raise ValueError(f'{label}: independent solver residual changed')
            regenerated[label][key] = archived[label][key]
    compare(regenerated, archived, product)


def benchmarks(output):
    def load(path):
        return json.loads((output/path).read_bytes())

    early_label = 't6.67000-p14-n4096-r91283-particle-pole'
    late_label = 't7.05000-p14-n8192-r91283-particle-pole'
    early = next(r for r in load('live/LIVE.json')['runs'] if r['label'] == early_label)
    late = next(r for r in load('live/CROSSINGS.json')['runs'] if r['label'] == late_label)
    cohorts = next(r for r in load('live/COHORTS.json')['runs'] if r['label'] == late_label)
    returned = next(r for r in load('return/EXTENDED.json')['cases'] if r['label'] == late_label)
    final = returned['iterations'][-1]['wall']
    directions = load('triad/TRIAD.json')
    gap = next(r['number_weighted']['any_boost_gap'] for r in directions['runs']
               if r['snapshot'] == early_label)
    shared = {r['flavor_mode']: r['number_weighted']['any_boost_gap']
              for r in directions['shared_runs'] if r['snapshot'] == early_label}
    regular, converted = sorted((r for r in load('flavor/EXTENDED.json')['runs']
                                if r['pairs'] == 16384), key=lambda row: row['flavor_change'])
    free = next(r for r in load('independent/LATE.json')['runs']
                if r['shared_mode'] == 'uniform' and r['pairs_per_birth_flavor'] == 4096)
    frozen = {}
    for row in load('freeze/FREEZE.json')['runs']:
        if row['pairs'] == 8192 and row['time'] == 6.67:
            frozen[row['flat_cross_section']] = next(x for x in row['checks'] if x['encounters'] == 8)
    return {
        'Axial field at shell time 6.67': {
            'Fluid reference': early['end_reference_E'][2],
            'Particle current': early['end_E'][2],
            'Particle first reversal': late['first_field_reversals'][0],
            'Particle second reversal': late['first_field_reversals'][1]},
        'Dedicated scattering at shell time 7.05': {
            'Encounters per particle': late['mean_encounters'], 'Fraction ever hit': late['hit'],
            'Largest pressure over energy': max(late['rest_pressures_over_energy']),
            'Common thermal histogram gap': late['common_gas_tv20'],
            'Encounters weighted by axial stress': cohorts['mean_contacts_axial_stress_weighted'],
            'Encounters with bulk return': final['mean_encounters'],
            'Largest pressure over energy with bulk return': max(final['rest_pressures_over_energy'])},
        'Three-flavor angular tests': {
            'Separate angular gap floor at 6.67': gap,
            'Shared flavor-preserving angular gap floor at 6.67': shared['elastic'],
            'Shared flavor-changing angular gap floor at 6.67': shared['convert'],
            'Shared flavor-preserving single-flavor largest pressure at 7.05': max(regular['rest_pressures_over_energy']),
            'Shared flavor-changing single-flavor largest pressure at 7.05': max(converted['rest_pressures_over_energy']),
            'Unrestricted flavor-changing combined largest pressure at 7.05': max(free['global_pressures_over_energy'])},
        'Eight additional encounters in the frozen 6.67 population': {
            'Resonant common gas histogram gap': frozen[False]['pooled_common_thermal_tv'],
            'Flat common gas histogram gap': frozen[True]['pooled_common_thermal_tv'],
            'Resonant contacts in highest energy tenth': frozen[False]['contacts_per_initial_top_energy_decile'],
            'Flat contacts in highest energy tenth': frozen[True]['contacts_per_initial_top_energy_decile']}}



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true', help='Regenerate all quantitative results')
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--jobs', type=int, default=2)
    args = parser.parse_args()
    if args.jobs < 1 or (args.output_dir and not args.verify):
        parser.error('--jobs must be positive; --output-dir requires --verify')
    authenticate()
    if args.verify:
        if args.output_dir:
            replay(args.output_dir, args.jobs)
        else:
            with tempfile.TemporaryDirectory(prefix='axial-current-results-') as folder:
                replay(Path(folder), args.jobs)


if __name__ == '__main__':
    main()
