"""Authenticate the scattering study and reconstruct its numerical results."""
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
PARENT = ROOT.parent / 'HDBLAST_POLARIZATION_20261121'
PARENT_DIGEST = '38c06ac62f996fbfbbd2e106e43c67de2896607500787630c617ef610e42916b'
PRODUCTS = (
    'bound/ELASTIC.json',
    'bound/CONTINUUM.json',
    'bound/ALL_FRAMES.json',
    'kinetic/KINETIC.json',
    'rate/RATE.json',
    'rate/ENCLOSURE.json',
    'survival/SURVIVAL.json',
    'live/LIVE.json',
)


def digest(path):
    with path.open('rb') as source:
        return hashlib.file_digest(source, 'sha256').hexdigest()


def authenticate(directory):
    manifest = directory / 'SHA256SUMS'
    available = {p.relative_to(directory).as_posix() for p in directory.rglob('*')
                 if p.is_file() and p != manifest and '__pycache__' not in p.parts}
    registered = set()
    for line in manifest.read_text().splitlines():
        checksum, name = line.split('  ', 1)
        if name not in available or name in registered:
            raise ValueError(f'Invalid manifest entry {directory.name}/{name}')
        file = directory / name
        if file.is_symlink() or digest(file) != checksum:
            raise ValueError(f'Input differs from manifest {directory.name}/{name}')
        registered.add(name)
    if registered != available:
        raise ValueError(f'Unregistered files in {directory.name}: {sorted(available-registered)}')
    return len(registered)


def run(script, *arguments):
    print(f'Replaying {script.parent.name}/{script.name}', flush=True)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='1')
    subprocess.run([sys.executable, '-I', '-B', str(script), *map(str, arguments)],
                   check=True, env=env)


def replay(dest, jobs):
    dest = dest.resolve()
    if any(dest.is_relative_to(study) for study in (ROOT, PARENT)):
        raise ValueError('Write replay outputs outside both input studies')
    if dest.exists() and (not dest.is_dir() or any(dest.iterdir())):
        raise ValueError('Choose a new or empty output directory')
    dest.mkdir(parents=True, exist_ok=True)
    for component in {Path(path).parent for path in PRODUCTS}:
        (dest / component).mkdir()
    compare = runpy.run_path(str(PARENT / 'parent/study.py'))['compare']
    with tempfile.TemporaryDirectory(prefix='hdblast-scattering-inputs-') as work:
        temp = Path(work)
        fields = temp / 'fields'
        particles = temp / 'particles'
        cases = ['dedicated_native', 'dedicated_half', 'shared_native', 'shared_half']
        run(PARENT/'fields/launch.py', '--parent', PARENT/'parent', '--history',
            PARENT/'parent/parent/history', '--cases', *cases, '--jobs', jobs, '--output-dir', fields)
        field_reference = json.loads((PARENT/'fields/FIELDS.json').read_bytes())
        field_results = json.loads((fields/'FIELDS.json').read_bytes())
        compare(field_results, {key: field_reference[key] for key in cases}, 'replayed pulse fields')
        traces = fields / 'traces'
        run(ROOT/'kinetic/particles.py', '--study', PARENT, '--fields', traces, '--out', particles,
            '--jobs', jobs, '--powers', 12, 13, 14)
        run(ROOT/'bound/elastic.py', '--snapshots', particles, '--output', dest/'bound/ELASTIC.json')
        run(ROOT/'bound/continuum.py', '--snapshots', particles, '--output', dest/'bound/CONTINUUM.json')
        run(ROOT/'bound/all_frames.py', '--snapshots', particles, '--output', dest/'bound/ALL_FRAMES.json')
        run(ROOT/'rate/exposure.py', '--parent', PARENT, '--fields', traces,
            '--output', dest/'rate/RATE.json', '--jobs', jobs)
        run(ROOT/'rate/enclosure.py', '--parent', PARENT, '--fields', traces,
            '--output', dest/'rate/ENCLOSURE.json', '--jobs', jobs)
        run(ROOT/'kinetic/scatter.py', '--particles', particles, '--output', dest/'kinetic/KINETIC.json',
            '--jobs', jobs)
        run(ROOT/'live/evolve.py', '--study', PARENT, '--fields', traces, '--output', dest/'live/LIVE.json',
            '--jobs', jobs)
        run(ROOT/'survival/cohorts.py', '--study', PARENT, '--particles', particles,
            '--output', dest/'survival/SURVIVAL.json')
    for name in PRODUCTS:
        rebuilt = json.loads((dest/name).read_bytes())
        original = json.loads((ROOT/name).read_bytes())
        compare(rebuilt, original, name)
        print('Reproduced', name, flush=True)
    for label, values in benchmarks(dest).items():
        payload = json.dumps(values, sort_keys=True, allow_nan=False)
        print(label, payload, flush=True)
        if os.environ.get('GITHUB_ACTIONS') == 'true':
            print(f'::notice file=study.py,line=1,title=HDblast {label}::{payload}', flush=True)


def benchmarks(output):
    def read(path):
        return json.loads((output/path).read_bytes())

    survival = read('survival/SURVIVAL.json')['runs']
    encounter = read('rate/RATE.json')['runs']
    turns = read('bound/CONTINUUM.json')['runs']
    kinetic = read('kinetic/KINETIC.json')['runs']
    live_result = read('live/LIVE.json')
    live = live_result['runs']
    result = {}
    for history, label in [('dedicated_half','Separate interaction'), ('shared_half','Shared interaction')]:
        regional = next(r for r in survival if r['history'] == history and r['case'] == 'small_pulse'
                        and r['birth_pairs'] == 16384)
        pair = next(r for r in encounter if r['history'] == history and r['case'] == 'early_late_4pct'
                    and r['pairs'] == 32768 and r['steps'] == 384)
        elastic = next(r for r in turns if r['snapshot'] == history+'_small_pulse_14.npz')
        box = next(r for r in kinetic if r['history'] == history and r['case'] == 'small_pulse'
                   and r['channels'] == 'mixed' and r['cone_degrees'] == 180 and r['particles_per_sign'] == 8192)
        late = next(s for s in box['samples'] if s['collisions_per_particle'] == 32)
        result[label] = {
            'Protected part of the 20-cell discrepancy at the specified rate cap, percentage points':
                100 * regional['tv_bound'],
            'Encounter count per particle on the unscattered reference paths':
                pair['free_pair_encounters_per_particle'],
            'Charge-resolved discrepancy floor across the tested continuum of axial scattering frames, percentage points':
                100 * elastic['minimum_charge_resolved_tv'],
            'Count discrepancy from the archived reference after 32 mixed-charge, full-angle encounters per particle in the box, percentage points':
                100 * late['archived_tv'],
        }
    for case, label in [('small_pulse','Small added pulse'), ('pressure_fit','Pressure-fitted pulse')]:
        rows = [r for r in live if r['source'] == case and r['birth_pairs'] == 16384 and r['steps'] == 512]
        off = next(r for r in rows if r['law'] == 'off')
        active = [r for r in rows if r['law'] == 'vector']
        if len(active) != 3:
            raise ValueError('Missing independent collision seeds')
        later = [r for r in live_result['other_work_markers'] if r['source'] == case and r['work_fraction'] == .05
                 and r['birth_pairs'] == 16384 and r['steps'] == 1024 and r['law'] == 'vector']
        if len(later) != 3:
            raise ValueError('Missing later-marker collision tests')
        result[label] = {
            'At 5% work, particles with at least one encounter, percent': 100 * sum(r['fraction_ever_hit'] for r in later)/len(later),
            'At 5% work, count discrepancy with physical collisions, percentage points': 100 * sum(r['tv20'] for r in later)/len(later),
            'At 1% work, particles with at least one encounter, percent': 100 * sum(r['fraction_ever_hit'] for r in active)/len(active),
            'Count discrepancy with collisions off, percentage points': 100 * off['tv20'],
            'Count discrepancy with the vector collision channel, percentage points': 100 * sum(r['tv20'] for r in active)/len(active),
            'Weak sideways pressure with collisions off, in units of the thermal pressure': off['sideways_pressure_over_target'][0],
            'Weak sideways pressure with vector collisions, in units of the thermal pressure': sum(r['sideways_pressure_over_target'][0] for r in active)/len(active),
        }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true', help='Regenerate all numerical products')
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--jobs', type=int, default=2)
    args = parser.parse_args()
    if args.jobs < 1 or (args.output_dir and not args.verify):
        parser.error('--jobs must be positive; --output-dir requires --verify')
    if digest(PARENT/'SHA256SUMS') != PARENT_DIGEST:
        raise ValueError('Predecessor differs from its tested public revision')
    print(f'Authenticated {authenticate(PARENT)} predecessor files and {authenticate(ROOT)} study files.',
          flush=True)
    if args.verify:
        if args.output_dir:
            replay(args.output_dir, args.jobs)
        else:
            with tempfile.TemporaryDirectory(prefix='scattering-results-') as temp:
                replay(Path(temp), args.jobs)


if __name__ == '__main__':
    main()
