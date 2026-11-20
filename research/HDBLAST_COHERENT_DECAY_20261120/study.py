"""Authenticate or reproduce the coupled neutral-parent conversion and heating study."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
PRODUCTS = ('coupled/JOINED', 'audit/EXPOSURE', 'quantum/HISTORY', 'quantum/NARROW_HISTORY',
            'quantum/FERMIONS', 'quantum/NARROW', 'kinetics/KINETICS', 'envelope/ENVELOPE', 'kinetics/MANY', 'kinetics/SWEEP')
PREVIOUS = {
    'PARENT.json': '17673e63c8e6a7393326f0f623738f4aa532263df0765839e7f177a196670acf',
    'production.py': '11f4db74efd44cbc0d8f30a2de8b14c928a96fc53e6391331c0d83a0c10c129d',
    'production_forward.cpp': 'dd305f17010ce54a59799db07e51060e0d299f880a6226ed4ad1f51e0f4fe600',
    'narrow_cases.json': '8953821a4793e0c04ef877530c9124940853cdfad6be1737de1bbf2d3bb182c4',
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def authenticate():
    manifest = ROOT/'SHA256SUMS'
    files = {p.relative_to(ROOT).as_posix(): p for p in ROOT.rglob('*')
             if p.is_file() and p != manifest and '__pycache__' not in p.parts}
    seen = set()
    for line in manifest.read_text().splitlines():
        digest, name = line.split('  ', 1)
        require(name in files and name not in seen, 'Invalid manifest entry '+name)
        require(not files[name].is_symlink(), 'Unexpected symbolic link '+name)
        with files[name].open('rb') as data:
            require(hashlib.file_digest(data, 'sha256').hexdigest() == digest, 'Checksum mismatch '+name)
        seen.add(name)
    require(seen == files.keys(), 'Unexpected files '+str(set(files)-seen))
    for name, digest in PREVIOUS.items():
        historical = ROOT/'prior'/name
        if name == 'PARENT.json':
            stream = gzip.open(historical.with_suffix('.json.gz'), 'rb')
        else:
            stream = historical.open('rb')
        with stream:
            require(hashlib.file_digest(stream, 'sha256').hexdigest() == digest, 'Previous study mismatch '+name)
    print(f'Authenticated {len(files)} files including four inputs from the preceding study.', flush=True)


def compare(actual, reference, path):
    if isinstance(reference, dict):
        require(isinstance(actual, dict) and actual.keys() == reference.keys(), 'Changed keys '+path)
        for key, value in reference.items():
            compare(actual[key], value, path+'/'+key)
    elif isinstance(reference, list):
        require(isinstance(actual, list) and len(actual) == len(reference), 'Changed length '+path)
        for i, (a, b) in enumerate(zip(actual, reference)):
            compare(a, b, path+'/'+str(i))
    elif type(reference) is float or (path.startswith('kinetics/KINETICS/') and
                                      path.endswith(('/first_ell','/last_ell','/count_all_six'))) or (
                                      path.startswith('kinetics/MANY/') and path.endswith('/label_upper')):
        # Mass and surface balance diagnostics subtract nearly equal terms. Their
        # physical ceilings and the conservation ledgers are checked independently.
        key = path.rsplit('/', 1)[-1]
        absolute = {'tube_mass_error': 1.e-6, 'sheet_balance_error': 4.e-9,
                    'relative_history_error': 1.5e-8}.get(key, 2.e-9)
        relative = 1.e-6
        if path.startswith('kinetics/KINETICS/') and '/delivery_windows/' in path:
            # Rounding a moving shell edge to the adjacent angular level changes
            # its degeneracy. The exact count at each run's own edges is checked below.
            if key in {'count_all_six', 'no_recycling_action_ceiling',
                       'clearance_limited_action_ceiling_at_fluid_rate',
                       'minimum_clearance_rate_at_pointwise_drag_floor',
                       'required_over_available_at_pointwise_drag_floor'}:
                relative = 3.e-4
            if key == 'max_step_change':
                absolute = 5.e-6
        require(type(actual) is type(reference) and math.isfinite(actual) and
                math.isclose(actual, reference, rel_tol=relative, abs_tol=absolute),
                f'{path}: {actual} versus {reference}')
    else:
        require(type(actual) is type(reference) and actual == reference, 'Changed value '+path)


def check(products):
    r = products['coupled/JOINED']['runs']
    p = r['primary']
    require(p['hot_intervals'] == 1 and 2.8 < math.exp(p['hot_efolds']) < 3.
            and r['finest']['hot_efolds'] > math.log(2.) and p['end_H'] < 0,
            'Primary sustained growth, refinement or contraction changed')
    require(p['gas_balance_error'] < 1.e-8 and p['sheet_balance_error'] < 1.e-7
            and p['tube_mass_error'] < 2.e-5, 'Primary conservation check failed')
    exposure = products['audit/EXPOSURE']
    hot = exposure['primary']['first_hot']
    require(2.4e6 < hot['min_Q_alpha_to_reproduce_actual_electric_work'] < 2.42e6
            and 25 < hot['integrated_excess_mass_phase'] < 25.2
            and exposure['primary']['electric_at_exit']['abs_eE_over_m2'] > 1.e5
            and exposure['historical']['first_hot']['integrated_excess_mass_phase'] < .19
            and 64.8 < exposure['primary']['trace_audit']['fixed_mass_three_flavor_test']['lower_over_upper'] < 65.,
            'Electric-work floor, strong-field scale or finite-phase control changed')
    q = products['quantum/FERMIONS']
    require(abs(q['static_control']['ratio']-1) < .001 and
            5.4 < q['cases']['stable']['windows']['through_turnaround']['ratio_to_velocity_damping'] < 5.7 and
            1.03 < q['cases']['heavy_stable']['windows']['through_turnaround']['ratio_to_velocity_damping'] < 1.05,
            'Static or delayed fermion production check failed')
    kinetics = products['kinetics/KINETICS']
    primary = kinetics['runs']['primary']
    for run in kinetics['runs'].values():
        for window in run['delivery_windows']:
            lo, hi = window['first_ell'], window['last_ell']
            require(0 < lo <= hi and window['count_all_six'] == 3*(hi*(hi+1)-(lo-1)*lo),
                    'The integer state count disagrees with its measured shell edges')
    delivery = next(w for w in primary['delivery_windows'] if w['fraction_bin'] == [.75, .9])
    require(kinetics['magnetic_spectrum_count_checks'] >= 25 and
            7.2e23 < primary['pointwise_swave_action_floor'] < 7.6e23 and
            2.4e23 < primary['hot_integrated_swave_action_floor'] < 2.7e23 and
            2.e19 < delivery['no_recycling_action_ceiling'] < 2.3e19 and
            delivery['minimum_clearance_rate_at_pointwise_drag_floor'] > 4.e4,
            'Compact-state inventory or elastic-scattering action bound changed')
    envelope = products['envelope/ENVELOPE']
    high = products['kinetics/MANY']['runs']['mass_4e10']
    highest_required_share = max(w['requested_over_capacity_at_screen_witness'] for w in high['windows'])
    require(abs(envelope['averaged_minus_resolved_relative']['hot_efolds']) < .01 and
            2.20 < high['hot_radius_factor'] < 2.30 and high['original_screen']['witness_action'] > 0
            and .105 < highest_required_share < .12,
            'High-mass averaged continuation or birth-state capacity changed')
    for label in ('mass_4e10', 'mass_4e12'):
        for mesh in ('', '_fine', '_half_step'):
            run = envelope['runs'][label + mesh]
            require(2.2 < math.exp(run['hot_efolds']) < 2.35 and run['stop_turnaround'] == 1,
                    'The heavy-mass growth interval or turnaround changed for '+label+mesh)
    inner = envelope['inner_source']
    require(abs(inner['half_resolution_relative']) < 1.e-8 and
            max(map(abs,inner['match_energy_relative'].values())) < 5.e-6,
            'Scaled drive energy disagrees with the full evolutions or its own refinement')
    sweep = products['kinetics/SWEEP']['runs']['mass_4e10']
    for row in sweep['runs'].values():
        require(row['floor_sweep_peak'] > 10. and row['main_sweep_peak'] > 200. and
                8.e-4 < row['main_sweep_deficit_fraction'] < 1.3e-3,
                'The 4e10 sharp-shell rate deficit or its integrated effect changed')
    high_sweep = products['kinetics/SWEEP']['runs']['mass_4e12']
    require(math.isclose(high_sweep['action'], 5*high_sweep['reference_drag_floor']),
            'The high-mass test action changed')
    for row in high_sweep['runs'].values():
        require(.82 < row['extended_main_sweep_peak'] < .83 and
                abs(row['extended_main_sweep_fraction_tested']-.99998) < 1e-8 and
                row['extended_main_sweep_deficit_fraction'] == 0. and
                row['extended_main_sweep_zero_crossings'] == 5.,
                'The high-mass sharp-shell rate check or delivery interval changed')
    narrow = products['quantum/NARROW']['matched_400']['windows']
    require(.95 < narrow['longer']['ratio_to_velocity_damping'] < 1.05 and
            narrow['crossing']['ratio_to_velocity_damping'] > 1., 'Matched finite-time calculation changed')


def replay(out):
    out = out.resolve()
    require(not out.is_relative_to(ROOT) and (not out.exists() or (out.is_dir() and not any(out.iterdir()))),
            'Choose a new or empty output directory outside the study.')
    out.mkdir(parents=True, exist_ok=True)
    targets = {n: out/(n+'.json') for n in PRODUCTS}
    for path in targets.values():
        path.parent.mkdir(exist_ok=True)
    historical = ROOT/'prior/PARENT.json.gz'
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')

    def run(script, *args):
        print('Replaying '+script, flush=True)
        subprocess.run([sys.executable, '-I', '-B', str(ROOT/script), *map(str,args)], check=True, env=env)

    def quantum():
        run('quantum/extract.py', historical, '--output', targets['quantum/HISTORY'])
        run('quantum/fermions.py','--data',targets['quantum/HISTORY'],'--output',targets['quantum/FERMIONS'])
        run('quantum/narrow_trace.py', ROOT/'prior', '--output', targets['quantum/NARROW_HISTORY'])
        run('quantum/narrow.py','--data',targets['quantum/NARROW_HISTORY'],'--output',targets['quantum/NARROW'])

    joined = targets['coupled/JOINED']
    with ThreadPoolExecutor(3) as pool:
        oscillators = pool.submit(quantum)
        run('coupled/joined.py','--output',joined,'--jobs','2')
        envelope = pool.submit(run,'envelope/envelope.py','--joined',joined,
                               '--output',targets['envelope/ENVELOPE'],'--jobs','2')
        run('coupled/verify.py',joined)
        run('audit/exposure.py','--source',joined,'--output',targets['audit/EXPOSURE'])
        run('kinetics/kinetics.py','--study',joined,'--output',targets['kinetics/KINETICS'])
        oscillators.result()
        envelope.result()
        many = pool.submit(run, 'kinetics/many.py', '--envelope',targets['envelope/ENVELOPE'],
                           '--output',targets['kinetics/MANY'])
        run('kinetics/sweep.py','--envelope',targets['envelope/ENVELOPE'],'--output',targets['kinetics/SWEEP'])
        many.result()
    products={n:json.loads(f.read_bytes()) for n,f in targets.items()}
    check(products)
    for name, value in products.items():
        compare(value, json.loads((ROOT/(name+'.json')).read_bytes()), name)
    print(f'All {len(products)} numerical products and physical checks passed.',flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    if args.output_dir and not args.verify:
        parser.error('--output-dir requires --verify')
    authenticate()
    if args.verify:
        if args.output_dir:
            replay(args.output_dir)
        else:
            with tempfile.TemporaryDirectory(prefix='coherent-decay-replay-') as tmp:
                replay(Path(tmp))
