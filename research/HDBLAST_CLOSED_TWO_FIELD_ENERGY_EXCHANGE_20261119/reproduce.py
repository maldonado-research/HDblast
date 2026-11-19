"""Check the archive and reproduce both the field evolutions and the energy balances."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
BASE = 'e00701af7b3d698c2a25eca6dcaa90ed3e84ddca'
PINS_SHA = '557ec21a42e2e325a1b02a69e46ef7fab6e3856d7f2e2ee049248b5b4dc5e61a'
PRODUCTS = {'coupled': 'coupled/COUPLED.json', 'work': 'work/WORK.json',
            'homogeneous': 'homogeneous/results.json', 'phase': 'checks/PHASE_WORK.json',
            'bounds': 'bounds/BOUNDS.json', 'return': 'coupled/RETURN.json'}
SCRIPTS = {'coupled': 'coupled/coupled.py', 'work': 'work/work.py',
           'homogeneous': 'homogeneous/exchange.py', 'phase': 'checks/phase_work.py',
           'bounds': 'bounds/envelope.py', 'return': 'coupled/return_pass.py'}


def require(value, message):
    if not value:
        raise ValueError(message)


def read(path):
    return json.loads(path.read_bytes())


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def verify_files():
    manifest = HERE/'SHA256SUMS'
    seen = set()
    for line in manifest.read_text().splitlines():
        checksum, name = line.split('  ', 1)
        path = HERE/name
        require(not Path(name).is_absolute() and Path(name).as_posix() == name
                and path == path.resolve() and path.is_relative_to(HERE)
                and path != manifest and path not in seen, 'Invalid inventory path')
        require(digest(path) == checksum, f'Hash differs for {name}')
        seen.add(path)
    actual = {p for p in HERE.rglob('*') if p.is_file() and p != manifest
              and not (p.parent.name == '__pycache__' and p.suffix == '.pyc')}
    require(actual == seen, 'Archive coverage differs')
    pin_file = HERE/'SOURCES.json'
    require(digest(pin_file) == PINS_SHA, 'Unexpected source manifest')
    pins = read(pin_file)
    require(pins['base'] == BASE and pins['previous_doi'] == '10.5281/zenodo.24379598',
            'Historical revision differs')
    prior = set()
    for record in pins['sources']:
        path = HERE/record['bundled']
        require(path in seen and path not in prior and path.stat().st_size == record['bytes']
                and digest(path) == record['sha256'], 'Historical input mismatch')
        prior.add(path)
    require(prior == {p for p in seen if p.is_relative_to(HERE/'prior')}, 'Historical coverage differs')
    print(f'Verified {len(seen)} files and {len(prior)} historical objects.', flush=True)


def compare(actual, expected, path=()):
    label = '/'.join(path)
    if isinstance(expected, dict):
        require(isinstance(actual,dict) and actual.keys()==expected.keys(), f'Different keys at {label}')
        for key in expected:
            compare(actual[key],expected[key],(*path,key))
    elif isinstance(expected, list):
        require(isinstance(actual,list) and len(actual)==len(expected), f'Different size at {label}')
        for i, (a,e) in enumerate(zip(actual,expected)):
            compare(a,e,(*path,str(i)))
    elif isinstance(expected, float):
        floor = 8e-11
        # Absolute energies for the weak classical daughter can be much smaller than unity.
        # Protect significant transfers and the directly computed work in those units.
        if path[0] in ('coupled','return') and path[-1] in {'transfer_into_daughter','daughter_pressure_work',
                                              'E_daughter','daughter_energy_change','E_daughter_kinetic',
                                              'E_daughter_gradient','initial_energy'}:
            floor = 5e-14
        require(type(actual) in (int,float) and math.isfinite(actual) and math.isfinite(expected)
                and math.isclose(actual,expected,rel_tol=2e-7,abs_tol=floor),
                f'{label}: expected {expected!r}, obtained {actual!r}')
    elif path and path[-1]=='ode_evaluations':
        require(type(actual) is int and actual>0 and abs(actual-expected)<=max(15,.05*expected),
                'Adaptive solver evaluation count differs')
    else:
        require(type(actual) is type(expected) and actual==expected, f'Scalar differs at {label}')


def physics(products):
    coupled = products['coupled']
    require(coupled['passed'] and products['phase']['checks']['passed'], 'Internal validations failed')
    old = read(HERE/'prior/passage/PASSAGE.json')['runs']
    for seed in ('0.001','-0.001'):
        unseeded = coupled['runs'][seed+'/velocity/0/matter']
        require(all(r['E_daughter']==0. for r in unseeded['snapshots'].values()),'Empty field moved')
        crossing = unseeded['snapshots']['crossing']['t']
        require(abs(crossing-old[seed]['crossing_x']) < 2e-8,'Historical crossing mismatch')
        bound = products['bounds'][seed]['fine']
        cap = bound['abs_coupling_work_over_E_begin_upper']
        require(abs(crossing-(bound['begin']+bound['end'])/2)<2e-8,'Different physical window')
        phase = products['phase']['runs'][seed]
        for form in phase['full_initial_bands'].values():
            low,high=(form['bracket_start_energy_ratio_'+end] for end in ('minimum','maximum'))
            require(-cap < low < 0 < high < cap, 'Optimized state outside continuum bound')
        for form in phase['fixed_n24_support_high_band_check']['propagated_bands'].values():
            lower = form['bracket_start_energy_ratio_minimum']
            upper = form['bracket_start_energy_ratio_maximum']
            witness = form['unchanged_n24_bracket_maximizing_witness']
            ratio = witness['coupling_work']/witness['bracket_start_energy']
            require(-cap < lower < 0 < ratio <= upper+1e-13 and upper < cap,
                    'High-band witnesses violate the energy inequality')
            require(math.isclose(ratio,witness['work_over_bracket_start_energy'],abs_tol=1e-12)
                    and ratio > .97*cap, 'Fixed initial wave differs from its predicted work')
        for mode in products['work']['runs'][seed]['kinetic']['matter']['modes'].values():
            src = mode['intervals']['bracket']
            energy0 = mode['samples']['before']['energy']
            energy1 = mode['samples']['after']['energy']
            require(abs(src['coupling_work']) <= cap*energy0, 'Excess-state work outside bound')
            require(bound['E_end_over_E_begin_lower']<=energy1/energy0<=bound['E_end_over_E_begin_upper'],
                    'Positive excess-state energy outside bound')
            require(abs(energy1-energy0-src['coupling_work']-src['metric_work'])<1e-10,
                    'Excess-state work balance failed')
    for key, run in coupled['runs'].items():
        for row in run['balances'].values():
            ws = row['transfer_into_daughter']
            require(abs(row['daughter_energy_change']+row['daughter_pressure_work']-ws)<1e-11
                    and abs(row['parent_energy_change']+row['parent_pressure_work']+ws)<1e-9,
                    f'Separate energy balance failed for {key}')
            if not run['matter_weight']:
                require(ws==0.,'Direct energy exchange with constant field metric')
    continuation = products['return']
    require(continuation['passed'], 'Return trajectory controls failed')
    events = read(HERE/'prior/rollout/ROLLOUT.json')['branches']['phi_m=-0.5']['runs']['0.001']['fine']['events']
    control = continuation['runs']['velocity/0/matter']['snapshots']
    for now, then in [('first_crossing','first_minimum'),('turning_point','first_turn_after_entry'),
                      ('return_crossing','return_through_minimum')]:
        require(abs(control[now]['t']-events[then]['x']) < 2e-8,'Historical polar event differs')
    for key, trajectory in continuation['runs'].items():
        primary = coupled['runs']['0.001/'+key]
        first = trajectory['balances']['launch_to_first']
        end = trajectory['balances']['launch_to_return']
        lo = primary['balances']['to_crossing']
        require(abs(first['transfer_into_daughter']-lo['transfer_into_daughter']) < 1e-12,
                'First crossing disagrees with independent return run')
        for field in ('transfer_into_daughter','parent_pressure_work','daughter_pressure_work'):
            partial = sum(trajectory['balances'][interval][field] for interval in
                          ('launch_to_first','first_to_turn','turn_to_return'))
            require(abs(partial-end[field])<2e-9, 'Return integral differs from segment sum')
        if key.startswith('displacement/') and key.endswith('/matter'):
            require(first['transfer_into_daughter'] < 0 < end['transfer_into_daughter']
                    and end['daughter_energy_change'] < 0, 'Cumulative displacement transfer failed')
    for branch in products['homogeneous']['branches'].values():
        for run in branch['runs'].values():
            ini = run['initial']
            for sample in ('probe','first_minimum_crossing','at_stop'):
                row = run[sample]
                if row is not None:
                    left = row['a']**6*row['Z']*row['r']
                    right = ini['a']**6*ini['Z']*ini['r']
                    require(math.isclose(left,right,rel_tol=2e-8,abs_tol=2e-12),
                            'Homogeneous conserved charge failed')


def regenerate(destination, jobs):
    destination = destination.resolve()
    require(not destination.is_relative_to(HERE),'Use a separate output directory')
    env = dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    def run(name):
        output = destination/PRODUCTS[name]
        output.parent.mkdir(parents=True,exist_ok=True)
        require(not output.exists(),f'Output already exists: {output}')
        command = [sys.executable,'-I','-B',str(HERE/SCRIPTS[name]),'--output',str(output)]
        print('Recomputing '+SCRIPTS[name],flush=True)
        proc = subprocess.run(command,cwd=HERE,env=env,stdout=subprocess.PIPE,
                              stderr=subprocess.STDOUT,text=True)
        require(proc.returncode==0, f"{name} exited {proc.returncode}\n{proc.stdout[-4000:]}")
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        list(pool.map(run,PRODUCTS))
    data = {name:read(destination/path) for name,path in PRODUCTS.items()}
    for name,path in PRODUCTS.items():
        compare(data[name],read(HERE/path),(name,))
    physics(data)
    print('Verified all six numerical products, independent balances and continuum bounds.',flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true',help='Repeat the numerical experiments')
    parser.add_argument('--jobs',type=int,default=2)
    parser.add_argument('--output-dir',type=Path,help='Keep regenerated files in this separate directory')
    args=parser.parse_args()
    require(args.jobs>=1 and (args.output_dir is None or args.verify),'Incompatible replay arguments')
    verify_files()
    if args.verify:
        if args.output_dir:
            regenerate(args.output_dir,args.jobs)
        else:
            with TemporaryDirectory(prefix='hdblast-energy-') as tmp:
                regenerate(Path(tmp),args.jobs)
