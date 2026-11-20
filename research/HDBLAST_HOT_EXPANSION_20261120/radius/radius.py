"""Replay finite-source radius controls and bisect the launch-time enclosed mass crossing."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / 'prior'
SOURCE_HASH = 'e3e20e80d7558895ec189f46e0989c577acccdd81b46828b39aeb463577be611'
spec = importlib.util.spec_from_file_location('prior_decay', PRIOR / 'decay.py')
prior = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prior)
NORM, STOP = prior.waveform_constants()


def run(exe, p):
    argv = (p['cells'], .002, p['dose'], p['end'], NORM, STOP, p['cfl'],
            p['late_cells'], .08, p['late_K'], 0, p['D'], p['Q'], p['nu'],
            p['psi_star'], .03, p['kappa'], .95, 1.05, p['a0'], p['height'], p['calibration'])
    proc = subprocess.run([str(exe), *map(str, argv)], capture_output=True, text=True, check=True)
    row = {key: float(value) for key, value in map(str.split, proc.stdout.splitlines())}
    if any(not math.isfinite(v) for v in row.values()):
        raise ValueError('Nonfinite wall diagnostic')
    row['hot_growth'] = math.exp(row.get('ninety_intrinsic_expanding_efolds', 0.))
    return row


def mass_bracket(exe, defaults, cells):
    params = dict(defaults, cells=cells, cfl=.08, end=.014)
    lo, hi = 2.02, 2.03
    def measure(radius):
        row = run(exe, dict(params, a0=radius))
        return {'a0':radius,'minimum_bulk_mass':row['min_bulk_mass'],
                'end_tau':row['end_tau'], 'excluded':bool(row.get('negative_mass_stop'))}
    lower, upper = measure(lo), measure(hi)
    if lower['excluded'] or not upper['excluded']:
        raise ValueError('Starting interval does not bracket the mass crossing')
    for _ in range(12):
        row = measure(.5*(lower['a0']+upper['a0']))
        if row['excluded']:
            upper = row
        else:
            lower = row
    return {'cells':cells, 'end':params['end'], 'lower':lower, 'upper':upper}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=3)
    parser.add_argument('--cases', type=Path, default=HERE/'cases.json')
    args=parser.parse_args()
    src = PRIOR / 'decay_forward.cpp'
    if hashlib.sha256(src.read_bytes()).hexdigest()!=SOURCE_HASH:
        raise ValueError('Prior charged-wall source differs from the published input')
    config=json.loads(args.cases.read_text())
    cases={k:{**config['defaults'],**overrides} for k,overrides in config['cases'].items()}
    with tempfile.TemporaryDirectory(prefix='hdblast-radius-') as tmp:
        exe=Path(tmp)/'solver'
        cpp=Path(tmp)/'radius_forward.cpp'
        subprocess.run(['patch','--batch','--silent','--output',str(cpp),str(src),str(HERE/'radius_forward.patch')],check=True)
        subprocess.run(['g++','-O3','-std=c++17',str(cpp),'-o',str(exe)],check=True)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            runs=dict(zip(cases,pool.map(lambda p: run(exe,p), cases.values())))
            brackets=list(pool.map(lambda cells: mass_bracket(exe,config['defaults'],cells),[16,24,32]))
    result={'cases':cases,'runs':runs,'launch_mass_brackets':brackets,'prior_cpp_sha256':SOURCE_HASH}
    args.output.write_text(json.dumps(result,allow_nan=False,indent=2,sort_keys=True)+'\n')


if __name__=='__main__':
    main()
