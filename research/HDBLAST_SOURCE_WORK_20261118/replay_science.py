"""Verify the published source-work diagnostics. Use --full to recompute the quantum fields."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
PREVIOUS = HERE.parent / 'HDBLAST_PULSE_PARTICLE_RESPONSE_20261118'
ARMS = {'quiet':'quiet', 'near_plus':'near_positive', 'near_minus':'near_negative', 'far_plus':'far_positive'}


def read(path):
    return json.loads(path.read_bytes())


def authenticate():
    sources = read(HERE / 'SOURCES.json')
    for entry in sources['sources']:
        raw = (HERE.parents[1] / entry['path']).read_bytes()
        if len(raw) != entry['bytes'] or hashlib.sha256(raw).hexdigest() != entry['sha256']:
            raise ValueError('Historical file differs: ' + entry['path'])
    manifest = HERE/'SHA256SUMS'
    if manifest.exists():
        for line in manifest.read_text().splitlines():
            digest, name = line.split(maxsplit=1)
            if hashlib.sha256((HERE/name.removeprefix('*')).read_bytes()).hexdigest() != digest:
                raise ValueError('Capsule file differs: ' + name)
    print(f"Verified {len(sources['sources'])} historical inputs at {sources['base'][:12]}", flush=True)


def run(script, *arguments, output=None):
    opts = ['-O'] if sys.flags.optimize else []
    command = [sys.executable, '-I', '-B', *opts, script, *arguments]
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    if output is None:
        subprocess.run(command, check=True, env=env)
    else:
        with output.open('wb') as handle:
            subprocess.run(command, stdout=handle, check=True, env=env)


def close(x, y, label, absolute, relative=0.):
    if not (math.isfinite(x) and math.isfinite(y) and math.isclose(x,y,rel_tol=relative,abs_tol=absolute)):
        raise ValueError(f'{label}: {x!r} differs from {y!r}')


def compare(x, y, label, absolute, relative):
    """Compare the complete structure and numerical outputs at declared tolerances."""
    if isinstance(x, dict) and isinstance(y, dict) and x.keys()==y.keys():
        for key in x:
            compare(x[key], y[key], f'{label}/{key}', absolute, relative)
    elif isinstance(x, list) and isinstance(y, list) and len(x)==len(y):
        for i,(value,expected) in enumerate(zip(x,y)):
            compare(value,expected,f'{label}/{i}',absolute,relative)
    elif type(x) is float and type(y) in (float,int):
        close(x,y,label,absolute,relative)
    elif type(x) is not type(y) or x!=y:
        raise ValueError(f'{label}: changed report field')


def check_pairings(reports):
    kernel = read(reports/'kernel/KERNEL.json')
    phase = read(reports/'subtraction/PHASE_WORK.json')
    reservoir = read(reports/'decay/RESERVOIR.json')
    for grid, filename, setting in [('coarse','COARSE','coarse_32'),('fine','LEDGER','fine_48')]:
        ledger = read(reports/f'ledger/{filename}.json')
        for name, phase_name in ARMS.items():
            raw = ledger['arms'][name]
            pair = phase[setting]['arms'][phase_name]
            predicted = kernel[grid]['quiet_N']
            if name != 'quiet':
                predicted += kernel[grid]['arms'][name]['delta_N']
            close(raw['number_second_order'],predicted,f'{grid}/{name}/late number',1e-6)
            close(raw['number_instantaneous'],pair['original_final_inst_number'],'instantaneous number',2e-6)
            for key in ('mass_work','expansion_work'):
                independent=pair['early'][key]+pair['original'][key]
                close(raw['work']['full'][key],independent,key,5e-5)
            for key in ('ward_residual','vacuum_residual','max_mode_norm_error'):
                close(raw[key],0.,key,5e-6)
            if name != 'quiet':
                delta=ledger['pulse_minus_quiet'][name]
                excess=delta['work']['full']['mass_work']
                close(excess,phase[setting]['pulse_minus_quiet'][phase_name]['mass_work'],'paired scalar work',4e-5)
                matter=reservoir[f'{grid}/y1']
                radiation=100*(matter[name]['2.1']['daughter_physical_density']/matter['quiet']['2.1']['daughter_physical_density']-1)
                expected_sign = -1 if name == 'near_plus' else 1
                if expected_sign * radiation <= 0:
                    raise ValueError(f'{grid}/{name}: unexpected sign of daughter radiation')
    print('Verified shared finite-band particle counts, scalar and pressure work, and daughter radiation contrasts',flush=True)


def quick():
    run(HERE/'kernel/audit_kernel.py')
    run(HERE/'subtraction/verify_phase.py')
    with tempfile.TemporaryDirectory() as temp:
        out=Path(temp)/'check.json'
        run(HERE/'decay/check.py',output=out)
        compare(read(out),read(HERE/'decay/CHECK.json'),'decay signatures',3e-9,2e-6)
    check_pairings(HERE)


def full(temp):
    for name in ('ledger','kernel','subtraction','decay'):
        (temp/name).mkdir()
    for name,args in [('COARSE',['--grid','coarse']),('LEDGER',[]),('REFINED',['--order','96','--steps','8192'])]:
        file=f'ledger/{name}.json'
        run(HERE/'ledger/ledger.py',PREVIOUS,*args,output=temp/file)
        compare(read(temp/file),read(HERE/file),file,3e-6,3e-6)
        print(f'Reproduced {file}',flush=True)

    for f in ('kernel/audit_kernel.py','decay/check.py'):
        shutil.copyfile(HERE/f,temp/f)
    csv_path=temp/'kernel/RESPONSE.csv'
    run(HERE/'kernel/kernel.py',PREVIOUS,'--write-kernel',csv_path,output=temp/'kernel/KERNEL.json')
    run(temp/'kernel/audit_kernel.py')
    actual, stored = [read(p/'kernel/KERNEL.json') for p in (temp,HERE)]
    for grid in stored:
        close(actual[grid]['quiet_N'],stored[grid]['quiet_N'],grid,2e-7)
        for arm, reference in stored[grid]['arms'].items():
            row=actual[grid]['arms'][arm]
            for field in ('delta_percent','linear_prediction_percent'):
                close(row[field],reference[field],f'{grid}/{arm}/{field}',2e-5)
            close(row['linear_to_finite_ratio'],reference['linear_to_finite_ratio'],'relative response',5e-4)
            compare(row['components'],reference['components'],'independent history components',2e-7,1e-5)
        design=stored[grid]['designed_fixed_geometry_controls']
        for key in ('lowest','highest'):
            compare(actual[grid]['designed_fixed_geometry_controls']['controls'][key]['direct_delta_N'],
                    design['controls'][key]['direct_delta_N'],'equal-area control',2e-7,1e-5)
    with csv_path.open() as f,(HERE/'kernel/RESPONSE.csv').open() as g:
        left,right=csv.reader(f),csv.reader(g)
        if next(left)!=next(right): raise ValueError('Kernel table header')
        for a,b in zip(left,right,strict=True):
            for v,w in zip(a,b,strict=True): close(float(v),float(w),'sampled quiet kernel',1e-6,2e-5)
    print('Recomputed the quiet retarded kernel, four rolling responses and designed controls',flush=True)
    out=temp/'subtraction/PHASE_WORK.json'
    run(HERE/'subtraction/phase_work.py',PREVIOUS,output=out)
    run(HERE/'subtraction/verify_phase.py',out)
    for program,file in [('reservoir.py','RESERVOIR.json'),('split.py','SPLIT.json')]:
        path='decay/'+file
        run(HERE/'decay'/program,PREVIOUS,output=temp/path)
        compare(read(temp/path),read(HERE/path),path,3e-5,3e-6)
    run(temp/'decay/check.py',output=temp/'decay/CHECK.json')
    for a,b in zip(read(temp/'decay/CHECK.json'),read(HERE/'decay/CHECK.json'),strict=True):
        if a['case']!=b['case']: raise ValueError('Reservoir case')
        close(a['change_percent'],b['change_percent'],'daughter change, percentage points',3e-5)
    print('Recomputed the damped field and independent split evolution',flush=True)
    check_pairings(temp)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--full',action='store_true',help='recompute the numerical results with NumPy and SciPy')
    options=parser.parse_args()
    authenticate()
    quick()
    if options.full:
        with tempfile.TemporaryDirectory() as directory:
            full(Path(directory))
