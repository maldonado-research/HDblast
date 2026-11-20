"""Run the local field-triggered mobility test on the coupled charged-wall solver."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'coupled'/'charged_forward.cpp'
SHA='3bd89d983f8243a3371099e132032ce19fe6975e7fde8b7f9a7581632bbce027'
default=dict(cells=24,late_cells=6,cfl=.08,depth=.002,dose=.123,end=16.,D=.0001,Q=1.e8,
             nu=6.,epsilon=.01,ramp=.5,scalar_floor=.38)
cases={
    'primary':{'cells':48},'coarse':{},'fine':{'cells':64},
    'time_space':{'late_cells':10,'cfl':.05},'faster':{'ramp':.25},
    'slower':{'ramp':1.},'less_selective':{'epsilon':.05},'broad':{'epsilon':.1},
    'strict_gate':{'epsilon':.0001},'no_trigger':{'epsilon':0.},'no_pulse':{'dose':0.},
    'early_window':{'scalar_floor':.004,'end':.12},
    'early_refined':{'scalar_floor':.004,'end':.12,'cells':48,'late_cells':10,'cfl':.05},
    'wrong_scalar':{'scalar_floor':1.5},'collisionless':{'nu':0.},'no_force':{'Q':0.}
}

def build(directory):
    if hashlib.sha256(SOURCE.read_bytes()).hexdigest()!=SHA:
        raise ValueError('Unexpected coupled solver revision')
    source=directory/'charged_forward.cpp'
    source.write_bytes(SOURCE.read_bytes())
    subprocess.run(['patch','-s','-p1','-i',str(HERE/'mobility.diff')],cwd=directory,check=True)
    executable=directory/'gate'
    subprocess.run(['g++','-O3','-std=c++17',str(source),'-o',str(executable)],check=True)
    return executable

def run(exe,params,waveform):
    norm,stop=waveform
    p=default|params
    args=[p['cells'],p['depth'],p['dose'],p['end'],norm,stop,p['cfl'],p['late_cells'],
          .08,13200,0,p['D'],p['Q'],p['nu'],p['epsilon'],p['ramp'],p['scalar_floor']]
    proc=subprocess.run([str(exe),*map(str,args)],text=True,capture_output=True,check=True)
    row={key:float(value) for key,value in (line.split() for line in proc.stdout.splitlines())}
    if not all(math.isfinite(v) for v in row.values()):
        raise ValueError('Nonfinite model results')
    for key,bound in [('friedmann_error',2e-6),('surface_balance_error',2e-6),
                      ('sheet_balance_error',2e-6),('gas_balance_error',2e-6)]:
        if row[key]>bound:raise ValueError((key,row[key],bound))
    if 'activation_radius' in row:
        if row['activation_charge_work']!=0. or row['activation_gas_temperature']!=0.:
            raise ValueError('Energy appeared before contact')
        if row['activation_psi']<p['scalar_floor'] or row['activation_magnetic_share']>=p['epsilon']:
            raise ValueError('Incorrect threshold')
    elif row['end_charge_work']!=0. or row['end_temperature']!=0.:
        raise ValueError('Nonzero work without contact')
    return row

def waveform():
    spec=importlib.util.spec_from_file_location('charged_waveform',SOURCE.with_name('charged.py'))
    source=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(source)
    return source.waveform_constants()

def calculate(jobs):
    pulse=waveform()
    with tempfile.TemporaryDirectory(prefix='hdblast-mobility-') as work:
        exe=build(Path(work))
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            rows=dict(zip(cases,pool.map(lambda k:run(exe,cases[k],pulse),cases),strict=True))
    for key in ('wrong_scalar','no_pulse','no_trigger'):
        if 'activation_tau' in rows[key]:raise ValueError('Unexpected trigger in '+key)
    if rows['collisionless']['collisionless_impulse_error']>2e-6 or rows['collisionless']['end_temperature']!=0.:
        raise ValueError('Collisionless gas failed impulse or adiabatic check')
    return dict(parent_solver_sha256=SHA,default_parameters=default,varied_parameters=cases,runs=rows)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--jobs',type=int,default=2)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.jobs),sort_keys=True,indent=2,allow_nan=False)+'\n')
