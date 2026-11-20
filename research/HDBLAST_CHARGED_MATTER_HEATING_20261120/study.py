"""Verify the frozen sources and rerun the closed-wall charged matter calculations."""
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
PRODUCTS = {'background':'TRACE','coupled':'CHARGED','kinetic':'KINETIC','theory':'SPECTRUM',
            'quantum':'PAIRS','mass':'MASS','gate':'GATE','synthesis':'BRIDGE'}
PINNED_SOLVER = '94d98015df5f38f98ccccbc3557020841a9dda31fcc3a28a328b9e2fdb8badd3'


def require(ok,message):
    if not ok:
        raise ValueError(message)


def checksum(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle,'sha256').hexdigest()


def verify_files():
    entries = {}
    for line in (ROOT/'SHA256SUMS').read_text().splitlines():
        digest, name = line.split('  ',1)
        path = Path(name)
        require(not path.is_absolute() and str(path)==name and '..' not in path.parts and
                not (ROOT/path).is_symlink() and name not in entries,'Invalid manifest path '+name)
        entries[name] = digest
    actual = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file()
              and '__pycache__' not in p.parts and p.name!='SHA256SUMS'}
    require(set(entries)==actual,'Files differ from manifest: '+str(set(entries)^actual))
    for name,digest in entries.items():
        require(checksum(ROOT/name)==digest,'Checksum mismatch: '+name)
    require(entries.get('background/forward.cpp')==PINNED_SOLVER,'Historical Maxwell solver differs')
    print(f'Verified {len(entries)} files, including the historical Maxwell solver.',flush=True)


def paths(base):
    return {section: base/('quantum' if section=='mass' else section)/(name+'.json') for section,name in PRODUCTS.items()}


def compare_values(a,b,path):
    loc='/'.join(map(str,path))
    if isinstance(b,dict):
        require(isinstance(a,dict) and a.keys()==b.keys(),'Different keys at '+loc)
        for k in b:
            compare_values(a[k],b[k],(*path,k))
    elif isinstance(b,list):
        require(isinstance(a,list) and len(a)==len(b),'Different list at '+loc)
        for i,(left,right) in enumerate(zip(a,b)):
            compare_values(left,right,(*path,i))
    elif type(b) in (int,float) and type(a) in (int,float):
        # Full machine runs agree far more closely on the dimensionless headline
        # observables. Diagnostics involving differences of small numbers need an
        # absolute floor. Separate physical invariants below protect the controls.
        absolute=2.e-10
        if path[-1] in ('equal_energy_shell2_relative_direct_error','charged_shell2_timestep_relative_discrepancy'):
            absolute=2.e-9
        if path[0]=='quantum' and ('occupation' in path or 'low_modes' in path):
            absolute=1.e-13
        require(math.isclose(a,b,rel_tol=5.e-7,abs_tol=absolute),f'{loc}: {a} differs from {b}')
    else:
        require(type(a) is type(b) and a==b,'Value differs at '+loc)


def verify_trace(got):
    import numpy as np
    original=ROOT/'background/trace.csv.gz'
    def load(path):
        with gzip.open(path,'rt') as raw:
            header=raw.readline()
            values=np.loadtxt(raw,delimiter=',')
        return header,values
    baseline,ref=load(original)
    header,data=load(got)
    require(header==baseline and data.shape==ref.shape,'Changed trace layout')
    require(np.all(np.diff(data[:,0])>0) and np.all(data[:,1]>0),'Nonpositive time step or radius')
    scale=np.maximum(1.e-8,np.max(abs(ref),axis=0))
    require(np.all(np.max(abs(ref-data),axis=0)<=2.e-8*scale+1.e-13),'Regenerated trace differs')
    names=header.strip().split(',')
    col={name:data[:,i] for i,name in enumerate(names)}
    residual=col['sheet']-.5*col['K']*(col['X']**2+col['Bmag']**2)
    require(np.max(abs(residual))<1.e-12 and np.all(col['rad']==0),'Neutral Maxwell budget or fluid changed')


def controls(results):
    k=results['kinetic']
    require(k['rate_scan'][0]['random_share']==0,'Collisionless cold beams gained random energy')
    require(0<k['setup_and_bounds']['all_rate_relativistic_kinetic_over_EM_bound']<2.e-7,
            'Work bound for dilute tracers failed')
    s=results['theory']['controlled_relaxation']
    require(math.isclose(s['angular_only']['current_relative_to_initial'],.01,abs_tol=1.e-11) and
            math.isclose(s['kinetic_only']['current_relative_to_initial'],.01,abs_tol=1.e-11),
            'Matched damping of the two currents failed')
    require(s['angular_only']['final_radial_number_TV'] > 90*s['kinetic_only']['final_radial_number_TV'],
            'Radial spectrum control failed')
    q=results['quantum']['protocols']['zero_crossing']
    shell=[r for r in q['q1']['low_modes'] if r['n']==2]
    require(len(shell)==3 and max(r['omega'] for r in shell)-min(r['omega'] for r in shell)<1.e-10 and
            max(r['occupation'] for r in shell)>4.e4*min(r['occupation'] for r in shell),
            'Equal-energy quantum occupation test failed')
    require(q['controls']['max_free_occupation']<1.e-22 and q['q1']['max_wronskian_error']<1.e-9,
            'Free quantum control or conserved mode norm failed')
    low=q['q1']['low_modes']
    state_mean=lambda n:sum(r['occupation'] for r in low if r['n']==n)/(n+1)
    require(590<state_mean(1)/state_mean(0)<610,'Mean population inversion changed')
    require(q['controls']['equal_energy_shell2_relative_direct_error']<2.e-6,
            'Independent equal-energy occupation integration failed')
    mass=results['mass']['cases']
    require(all(r['charged_shell2_timestep_relative_discrepancy']<3.e-5 for r in mass) and
            mass[0]['charged_shell2_max_to_min']>5.e4 and mass[-1]['charged_shell2_max_to_min']<1.1,
            'Mass sensitivity or time refinement check failed')
    gate=results['gate']['runs']
    require(gate['fine']['maximum_intrinsic_fraction']>.9 and
            gate['broad']['maximum_intrinsic_fraction']<.9 and
            gate['no_force']['end_temperature']==gate['collisionless']['end_temperature']==0,
            'Locally regulated onset or zero-heating controls failed')
    for name in ('taper_short','taper_long'):
        data=results['quantum']['protocols'][name]['q1']
        require(0<data['probability_of_any_pair_finite_band']<=data['total_quanta_finite_band']/2<.002,
                'Finite-band probability control failed')
    c=results['coupled']['runs']
    primary=c['primary']
    require(primary['bounce_radius']>0 and primary['end_H']<0 and primary['exit_complete_tau']>primary['bounce_tau'],
            'The charged wall did not complete the reported bounce and exit')
    require(primary['maximum_intrinsic_fraction']>.9 and .010 < primary['ninety_intrinsic_expanding_efolds'] < .012 and primary['end_radiation']==0,
            'The selected gas did not acquire a majority intrinsic share with the separate fluid disabled')
    require(c['activation_4p7']['maximum_intrinsic_fraction'] < .02, 'Early mobility did not suppress heating')
    stored=primary['activation_radius']**4*primary['activation_sheet_density']
    inflow=primary['comoving_wave_to_sheet']-primary['activation_bulk_sheet_flux']
    release=primary['comoving_scalar_to_sheet']-primary['activation_scalar_sheet_flux']
    budget=stored+inflow+release-primary['end_charge_work_a4']-primary['end_comoving_sheet_energy']
    require(abs(budget)<1.e-6 and inflow<0 and abs(release)<1.e-8, 'Post-activation Maxwell work budget failed')
    require(min(primary['min_bulk_mass'],primary['min_lapse'],primary['min_radial_gradient_norm'])>0,
            'Charged wall failed the positive-mass branch test')
    for name,model in results['coupled']['parameters'].items():
        if model['nu']==0:
            require(abs(c[name]['end_intrinsic_random_density'])<1.e-12,'Adiabatic cold control heated: '+name)
    require(results['synthesis']['matching_example']['total_seed_carriers_per_unit_S0']>60000,
            'Particle count conversion changed')


def reproduce(out):
    out=out.resolve()
    require(not out.is_relative_to(ROOT) and not any(out.glob('*')),'Use an empty output folder outside the study')
    out.mkdir(parents=True,exist_ok=True)
    targets=paths(out)
    for path in targets.values():
        path.parent.mkdir(exist_ok=True)
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    def run(path,*args):
        print('Recomputing '+path,flush=True)
        subprocess.run([sys.executable,'-I','-B',str(ROOT/path),*map(str,args)],env=env,check=True)
    def wave_and_probes():
        run('background/generate.py','--output-dir',out/'background')
        trace=out/'background/trace.csv.gz'
        verify_trace(trace)
        for area,prog,flag in [('kinetic','kinetic','--input'),('theory','spectrum','--input'),
                               ('quantum','pairs','--trace')]:
            run(area+'/'+prog+'.py',flag,trace,'--output',targets[area])
        run('quantum/mass.py','--trace',trace,'--output',targets['mass'])
    with ThreadPoolExecutor(max_workers=2) as pool:
        background=pool.submit(wave_and_probes)
        run('coupled/charged.py','--output',targets['coupled'],'--jobs','2')
        run('gate/gate.py','--output',targets['gate'],'--jobs','2')
        background.result()
    run('coupled/phase.py',targets['coupled'])
    run('synthesis/bridge.py','--inputs',out,'--output',targets['synthesis'])
    results={name:json.loads(p.read_bytes()) for name,p in targets.items()}
    controls(results)
    for name,p in paths(ROOT).items():
        compare_values(results[name],json.loads(p.read_bytes()),(name,))
    print('All eight products, the dense field history and independent physical controls passed.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify',action='store_true')
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args()
    verify_files()
    if args.output_dir and not args.verify:
        parser.error('--output-dir requires --verify')
    if args.verify:
        if args.output_dir:
            reproduce(args.output_dir)
        else:
            with tempfile.TemporaryDirectory(prefix='hdblast-charged-replay-') as tmp:
                reproduce(Path(tmp))
