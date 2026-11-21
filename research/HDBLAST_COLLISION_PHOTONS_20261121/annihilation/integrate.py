"""Two-photon annihilation on the previously evolved homogeneous dedicated-shell histories."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import subprocess
import tempfile

import numpy as np


def build(study, polar, scatter, work):
    prefix, sep, _ = (polar/'feedback/homogeneous.cpp').read_text().partition('int main(int argc,char** argv)')
    if not sep:
        raise ValueError('Unexpected parent particle representation')
    src=(study/'live/solver.cpp').read_text()
    for old,new in (
            ('if(argc!=12)', 'if(argc!=13)'),
            ('    PoleCollisions collisions', '    Annihilation annihilation(particles.size());\n    PoleCollisions collisions'),
            ('        for(size_t i=0;i<n;++i){\n            auto& v=particles[i];\n            double start=',
             '        annihilation.begin(particles,n);\n        for(size_t i=0;i<n;++i){\n            auto& v=particles[i];\n            double start='),
            ('        collisions.step(particles,','        annihilation.add(particles,n,weight,l[0],r[0],l[1],r[1]);\n        collisions.step(particles,'),
            ('    collisions.save(argv[8]);','    annihilation.save(argv[12]);\n    collisions.save(argv[8]);')):
        if src.count(old)!=1: raise ValueError('Unexpected parent solver marker '+old)
        src=src.replace(old,new)
    cpp=work/'instrument.cpp'
    cpp.write_text(prefix+'\n#include "pairs.hpp"\n'+src)
    program=work/'instrument'
    subprocess.run(['g++','-O3','-ffp-contract=off','-std=c++17','-I',str(Path(__file__).resolve().parent),'-I',str(scatter/'live'),str(cpp),'-o',str(program)],check=True)
    return program


def case(task, state):
    time,power,steps,seed,pole=task
    hist,orbits,table_builder,program,work,alpha,refs=state
    mark=dict(time=time,**{k:float(np.interp(time,hist['t'],hist[k])) for k in ('D','a','theta','R','X','B')})
    table=table_builder(hist,time,steps,None)
    birth_seed=orbits['SCRAMBLE']+{12:1,13:71,14:121}[power]
    tags=orbits['seeds'](hist,mark,power,birth_seed)
    tags[:,1]=np.interp(tags[:,0],hist['t'],hist['a'])
    mu=float(hist['daughter_mass_L']); action=float(hist['action']); density=action/(6*mu**3)
    with tempfile.TemporaryDirectory(dir=work) as td:
        d=Path(td)
        table.tofile(d/'input');tags.tofile(d/'births')
        subprocess.run([str(program),str(d/'input'),str(d/'births'),str(mark['D']),'2','1',str(density if pole else 0.),str(seed),
                        str(d/'counts'),str(d/'particles'),str(d/'trace'),str(d/'ledger'),str(d/'anns')],check=True)
        events=float(np.loadtxt(d/'counts')[1])
        gas=float((d/'ledger').read_text().splitlines()[-1].split()[0])
        vals=np.loadtxt(d/'anns')*(np.pi*alpha**2*density)
    ref=next(r for r in refs if r['mode']=='particle' and r['time']==time and r['pole']==pole and r['steps']==steps and r['seed']==seed and r['pairs']==2**power)
    if events/2**power!=ref['mean_encounters'] or not np.isclose(gas,ref['gas'],rtol=2e-6): raise ValueError('Parent trajectory mismatch')
    count=vals[:-1,0]/mark['D']; energy=vals[:-1,1]/gas
    mean_count,mean_energy=count.mean(),energy.mean()
    # Geometry-only maximum hazard on this interval, ignoring dilution and all pair momenta.
    amin=float(np.min(table[:,1]))
    any_daughter_bound=17./8*np.pi*alpha**2*density*mark['D']*(time-table[0,0])/amin**3
    if not (0<count.max()<=vals[-1,0]/mark['D']<=any_daughter_bound*(1+1.e-12) and
            0<energy.max()<=vals[-1,1]/gas):
        raise ValueError('Positive estimator, density ceiling, or energy ceiling failed')
    return dict(time=time,log2_pairs=power,steps=steps,seed=seed,pole=pole,charged_rest=mark['D'],gas_conformal=gas,
                births_conformal_rest_at_final_radius=mark['a']*mark['D'], mean_encounters=events/2**power,
                photons_per_final_daughter=float(mean_count), photon_energy_over_final_gas=float(mean_energy),
                photon_energy_conformal=float(vals[:-1,1].mean()),
                photon_energy_over_final_charged_rest_conformal=float(vals[:-1,1].mean()/(mark['a']*mark['D'])),
                expected_photons_whole_unit_round_wall=float(2*np.pi**2*action/mu*mark['D']*mean_count),
                any_daughter_geometry_only_optical_depth_bound=float(any_daughter_bound),
                count_bound_per_final_daughter=float(vals[-1,0]/mark['D']),energy_bound_over_final_gas=float(vals[-1,1]/gas),
                quadrature_count_per_daughter=count.tolist(),quadrature_energy_over_gas=energy.tolist())


def main(args):
    hist=dict(np.load(args.history))
    if not np.all(hist['K']==2.54) or float(hist['channel_fraction'])!=1. or hist['t'][-1]<7.05:
        raise ValueError('Expected extended dedicated history')
    alpha=float(3*hist['charge_to_mass']**2*hist['daughter_mass_L']**2/(4*np.pi*hist['action']*2.54))
    orbits=runpy.run_path(str(args.polarization/'orbits/pairs.py'))
    table_builder=runpy.run_path(str(args.study/'live/evolve.py'))['field_table']
    refs=sum((json.loads((args.study/f).read_bytes())['runs'] for f in ('live/LIVE.json','live/CROSSINGS.json')),[])
    tasks=[(6.45,12,1024,18431,True),(6.67,12,2048,18431,True),(6.67,14,4096,91283,True),
           (7.05,12,4096,18431,True),(7.05,13,4096,70163,True),(7.05,13,8192,70163,True),
           (7.05,13,4096,70163,False),(7.05,14,8192,91283,True)]
    with tempfile.TemporaryDirectory(prefix='dirac-annihilation-') as td:
        work=Path(td); program=build(args.study,args.polarization,args.scattering,work)
        state=(hist,orbits,table_builder,program,work,alpha,refs)
        with ThreadPoolExecutor(args.jobs) as pool: data=list(pool.map(lambda t:case(t,state),tasks))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(dict(alpha=alpha,normalized_kernel_bound=17./8.,samples_per_batch_per_step=256,runs=data),
                                     indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ('study','polarization','scattering','history','output'): p.add_argument('--'+arg,type=Path,required=True)
    p.add_argument('--jobs',type=int,default=2)
    main(p.parse_args())
