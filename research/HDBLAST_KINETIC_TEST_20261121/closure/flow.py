"""Evolve the wall with mutually consistent rest-frame stresses and collision distributions."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import runpy
import shlex
import subprocess
import tempfile

import numpy as np
from scipy.special import kve

HERE=Path(__file__).resolve().parent


def eos_table(path):
    # In the table e is specific thermal energy, excluding rest mass, divided by Theta.
    # A centered derivative avoids cancellation between large Bessel-ratio terms.
    x=np.linspace(-4.,4.,401)*np.log(10.)
    def coefficient(a):
        t=np.exp(a)
        return (kve(1,1/t)/kve(2,1/t)-1.)/t+3.
    dx=1.e-3
    derivative=(-coefficient(x+2*dx)+8*coefficient(x+dx)-8*coefficient(x-dx)+coefficient(x-2*dx))/(12*dx)
    out=''
    for name, data in [('synge_nodes',coefficient(x)),('synge_slopes',derivative*np.log(10.)/50)]:
        out+=('inline constexpr double '+name+'[] = {\n'+',\n'.join(format(v,'.17g') for v in data)+'\n};\n')
    path.write_text(out)


def source(prior):
    base=prior/'prior'
    fb=runpy.run_path(str(prior/'feedback/feedback.py'))
    alter=fb['change']
    transform=runpy.run_path(str(base/'envelope/build.py'))['transform']
    src=transform((base/'coupled/joined_forward.cpp').read_text(),centered=True)
    make=runpy.run_path(str(prior/'boost/evolve.py'))['make_source']
    src=make(src,fb,(prior/'feedback/angular.hpp').read_text())
    src=alter(src,'if (argc != 34)', 'if (argc != 36)')
    src=alter(src,'struct Gas {', 'int kinetic_closure=0;\ndouble flavor_weight=1.;\nstruct Gas {')
    src=alter(src,'collision_mode=std::stoi(argv[33]);',
              'collision_mode=std::stoi(argv[33]); kinetic_closure=std::stoi(argv[34]); flavor_weight=std::stod(argv[35]);')
    src=alter(src,'    g.energy=g.dilution*E;', '''    if(kinetic_closure==1) {
        const double ratio=moment/g.comoving_rest, surplus=extra/g.comoving_rest;
        g.theta=synge::invert(ratio,surplus,g.theta);
        thermal=synge::thermal(g.theta).specific;
        const double enthalpy=1.+thermal+g.theta;
        g.u=ratio/enthalpy;
        g.gamma=std::hypot(1.,g.u); g.velocity=g.u/g.gamma;
        g.lab_pressure=g.comoving_rest*g.theta/g.gamma;
    }
    g.energy=g.dilution*E;''')
    src=alter(src,'const double enthalpy=2.5*previous_gas.theta+std::hypot(1.,1.5*previous_gas.theta);',
              '''const double tgas=previous_gas.theta;
            const double enthalpy=kinetic_closure==1 ? 1.+tgas+synge::thermal(tgas).specific :
                2.5*tgas+std::hypot(1.,1.5*tgas);''')
    src=alter(src,'return proper_density*boosted_rate(gas.theta,gas.u)/(daughter_mass*daughter_mass);',
              '''double c=boosted_rate(gas.theta,gas.u);
    if(kinetic_closure==1) {
        // The archived Juttner force had already been divided by h_TM.
        const double h_tm=2.5*gas.theta+std::hypot(1.,1.5*gas.theta);
        c*=h_tm/(1.+gas.theta+synge::thermal(gas.theta).specific);
    }
    return flavor_weight*proper_density*c/(daughter_mass*daughter_mass);''')
    return '#include "synge.hpp"\n'+src


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--jobs',type=int,default=2)
    parser.add_argument('--spatial',action='store_true',help='Vary early and late radial element counts for the shared pole')
    parser.add_argument('--convergence',action='store_true',help='Resolve early collar and remeshing sensitivity at four additional element counts')
    args=parser.parse_args()
    parent=args.prior.resolve()
    base=parent/'prior'
    ref=json.loads((parent/'boost/EVOLUTION.json').read_text())
    cfg=json.loads((base/'envelope/ENVELOPE.json').read_text())['parameters']['mass_4e12']
    run=runpy.run_path(str(base/'coupled/joined.py'))['run']
    fields=('maximum_radius maximum_tau stop_turnaround first_hot_time first_hot_aH first_hot_radius '
            'peak_hot_fraction doubling_hot_fraction doubling_aH doubling_theta max_gamma max_nu '
            'end_electric_work abs_electric first_hot_integral_nu end_rest_drag_work gas_balance_error '
            'gas_primitive_error surface_balance_error friedmann_error end_rest end_decay_energy').split()
    with tempfile.TemporaryDirectory(prefix='hdblast-closure-') as temp:
        d=Path(temp)
        eos_table(d/'synge_table.inc')
        lookup=(parent/'boost/lookup.hpp').read_text()
        lookup=lookup.replace('if (theta == 0.) return cold_rate(theta,u);',
                 'if (theta == 0.) return shell_mode ? shellbath::rate(theta,u) : cold_rate(theta,u);')
        lookup=lookup.replace('if (it < 0.) return cold_rate(theta,u);',
                 'if (it < 0.) return shell_mode ? shellbath::rate(theta,u) : cold_rate(theta,u);')
        (d/'lookup.hpp').write_text('#include "shell.hpp"\ninline bool shell_mode=false;\n'+lookup)
        src=source(parent).replace('flavor_weight=std::stod(argv[35]);','flavor_weight=std::stod(argv[35]); shell_mode=kinetic_closure==2;')
        (d/'evolve.cpp').write_text(src)
        for fn,exe in [(HERE/'grid.cpp',d/'grid'),(d/'evolve.cpp',d/'evolve')]:
            subprocess.run(['g++','-O3','-std=c++17','-include','algorithm','-I',str(d),'-I',str(HERE),str(fn),'-o',str(exe)],check=True)
        for resolution in (1,2):
            with (d/('shell'+str(resolution)+'.dat')).open('wb') as out:
                subprocess.run([str(d/'grid'),str(resolution)],stdout=out,check=True)
        jobs=[(name,index,flavor,step) for name,index in [('thermal_synge',1),('sphere_tm',2)]
              for flavor in ['dedicated','shared'] for step in ['native','refined']]
        jobs.append(('hybrid_reference',0,'dedicated','native'))
        if args.spatial:
            jobs=[(name,index,'shared',mesh) for name,index in [('thermal_synge',1),('sphere_tm',2)]
                  for mesh in ['original','early_double','late_double','both_double']]
        if args.convergence:
            jobs=[(name,index,'shared',str(elements)) for name,index in [('thermal_synge',1),('sphere_tm',2)]
                  for elements in (36,72,96,144)]
        def evaluate(job):
            name,index,flavor,step=job
            key='_'.join([name,flavor,step])
            table=d/('shell'+str(2 if step=='refined' or args.spatial or args.convergence else 1)+'.dat') if index==2 else parent/'boost/BOOST_TABLE.dat'
            wrapper=d/key
            tail=['2',repr(ref['action']),str(table),'1',str(index),repr(1. if flavor=='dedicated' else 1./3)]
            wrapper.write_text('#!/bin/sh\nexec '+shlex.join([str(d/'evolve')])+' "$@" '+shlex.join(tail)+'\n')
            wrapper.chmod(0o755)
            p=dict(cfg,sample=0)
            if step=='refined': p['cfl']*=.5
            if args.convergence: p['early_cells']=int(step)
            if args.spatial:
                if step in ('early_double','both_double'): p['early_cells']*=2
                if step in ('late_double','both_double'): p['late_cells']*=2
            try: r=run(wrapper,p)
            except subprocess.CalledProcessError as e: raise RuntimeError(key+': '+str(e.stderr)) from e
            row={k:r[k] for k in fields if k in r}
            row['hot_growth']=math.exp(r.get('hot_efolds',0.))
            if not row['stop_turnaround'] or max(row[k] for k in ['gas_balance_error','surface_balance_error','gas_primitive_error'])>1e-5:
                raise ValueError(f'Invalid evolution {key}: {row}')
            print(key,row['hot_growth'],row['peak_hot_fraction'],row['max_nu'],flush=True)
            return key,row
        with ThreadPoolExecutor(args.jobs) as pool:
            rows=dict(pool.map(evaluate,jobs))
    if not(args.spatial or args.convergence) and abs(rows['hybrid_reference_dedicated_native']['hot_growth']-ref['runs']['pole_cap']['hot_radius_factor'])>1e-8:
        raise ValueError('Historical control failed')
    args.output.write_text(json.dumps(rows,indent=2,sort_keys=True,allow_nan=False)+'\n')

if __name__=='__main__':
    main()
