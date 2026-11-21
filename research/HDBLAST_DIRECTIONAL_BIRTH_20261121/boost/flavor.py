"""Test the change in the charge-odd drag moment when three equal flavors share one pole."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import runpy
import shlex
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--base',type=Path,default=HERE.parent/'prior')
    cli.add_argument('--feedback',type=Path,default=HERE.parent/'feedback')
    cli.add_argument('--table',type=Path,default=HERE/'BOOST_TABLE.dat')
    cli.add_argument('--output',type=Path,required=True)
    cli.add_argument('--jobs',type=int,default=2)
    args=cli.parse_args()
    env=json.loads((args.base/'envelope/ENVELOPE.json').read_text())
    feedback=runpy.run_path(str(args.feedback/'feedback.py'))
    transform=runpy.run_path(str(args.base/'envelope/build.py'))['transform']
    run=runpy.run_path(str(args.base/'coupled/joined.py'))['run']
    src=transform((args.base/'coupled/joined_forward.cpp').read_text(),centered=True)
    src=runpy.run_path(str(HERE/'evolve.py'))['make_source'](src,feedback,(args.feedback/'angular.hpp').read_text())
    change=feedback['change']
    src=change(src,'int collision_mode = 0;', 'double branching=1.;\nint collision_mode = 0;')
    src=change(src,'if (argc != 34)','if (argc != 35)')
    src=change(src,'collision_mode=std::stoi(argv[33]);','collision_mode=std::stoi(argv[33]); branching=std::stod(argv[34]);')
    src=change(src,'return proper_density*boosted_rate(', 'return branching*proper_density*boosted_rate(')
    action=5*env['screen']['mass_4e12']['elastic_drag_action_floor']
    jobs=[(name,fraction,mesh) for name,fraction in [('inclusive',1/3),('elastic_component',1/9)]
          for mesh in ('native','half')]
    fields='peak_hot_fraction first_hot_time first_hot_radius first_hot_aH doubling_hot_fraction doubling_aH maximum_radius maximum_tau max_gamma max_nu end_electric_work abs_electric integral_nu end_rest_drag_work stop_turnaround surface_balance_error gas_balance_error friedmann_error'.split()
    with tempfile.TemporaryDirectory(prefix='hdblast-shared-pole-') as temp:
        dest=Path(temp)
        cpp=dest/'flavor.cpp'
        cpp.write_text(src)
        exe=dest/'flavor'
        subprocess.run(['g++','-std=c++17','-O3','-I',str(HERE),str(cpp),'-o',str(exe)],check=True)
        def evaluate(job):
            name,branch,mesh=job
            wrapper=dest/(name+'_'+mesh)
            cmd=shlex.join(['2',repr(action),str(args.table.resolve()),'1',repr(branch)])
            wrapper.write_text('#!/bin/sh\nexec '+shlex.quote(str(exe))+' "$@" '+cmd+'\n')
            wrapper.chmod(0o755)
            p=dict(env['parameters']['mass_4e12'],sample=0)
            if mesh=='half':
                p['cfl']*=.5
            r=run(wrapper,p)
            if not r['stop_turnaround'] or max(r['gas_balance_error'], r['surface_balance_error'])>1e-5:
                raise ValueError('Invalid evolution '+str(job))
            row={k:r[k] for k in fields if k in r}
            row.update(moment_weight=branch,hot_radius_factor=math.exp(r.get('hot_efolds',0.)))
            print(wrapper.name,row['peak_hot_fraction'],row['hot_radius_factor'],flush=True)
            return wrapper.name,row
        with ThreadPoolExecutor(args.jobs) as pool:
            rows=dict(pool.map(evaluate,jobs))
    args.output.write_text(json.dumps(rows,indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
