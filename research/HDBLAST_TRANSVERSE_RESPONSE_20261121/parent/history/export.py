"""Capture each completed step of the pinned directional-birth wall evolution."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import runpy
import shlex
import subprocess
import tempfile


COLUMNS = 't a H X D theta u R tried_rest e intrinsic c d target'.split()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior', type=Path, required=True)
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--jobs', type=int, default=2)
    args = parser.parse_args()
    source = args.prior.resolve()
    dest = args.outdir.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    feedback = runpy.run_path(str(source/'feedback/feedback.py'))
    change = feedback['change']
    base = source/'prior'
    envelope = json.loads((base/'envelope/ENVELOPE.json').read_bytes())
    reference = json.loads((source/'boost/EVOLUTION.json').read_bytes())
    shared_ref = json.loads((source/'boost/FLAVOR.json').read_bytes())
    config = dict(envelope['parameters']['mass_4e12'], sample=0)
    action = reference['action']
    transform = runpy.run_path(str(base/'envelope/build.py'))['transform']
    src = transform((base/'coupled/joined_forward.cpp').read_text(), centered=True)
    src = runpy.run_path(str(source/'boost/evolve.py'))['make_source'](
        src, feedback, (source/'feedback/angular.hpp').read_text())
    src = '#include <fstream>\n' + src
    src = change(src, 'if (argc != 34)', 'if (argc != 36)')
    src = change(src, 'int collision_mode = 0;', 'double branching=1.;\nint collision_mode = 0;')
    src = change(src, 'collision_mode=std::stoi(argv[33]);',
                 'collision_mode=std::stoi(argv[33]); branching=std::stod(argv[34]);\n'
                 '    std::ofstream trace(argv[35], std::ios::binary);\n'
                 '    trace.exceptions(std::ios::badbit | std::ios::failbit);')
    src = change(src, 'return proper_density*boosted_rate(', 'return branching*proper_density*boosted_rate(')
    anchor = '        const double margin=gas.intrinsic_random-.9*wnew.e;'
    src = change(src, anchor, '''        {
            const auto ch = oscillator(y);
            const double R = ch.mass_ratio;
            const double z = R>1. ? std::sqrt((R-1.)*(R+1.)) : 0.;
            const double c = z>0. ? h*z+R*ch.mass_slope*y[PI]/(2.*daughter_mass*z) : 0.;
            const double d = charge_mass*wnew.x;
            const double target = envelope_active && ch.power>0. && z>0. ?
                std::pow(std::acos(-1.),2)/6.*feedback_action/std::pow(daughter_mass,4)*ch.power/(R*z*z) : 0.;
            const double row[] = {t+dt,a,h,wnew.x,gas.comoving_rest,gas.theta,gas.u,R,
                                  y[TRIED_REST],wnew.e,gas.intrinsic_random,c,d,target};
            trace.write(reinterpret_cast<const char*>(row), sizeof(row));
        }
''' + anchor)
    run = runpy.run_path(str(base/'coupled/joined.py'))['run']
    cases = [(a,b,s,step) for a,b in [('dedicated',1.),('shared',1/3)]
             for s,step in [('native',1.),('half',.5)]]
    with tempfile.TemporaryDirectory(prefix='hdblast-event-clock-') as t:
        tmp = Path(t)
        (tmp/'capture.cpp').write_text(src)
        exe = tmp/'capture'
        subprocess.run(['g++','-std=c++17','-O3','-I',str(source/'boost'),
                        str(tmp/'capture.cpp'),'-o',str(exe)],check=True)
        def evaluate(case):
            channel, weight, resolution, step = case
            name = channel + '_' + resolution
            path = dest/(name+'.bin')
            wrapper = tmp/name
            tail = ['2',repr(action),str(source/'boost/BOOST_TABLE.dat'),'1',repr(weight),str(path)]
            wrapper.write_text('#!/bin/sh\nexec '+shlex.quote(str(exe))+' "$@" '+shlex.join(tail)+'\n')
            wrapper.chmod(0o755)
            row = run(wrapper,dict(config, cfl=step*config['cfl']))
            if channel=='dedicated':
                ref = reference['runs']['pole_cap' + ('_half' if step==.5 else '')]
            else:
                ref = shared_ref['inclusive_'+resolution]
            for field in ('maximum_radius','maximum_tau','first_hot_time','peak_hot_fraction','doubling_aH'):
                if not math.isclose(row[field],ref[field],rel_tol=3e-5,abs_tol=2e-8):
                    raise ValueError('Native replay differs from reference: '+name+' '+field)
            fields = ['first_hot_time','doubling_time','maximum_tau','opening_time','peak_hot_time',
                      'end_tau','end_rest','maximum_radius','first_hot_aH','doubling_aH',
                      'end_electric_work','abs_electric','tried_rest','lost_rest']
            out = {field:row[field] for field in fields}
            out.update(channel_fraction=weight, rows=path.stat().st_size//(8*len(COLUMNS)))
            print(name, out['rows'], 'first hot', out['first_hot_time'], flush=True)
            return name,out
        with ThreadPoolExecutor(args.jobs) as pool:
            runs = dict(pool.map(evaluate,cases))
    metadata = dict(columns=COLUMNS, format='native-endian contiguous IEEE-754 binary64 rows',
                    daughter_mass_L=config['daughter_mass'],charge_to_mass=config['charge_to_mass'],
                    action=action, runs=runs)
    (dest/'metadata.json').write_text(json.dumps(metadata,indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
