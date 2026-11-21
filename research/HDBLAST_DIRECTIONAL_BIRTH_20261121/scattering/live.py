"""Restart the transport closure at the first expanding hot, low-drift state and evolve gravity."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import runpy
import subprocess
import tempfile

import numpy as np

HERE=Path(__file__).resolve().parent
R=runpy.run_path(str(HERE/'resonance.py'))


def table_source(count=1024):
    # Direct quadrature of the exact positive distribution at geometric temperature nodes.
    z,w,_=R["distribution"](8192)
    mid=np.sqrt(z[:-1]*z[1:])
    th=np.geomspace(100.,1.e5,count+1)
    table=[]
    for temp in th:
        r=R['VECTOR_MASS']/temp; b=R['VECTOR_WIDTH']
        table.append(.75*float(w @ (1/(1+((r/mid)**2-1)**2/b**2))))
    values=','.join(repr(v) for v in table)
    return f"""
bool pole_active=false;
int pole_mode=0;
const std::array<double,{count+1}> thermal_table{{{values}}};
double cubic_weight(double u) {{
    u=std::clamp(u,0.,1.); return u*u*(3.-2.*u);
}}
double pole_weight(const Gas& gas) {{
    if(pole_mode<3) return pole_active && pole_mode==2 ? 1. : 0.;
    return cubic_weight((gas.theta-100.)/100.)*cubic_weight((.1-std::abs(gas.velocity))/.05);
}}
double thermal_part(const Gas& gas,double weight) {{
    if(weight==0.) return 0.;
    if(gas.theta<100. || gas.theta>1.e5 || std::abs(gas.velocity)>.101)
        throw std::runtime_error("Thermal continuation outside its low-drift domain");
    const double u=std::log(gas.theta/100.)*({count}./std::log(1000.));
    const int j=std::clamp(static_cast<int>(u),0,{count-1});
    const double value=thermal_table[j]+(u-j)*(thermal_table[j+1]-thermal_table[j]);
    return weight*value*std::acos(-1.)*feedback_action*gas.rest/
           (6.*gas.gamma*std::pow(daughter_mass,3)*gas.theta*gas.theta);
}}
double remaining_drag(double weight) {{
    if(pole_mode==3 || (pole_mode==1 && pole_active)) return 0.;
    return (1.-weight)*drag_rate;
}}
double resonance_drag(const Gas& gas) {{
    const double weight=pole_weight(gas);
    return thermal_part(gas,weight)+remaining_drag(weight);
}}
"""


def with_pole(src,change,count):
    src=change(src,'struct Solver {',table_source(count)+'\nstruct Solver {')
    src=change(src,'ABS_BIRTH_KICK,','ABS_BIRTH_KICK, ABS_EWORK, WEIGHTED_ABS_EWORK, MICRO_HEAT_WORK,')
    src=change(src,'dy[REST_WORK_DRAG]=drag_rate*gas.u*moment;',
               '''dy[REST_WORK_DRAG]=drag_rate*gas.u*moment;
        dy[ABS_EWORK]=std::abs(dy[ELECTRIC_WORK]);
        dy[WEIGHTED_ABS_EWORK]=pole_weight(gas)*std::abs(dy[ELECTRIC_WORK]);
        dy[MICRO_HEAT_WORK]=thermal_part(gas,pole_weight(gas))*gas.u*moment;''')
    src=change(src,'    out["end_tau"] = t;', '''    out["absolute_electric_work"]=y[ABS_EWORK];
    out["weighted_absolute_electric_work"]=y[WEIGHTED_ABS_EWORK];
    out["thermal_drag_work"]=y[MICRO_HEAT_WORK];
    out["all_drag_work"]=y[REST_WORK_DRAG];
    out["end_tau"] = t;''')
    src=change(src,'        const double moment=y[DIRECTED];',
               '        const double moment=y[DIRECTED];\n        const double local_drag=resonance_drag(gas);')
    for old in ('(h+drag_rate)*moment', '+drag_rate*moment*vc', '=drag_rate*gas.u*moment'):
        src=change(src,old,old.replace('drag_rate','local_drag'))
    return change(src,'        const double margin=gas.intrinsic_random-.9*wnew.e;',
    '''        if(!pole_active && h>0. && y[PHI]>=fall_phi+fall_width && gas.intrinsic_random>.9*wnew.e &&
           gas.theta>=100. && std::abs(gas.velocity)<=.1) {
            pole_active=true;
            out["pole_start_time"]=t+dt;
            out["pole_start_radius"]=a;
            out["pole_start_theta"]=gas.theta;
            out["pole_start_electric_work"]=y[ELECTRIC_WORK];
            out["pole_start_input_at_birth"]=y[DECAY_ENERGY];
            out["pole_start_rest"]=y[DAUGHTER_REST];
        }
        if(pole_active && h>0.) {
            const double nu=resonance_drag(gas);
            if(!out.count("pole_drag_min")) {
                out["pole_drag_min"]=nu;
                out["pole_drag_max"]=nu;
            }
            out["pole_drag_min"]=std::min(out["pole_drag_min"],nu);
            out["pole_drag_max"]=std::max(out["pole_drag_max"],nu);
            out["pole_max_abs_drift"]=std::max(out["pole_max_abs_drift"],std::abs(gas.velocity));
            out["pole_drag_integral"]+=dt*nu;
        }
        const double margin=gas.intrinsic_random-.9*wnew.e;''')


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--base",type=Path,required=True)
    cli.add_argument("--output",type=Path,required=True)
    cli.add_argument("--jobs",type=int,default=2)
    cli.add_argument("--only",nargs='+')
    args=cli.parse_args()
    prior=args.base
    env=json.loads((prior/'envelope/ENVELOPE.json').read_text())
    src=runpy.run_path(str(prior/'envelope/build.py'))['transform']((prior/'coupled/joined_forward.cpp').read_text(),centered=True)
    if hashlib.sha256(src.encode()).hexdigest()!=env['centered_solver_sha256']:
        raise ValueError('Parent solver changed')
    feedback=runpy.run_path(str(HERE.parent/'feedback/feedback.py'))
    src=feedback['with_feedback'](src,(HERE.parent/'feedback/angular.hpp').read_text())
    change=feedback['change']
    src=change(src,'if (argc != 32)','if (argc != 33)')
    src=change(src,'feedback_law=FEEDBACK_LAW; feedback_action=FEEDBACK_ACTION;',
               'feedback_law=std::stoi(argv[30]); feedback_action=std::stod(argv[31]); pole_mode=std::stoi(argv[32]);')
    p=env['parameters']['mass_4e12']
    action=5*env['screen']['mass_4e12']['elastic_drag_action_floor']
    runner=runpy.run_path(str(prior/'coupled/joined.py'))['run']
    results={}
    jobs=[]
    with tempfile.TemporaryDirectory(prefix='hdblast-live-resonance-') as directory:
        directory=Path(directory)
        for count,combos in [(1024,[(0,2,'vector'),(0,1,'zero'),(2,2,'capped_vector'),
                                   (3,3,'early_pole_alone'),(3,4,'early_pole_matched'),
                                   (3,2,'vacancy_vector'),(3,1,'vacancy_zero')]),
                             (2048,[(3,2,'vacancy_vector_finer'),(3,3,'early_pole_alone_finer')])]:
            exe=directory/('pole'+str(count)); cpp=exe.with_suffix('.cpp')
            cpp.write_text(with_pole(src,change,count))
            subprocess.run(['g++','-std=c++17','-O3',str(cpp),'-o',str(exe)],check=True)
            for law,mode,name in combos:
                if args.only and name not in args.only: continue
                fn=directory/name
                fn.write_text('#!/bin/sh\nexec "$(dirname "$0")/'+exe.name+'" "$@" '+str(law)+' '+repr(action)+' '+str(mode)+'\n')
                fn.chmod(0o755)
                jobs.append((name,fn,dict(p,sample=0,cfl=p['cfl']*(.5 if count>1024 else 1.))))
        def run(job):
            name,fn,pr=job; m=runner(fn,pr)
            keys=('hot_efolds','first_hot_time','doubling_hot_fraction','maximum_radius','maximum_tau',
                  'end_rest','end_decay_energy','end_electric_work','peak_hot_fraction','stop_turnaround',
                  'gas_balance_error','oscillator_balance_error','surface_balance_error','friedmann_error',
                  'lost_rest','tried_rest','birth_kick','absolute_electric_work','weighted_absolute_electric_work',
                  'thermal_drag_work','all_drag_work')
            row={k:m[k] for k in keys if k in m}
            row.update({k:v for k,v in m.items() if k.startswith('pole_')})
            row['hot_radius_factor']=math.exp(m.get('hot_efolds',0))
            if max(m[k] for k in ('gas_balance_error','oscillator_balance_error','surface_balance_error'))>1e-5:
                raise ValueError('Live conservation failed')
            if abs(m['end_rest']+m['lost_rest']-m['tried_rest'])>1e-10:
                raise ValueError('Live source accounting failed')
            if not m.get('stop_turnaround'):
                raise ValueError('Expected turnaround')
            return name,row
        with ThreadPoolExecutor(args.jobs) as pool:
            results=dict(pool.map(run,jobs))
    args.output.write_text(json.dumps(dict(trial_action=action,runs=results),indent=2,sort_keys=True)+'\n')


if __name__=='__main__':
    main()
