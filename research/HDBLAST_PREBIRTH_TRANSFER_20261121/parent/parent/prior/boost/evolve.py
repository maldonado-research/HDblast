"""Evolve the entire closed-wall history using one fixed elastic vector pole."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import runpy
import shlex
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent


def make_source(src, feedback, header):
    change = feedback['change']
    src = feedback['with_feedback'](src, header)
    src = change(src, 'if (argc != 32)', 'if (argc != 34)')
    src = change(src, 'feedback_law=FEEDBACK_LAW; feedback_action=FEEDBACK_ACTION;',
                 'feedback_law=std::stoi(argv[30]); feedback_action=std::stod(argv[31]);\n'
                 '    load_thermal_table(argv[32]); collision_mode=std::stoi(argv[33]);')
    src = '#include "lookup.hpp"\n' + src
    anchor = 'Values wall(const vector<double> &y, double i) {'
    src = change(src, anchor, '''int collision_mode = 0;
double collision_rate(const Gas& gas) {
    if (collision_mode == 0 || gas.rest <= 0.) return 0.;
    if (collision_mode == 2) return drag_rate;
    // gas.rest = D/a^3 sums charge signs and spatial axes, in tension units.
    const double proper_density = feedback_action*gas.rest/(6.*gas.gamma*daughter_mass);
    return proper_density*boosted_rate(gas.theta,gas.u)/(daughter_mass*daughter_mass);
}
''' + anchor)
    src = change(src, '        const double moment=y[DIRECTED];',
                 '        const double nu = collision_rate(gas);\n        const double moment=y[DIRECTED];')
    for old, new in [('(h+drag_rate)*moment','(h+nu)*moment'), ('+drag_rate*moment*vc','+nu*moment*vc'),
                     ('=drag_rate*gas.u*moment','=nu*gas.u*moment')]:
        src = change(src, old, new)
    src = change(src, 'ABS_BIRTH_KICK,', 'ABS_BIRTH_KICK, ABS_ELECTRIC, FAST_ELECTRIC, INTEGRAL_NU,')
    src = change(src, '        dy[ELECTRIC_WORK_A4]', '''        dy[ABS_ELECTRIC]=std::abs(dy[ELECTRIC_WORK]);
        dy[FAST_ELECTRIC]=(gas.theta<100. || std::abs(gas.velocity)>.1) ? dy[ABS_ELECTRIC] : 0.;
        dy[INTEGRAL_NU]=nu;
        dy[ELECTRIC_WORK_A4]''')
    src = change(src, '        if(charge_mass>0. && previous_gas.comoving_rest>0.) {',
                 '''        if(collision_mode==1 && previous_gas.comoving_rest>0.) {
            const double nu=collision_rate(previous_gas);
            if(nu>0.) dt=std::min(dt,.15/nu);
        }
        if(charge_mass>0. && previous_gas.comoving_rest>0.) {''')
    src = change(src, '        const double margin=gas.intrinsic_random-.9*wnew.e;',
                 '''        if(gas.rest>0.) {
            out["max_nu"]=std::max(out["max_nu"],collision_rate(gas));
            out["max_theta"]=std::max(out["max_theta"],gas.theta);
            out["max_gamma"]=std::max(out["max_gamma"],gas.gamma);
            if(gas.theta>0. && gas.theta<1e-5) out["cold_continuation_used"]=1.;
        }
        const double margin=gas.intrinsic_random-.9*wnew.e;''')
    src = change(src, '        const double old_wall_mass=y[MASS];',
                 '        const double old_integral_nu=y[INTEGRAL_NU], old_abs_electric=y[ABS_ELECTRIC];\n'
                 '        const double old_wall_mass=y[MASS];')
    for stage in ('first_hot','doubling'):
        anchor = 'store_budget("'+stage+'",f);'
        src = change(src, anchor, anchor+'\n'
                     '                    out["'+stage+'_integral_nu"]=value(old_integral_nu,y[INTEGRAL_NU],f);\n'
                     '                    out["'+stage+'_abs_electric"]=value(old_abs_electric,y[ABS_ELECTRIC],f);')
    src = change(src, '    out["end_tau"] = t;', '''    out["abs_electric"]=y[ABS_ELECTRIC];
    out["fast_electric"]=y[FAST_ELECTRIC]; out["integral_nu"]=y[INTEGRAL_NU];
    out["end_rest_drag_work"]=y[REST_WORK_DRAG]; out["end_rest_mix_work"]=y[REST_WORK_MIX];
    out["end_tau"] = t;''')
    return src


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--base', type=Path, default=HERE.parent/'prior')
    cli.add_argument('--feedback', type=Path, default=HERE.parent/'feedback')
    cli.add_argument('--table', type=Path, default=HERE/'BOOST_TABLE.dat')
    cli.add_argument('--output', type=Path, required=True)
    cli.add_argument('--jobs', type=int, default=2)
    args = cli.parse_args()
    env = json.loads((args.base/'envelope/ENVELOPE.json').read_text())
    data = runpy.run_path(str(args.feedback/'feedback.py'))
    transform = runpy.run_path(str(args.base/'envelope/build.py'))['transform']
    run = runpy.run_path(str(args.base/'coupled/joined.py'))['run']
    src = transform((args.base/'coupled/joined_forward.cpp').read_text(), centered=True)
    src = make_source(src, data, (args.feedback/'angular.hpp').read_text())
    action = 5*env['screen']['mass_4e12']['elastic_drag_action_floor']
    fields = set(data['FIELDS']) | set('max_nu max_theta max_gamma cold_continuation_used first_hot_aH first_hot_integral_nu first_hot_abs_electric doubling_integral_nu doubling_abs_electric doubling_time opening_time peak_hot_time abs_electric fast_electric integral_nu end_gas_density end_gas_random end_gas_pressure_work end_rest_drag_work end_rest_mix_work'.split())
    jobs = [('pole_cap',2,1,1),('pole_cap_half',2,1,.5),('pole_vacancy',3,1,1),
            ('pole_vacancy_half',3,1,.5),('zero_drag',2,0,1),('prescribed_drag',2,2,1)]
    with tempfile.TemporaryDirectory(prefix='hdblast-relativistic-drag-') as temp:
        dest = Path(temp)
        (dest/'evolve.cpp').write_text(src)
        exe = dest/'evolve'
        subprocess.run(['g++','-std=c++17','-O3','-I',str(HERE),str(dest/'evolve.cpp'),'-o',str(exe)],check=True)
        def evaluate(job):
            name, law, mode, step = job
            wrapper = dest/name
            cmd = shlex.join([str(law),repr(action),str(args.table.resolve()),str(mode)])
            wrapper.write_text('#!/bin/sh\nexec '+shlex.quote(str(exe))+' "$@" '+cmd+'\n')
            wrapper.chmod(0o755)
            p = dict(env['parameters']['mass_4e12'], sample=0)
            p['cfl'] *= step
            try:
                raw = run(wrapper,p)
            except subprocess.CalledProcessError as err:
                raise RuntimeError(f'{name} failed: {err.stderr}') from err
            r = {k:raw[k] for k in sorted(fields) if k in raw}
            if not r['stop_turnaround'] or max(r[k] for k in ('surface_balance_error','gas_balance_error','oscillator_balance_error'))>1e-5:
                raise ValueError(f'{name} failed evolution checks')
            if abs(r['tried_rest']-r['lost_rest']-r['end_rest'])>1e-10:
                raise ValueError(f'{name} failed rest accounting')
            r['cool_or_fast_absolute_electrical_work_fraction']=r['fast_electric']/r['abs_electric']
            r['hot_radius_factor']=math.exp(r.get('hot_efolds',0.))
            print(name,r['peak_hot_fraction'],r['hot_radius_factor'],r.get('max_nu',0.),flush=True)
            return name,r
        with ThreadPoolExecutor(args.jobs) as pool:
            rows=dict(pool.map(evaluate,jobs))
    result={'daughter_mass_L':4e12,'pole_mass_over_daughter_mass':8300.,'width_over_mass':.30,
            'action':action,'runs':rows}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')

if __name__=='__main__':
    main()
