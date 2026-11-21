"""Evolve parent depletion, the daughter current, and gravity with a local angular source."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import runpy
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
LAWS=("unlimited","angle_average","directional_cap","joint_vacancy")

def change(src,old,new):
    if src.count(old)!=1:
        raise ValueError("Unexpected parent solver near "+old[:70])
    return src.replace(old,new)

def with_feedback(src, header):
    src=change(src,'if (argc != 30)', 'if (argc != 32)')
    src=change(src,'struct Oscillator {',header+'\nint feedback_law=0;\ndouble feedback_action=0.;\nstruct Oscillator {')
    src=change(src,'COMOVING_PARENT_ACTION,','COMOVING_PARENT_ACTION, TRIED_REST, LOST_REST, BIRTH_KICK, ABS_BIRTH_KICK,')
    src=change(src,'        const auto chi=oscillator(y);\n        dy[COMOVING_PARENT_ACTION]',
    '''        auto chi=oscillator(y);
        const double trial=chi.power;
        BirthMoments angular;
        double injection_momentum=0.;
        if(envelope_active && trial>0. && chi.mass_ratio>1.) {
            const double R=chi.mass_ratio;
            const double s=std::sqrt((R-1.)*(R+1.));
            const double c=h*s+R*(chi.mass_slope*pi/(2.*daughter_mass))/s;
            const double d=charge_mass*w.x;
            // Divide both normal speeds and requested flux per area by a*mu.
            const double target=(std::pow(std::acos(-1.),2)/6.)*
                (feedback_action/std::pow(daughter_mass,4))*trial/(R*s*s);
            angular=birth_moments(feedback_law,c,d,target);
            chi.width*=angular.accept;
            chi.power*=angular.accept;
            injection_momentum=a3*trial/R*s*angular.first;
        }
        dy[TRIED_REST]=a3*trial/chi.mass_ratio;
        dy[LOST_REST]=dy[TRIED_REST]*(1.-angular.accept);
        dy[BIRTH_KICK]=injection_momentum;
        dy[ABS_BIRTH_KICK]=std::abs(injection_momentum);
        dy[COMOVING_PARENT_ACTION]''')
    src=change(src,'dy[DIRECTED]=charge_mass*gas.comoving_rest*w.x-(h+drag_rate)*moment;',
               'dy[DIRECTED]=charge_mass*gas.comoving_rest*w.x-(h+drag_rate)*moment+injection_momentum;')
    src=change(src,'+(chi.mass_ratio-inv_gamma_c)*dy[DAUGHTER_REST];',
              '+(chi.mass_ratio-inv_gamma_c)*dy[DAUGHTER_REST]-vc*injection_momentum;')
    src=change(src,'dy[REST_WORK_MIX]=(gas.gamma-1.)*dy[DECAY_ENERGY];',
              'dy[REST_WORK_MIX]=(gas.gamma-1.)*dy[DECAY_ENERGY]-gas.u*injection_momentum;')
    src=change(src,'        double olda = y[A_W], oldh', '        const double old_wall_mass=y[MASS];\n        double olda = y[A_W], oldh')
    src=change(src,'        if (oldh > 0 && h <= 0 && t > .001)','''        if(oldh>0 && h<=0 && t>.001 && !out.count("turn_wall_mass")) {
            const double z=oldh/(oldh-h);
            out["turn_wall_mass"]=old_wall_mass+z*(y[MASS]-old_wall_mass);
            out["turn_wall_e"]=previous_wall.e+z*(wnew.e-previous_wall.e);
        }
        if (oldh > 0 && h <= 0 && t > .001)''')
    src=change(src,'        const double margin=gas.intrinsic_random-.9*wnew.e;','''        if(!out.count("maximum_tau")) {
            const double stored=std::pow(a,4)*wnew.e;
            if(stored>out["max_stored_a4e"]) {
                out["max_stored_a4e"]=stored;
                out["max_stored_time"]=t+dt;
            }
            out["max_curvature_budget"]=std::max(out["max_curvature_budget"],a*a*(1.+std::pow(a*h,2)));
            out["min_wall_mass_to_turn"]=out.count("min_wall_mass_to_turn") ?
                std::min(out["min_wall_mass_to_turn"],y[MASS]) : y[MASS];
        }
        const double margin=gas.intrinsic_random-.9*wnew.e;''')
    src=change(src,'    out["end_tau"] = t;', '''    out["tried_rest"]=y[TRIED_REST];
    out["lost_rest"]=y[LOST_REST];
    out["birth_kick"]=y[BIRTH_KICK];
    out["absolute_birth_kick"]=y[ABS_BIRTH_KICK];
    out["end_tau"] = t;''')
    src=change(src,'width_scale=std::stod(argv[27]);','feedback_law=FEEDBACK_LAW; feedback_action=FEEDBACK_ACTION;\n    width_scale=std::stod(argv[27]);')
    return src

FIELDS="""tried_rest lost_rest birth_kick absolute_birth_kick end_rest end_decay_energy end_osc_action
end_electric_work end_birth_kinetic end_pair_count first_hot_time first_hot_radius first_hot_rest
hot_efolds doubling_hot_fraction peak_hot_fraction maximum_radius maximum_tau stop_turnaround
end_tau opening_time gas_balance_error oscillator_balance_error friedmann_error sheet_balance_error
surface_balance_error pair_number_error low_drift_drag_floor turn_wall_e turn_wall_mass max_stored_a4e max_stored_time max_curvature_budget min_wall_mass_to_turn doubling_rest doubling_radius doubling_theta doubling_aH""".split()

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base',type=Path,default=HERE.parent/'prior')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--jobs',type=int,default=2)
    args=parser.parse_args()
    data=json.loads((args.base/'envelope/ENVELOPE.json').read_text())
    native=(args.base/'coupled/joined_forward.cpp').read_text()
    transform=runpy.run_path(str(args.base/'envelope/build.py'))['transform']
    run=runpy.run_path(str(args.base/'coupled/joined.py'))['run']
    header=(HERE/'angular.hpp').read_text()
    jobs=[]
    generated={}
    with tempfile.TemporaryDirectory(prefix="hdblast-directional-flow-") as folder:
        folder=Path(folder)
        for case in ('mass_4e10','mass_4e12'):
            high=case=='mass_4e12'
            src=transform(native,centered=high)
            expected=data['centered_solver_sha256' if high else 'transformed_solver_sha256']
            if hashlib.sha256(src.encode()).hexdigest()!=expected:raise ValueError("Source mismatch")
            scale=data['screen'][case]
            action=5*scale['elastic_drag_action_floor'] if high else scale['witness_action']
            generated[case]=with_feedback(src,header)
            cpp=folder/(case+'.cpp')
            # Native runner has a stable CLI. Configure the local source through two extra arguments.
            source=generated[case].replace('feedback_law=FEEDBACK_LAW; feedback_action=FEEDBACK_ACTION;',
                      'feedback_law=std::stoi(argv[30]); feedback_action=std::stod(argv[31]);')
            cpp.write_text(source)
            exe=cpp.with_suffix('')
            subprocess.run(['g++','-O3','-std=c++17',str(cpp),'-o',str(exe)],check=True)
            for law,name in enumerate(LAWS):
                wrapper=folder/(case+'_'+name)
                wrapper.write_text('#!/bin/sh\nexec "$(dirname "$0")/'+case+'" "$@" '+str(law)+' '+repr(action)+'\n')
                wrapper.chmod(0o755)
                for resolution in (('native','half') if law>=2 else ('native',)):
                    p=dict(data['parameters'][case],sample=0)
                    if resolution=='half':p['cfl']*=.5
                    jobs.append((case,name,resolution,wrapper,p))
        def compute(job):
            case,name,resolution,exe,p=job
            r=run(exe,p)
            row={k:r[k] for k in FIELDS if k in r}
            if abs(row['tried_rest']-row['lost_rest']-row['end_rest'])>1e-10:
                raise ValueError('Trial and accepted rest accounting failed')
            if max(row[k] for k in ['gas_balance_error','oscillator_balance_error','surface_balance_error'])>1e-5:
                raise ValueError('Live energy ledger failed')
            row['refused_fraction_of_attempts']=row['lost_rest']/row['tried_rest']
            row['hot_radius_factor']=math.exp(row.get('hot_efolds',0.))
            return case,name,resolution,row
        results={}
        with ThreadPoolExecutor(args.jobs) as pool:
            for case,name,resolution,row in pool.map(compute,jobs):
                results.setdefault(case,{}).setdefault(name,{})[resolution]=row
    for case,laws in results.items():
        ref=data['runs'][case]
        if abs(laws['unlimited']['native']['hot_efolds']-ref['hot_efolds'])>2e-7:
            raise ValueError('Unlimited source differs from archived reference')
    args.output.write_text(json.dumps(results,sort_keys=True,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':
    main()
