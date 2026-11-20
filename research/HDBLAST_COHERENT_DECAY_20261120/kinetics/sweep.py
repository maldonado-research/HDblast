"""Sharp-shell semiclassical fermion-state sweep on the large-mass averaged preparations."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
SOURCE_STEP=r"""
        if(gas.comoving_rest > TARGET_REST*.01 && previous_gas.comoving_rest < TARGET_REST*.99) {
            const double m=daughter_mass;
            const auto kbirth=[&](double rad, double phi) {
                const double R=mass_low+(mass_high-mass_low)*step((phi-threshold_phi)/field_width+.5)[0];
                return m*rad*std::sqrt(std::max(0.,(R-1.)*(R+1.)));
            };
            const double fr0=previous_gas.comoving_rest/TARGET_REST,fr1=gas.comoving_rest/TARGET_REST;
            const double u0=std::clamp((.01-fr0)/(fr1-fr0),0.,1.),u1=std::clamp((.99-fr0)/(fr1-fr0),0.,1.);
            if(u1>u0 && fr1>fr0) {
                const double k0=kbirth(olda,oldphi),k1=kbirth(a,y[PHI]),dk=std::abs((k1-k0)/dt);
                const double k=k0+.5*(u0+u1)*(k1-k0);
                const double x0=m*charge_mass*olda*previous_wall.x,x1=m*charge_mass*a*wnew.x;
                const double rest=(fr1-fr0)*TARGET_REST;
                const double injection=2.*std::pow(std::acos(-1.),2)*rest/(dt*m);
                const double subdt=dt*(u1-u0),lo=x0+u0*(x1-x0),hi=x0+u1*(x1-x0);
                main_sweep.collect(t+u0*dt,subdt,lo,hi,dk,k,ACTION*injection,rest*(u1-u0),fr0+u0*(fr1-fr0),fr0+u1*(fr1-fr0));
                floor_sweep.collect(t+u0*dt,subdt,lo,hi,dk,k,FLOOR*injection,rest*(u1-u0),fr0+u0*(fr1-fr0),fr0+u1*(fr1-fr0));
            }
        }
"""

def make_source(source,rest,witness,floor):
    def change(s,a,b):
        if s.count(a)!=1:raise ValueError('Missing or ambiguous diagnostic anchor '+a)
        return s.replace(a,b)
    s=change(source,'int main(int argc, char **argv) {',(HERE/'sweep.hpp').read_text()+'\nint main(int argc, char **argv) {\n    check_sweep_integral(); SweepCollector main_sweep,floor_sweep,extended_main_sweep,extended_floor_sweep;')
    anchor='        const double margin=gas.intrinsic_random-.9*wnew.e;'
    part=SOURCE_STEP.replace('TARGET_REST',repr(rest)).replace('ACTION',repr(witness)).replace('FLOOR',repr(floor))
    extended=part.replace('*.01','*.00001').replace('*.99','*.99999').replace('(.01-fr0)','(.00001-fr0)').replace('(.99-fr0)','(.99999-fr0)')
    extended=extended.replace('main_sweep.collect','extended_main_sweep.collect').replace('floor_sweep.collect','extended_floor_sweep.collect')
    s=change(s,anchor,part+extended+anchor)
    anchor='    out["end_tau"] = t;'
    saves=''.join('\n    '+key+'.save(out,"'+key+'_",'+repr(rest)+');'
                  for key in ('main_sweep','floor_sweep','extended_main_sweep','extended_floor_sweep'))
    return change(s,anchor,anchor+saves)


def run_one(source,run,params,reference,screen,witness_name):
    rest=reference['end_rest'];floor=screen['elastic_drag_action_floor']
    action=screen['witness_action'] if witness_name=='screen' else 5*floor
    with tempfile.TemporaryDirectory(prefix='hdblast-sharp-') as tmp:
        cpp=Path(tmp)/'measure.cpp';exe=Path(tmp)/'measure'
        cpp.write_text(make_source(source,rest,action,floor))
        subprocess.run(['g++','-O3','-std=c++17',str(cpp),'-o',str(exe)],check=True)
        rows={}
        for label,scale in [('native',1.),('half_step',.5)]:
            z=run(exe,dict(params,sample=0,cfl=params['cfl']*scale))
            row={k:v for k,v in z.items() if 'sweep_' in k}
            row['end_rest']=z['end_rest'];row['hot_factor']=math.exp(z['hot_efolds'])
            if not math.isclose(rest,z['end_rest'],rel_tol=3e-8 if scale==1 else .0003):
                raise ValueError('Reference trajectory mismatch')
            if abs(row['main_sweep_fraction_tested']-.98)>1e-8:raise ValueError('Incomplete delivery interval')
            row['continuum_pointwise_S0_upper']=action/row['main_sweep_peak']
            row['min_fraction_of_drag_unitarity_ceiling_at_pointwise_upper']=floor/row['continuum_pointwise_S0_upper']
            if abs(row['extended_main_sweep_fraction_tested']-.99998)>1e-8:raise ValueError('Incomplete extended delivery interval')
            rows[label]=row
    return dict(fraction_windows=dict(main_sweep=[.01,.99],floor_sweep=[.01,.99],extended_main_sweep=[.00001,.99999],extended_floor_sweep=[.00001,.99999]),mu=params['daughter_mass'],charge_to_mass=params['charge_to_mass'],reference_drag_floor=floor,
                action=action,action_choice=witness_name,runs=rows)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--envelope',type=Path,default=HERE.parent/'envelope/ENVELOPE.json')
    parser.add_argument('--cases',nargs='+',default=['mass_4e10','mass_4e12'])
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    env=json.loads(args.envelope.read_text())
    coupled=HERE.parent/'coupled'
    run=runpy.run_path(str(coupled/'joined.py'))['run']
    transform=runpy.run_path(str(HERE.parent/'envelope/build.py'))['transform']
    native_source=(coupled/'joined_forward.cpp').read_text()
    hashes={}
    results={}
    for n in args.cases:
        centered=env['parameters'][n]['daughter_mass']>=1e11
        source=transform(native_source,centered=centered)
        expected=env['centered_solver_sha256' if centered else 'transformed_solver_sha256']
        if hashlib.sha256(source.encode()).hexdigest()!=expected:
            raise ValueError('Averaged solver source mismatch')
        hashes[n]=expected
        results[n]=run_one(source,run,env['parameters'][n],env['runs'][n],env['screen'][n],
                           'screen' if n=='mass_4e10' else 'five_times_drag_floor')
    args.output.write_text(json.dumps({'per_case_solver_sha256':hashes,'runs':results},
                                     sort_keys=True,indent=2,allow_nan=False)+'\n')
