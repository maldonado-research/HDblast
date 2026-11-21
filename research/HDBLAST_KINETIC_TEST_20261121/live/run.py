"""Re-evolve the full finite collar using a momentum-particle current and stress."""
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


def make_source(source,driven):
    fb=runpy.run_path(str(source/'feedback/feedback.py'))
    change=fb['change']
    base=source/'prior'
    original=(base/'coupled/joined_forward.cpp').read_text()
    src=runpy.run_path(str(base/'envelope/build.py'))['transform'](original,centered=True)
    src=runpy.run_path(str(source/'boost/evolve.py'))['make_source'](src,fb,(source/'feedback/angular.hpp').read_text())
    ensemble=(driven/'ensemble.cpp').read_text()
    start=ensemble.index('constexpr double pi=')
    end=ensemble.index('int main(int argc')
    namespace=ensemble[:start]+'\nnamespace kin {\n'+ensemble[start:end]+'\n}\n'
    src=namespace+src
    src=change(src,'struct Gas {',(HERE/'coupling.hpp').read_text()+'\nstruct Gas {')
    src=change(src,'Gas charged_gas(const vector<double> &y) {','Gas charged_gas(const vector<double> &y,bool detailed=false) {')
    src=change(src,'const auto previous_gas=charged_gas(y);','const auto previous_gas=charged_gas(y,true);')
    marker='        const auto gas=charged_gas(y);'
    ep=src.index('            out["prepare_radius"]=a;')
    idx=src.index(marker,ep)
    src=src[:idx]+src[idx:].replace(marker,'        const auto gas=charged_gas(y,true);',1)
    src=change(src,'    const auto osc=oscillator(y);\n    const auto gas=charged_gas(y);','    const auto osc=oscillator(y);\n    const auto gas=charged_gas(y,true);')
    src=change(src,'    g.gamma=1.;','''    g.gamma=1.;
    if(cloud) {
        g.dilution=1./std::pow(y[A_W],3);
        const auto kin=stress(y[A_W],y[DIRECTED],detailed);
        g.rest=g.comoving_rest*g.dilution;
        g.current=kin.current; g.velocity=kin.average_v; g.gas_momentum=0.;
        g.energy=(y[LAB_EXCESS]+g.comoving_rest)*g.dilution;
        g.pressure=kin.pressure;g.intrinsic_random=kin.intrinsic;return g;
    }''')
    src=change(src,'if (argc != 34)','if (argc != 38)')
    src=change(src,'collision_mode=std::stoi(argv[33]);','''collision_mode=std::stoi(argv[33]);
    kinetic_particles=std::stoi(argv[34]);
    if(kinetic_particles>0) {
        const double scale=std::stod(argv[37])/kinetic_particles;
        const double branch=std::stod(argv[36]);
        cloud=new kin::Wall{scale,feedback_action/(6.*std::pow(daughter_mass,3)),branch,branch,true,{},kin::RNG(std::stoul(argv[35]))};
        cloud->ps.reserve(kinetic_particles+5);
    }''')
    src=change(src,'if (collision_mode == 0 || gas.rest <= 0.) return 0.;',
               'if (cloud || collision_mode == 0 || gas.rest <= 0.) return 0.;')
    src=change(src,'        if(charge_mass>0. && previous_gas.comoving_rest>0.) {',
               '        if(!cloud && charge_mass>0. && previous_gas.comoving_rest>0.) {')
    src=change(src,'        dy[GAS_PRESS_WORK]=3.*a3*h*gas.pressure;', '''        dy[GAS_PRESS_WORK]=3.*a3*h*gas.pressure;
        if(cloud) {
            // Re-purpose DIRECTED as the signed comoving canonical impulse and
            // LAB_EXCESS as E_gas,comoving minus the accepted daughter rest.
            dy[DIRECTED]=a*charge_mass*w.x;
            dy[LAB_EXCESS]=dy[DECAY_ENERGY]-dy[DAUGHTER_REST]+dy[ELECTRIC_WORK]-dy[GAS_PRESS_WORK];
        }''')
    src=change(src,'        const double old_wall_mass=y[MASS];',
        '        const double old_impulse=y[DIRECTED];\n'
        '        const auto old_stress=cloud ? stress(y[A_W],old_impulse,true) : ParticleStress{};\n'
        '        const double old_wall_mass=y[MASS];')
    src=change(src,'        steps++;','''        if(cloud) {
            const double a=y[A_W], I=y[DIRECTED], D=y[DAUGHTER_REST];
            const auto os=oscillator(y);
            const auto local=wall(y,y[s.index(0,I /* placeholder */)]);
            const double z=std::sqrt(std::max(0.,(os.mass_ratio-1.)*(os.mass_ratio+1.)));
            const double cm=z>0 ? y[H]*z+os.mass_ratio*os.mass_slope*y[PI]/(2.*daughter_mass*z) : 0.;
            const double dm=charge_mass*local.x;
            const double target=z>0 ? std::pow(std::acos(-1.),2)/6.*feedback_action/std::pow(daughter_mass,4)*os.power/(os.mass_ratio*z*z) : 0.;
            while(D>=(cloud->ps.size()+.5)*cloud->w) {
                const double q=(cloud->ps.size()+.5)*cloud->w;
                const double frac=std::clamp((q-old_sources[0])/(D-old_sources[0]),0.,1.);
                const double ib=old_impulse+frac*(I-old_impulse), ab=olda+frac*(a-olda);
                cloud->birth({t+dt*frac,ab,z,cm,dm,target},ab,ib);
                ++generation;
            }
            const auto before=cloud->hits;
            cloud->I=.5*(old_impulse+I);
            cloud->collide(.5*(olda+a),dt);
            if(cloud->hits!=before) ++generation;
            cloud->I=I;
        }
        steps++;'''.replace('s.index(0,I /* placeholder */)','s.index(0,::I)'))
    src=change(src,'        const double margin=gas.intrinsic_random-.9*wnew.e;', '''        if(cloud && !cloud->ps.empty()) {
            const auto part=stress(a,y[DIRECTED],true);
            auto record=[&](Interval& gate,double f0,double f1,bool use_angle) {
                gate.update(dt,olda,a,oldh,h,oldphi-fall_phi-fall_width,y[PHI]-fall_phi-fall_width,
                            f0,f1,use_angle ? old_stress.longitudinal_share : 1./3.,
                            use_angle ? part.longitudinal_share : 1./3.);
            };
            record(landau_angular,previous_gas.intrinsic_random/previous_wall.e,gas.intrinsic_random/wnew.e,true);
            record(lab_hot,old_stress.lab_dispersion/previous_wall.e,part.lab_dispersion/wnew.e,false);
            record(lab_angular,old_stress.lab_dispersion/previous_wall.e,part.lab_dispersion/wnew.e,true);
            if(kinetic_particles>=2048) {
                record(qualifying_landau,previous_gas.intrinsic_random/previous_wall.e,gas.intrinsic_random/wnew.e,false);
                record(post_release,1.,1.,false);
            }
            if(h>0. && y[PHI]>=fall_phi+fall_width)
                out["peak_lab_dispersion"]=std::max(out["peak_lab_dispersion"],part.lab_dispersion/wnew.e);
            out["max_dt_omega_kinetic"]=std::max(out["max_dt_omega_kinetic"],dt*charge_mass*std::sqrt(part.susceptibility/wnew.k));
            out["max_dt_omega_cold_bound"]=std::max(out["max_dt_omega_cold_bound"],
                        dt*charge_mass*std::sqrt(cloud->w*cloud->ps.size()/(a*a*a*wnew.k)));
            if(oldh>0. && h<=0. && t>.001 && !out.count("turn_shape_ks")) snapshot(out,"turn",a,y[DIRECTED],wnew.e);
            const double mismatch=std::abs(part.energy-gas.energy)*a*a*a;
            out["max_particle_energy_error"]=std::max(out["max_particle_energy_error"],mismatch);
            out["max_comoving_particle_energy"]=std::max(out["max_comoving_particle_energy"],part.energy*a*a*a);
            if(h>0. && y[PHI]>=fall_phi+fall_width) {
                out["max_eckart_hot"]=std::max(out["max_eckart_hot"],part.eckart/wnew.e);
            }
        }
        const double margin=gas.intrinsic_random-.9*wnew.e;''')
    for k in ('first_hot','doubling'):
        anchor='store_budget("'+k+'",f);'
        src=change(src,anchor,anchor+' if(cloud) snapshot(out,"'+k+'",a,y[DIRECTED],wnew.e);')
    src=change(src,'    out["end_tau"] = t;', '''    if(cloud) {
        out["particle_count"]=cloud->ps.size();out["scatter_events"]=cloud->hits;
        out["max_event_conservation_error"]=cloud->event_residual;
        out["particle_energy_error_over_peak"]=out["max_particle_energy_error"]/out["max_comoving_particle_energy"];
        out["end_sampled_rest"]=cloud->w*cloud->ps.size();
        out["integrated_birth_energy"]=y[DECAY_ENERGY];
        out["particle_birth_energy"]=cloud->injection;
        out["particle_birth_energy_relative_error"]=std::abs(cloud->injection/y[DECAY_ENERGY]-1.);
        landau_angular.write(out,"joint_hot_angular");lab_hot.write(out,"lab_hot");lab_angular.write(out,"lab_hot_angular");
        if(kinetic_particles>=2048) {qualifying_landau.write(out,"landau_index");post_release.write(out,"post_birth_expanding");}
    }
    out["end_tau"] = t;''')
    return src


def main():
    pa=argparse.ArgumentParser(description=__doc__)
    pa.add_argument('--prior',type=Path,required=True)
    pa.add_argument('--driven',type=Path,default=HERE.parent/'driven')
    pa.add_argument('--output',type=Path,required=True)
    pa.add_argument('--jobs',type=int,default=2)
    pa.add_argument('--pilot',action='store_true')
    pa.add_argument('--extended',action='store_true',help='Run the dedicated 2048-particle resolution pair')
    pa.add_argument('--only',nargs='+')
    args=pa.parse_args()
    source=args.prior.resolve()
    ref=json.loads((source/'boost/EVOLUTION.json').read_text())
    env=json.loads((source/'prior/envelope/ENVELOPE.json').read_text())
    run=runpy.run_path(str(source/'prior/coupled/joined.py'))['run']
    src=make_source(source,args.driven)
    controls=[('control',0,1.,.08,1717),('no_scattering_512',512,0.,.08,1717)]
    dedicated=[(256,.08,1717),(256,.08,7127),(512,.08,1717),(512,.08,7127),
               (512,.04,1717),(1024,.08,1717)]
    configs=controls+[(f'dedicated_{n}_{seed}_{cfl}',n,1.,cfl,seed) for n,cfl,seed in dedicated]
    if args.pilot:
        configs=[controls[0],('dedicated128',128,1.,.08,1717)]
    if args.extended:
        configs=[(f'dedicated_2048_{seed}_0.08',2048,1.,.08,seed) for seed in (1717,7127)]
    if args.only:
        targets=set(args.only)
        configs=[q for q in configs if q[0] in targets]
        if set(q[0] for q in configs)!=targets:raise ValueError('Unknown requested live case')
    keys='steps maximum_radius maximum_tau min_lapse min_areal_radius min_bulk_mass friedmann_error surface_balance_error gas_balance_error oscillator_balance_error end_tau end_rest end_electric_work end_gas_density end_gas_random end_gas_pressure_work first_hot_time first_hot_radius hot_intervals peak_hot_time peak_hot_fraction peak_hot_radius doubling_time doubling_hot_fraction stop_turnaround particle_count scatter_events max_event_conservation_error max_particle_energy_error max_comoving_particle_energy particle_energy_error_over_peak end_sampled_rest max_eckart_hot'.split()
    with tempfile.TemporaryDirectory(prefix='hdblast-live-cloud-') as folder:
        tmp=Path(folder); cpp=tmp/'live.cpp';ex=tmp/'live'; cpp.write_text(src)
        subprocess.run(['g++','-O3','-std=c++17','-I',str(source/'boost'),'-I',str(args.driven.resolve()),str(cpp),'-o',str(ex)],check=True)
        def one(cfg):
            name,n,br,cfl,seed=cfg
            past=ref['runs']['pole_cap']
            opts=['2',str(ref['action']),str(source/'boost/BOOST_TABLE.dat'),'1',str(n),str(seed),str(br),str(past['end_rest'])]
            w=tmp/name; w.write_text('#!/bin/sh\nexec '+shlex.quote(str(ex))+' "$@" '+shlex.join(opts)+'\n');w.chmod(0o755)
            try: raw=run(w,dict(env['parameters']['mass_4e12'],cfl=cfl,sample=0))
            except subprocess.CalledProcessError as e:raise RuntimeError(f'Failed {name} {e.stderr}') from e
            row={k:raw[k] for k in keys if k in raw}
            row.update({k:v for k,v in raw.items() if k.startswith(('joint_','lab_hot','peak_lab_dispersion','max_dt_omega','integrated_birth_energy','particle_birth_energy','landau_index_','post_birth_expanding_')) or k.endswith(('_lab_dispersion_fraction','_shape_ks','_angular_trace_share','_fitted_temperature_over_m','_empirical_tube_share','_shape_pressure_check'))})
            if not row['stop_turnaround'] or row['min_bulk_mass']<0 or row['friedmann_error']>3e-7 or row['surface_balance_error']>3e-6:
                raise ValueError('Full collar or conservation check failed: '+name)
            if n and (row['particle_energy_error_over_peak']>.005 or row['max_event_conservation_error']>1e-10):
                raise ValueError('Direct particle energy consistency check failed: '+name)
            if not n and not all(math.isclose(row[k],past[k],rel_tol=2e-5,abs_tol=1e-8)
                                 for k in ('peak_hot_fraction','maximum_radius','maximum_tau','first_hot_time')):
                raise ValueError('Thermal control differs from its reference')
            row.update(hot_radius_growth=math.exp(raw.get('hot_efolds',0.)),particles=n,seed=seed,branch=br,cfl=cfl)
            print(name,row,flush=True)
            return name,row
        with ThreadPoolExecutor(args.jobs) as pool:
            rows=dict(pool.map(one,configs))
    args.output.write_text(json.dumps(rows,indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
