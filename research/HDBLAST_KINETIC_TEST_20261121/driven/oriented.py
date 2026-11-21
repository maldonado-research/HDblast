"""Resolve the fixed triad and magnetic azimuthal motion in the driven particle replay."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import runpy
import subprocess
import tempfile

import numpy as np

HERE=Path(__file__).resolve().parent
FIELDS_ADDED=("fourth_cosine","fourth_sine","max_transverse_norm_residual","max_phase_step",
              "rotations_over_one_radian_fraction")


def transform(src):
    def change(old,new):
        nonlocal src
        if src.count(old)!=1:
            raise ValueError("Unexpected driven source near "+old[:110])
        src=src.replace(old,new)
    change("struct Particle {double qperp,qz; int hits;};",
           'struct Particle {double qperp,qz; int hits; double qx,qy;};\n#include "oriented.hpp"')
    change("    long long hits=0, proposed=0;", """    long long hits=0, proposed=0, phase_steps=0, large_steps=0;
    double max_phase=0.;
    bool magnetic=true, event_phase=true;""")
    change("        for(const auto& p:ps) ledger(p,left.a,right.a,I,I+dI);","""        const double magnetic_column=Q*h*(mid.Bref+2.*z_mid/(mid.a*mid.a));
        for(auto& p:ps) {
            ledger(p,left.a,right.a,I,I+dI);
            if(magnetic) {
                const double phase=magnetic_column*average_inverse_energy(p.qperp,p.qz+I,p.qz+I+dI,mid.a);
                rotate_gyro(p,phase);
                max_phase=std::max(max_phase,std::abs(phase));
                ++phase_steps;
                if(std::abs(phase)>1.) ++large_steps;
            }
        }""")
    change("        ps.push_back({a*b.pstar*std::sqrt(1.-u*u),a*b.pstar*u-impulse,0});",
           """        const double pt=b.pstar*std::sqrt(1.-u*u), az=2*pi*random.unit();
        ps.push_back(canonical_particle({pt*std::cos(az),pt*std::sin(az),b.pstar*u},a,impulse,0));""")
    change("    void collide(double a,double dt) {","    void collide(double a,double dt,double magnetic_field) {")
    change("wall.collide(mid.a,rows[i+4].t-rows[i].t);", "wall.collide(mid.a,rows[i+4].t-rows[i].t,mid.Bref+2.*wall.Z/(mid.a*mid.a));")
    change("            const auto pa=ps[i],pb=ps[j];","""            auto pa=ps[i],pb=ps[j];
            // Stored azimuths refer to the midpoint. Rephase each candidate to its actual
            // event time under the frozen collision operator, then undo it for the products.
            const double orbit_offset=event_phase ? Q*magnetic_field*(now-.5*dt) : 0.;
            if(event_phase) {
                rotate_gyro(pa,orbit_offset/mid_step_energy(pa,a,I));
                rotate_gyro(pb,orbit_offset/mid_step_energy(pb,a,I));
            }
""")
    start="            const double angle=2.*pi*random.unit();"
    end="            const double e=std::sqrt(1.+dot(p,p)), f=std::sqrt(1.+dot(q,q));"
    old=src[src.index(start):src.index(end)]
    change(old,"""            const Vector p=quarter_turn({pa.qx/a,pa.qy/a,(pa.qz+I)/a},random.index(4));
            const Vector qp=quarter_turn({pb.qx/a,pb.qy/a,(pb.qz+I)/a},random.index(4));
            const Vector q=charge_partner(qp);
""")
    start="            Vector dest={0,0,1};"
    end="            ++hits;"
    change(src[src.index(start):src.index(end)],"""            const int dest=weight<1. ? random.index(3) : 2;
            const Vector positive=final_basis(pnew,dest);
            const Vector negative=final_basis(qnew,dest);
            ps[i]=canonical_particle(positive,a,I,pa.hits+1);
            ps[j]=canonical_particle(charge_partner(negative),a,I,pb.hits+1);
            if(event_phase) {
                rotate_gyro(ps[i],-orbit_offset/mid_step_energy(ps[i],a,I));
                rotate_gyro(ps[j],-orbit_offset/mid_step_energy(ps[j],a,I));
            }
""")
    change("        std::vector<distribution::Observation> list;",
           """        double c4=0.,s4=0.,rmax=0.;
        for(const auto& particle:ps) {
            const double x=particle.qx,y=particle.qy,z=particle.qz+I, norm2=x*x+y*y+z*z;
            if(norm2>0.) {
                c4+=(x*x*x*x-6.*x*x*y*y+y*y*y*y)/(norm2*norm2);
                s4+=4.*x*y*(x*x-y*y)/(norm2*norm2);
            }
            rmax=std::max(rmax,std::abs(std::hypot(x,y)-particle.qperp)/std::max(1.,particle.qperp));
        }
        // The four rotations about z preserve both fourth harmonics.
        std::vector<distribution::Observation> list;""")
    change("<<(phase_time>0.?phase_pressure/phase_time:pz/tr)<<'\\n';",
           "<<(phase_time>0.?phase_pressure/phase_time:pz/tr)<<' '<<c4/n<<' '<<s4/n<<' '<<rmax<<' '<<max_phase<<' '<<double(large_steps)/std::max(1LL,phase_steps)<<'\\n';")
    change("if(argc!=6)", "if(argc!=7)")
    change("    wall.I=rows[0].Iref;", "    wall.magnetic=std::stoi(argv[6])!=0; wall.event_phase=std::stoi(argv[6])==1;\n    wall.I=rows[0].Iref;")
    return src


def run(args):
    baseline=runpy.run_path(str(HERE/"run.py"))
    fields=baseline["FIELDS"]+list(FIELDS_ADDED)
    keys=json.loads((args.native/"metadata.json").read_text())
    gyro=json.loads(args.gyro.read_text()) if args.gyro else {"groups":{}}
    seeds=baseline["SEEDS"][:1] if args.pilot else baseline["SEEDS"][:2]
    cfgs=[("baseline",2048,"native",.001),("finer_time",2048,"native",.0005),("more_particles",4096,"half",.0005)]
    if args.pilot: cfgs=[("pilot",1024,"native",.0015)]
    with tempfile.TemporaryDirectory(prefix="hdblast-oriented-") as tmp:
        dest=Path(tmp);source=dest/"oriented.cpp";exe=dest/"oriented"
        source.write_text(transform((HERE/"ensemble.cpp").read_text()))
        subprocess.run(["g++","-O3","-std=c++17","-I",str(HERE),str(source),"-o",str(exe)],check=True)
        tasks=[];markers={}
        for cfg,n,mesh,step in cfgs:
            for channel in ("dedicated","shared"):
                trace=dest/(cfg+"_"+channel)
                markers[channel+"_"+mesh]=baseline["prepare"](args.native,keys,channel,n,mesh,trace,step)
                modes=("magnetic","split_gyro","no_gyromotion") if cfg=="baseline" or args.pilot else ("magnetic",)
                tasks.extend((cfg,channel,mode,seed,trace) for mode in modes for seed in seeds)
        def evaluate(task):
            cfg,channel,mode,seed,trace=task
            cmd=[str(exe),str(trace),str(seed),"1","1",str(1. if channel=="dedicated" else 1./3),
                 str({"magnetic":1,"split_gyro":2,"no_gyromotion":0}[mode])]
            raw=subprocess.run(cmd,check=True,capture_output=True,text=True).stdout
            arr=np.fromstring(raw,sep=" ").reshape(len(baseline["MARKERS"]),len(fields)+1)
            if not np.isfinite(arr).all():raise ValueError("Nonfinite oriented run")
            result={}
            for i,tag in enumerate(baseline["MARKERS"]):
                r=dict(zip(fields,map(float,arr[i,1:])))
                r["energy_residual"]=abs(r["gas_comoving"]-r["injection"]-r["electric_work"]+r["redshift_loss"])/r["gas_comoving"]
                for frame in ("landau","eckart"):
                    r[frame+"_fraction"]=r[frame+"_density"]/r["replaced_total_density"]
                if max(r["energy_residual"],r["pair_error"],r["max_transverse_norm_residual"])>2e-6:
                    raise ValueError((tag,r))
                result[tag]=r
            print(cfg,channel,mode,seed,"first hot",result["reference_first_hot"]["landau_fraction"],
                  "doubling", result["radius_doubling"]["landau_fraction"],flush=True)
            return (cfg,channel,mode),result
        groups={}
        with ThreadPoolExecutor(args.jobs) as pool:
            for k,r in pool.map(evaluate,tasks):
                groups.setdefault(k,[]).append(r)
    out={}
    for (cfg,channel,mode),rs in groups.items():
        name="_".join((cfg,channel,mode));out[name]={}
        for tag in baseline["MARKERS"]:
            mean={key:float(np.mean([s[tag][key] for s in rs])) for key in rs[0][tag]}
            extremes={key:[min(s[tag][key] for s in rs),max(s[tag][key] for s in rs)]
                      for key in mean}
            rec=dict(mean=mean,seed_range=extremes)
            key=cfg+"_"+channel+"_kinetic"
            if key in gyro["groups"]:
                ref=gyro["groups"][key][tag]["mean"]
                rec["gyrotropic_landau_fraction"]=ref["landau_fraction"]
            out[name][tag]=rec
    return dict(configurations=cfgs,seeds=seeds,markers=markers,groups=out)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--native",type=Path,required=True)
    p.add_argument("--gyro",type=Path,help="Gyrotropic reference JSON for side-by-side comparison")
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--jobs",type=int,default=2)
    p.add_argument("--pilot",action="store_true")
    args=p.parse_args()
    args.output.write_text(json.dumps(run(args),indent=2,sort_keys=True,allow_nan=False)+"\n")
