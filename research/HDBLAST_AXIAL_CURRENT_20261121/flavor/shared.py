"""Evolve shared vector-pole particle populations with cyclic flavor exchange."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import shlex
import subprocess
import tempfile

import numpy as np


HERE = Path(__file__).resolve().parent
WALL = HERE.parent/"live"


def once(source, original, updated):
    if source.count(original) != 1:
        raise ValueError(f"Unexpected input implementation near {original}")
    return source.replace(original, updated)


def build(study, scattering, work):
    header = (scattering/"live/collide.hpp").read_text()
    header = once(header, "class PoleCollisions {", '#include "conversion.hpp"\nclass PoleCollisions {\n    FlavorJump flavors;')
    header = once(header, "bool isotropic):", "bool isotropic, bool convert):")
    header = once(header, "rng(seed),contacts", "flavors(seed, convert),rng(seed),contacts")
    header = once(header, "second.u=radius*other;", "second.u=radius*other;\n        flavors.project(first.u,second.u,P,radius,total);")
    header = once(header, "public:\n", "public:\n    std::array<uint64_t,3> flavor_counts() const {return flavors.counts();}\n    void flavor_report(const char* name) const {flavors.save(name);}\n")
    (work/"collide.hpp").write_text(header)
    prefix, separator, _ = (study/"feedback/homogeneous.cpp").read_text().partition("int main(int argc,char** argv)")
    if not separator:
        raise ValueError("Pinned electromagnetic types changed")
    solver = (WALL/"solver.cpp").read_text()
    solver = once(solver, 'if(!trace||!output||!ledger)throw std::runtime_error("Incomplete output");', 'if(!trace||!output||!ledger)throw std::runtime_error("Incomplete output");\n    return 0;')
    solver = once(solver, "if(argc!=12)", "if(argc!=14)")
    solver = once(solver, "collisions(particles.size(),density,seed,false)", "collisions(particles.size(),density/3.,seed,false,std::stoi(argv[12])!=0)")
    solver = once(solver, "collisions.save(argv[8]);", "collisions.save(argv[8]);\n    collisions.flavor_report(argv[13]);")
    binary = work/"shared"
    (work/"shared.cpp").write_text(prefix+solver)
    subprocess.run(["g++", "-O3", "-ffp-contract=off", "-std=c++17", "-I", str(HERE), str(work/"shared.cpp"), "-o", str(binary)],check=True)
    return binary


def run(opts):
    with np.load(opts.history) as file:
        history=dict(file)
    if not np.allclose(history["K"], 2.54, atol=2e-12,rtol=0) or float(history["channel_fraction"])!=1./3.:
        raise ValueError("Expected a shared interaction with constant gauge coefficient")
    evolution=runpy.run_path(str(WALL/"evolve.py"))
    closure=runpy.run_path(str(opts.scattering/"live/evolve.py"))
    coarse=runpy.run_path(str(opts.polarization/"orbits/coarse.py"))
    coarse["thermal_cells"].__globals__["quad"]=closure["thermal_quadrature"]
    orbits=runpy.run_path(str(opts.polarization/"orbits/pairs.py"))
    tasks=[]
    for t, step in ((6.30944,512),(6.45,1024),(6.67,2048)):
        for power,n,seed in ((12,step,18431),(13,step,70163),(13,2*step,70163),(14,2*step,91283)):
            for convert in (False,True):
                tasks.append((t,power,n,seed,"particle",True,convert))
    if opts.endpoint is not None:
        t=opts.endpoint
        tasks=[(t,power,steps,seed,"particle",True,True) for power,steps,seed in
               ((12,4096,18431),(13,4096,70163),(13,8192,70163),(14,8192,91283))]
        tasks.append((t,14,8192,91283,"particle",True,False))
    if opts.single:
        tasks=[(float(opts.single[0]),int(opts.single[1]),int(opts.single[2]),int(opts.single[3]),"particle",True,bool(int(opts.single[4])))]
    if any(t[0]<history["t"][0] or t[0]>history["t"][-1] for t in tasks):
        raise ValueError("The recorded shared history must span every requested time")
    snapshots=opts.output.parent/"snapshots"
    for name in ("convert","elastic"):
        (snapshots/name).mkdir(exist_ok=True,parents=True)
    with tempfile.TemporaryDirectory(prefix="shared-vector-") as folder:
        work=Path(folder)
        binary=build(opts.polarization,opts.scattering,work)
        control_binary=work/"flavor_control"
        subprocess.run(["g++","-O3","-ffp-contract=off","-std=c++17","-I",str(work),"-I",str(HERE),
                        str(HERE/"control.cpp"),"-o",str(control_binary)],check=True)
        output=subprocess.check_output([str(control_binary)],text=True)
        values=np.fromstring(output,sep=" ")
        if len(values)!=9:
            raise ValueError("Incomplete separate flavor control")
        trials,events,hits,n0,n1,n2,vector_error,charge_error,control_chisquare=values
        mu=.8
        prob=-np.expm1(-mu)
        if (abs(events/trials-mu)>5*np.sqrt(mu/trials) or abs(hits/trials-prob)>5*np.sqrt(prob*(1-prob)/trials)
            or n0+n1+n2!=events or control_chisquare>25 or vector_error>1e-14 or charge_error!=0):
            raise ValueError(("Separate flavor control failed",values))
        control=dict(trials=int(trials),events=int(events),ever_hit=int(hits),destinations=[int(n0),int(n1),int(n2)],
                     predicted_events_per_trial=mu,predicted_hit_fraction=float(prob),flavor_chisquare=control_chisquare,
                     global_angular_relabel_error=vector_error,maximum_per_flavor_net_charge_change=charge_error)

        def one(job):
            *task,convert=job
            with tempfile.TemporaryDirectory(dir=work) as scratch:
                d=Path(scratch)
                report=d/"jump_counts"
                runner=d/"run"
                runner.write_text('#!/bin/sh\nexec '+shlex.quote(str(binary))+' "$@" '+shlex.join([str(int(convert)),str(report)])+'\n')
                runner.chmod(0o755)
                context=(history,orbits,coarse,closure,runner)
                destination="convert" if convert else "elastic"
                row=evolution["run_case"](context,task,d,snapshots/destination)
                n0,n1,n2,dp,de=np.loadtxt(report)
                if n0+n1+n2 != row["events"] or (not convert and n1+n2) or max(dp,de)>1e-11:
                    raise ValueError("Population counts or per-event invariants violated")
                # Include the summed stress of the three signed-charge populations.
                with np.load(snapshots/destination/(row["label"]+".npz")) as snap:
                    p=snap["p"]
                    states=np.concatenate((p,p[:,[2,0,1]],p[:,[1,2,0]]))
                    pressures=closure["combined_stress"](states)
                row.update(flavor_change=convert,destinations=list(map(int,(n0,n1,n2))),
                    converted_encounters_per_particle=(n1+n2)/row["pairs"],
                    physical_momentum_residual_after_conversion=dp,energy_residual_after_conversion=de,
                    all_flavor_pressures_over_rest_energy=pressures)
                return row

        with ThreadPoolExecutor(opts.jobs) as pool:
            result=list(pool.map(one,tasks))
    for convert,destination in ((True,"convert"),(False,"elastic")):
        rows=[row for row in result if row["flavor_change"]==convert]
        evolution["add_common_target"](rows, snapshots/destination, opts.scattering)
    return {"runs":result,"resonant_pair_control":control}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("history","scattering","polarization","output"):
        parser.add_argument("--"+name,type=Path,required=True)
    parser.add_argument("--jobs",type=int,default=2)
    schedule=parser.add_mutually_exclusive_group()
    schedule.add_argument("--endpoint",type=float,help="Run a separate long-time marker with three seeds and a matched control")
    schedule.add_argument("--single",nargs=5,metavar=("TIME","POWER","STEPS","SEED","CONVERT"))
    args=parser.parse_args()
    args.output.parent.mkdir(exist_ok=True,parents=True)
    out=run(args)
    args.output.write_text(json.dumps(out,sort_keys=True,indent=2,allow_nan=False)+"\n")
