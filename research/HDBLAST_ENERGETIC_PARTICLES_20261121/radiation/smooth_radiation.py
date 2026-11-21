"""Postprocess classical coherent-field cooling on the archived HDblast dedicated wall trajectories."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import subprocess
import tempfile

import numpy as np


# Physical shell-frame momenta are in units of the daughter mass.
CPP = r"""
double proper_force_squared(V p,V electric,V magnetic) {
    double magnitude=length(p);
    if(magnitude==0.)return dot(electric,electric);
    const V axis=(1./magnitude)*p;
    const double parallel=dot(electric,axis);
    const V transverse=cross(axis,cross(electric,axis));
    const V boost=energy(1.,p)*transverse+cross(p,magnetic);
    return parallel*parallel+dot(boost,boost);
}
void radiation_controls(){
    auto check=[](double observed,double expected){
        if(std::abs(observed-expected)>3e-14*std::max(1.,std::abs(expected)))
            throw std::runtime_error("Analytic radiation check");
    };
    // Parallel electric acceleration is independent of particle energy; a transverse magnetic field scales as p².
    check(proper_force_squared({0.,0.,1e8},{0.,0.,2.5},{0.,0.,3.}),6.25);
    check(proper_force_squared({0.,0.,1e8},{},{0.,3.,0.}),9e16);
    check(proper_force_squared({0.,0.,0.},{1.,2.,3.},{2.,4.,6.}),14.);
    V p={1.5,2.,-3.5},e={1.2,-2.3,.4},b={-.1,1.5,.7};
    V k=energy(1.,p)*e+cross(p,b);
    check(proper_force_squared(p,e,b),dot(k,k)-dot(p,e)*dot(p,e));
}
"""
START = r"""
    radiation_controls();
    double largest_chi_squared=0., envelope_chi_squared=0.;
    const double mu=std::stod(argv[12]), alpha=std::stod(argv[13]);
    const double critical=Q/mu, radiation_factor=(2./3.)*alpha*mu;
    std::vector<double> radiation(2*births.size(),0.), upper(2*births.size(),0.), physical(2*births.size(),0.);
"""
FIELDS = r"""
        const V E=(1./(a*a))*avg,B=(1./(a*a))*wave.b;
        const double envelope=length(E)+length(B);
"""
COLLECT = r"""
            // Before the second magnetic half-rotation. Tailor the duration to each birth.
            const double age=r[0]-std::max(v.birth,l[0]);
            const V midpoint=(.5/a)*(v.u+next);
            const double chi2=critical*critical*proper_force_squared(midpoint,E,B);
            const double max_chi=critical*std::max(v.g0,energy(a,next))/a*envelope;
            if(chi2<0. || !std::isfinite(chi2) || chi2>max_chi*max_chi*(1.+1e-12))
                throw std::runtime_error("Force invariant or bound failed");
            const double loss=age*radiation_factor*chi2;
            radiation[i]+=a*loss; physical[i]+=loss;
            upper[i]+=a*age*radiation_factor*max_chi*max_chi;
            largest_chi_squared=std::max(largest_chi_squared,chi2);
            envelope_chi_squared=std::max(envelope_chi_squared,max_chi*max_chi);
"""
SAVE = r"""
    std::ofstream losses(argv[14],std::ios::binary);
    for(size_t i=0;i<particles.size();i++){
        const std::array<double,4> row={radiation[i],upper[i],physical[i],energy(final_a,particles[i].u)};
        losses.write(reinterpret_cast<const char*>(row.data()),sizeof(row));
    }
    std::ofstream peak(argv[15]);
    peak<<std::setprecision(17)<<std::sqrt(largest_chi_squared)<<" "<<std::sqrt(envelope_chi_squared)<<"\n";
    if(!peak||!losses)throw std::runtime_error("Radiation output");
"""


def replace_once(src, marker, value):
    if src.count(marker)!=1:
        raise ValueError(f"Upstream source changed near {marker}")
    return src.replace(marker,value)


def build(study,polarization,scattering,work):
    old=(polarization/"feedback/homogeneous.cpp").read_text()
    prefix,sep,_=old.partition("int main(int argc,char** argv)")
    if not sep:
        raise ValueError("Unexpected momentum representation")
    src=(study/"live/solver.cpp").read_text()
    src=replace_once(src,"if(argc!=12)","if(argc!=16)")
    for marker,addition in (("const double weight=",START),
            ("current={};\n        for(size_t i=0;i<n;++i){\n            auto& v=",FIELDS),
            ("const V flux=(weight*v.factor",COLLECT),
            ("collisions.save(argv[8]);",SAVE)):
        src=replace_once(src,marker,addition+"\n"+marker)
    path=work/"radiation.cpp"
    path.write_text(prefix+CPP+src)
    binary=work/"radiation"
    subprocess.run(["g++","-O3","-ffp-contract=off","-std=c++17","-I",str(scattering/"live"),str(path),"-o",str(binary)],check=True)
    return binary


def run(context, task):
    hist,orbits,table_builder,binary,work,alpha= context
    time,power,steps,seed,pole=task
    with tempfile.TemporaryDirectory(dir=work) as td:
        d=Path(td)
        mark=dict(time=time,**{k:float(np.interp(time,hist["t"],hist[k])) for k in ("D","a","theta","R","X","B")})
        table=table_builder(hist,time,steps,None)
        table.tofile(d/"inputs")
        birth_seed=orbits["SCRAMBLE"]+{11:19,12:1,13:71,14:121}[power]
        tags=orbits["seeds"](hist,mark,power,birth_seed)
        tags[:,1]=np.interp(tags[:,0],hist["t"],hist["a"])
        tags.tofile(d/"births")
        mu=float(hist["daughter_mass_L"])
        density=float(hist["action"])/(6*mu**3) if pole else 0.
        subprocess.run([str(binary),str(d/"inputs"),str(d/"births"),str(mark["D"]),"2","1",str(density),str(seed),
            str(d/"counts"),str(d/"particles"),str(d/"trace"),str(d/"ledger"),str(mu),str(alpha),str(d/"loss"),str(d/"peak")],check=True)
        loss,upper,local,final=np.fromfile(d/"loss",dtype="=f8").reshape(-1,4).T
        peak,envelope=np.loadtxt(d/"peak")
        events=np.loadtxt(d/"counts")[1]
        gas=float((d/"ledger").read_text().splitlines()[-1].split()[0])
        w=mark["D"]/len(loss)
        if not np.allclose(w*sum(final),gas,rtol=3e-13) or len(loss)!=2**(power+1) or not (upper>=loss).all():
            raise ValueError("Particle count, energy, or radiation bound failed")
        # Separate endpoint calculation using the direct four-vector invariant in extended precision.
        p=np.fromfile(d/"particles",dtype="=f8").reshape(-1,4)[:,1:].astype(np.longdouble)
        trace=np.fromfile(d/"trace",dtype="=f8").reshape(-1,13)[-1].astype(np.longdouble)
        electric,magnetic=trace[1:4],trace[4:7]
        gamma=np.sqrt(1.+(p*p).sum(axis=1))
        force=gamma[:,None]*electric+np.cross(p,magnetic)
        chi2=((float(hist["charge_to_mass"])/mu)**2)*((force*force).sum(axis=1)-(p@electric)**2)
        if np.any(chi2<0.):
            raise ValueError("Unresolved endpoint force invariant")
        endpoint_rate=float((2./3.)*alpha*mu*chi2.sum()/gamma.sum())
        fast=final>=np.quantile(final,.9)
        return dict(time=time,log2_pairs=power,steps=steps,collision_seed=seed,birth_seed=birth_seed,pole=pole,
            reference_work_fraction=float(table[-1,-1]),mean_collisions=float(events/(2**power)),gas_conformal=gas,
            cumulative_emitted_conformal=float(w*sum(loss)),cumulative_emitted_local=float(w*sum(local)),
            emitted_over_final_gas=float(sum(loss)/sum(final)),force_envelope_over_final_gas=float(sum(upper)/sum(final)),
            maximum_particle_loss_over_its_final_conformal_energy=float(np.max(loss/final)),
            fraction_of_emission_from_final_top_energy_decile=float(sum(loss[fast])/sum(loss)),
            emitted_over_final_top_decile_energy=float(sum(loss[fast])/sum(final[fast])),max_chi=float(peak),
            max_chi_envelope=float(envelope),endpoint_max_chi=float(np.sqrt(chi2.max())),
            endpoint_emission_over_gas_per_wall_time=endpoint_rate)


def main(args):
    history=dict(np.load(args.history))
    if np.any(history["K"]!=2.54) or float(history["channel_fraction"])!=1. or history["t"][-1]<7.05:
        raise ValueError("Requires the original extended dedicated trajectory on the gauge plateau")
    q=float(history["charge_to_mass"]);mu=float(history["daughter_mass_L"]);act=float(history["action"])
    alpha=float(3*q*q*mu*mu/(4*np.pi*act*2.54))
    if q!=1e7 or mu!=4e12:
        raise ValueError("Unrecognized inherited charge or daughter mass")
    orbits=runpy.run_path(str(args.polarization/"orbits/pairs.py"))
    table=runpy.run_path(str(args.study/"live/evolve.py"))["field_table"]
    tasks=[(6.45,12,1024,18431,True),(6.67,12,2048,18431,True),(6.67,14,4096,91283,True),(7.05,12,4096,18431,True),
           (7.05,13,4096,70163,True),(7.05,13,8192,70163,True),(7.05,13,4096,70163,False),
           (7.05,14,8192,91283,True)]
    with tempfile.TemporaryDirectory(prefix="wall-coherent-radiation-") as temp:
        td=Path(temp)
        context=(history,orbits,table,build(args.study,args.polarization,args.scattering,td),td,alpha)
        with ThreadPoolExecutor(args.jobs) as pool:
            data=list(pool.map(lambda task:run(context,task),tasks))
    parent=sum((json.loads((args.study/name).read_bytes())["runs"] for name in ("live/LIVE.json","live/CROSSINGS.json")),[])
    for r in data:
        ref=next(x for x in parent if x["time"]==r["time"] and x["mode"]=="particle" and x["seed"]==r["collision_seed"]
                 and x["steps"]==r["steps"] and x["pole"]==r["pole"] and x["pairs"]==2**r["log2_pairs"])
        if not np.isclose(r["gas_conformal"],ref["gas"],rtol=2.e-6) or r["mean_collisions"]!=ref["mean_encounters"]:
            raise ValueError("Instrumented late trajectory differs from archived trajectory")
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(dict(alpha=alpha,Q=q,mu=mu,K=2.54,action=act,runs=data),indent=2,sort_keys=True,allow_nan=False)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for option in ("history","study","polarization","scattering","output"):
        p.add_argument("--"+option,type=Path,required=True)
    p.add_argument("--jobs",type=int,default=2)
    main(p.parse_args())
