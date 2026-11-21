"""Record the incoming and outgoing momenta of each accepted elastic wall collision."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import tempfile

import numpy as np


def insert_once(source, key, replacement):
    if source.count(key) != 1:
        raise ValueError(f"Unrecognized parent source anchor {key!r}")
    return source.replace(key, replacement)


def build_recording_solver(options, wall, work):
    header=(options.scattering/"live/collide.hpp").read_text()
    header=insert_once(header,"class PoleCollisions {","class PoleCollisions {\n    std::ofstream tape;")
    header=insert_once(header,"unsigned i, unsigned j) {","unsigned i, unsigned j, double time) {")
    header=insert_once(header,"encounter(p[i],p[j],radius,i,j);","encounter(p[i],p[j],radius,i,j,left+offset);")
    header=insert_once(header,"public:",'''public:
    void set_tape(const std::string& file) {
        tape.open(file,std::ios::binary);
        if (!tape) throw std::runtime_error("Cannot open collision ledger");
    }
''')
    anchor="first.u=radius*next;"
    report='''const std::array<double,18> event{
            time,radius,static_cast<double>(i),static_cast<double>(j),s,mu,
            p.x,p.y,p.z,q.x,q.y,q.z,next.x,next.y,next.z,other.x,other.y,other.z};
        tape.write(reinterpret_cast<const char*>(event.data()),sizeof(event));
        if (!tape) throw std::runtime_error("Cannot write collision ledger");
        '''
    header=insert_once(header,anchor,report+anchor)
    (work/"collide.hpp").write_text(header)
    origin=(options.axial/"live/solver.cpp").read_text()
    anchor="PoleCollisions collisions(particles.size(),density,seed,false);"
    (work/"solver.cpp").write_text(insert_once(origin,anchor,
        anchor+'\n    collisions.set_tape(std::string(argv[0])+".events");'))
    # The parent's build prepends its pinned vector, wave and integrator definitions.
    wall["build"].__globals__["HERE"]=work
    return wall["build"](options.polarization,options.scattering,work)


def minkowski(p,q):
    """Stable equal-unit-mass four-vector product for nearly parallel momenta."""
    lengthp=np.linalg.norm(p,axis=1)
    lengthq=np.linalg.norm(q,axis=1)
    prod=lengthp*lengthq
    rest=(1+lengthp**2+lengthq**2)/(np.hypot(1.,lengthp)*np.hypot(1.,lengthq)+prod)
    pdir=np.divide(p,lengthp[:,None],out=np.zeros_like(p),where=lengthp[:,None]>0.)
    qdir=np.divide(q,lengthq[:,None],out=np.zeros_like(q),where=lengthq[:,None]>0.)
    return rest+.5*prod*np.sum((pdir-qdir)**2,axis=1)


def verify_and_save(history, orbits, reference, row, tape, snapshot, target):
    final=snapshot["p"]
    raw=np.fromfile(tape,dtype="=f8").reshape(-1,18)
    time,a,raw_i,raw_j,s,mu=raw[:,:6].T
    i,j=raw_i.astype(np.int32),raw_j.astype(np.int32)
    pi=raw[:,6:12].reshape(-1,2,3)
    po=raw[:,12:18].reshape(-1,2,3)
    ei=np.hypot(1.,np.linalg.norm(pi,axis=2))
    eo=np.hypot(1.,np.linalg.norm(po,axis=2))
    n=2*row["pairs"]
    tags=orbits["seeds"](history,{"D":float(np.interp(row["time"],history["t"],history["D"]))},
                          row["pairs"].bit_length()-1,row["birth_seed"])
    birth=np.repeat(np.sort(tags[:,0],kind="stable"),2)
    di=np.maximum(abs(2*(1+minkowski(pi[:,0],pi[:,1]))-s),abs(2*(1+minkowski(po[:,0],po[:,1]))-s))/s
    angular=abs((1+4*(1-minkowski(pi[:,0],po[:,0]))/(s-4))-mu)
    de=abs(ei.sum(axis=1)-eo.sum(axis=1))/ei.sum(axis=1)
    dp=np.linalg.norm(pi.sum(axis=1)-po.sum(axis=1),axis=1)/ei.sum(axis=1)
    pole_s=8300.**2
    width=np.sqrt(1.-4./s)*.30*(s+2.)
    detuning=s-pole_s
    bands={"below":detuning < -width, "within":abs(detuning)<=width, "above":detuning > width}
    event_energy=ei.sum(axis=1)
    bands_report={name:{"event_fraction":float(mask.mean()),
                        "pair_energy_fraction":float(event_energy[mask].sum()/event_energy.sum())}
                  for name,mask in bands.items()}
    endpoint=row["time"]
    weight=np.interp(endpoint,history["t"],history["D"])/n
    gas=float(weight*np.interp(endpoint,history["t"],history["a"])*np.hypot(1.,np.linalg.norm(final,axis=1)).sum())
    counted=np.bincount(np.r_[i,j],minlength=n)
    if (len(raw)!=reference["events"] or row["events"]!=reference["events"] or np.any(raw_i!=i) or
            np.any(raw_j!=j) or np.any((i<0)|(j<0)|(i>=n)|(j>=n)|(i%2!=0)|(j%2!=1)) or
            not np.isfinite(raw).all() or np.any(s<=4) or np.max(abs(mu))>1 or
            np.any(np.diff(time)<0.) or np.any(time<np.maximum(birth[i],birth[j])) or np.any(time>endpoint) or
            np.sum(counted)!=2*reference["events"] or np.mean(counted>0)!=reference["hit"] or
            not np.isclose(gas,reference["gas"],rtol=2.e-8) or max(de.max(),dp.max())>1.e-12 or
            di.max()>1.e-9 or angular.max()>1.e-9):
        raise ValueError(f'Collision ledger mismatch in {row["label"]} (conservation {de.max()}, {dp.max()}, kinematics {di.max()}, {angular.max()})')
    pair_code=i.astype(np.int64)*n+j
    distinct=np.unique(pair_code).size
    metadata={k:float(history[k]) for k in ("action","charge_to_mass","daughter_mass_L","channel_fraction")}
    metadata.update(weight=float(weight),gas=gas,K=2.54,endpoint_a=float(np.interp(endpoint,history["t"],history["a"])),
                    D=float(n*weight),time=endpoint,seed=row["seed"],pairs=row["pairs"],steps=row["steps"],label=row["label"])
    # Publish the NPZ atomically. The analysis workers can read individual cases as soon as they appear.
    temp=target.with_suffix(".partial.npz")
    np.savez(temp,t=time,a=a,i=i,j=j,s=s,mu=mu,p_in=pi,p_out=po,E_in=ei,E_out=eo,
             birth=birth,p_final=final,wall_t=snapshot["t"],wall_E=snapshot["E"],wall_B=snapshot["B"],**metadata)
    temp.replace(target)
    def quant(x):
        return np.quantile(x,[.01,.5,.99],method="linear").tolist()
    return dict(label=row["label"],time=endpoint,pairs=row["pairs"],steps=row["steps"],seed=row["seed"],events=len(raw),
                hit=float(np.mean(counted>0)),mean_encounters=float(np.mean(counted)),weight=float(weight),gas=gas,
                gas_relative_replay_error=float(abs(gas-reference["gas"])/reference["gas"]),
                max_event_relative_energy_error=float(de.max()),max_event_scaled_momentum_error=float(dp.max()),
                max_relative_pair_s_error=float(di.max()),max_cm_cosine_error=float(angular.max()),
                fraction_events_reusing_pair=float(1-distinct/len(raw)),original_birth_pair_events=int(np.count_nonzero(j==i+1)),
                distinct_event_pairs=int(distinct),unique_positive_tracks_hit=int(np.count_nonzero(counted[::2])),
                unique_negative_tracks_hit=int(np.count_nonzero(counted[1::2])),first_event_time=float(time[0]),
                last_event_time=float(time[-1]),sqrt_s_q01_q50_q99=quant(np.sqrt(s)),
                pair_gamma_cm_q01_q50_q99=quant(ei.sum(axis=1)/np.sqrt(s)),
                cm_cosine_q01_q50_q99=quant(mu),
                running_width_bands=bands_report,
                mean_cm_cosine=float(np.mean(mu)),mean_squared_cm_cosine=float(np.mean(mu**2)))


def run(opts):
    opts.artifacts.mkdir(exist_ok=True,parents=True)
    opts.output.parent.mkdir(exist_ok=True,parents=True)
    with np.load(opts.history) as src:
        h=dict(src)
    if h["t"][-1]<7.05 or not np.allclose(h["K"],2.54,rtol=0.,atol=2.e-12):
        raise ValueError("Need the full dedicated branch through wall time 7.05")
    if float(h["channel_fraction"])!=1.:
        raise ValueError("Need the dedicated-interaction normalization")
    wall=runpy.run_path(str(opts.axial/"live/evolve.py"))
    prior=runpy.run_path(str(opts.scattering/"live/evolve.py"))
    orbits=runpy.run_path(str(opts.polarization/"orbits/pairs.py"))
    coarse=runpy.run_path(str(opts.polarization/"orbits/coarse.py"))
    coarse["thermal_cells"].__globals__["quad"]=prior["thermal_quadrature"]
    refs=[r for fname in ("CROSSINGS","LIVE") for r in json.loads((opts.axial/f"live/{fname}.json").read_text())["runs"]
          if r["time"] in (6.67,7.05) and r["mode"]=="particle" and r["pole"]]
    if len(refs)!=8:
        raise ValueError("Expected four archived refinements at each endpoint")
    with tempfile.TemporaryDirectory(prefix="wall-scatter-event-") as wd:
        work=Path(wd)
        binary=build_recording_solver(opts,wall,work)
        def one(ref):
            task=(ref["time"],ref["pairs"].bit_length()-1,ref["steps"],ref["seed"],ref["mode"],ref["pole"])
            label=ref["label"]
            alias=work/label
            alias.symlink_to(binary)
            context=h,orbits,coarse,prior,alias
            row=wall["run_case"](context,task,work,work)
            with np.load(work/(label+".npz")) as sample:
                if not np.all(sample["sign"][::2]==1) or not np.all(sample["sign"][1::2]==-1):
                    raise ValueError("Unexpected charge order in end-state tracks")
                return verify_and_save(h,orbits,ref,row,Path(str(alias)+".events"),sample,opts.artifacts/(label+".npz"))
        with ThreadPoolExecutor(opts.jobs) as pool:
            rows=list(pool.map(one,refs))
    opts.output.write_text(json.dumps({"runs":rows},sort_keys=True,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for argument in ("axial","polarization","scattering","history","output","artifacts"):
        parser.add_argument("--"+argument,type=Path,required=True)
    parser.add_argument("--jobs",type=int,default=2)
    run(parser.parse_args())
