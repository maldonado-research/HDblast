"""Isolate the role of the cross section at identical numbers of broad elastic encounters."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import subprocess
import tempfile

import numpy as np


HERE=Path(__file__).resolve().parent
CHECKPOINTS=(0,1,2,4,8,16)


def build(polarization, scattering, tmp):
    source=(polarization/"feedback/homogeneous.cpp").read_text().split("int main(")[0]
    (tmp/"vector.hpp").write_text(source)
    header=(scattering/"live/collide.hpp").read_text()
    select="const double accept=kernel/kernel_max;"
    anchor="public:\n    PoleCollisions"
    if header.count(select)!=1 or header.count(anchor)!=1:
        raise ValueError("Unexpected vector-pole implementation")
    header=header.replace(select,"const double accept=flat ? beta*s/(4.*ep*eq) : kernel/kernel_max;")
    header=header.replace(anchor,"""public:
    bool flat=false;
    double attempt(std::vector<Particle>&);
    unsigned hits(size_t i) const {return contacts.at(i);}
    PoleCollisions""")
    (tmp/"pole.hpp").write_text(header)
    exe=tmp/"box"
    subprocess.run(["g++","-O3","-ffp-contract=off","-std=c++17","-I",str(tmp),str(HERE/"box.cpp"),"-o",str(exe)],check=True)
    return exe


def select_snapshot(folder, time):
    files=sorted(folder.glob("t*-p*-particle-pole.npz"))
    candidates=[]
    for f in files:
        tm=float(f.name.split("-")[0][1:])
        if abs(tm-time)<1e-4:
            parts=f.stem.split("-")
            candidates.append((int(parts[1][1:]),int(parts[2][1:]),f))
    if not candidates:
        raise ValueError(f"Missing particle-field snapshot near {time:g}")
    return max(candidates)[-1]


def run(opts):
    kinetic=runpy.run_path(str(opts.scattering/"kinetic/scatter.py"))
    live=runpy.run_path(str(opts.scattering/"live/evolve.py"))
    kinetic["target_cells"].__globals__["quad"]=live["thermal_quadrature"]
    axis=np.array([0.,0.,1.])
    jobs=[]
    for tm in (.0,6.67,7.05):
        folder=opts.extended if tm==7.05 else opts.snapshots
        file=select_snapshot(folder,6.278545 if tm==0 else tm)
        for power,seed in ((12,1583),(13,62873)):
            for flat in (False,True):
                jobs.append((file,power,seed,flat))
    with tempfile.TemporaryDirectory(prefix="vector-box-") as tmpdir:
        tmp=Path(tmpdir)
        exe=build(opts.polarization,opts.scattering,tmp)
        def one(job):
            path,power,seed,flat=job
            with np.load(path) as s:
                p,sign,time=s["p"],s["sign"],float(s["time"])
            if not (np.all(sign[::2]==1) and np.all(sign[1::2]==-1)):
                raise ValueError("Particles are not ordered by charge-paired births")
            n=2**power
            rng=np.random.default_rng(seed)
            chosen=rng.choice(len(p)//2,n,replace=False)
            idx=np.stack((2*chosen,2*chosen+1),axis=1).ravel()
            p=p[idx];sign=sign[idx]
            theta,drift=kinetic["equilibrium"](p)
            target=kinetic["target_cells"](theta,drift,axis)
            g0=kinetic["gamma"](p)
            high=g0>=np.quantile(g0,.90)
            with tempfile.TemporaryDirectory(dir=tmp) as temp:
                work=Path(temp)
                np.column_stack((sign,p)).tofile(work/"in")
                subprocess.run([str(exe),str(work/"in"),str(work/"snapshots"),str(seed),str(int(flat)),
                                str(CHECKPOINTS[-1]),str(work/"diag")],check=True)
                out=np.fromfile(work/"snapshots",dtype="=f8").reshape(len(CHECKPOINTS),1+8*n)
                candidates,events,_,de,dp,_,_,_=np.loadtxt(work/"diag")
            rows=[]
            for count,snap in zip(CHECKPOINTS,out,strict=True):
                vals=snap[1:].reshape(2*n,4)
                q,hits=vals[:,:3],vals[:,3]
                g=kinetic["gamma"](q)
                if hits.sum()!=count*2*n or np.linalg.norm(q.sum(axis=0)-p.sum(axis=0))/g0.sum()>1.e-10 or abs(g.sum()/g0.sum()-1)>1.e-10:
                    raise ValueError("Count or four-momentum violation")
                per_sign=[float(.5*np.abs(kinetic["histogram"](q[sign==s],theta,drift,axis)-target).sum()) for s in (1.,-1.)]
                per_sign_speed=[float(abs(np.mean(q[sign==s,2]/g[sign==s]))) for s in (1.,-1.)]
                pooled=kinetic["histogram"](q,theta,drift,axis)
                rows.append(dict(encounters=count,normalized_elapsed=float(snap[0]),hit=float(np.mean(hits>0)),
                    common_thermal_tv_per_sign=per_sign,pooled_common_thermal_tv=float(.5*np.abs(pooled-target).sum()),
                    energy_only_tv=float(.5*np.abs(pooled.sum(axis=1)-target.sum(axis=1)).sum()),
                    mean_axial_speed_per_sign=per_sign_speed,pressure_over_rest_energy=live["combined_stress"](q),
                    contacts_per_initial_top_energy_decile=float(np.mean(hits[high])),
                    contacts_per_initial_bottom_nine_deciles=float(np.mean(hits[~high])),
                    fraction_initial_energy_in_unhit=float(g0[hits==0].sum()/g0.sum())))
            return dict(time=time,pairs=n,seed=seed,flat_cross_section=flat,theta=float(theta),
                        equilibrium_four_velocity=drift.tolist(),checks=rows,proposals=int(candidates),
                        events=int(events),max_event_energy_residual=float(de),max_event_momentum_residual=float(dp))
        with ThreadPoolExecutor(opts.jobs) as pool:
            return list(pool.map(one,jobs))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("polarization","scattering","snapshots","extended","output"):
        parser.add_argument("--"+key,type=Path,required=True)
    parser.add_argument("--jobs",type=int,default=2)
    options=parser.parse_args()
    result=run(options)
    options.output.parent.mkdir(parents=True,exist_ok=True)
    options.output.write_text(json.dumps(dict(runs=result),sort_keys=True,allow_nan=False,indent=2)+"\n")
