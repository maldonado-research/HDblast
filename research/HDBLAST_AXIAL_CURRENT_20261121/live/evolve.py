"""Advance the three wall Maxwell fields and a homogeneous charged-particle population."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import subprocess
import tempfile

import numpy as np


HERE = Path(__file__).resolve().parent
AMPLITUDE = 0.04 / np.sqrt(2)
MODES = {"reference": 0, "fluid": 1, "particle": 2}


def build(polarization, scattering, work, record_contacts=False):
    original = (polarization/"feedback/homogeneous.cpp").read_text()
    prefix, separator, _ = original.partition("int main(int argc,char** argv)")
    if not separator:
        raise ValueError("Unexpected pinned Maxwell solver")
    if record_contacts:
        # Expose the counters already used to compute the pinned kernel's "hit" fraction.
        header = (scattering/"live/collide.hpp").read_text()
        anchor = 'if (!out) throw std::runtime_error("Missing collision diagnostics");'
        if header.count(anchor) != 1:
            raise ValueError("Unexpected pinned collision reporter")
        reporter = r'''std::ofstream tally(std::string(file)+".per_particle",std::ios::binary);
        for (unsigned count : contacts) {
            uint64_t value = count;
            tally.write(reinterpret_cast<const char*>(&value), sizeof(value));
        }
        if (!tally) throw std::runtime_error("Missing contact tally");'''
        (work/"collide.hpp").write_text(header.replace(anchor, anchor+"\n"+reporter))
    source = work/"solver.cpp"
    source.write_text(prefix + (HERE/"solver.cpp").read_text())
    binary = work/"solver"
    subprocess.run(["g++", "-O3", "-ffp-contract=off", "-std=c++17", "-I", str(scattering/"live"), str(source), "-o", str(binary)], check=True)
    return binary


def zeroes(t, y):
    """Times of sign change on a sampled trace. Discard the initial exactly-zero state."""
    nonzero = np.flatnonzero(y != 0.)
    if not len(nonzero):
        return []
    left = nonzero[:-1][np.signbit(y[nonzero[:-1]]) != np.signbit(y[nonzero[1:]])]
    crossings = []
    for i in left:
        j = nonzero[np.searchsorted(nonzero, i)+1]
        crossings.append(float(t[i] - y[i]*(t[j]-t[i])/(y[j]-y[i])))
    return crossings


def field_table(history, time, steps, extra):
    t = np.linspace(history["t"][0], time, steps+1)
    def at(key):
        return np.interp(t, history["t"], history[key])
    columns = [t, at("a")]
    for name in ("E", "B", "I"):
        third = {"E":"X", "B":"B", "I":"incoming"}[name]
        columns.extend([AMPLITUDE*at(name+"1"), AMPLITUDE*at(name+"2"), at(third)])
    columns.extend([at("J"), at("D"), at("work_fraction")])
    table = np.column_stack(columns)
    if extra is not None:
        if table[0,0] < extra["t"][0]-1e-12 or time > extra["t"][-1]+1e-12:
            raise ValueError("The return must cover this whole time interval")
        table[:,8:11] += np.column_stack([np.interp(t, extra["t"], extra["delta_incoming"][:,j]) for j in range(3)])
    return table


def run_case(context, task, work, output_dir, extra=None, record_contacts=False):
    history, orbits, coarse, live, binary = context
    time, power, steps, seed, mode, pole = task
    count = 2**power
    label = f"t{time:.5f}-p{power}-n{steps}-r{seed}-{mode}-{'pole' if pole else 'off'}"
    mark = {"time":time, **{key:float(np.interp(time, history["t"], history[key]))
                            for key in ("D","a","theta","R","X","B")}}
    table = field_table(history, time, steps, extra)
    # Vary the birth scramble with sample resolution. At fixed power the time-step refinements use identical tags.
    birth_seed = orbits["SCRAMBLE"] + {11:19, 12:1, 13:71, 14:121}.get(power, power)
    tags = orbits["seeds"](history, mark, power, birth_seed)
    tags[:,1] = np.interp(tags[:,0], history["t"], history["a"])
    density = float(history["action"])/(6.*float(history["daughter_mass_L"])**3)
    with tempfile.TemporaryDirectory(dir=work) as folder:
        folder=Path(folder)
        table.tofile(folder/"fields"); tags.tofile(folder/"births")
        subprocess.run([str(binary),str(folder/"fields"),str(folder/"births"),str(mark["D"]),str(MODES[mode]),"1",
                        str(density if pole else 0),str(seed),str(folder/"counts"),str(folder/"particles"),
                        str(folder/"trace"),str(folder/"ledger")],check=True)
        tally={}
        if record_contacts:
            count_per_particle=np.fromfile(folder/"counts.per_particle",dtype="=u8")
            if len(count_per_particle)!=2*count:
                raise ValueError("Incomplete event counts")
            tally={"contacts":count_per_particle, "birth":np.repeat(np.sort(tags[:,0],kind="stable"),2)}
        raw=np.fromfile(folder/"particles",dtype="=f8").reshape(-1,4)
        trace=np.fromfile(folder/"trace",dtype="=f8").reshape(-1,13)
        budget=[np.fromstring(line,sep=" ") for line in (folder/"ledger").read_text().splitlines()]
        _,events,hit,dE,dp,kmax,_,_=np.loadtxt(folder/"counts")
    if len(raw)!=2*count or len(trace)!=steps+1 or not np.isfinite(trace).all() or not np.isfinite(raw).all():
        raise ValueError("Incomplete simulation")
    initial,final,flux_in,flux_out,source_work,particle_work = budget[:6]
    gas,birth,expansion,max_iterations,iteration_error=budget[6]
    scale=max(gas, float(np.sum(initial)), float(np.sum(flux_in)),1.)
    field_residual=final-initial-flux_in+flux_out+source_work
    matter_residual=gas-birth-expansion-particle_work.sum()
    if np.max(abs(field_residual))/scale>2e-10 or abs(matter_residual)/scale>2e-10:
        raise ValueError(("Energy budget",field_residual/scale,matter_residual/scale))
    result_mat=np.zeros((len(raw),8))
    result_mat[:,0]=raw[:,0];result_mat[:,1]=1.;result_mat[:,5:]=raw[:,1:]
    c=coarse["compare"](history,mark,result_mat)
    target_pressure=mark["theta"]/np.hypot(1.,np.interp(time,history["t"],history["u"]))
    energies=np.hypot(1.,np.linalg.norm(raw[:,1:],axis=1))
    local=[]
    for sign in (1.,-1.):
        subset=raw[:,0]==sign
        p=raw[subset,1:3]
        local.append(np.linalg.eigvalsh((p.T/energies[subset]) @ p / subset.sum())/target_pressure)
    e=trace[:,1:4];b=trace[:,4:7];current=trace[:,7:10]
    combined=live["combined_stress"](raw[:,1:])
    reference=table[:,2:5]
    delta=e-reference
    midt=.5*(trace[1:,0]+trace[:-1,0])
    end_current=np.array([np.mean(current[max(1, len(current)-4):,j]) for j in range(3)])
    np.savez(output_dir/(label+".npz"), t=trace[:,0],E=e,B=b,J=current,J_used=trace[:,10:],reference_E=reference,
        delta_E=delta,reference_B=table[:,5:8], incoming=table[:,8:11],p=raw[:,1:],sign=raw[:,0],time=time,
        mode=mode, pole=pole,power=power,steps=steps,seed=seed,**tally)
    return dict(label=label,time=time,work_fraction=float(table[-1,-1]),mode=mode,pole=bool(pole),pairs=count,steps=steps,
                seed=seed,birth_seed=birth_seed,events=int(events),hit=float(hit),mean_encounters=float(events/count),
                tv20=c["coarse_total_variation"],coarse_counts=c["observed"],target_counts=c["thermal"],
                pressures_over_stream_target=np.mean(local,axis=0).tolist(),rest_pressures_over_energy=combined,
                end_E=e[-1].tolist(),end_B=b[-1].tolist(),end_reference_E=reference[-1].tolist(),
                peak_E=np.max(abs(e),axis=0).tolist(),peak_delta_E=np.max(abs(delta),axis=0).tolist(),
                peak_delta_B=np.max(abs(b-table[:,5:8]),axis=0).tolist(),
                first_field_reversals=zeroes(trace[:,0],e[:,2]),first_reference_field_reversals=zeroes(trace[:,0],reference[:,2]),
                first_current_reversals=zeroes(midt,current[1:,2]),first_fluid_current_reversals=zeroes(table[:,0],table[:,11]),
                end_current=end_current.tolist(), end_reference_current=float(table[-1,11]),
                max_field_residual_over_scale=float(np.max(abs(field_residual))/scale),
                matter_residual_over_scale=float(matter_residual/scale), collision_energy_residual=dE,
                collision_momentum_residual=dp,max_kernel_ratio=kmax, implicit_iterations=int(max_iterations),
                implicit_max_residual=iteration_error,gas=float(gas), birth=float(birth),expansion=float(expansion),
                initial_field_energy=initial.tolist(), final_field_energy=final.tolist(), reservoir_in=flux_in.tolist(),
                radiated_out=flux_out.tolist(), source_work=source_work.tolist(), particle_work=particle_work.tolist(),
                imposed_current_work_difference=(particle_work-source_work).tolist())


def add_common_target(rows, snapshots, scattering):
    """Compare with a common thermal gas at each measured energy and momentum."""
    kinetics = runpy.run_path(str(scattering/"kinetic/scatter.py"))
    beam = np.array([0., 0., 1.])
    for row in rows:
        with np.load(snapshots/(row["label"]+".npz")) as sample:
            p = sample["p"]
        theta, u = kinetics["equilibrium"](p)
        target = kinetics["target_cells"](theta, u, beam)
        observed = kinetics["histogram"](p, theta, u, beam)
        row["common_gas_theta"] = float(theta)
        row["common_gas_four_velocity"] = u.tolist()
        row["common_gas_tv20"] = float(0.5*np.abs(observed-target).sum())


def simulate(opts):
    with np.load(opts.history) as data:
        history=dict(data)
    if not np.allclose(history["K"],2.54,atol=2e-12,rtol=0):
        raise ValueError("Gauge coefficient must be constant over this interval")
    live=runpy.run_path(str(opts.scattering/"live/evolve.py"))
    orbits=runpy.run_path(str(opts.polarization/"orbits/pairs.py"))
    coarse=runpy.run_path(str(opts.polarization/"orbits/coarse.py"))
    coarse["thermal_cells"].__globals__["quad"]=live["thermal_quadrature"]
    if opts.return_input:
        with np.load(opts.return_input) as src:
            extra=dict(src)
    else:extra=None
    if opts.single:
        tm,power,n,seed,mode,pole=opts.single
        tasks=[(float(tm),int(power),int(n),int(seed),mode,bool(int(pole)))]
    else:
        times=[float(np.interp(f,history["work_fraction"],history["t"])) for f in (.01,.05)]+[6.45,6.67]
        tasks=[]
        for time in times:
            base=512 if time<6.4 else (1024 if time<6.5 else 2048)
            for power,n,seed in ((12,base,18431),(13,2*base,70163),(13,base,70163),(14,2*base,91283)):
                for mode,pole in (("fluid",True),("particle",True)):
                    tasks.append((time,power,n,seed,mode,pole))
                if power==12 or power==14:
                    tasks.extend((time,power,n,seed,mode,False) for mode in ("reference","fluid","particle"))
    if any(tm < history["t"][0] or tm > history["t"][-1] for tm,*_ in tasks):
        raise ValueError("Requested time outside the actual input field history")
    snapshots=opts.output.parent/"snapshots"
    snapshots.mkdir(exist_ok=True,parents=True)
    with tempfile.TemporaryDirectory(prefix="three-component-wall-") as path:
        work=Path(path)
        context=history,orbits,coarse,live,build(opts.polarization,opts.scattering,work)
        with ThreadPoolExecutor(opts.jobs) as pool:
            results=list(pool.map(lambda task:run_case(context,task,work,snapshots,extra),tasks))
    add_common_target(results, snapshots, opts.scattering)
    opts.output.write_text(json.dumps({"runs":results},indent=2,sort_keys=True,allow_nan=False)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("history","scattering","polarization","output"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--jobs",type=int,default=2)
    p.add_argument("--single",nargs=6,metavar=("TIME","LOG2_PAIRS","STEPS","SEED","MODE","POLE"))
    p.add_argument("--return-input",type=Path,help="NPZ with t and delta_incoming in physical field units")
    simulate(p.parse_args())
