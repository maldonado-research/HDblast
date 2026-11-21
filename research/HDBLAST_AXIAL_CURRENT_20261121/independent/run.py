"""Test the cyclic flavor reduction against a finite particle gas with independently evolving flavor occupancies."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import subprocess
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
WALL = HERE.parent/"live"


def prefix(source, marker):
    if source.count(marker) != 1:
        raise ValueError(f"Unexpected upstream source near {marker}")
    return source.partition(marker)[0]


def build(polarization, scattering, work):
    # Reuse just the verified CM kernel. The varying-density candidate clock is implemented in ensemble.cpp.
    source = (scattering/"live/collide.hpp").read_text()
    marker = "    void encounter(Particle&"
    if source.count(marker) != 1 or "kernel_max=1.118569e-6" not in source:
        raise ValueError("The inherited collision kernel or its known majorant changed")
    (work/"collide.hpp").write_text(source.replace(marker, "public:\n"+marker))
    vectors = prefix((polarization/"feedback/homogeneous.cpp").read_text(), "int main(int argc,char** argv)")
    wall = prefix((WALL/"solver.cpp").read_text(), "V source_current(")
    (work/"gas.cpp").write_text(vectors + wall + (HERE/"ensemble.cpp").read_text())
    binary = work/"gas"
    subprocess.run(["g++","-std=c++17","-O3","-ffp-contract=off",str(work/"gas.cpp"),"-o",str(binary)], check=True)
    return binary


def axisymmetric_gap(p, weights, fourth_rank):
    """Lower-bound angular TV against every axisymmetric target using four nonoverlapping caps."""
    weights=weights/weights.sum()
    lengths=np.linalg.norm(p,axis=1)
    directions=np.divide(p,lengths[:,None],out=np.zeros_like(p),where=lengths[:,None]>0)
    best=dict(lower_bound=0.,half_angle_deg=0.,four_cap_mass_floor=0.,axisymmetric_ceiling=0.)
    for angle in (2.,3.,5.,8.,10.,15.,20.,25.,30.):
        radians=np.radians(angle)
        # Whatever the target axis, at least two coordinate axes are >=45 degrees
        # away from both poles. A latitude hits each of those four caps over a
        # fraction <= asin(sqrt(2) sin(alpha))/pi, independently of the latitude.
        masses=(abs(directions)>=np.cos(radians)).T@weights
        floor=float(masses.sum()-masses.max())
        ceiling=float(4*np.arcsin(np.sqrt(2)*np.sin(radians))/np.pi)
        gap=floor-ceiling
        if gap>best["lower_bound"]:
            best=dict(lower_bound=float(gap),half_angle_deg=angle,four_cap_mass_floor=floor,
                      axisymmetric_ceiling=ceiling)
    cubic=float(2.5*(weights@(directions**4).sum(axis=1)-.6))
    if np.any(lengths==0.):
        raise ValueError("Angular moments require resolved nonzero momenta")
    moments=fourth_rank(directions,weights,cubic,cyclic=False)
    best["fourth_harmonic_gap"]=moments["natural_cyclic_axial_gap"]
    best["cubic_harmonic"]=cubic
    return best


def simulate(context, task):
    history, orbits, coarse, inherited, kinetics, table_builder, fourth_rank, binary, directory = context
    time, power, steps, seed, mode, identical = task
    samples = 2**power
    mark = dict(time=time)
    for key in ("D","a","theta","R","X","B"):
        mark[key] = float(np.interp(time, history["t"], history[key]))
    table = table_builder(history, time, steps, None)
    base_seed = orbits["SCRAMBLE"] + {11:19,12:1,13:71}[power]
    birth_sets = []
    for flavor in range(3):
        seed_offset = 0 if identical else flavor*1181
        tags = orbits["seeds"](history, mark, power, base_seed+seed_offset)
        tags[:,1] = np.interp(tags[:,0], history["t"], history["a"])
        birth_sets.append(np.column_stack((tags, np.full(samples,flavor))))
    density = float(history["action"])/(18.*float(history["daughter_mass_L"])**3)
    with tempfile.TemporaryDirectory(dir=directory) as scratch:
        scratch=Path(scratch)
        table.tofile(scratch/"table")
        np.concatenate(birth_sets).tofile(scratch/"tags")
        subprocess.run([str(binary),str(scratch/"table"),str(scratch/"tags"),str(mark["D"]),str(density),str(seed),
                        str(mode),str(scratch/"raw"),str(scratch/"ledger"),str(scratch/"counts")],check=True)
        raw=np.fromfile(scratch/"raw",dtype="=f8").reshape(-1,6)
        rows=[np.fromstring(line,sep=" ") for line in (scratch/"ledger").read_text().splitlines()]
        _,events,hit,delta_e,delta_p,max_kernel,_,_=np.loadtxt(scratch/"counts")
    if raw.shape!=(6*samples,6) or not np.isfinite(raw).all():
        raise ValueError("Particle history is incomplete")
    fields=np.array(rows[:3])
    initial,final,incident,escaped,work,birth,expansion,gas=rows[3]
    transitions=np.array(rows[4:7],dtype=int)
    peak_count_range=int(rows[7][0])
    scale=max(1.,gas,initial)
    field_error=(final-initial-incident+escaped+work)/scale
    gas_error=(gas-birth-expansion-work)/scale
    flavors,signs,lab,hits=raw[:,0].astype(int),raw[:,1],raw[:,2:5],raw[:,5]
    plus=np.bincount(flavors[signs>0],minlength=3)
    minus=np.bincount(flavors[signs<0],minlength=3)
    if (abs(field_error)>2e-10 or abs(gas_error)>2e-10 or delta_e>1e-12 or delta_p>1e-12
            or sum(hits)!=2*events or transitions.sum()!=events or not np.array_equal(plus,minus)
            or not np.array_equal(plus, samples+transitions.sum(axis=0)-transitions.sum(axis=1))):
        raise ValueError(("Conservation failure",field_error,gas_error,int(events),plus,minus,transitions))
    if mode<2 and np.any(plus!=samples):
        raise ValueError("Unexpected flavor transfer")
    if identical and np.max(np.ptp(fields,axis=0))>1e-12:
        raise ValueError(("Broken cyclic symmetry in collisionless control",fields))
    local_raw=np.zeros((6*samples,8))
    local_raw[:,0]=signs
    local_raw[:,1]=1.
    tvs=[]
    pressures=[]
    for f in range(3):
        which=flavors==f
        local_raw[which,5:]=np.roll(lab[which],-f,axis=1)
        pressures.append(inherited["combined_stress"](lab[which]))
        tvs.append(coarse["compare"](history,mark,local_raw[which])["coarse_total_variation"])
    tv=coarse["compare"](history,mark,local_raw)["coarse_total_variation"]
    beam=np.array([0.,0.,1.])
    theta,u=kinetics["equilibrium"](lab)
    observed=kinetics["histogram"](lab, theta, u, beam)
    target=kinetics["target_cells"](theta, u, beam)
    return dict(time=time,pairs_per_birth_flavor=samples,steps=steps,seed=seed,birth_seed_first=base_seed,
                shared_mode=("off","retain","uniform")[mode],identical_births=identical,
                events=int(events),mean_encounters=events/(3*samples),ever_hit=float(hit),
                positive_per_flavor=plus.tolist(),peak_pair_count_range=peak_count_range,transition_matrix=transitions.tolist(),
                end_local_electric=fields.tolist(),mean_end_local_electric=fields.mean(axis=0).tolist(),
                tv20_against_stream_target=float(tv),tv20_by_flavor=tvs,
                global_pressures_over_energy=inherited["combined_stress"](lab),per_flavor_pressures_over_energy=pressures,
                number_axisymmetric_gap=axisymmetric_gap(lab,np.ones(len(lab)),fourth_rank),
                energy_axisymmetric_gap=axisymmetric_gap(lab,np.hypot(1.,np.linalg.norm(lab,axis=1)),fourth_rank),
                global_tv20_to_common_fit=float(.5*np.abs(observed-target).sum()),
                field_balance_over_scale=float(field_error),matter_balance_over_scale=float(gas_error),
                max_collision_energy_residual=float(delta_e),max_collision_momentum_residual=float(delta_p),
                max_kernel_ratio=float(max_kernel),field_initial=float(initial),field_final=float(final),
                flux_in=float(incident),flux_out=float(escaped),work_to_particles=float(work),gas_energy=float(gas),
                injected_energy=float(birth),expansion_energy=float(expansion))


def run(args):
    with np.load(args.history) as data:
        history=dict(data)
    if not np.all(history["K"]==2.54) or float(history["channel_fraction"])!=1./3.:
        raise ValueError("Expected shared gauge history")
    inherited=runpy.run_path(str(args.scattering/"live/evolve.py"))
    kinetics=runpy.run_path(str(args.scattering/"kinetic/scatter.py"))
    coarse=runpy.run_path(str(args.polarization/"orbits/coarse.py"))
    coarse["thermal_cells"].__globals__["quad"]=inherited["thermal_quadrature"]
    orbits=runpy.run_path(str(args.polarization/"orbits/pairs.py"))
    table=runpy.run_path(str(WALL/"evolve.py"))["field_table"]
    fourth_rank=runpy.run_path(str(HERE.parent/"triad/diagnose.py"))["fourth_rank_test"]
    tasks=[]
    for time,steps in ((6.30944,384),(6.45,640),(6.67,1024)):
        for power,seed,multiplier in ((11,18431,1),(11,70163,1),(12,91283,2)):
            tasks.append((time,power,steps*multiplier,seed,2,False))
        tasks.extend((time,11,steps,18431,mode,clone) for mode,clone in ((1,False),(0,False),(0,True)))
    if args.extended:
        tasks=[(7.05,11,2048,18431,2,False),(7.05,11,2048,70163,2,False),
               (7.05,12,4096,91283,2,False),(7.05,11,2048,18431,1,False),
               (7.05,11,2048,18431,0,False)]
    if args.single:
        t,p,n,s,m,c = args.single
        tasks=[(float(t),int(p),int(n),int(s),int(m),bool(int(c)))]
    if any(t<history["t"][0] or t>history["t"][-1] for t,*rest in tasks):
        raise ValueError("Requested marker outside verified input trajectory")
    with tempfile.TemporaryDirectory(prefix="shared-ensemble-") as td:
        work=Path(td)
        context=(history,orbits,coarse,inherited,kinetics,table,fourth_rank,build(args.polarization,args.scattering,work),work)
        with ThreadPoolExecutor(args.jobs) as pool:
            results=list(pool.map(lambda job:simulate(context,job),tasks))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps({"runs":results},sort_keys=True,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("history","polarization","scattering","output"):
        parser.add_argument("--"+name,type=Path,required=True)
    parser.add_argument("--jobs",type=int,default=2)
    parser.add_argument("--extended",action="store_true",help="Run the late 7.05 marker suite")
    parser.add_argument("--single",nargs=6,metavar=("TIME","POWER","STEPS","SEED","MODE","IDENTICAL"))
    run(parser.parse_args())
