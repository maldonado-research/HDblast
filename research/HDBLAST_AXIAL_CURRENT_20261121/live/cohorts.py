"""Measure the accepted contact distribution and its partition of wall-frame stress."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import tempfile

import numpy as np
from scipy.stats import spearmanr


HERE = Path(__file__).resolve().parent
wall = runpy.run_path(str(HERE/"evolve.py"))


def partition(groups, contacts, energy, axial_stress, birth):
    result = []
    for name, mask in groups:
        n = int(mask.sum())
        result.append(dict(group=name, count=n, particle_fraction=float(n/len(mask)),
                           energy_fraction=float(energy[mask].sum()/energy.sum()),
                           axial_stress_fraction=float(axial_stress[mask].sum()/axial_stress.sum()),
                           mean_contacts=float(contacts[mask].mean()) if n else None,
                           zero_fraction=float(np.mean(contacts[mask]==0)) if n else None,
                           fraction_below_eight=float(np.mean(contacts[mask]<8)) if n else None,
                           mean_final_gamma=float(energy[mask].mean()) if n else None,
                           mean_birth_time=float(birth[mask].mean()) if n else None))
    return result


def summarize(row, sample):
    contacts=sample["contacts"]
    birth=sample["birth"]
    p=sample["p"]
    energy=np.hypot(1., np.linalg.norm(p,axis=1))
    axial=p[:,2]**2/energy
    if (int(contacts.sum())!=2*row["events"] or
            not np.isclose(np.mean(contacts>0),row["hit"],rtol=0.,atol=1e-16)):
        raise ValueError("Particle counters disagree with the wall ledger")
    ranges=[("0",0,1),("1-3",1,4),("4-7",4,8),("8-15",8,16),("16+",16,np.inf)]
    subsets=[(name,(contacts>=left)&(contacts<right)) for name,left,right in ranges]
    def quartiles(x):
        cuts=np.quantile(x,[.25,.5,.75])
        entries=np.searchsorted(cuts,x,side="right")
        return [(str(i+1),entries==i) for i in range(4)]
    return dict(label=row["label"], time=row["time"],pairs=row["pairs"],steps=row["steps"],seed=row["seed"],
                events=row["events"],mean_encounters=row["mean_encounters"],hit=row["hit"],
                mean_contacts_energy_weighted=float(np.average(contacts,weights=energy)),
                mean_contacts_axial_stress_weighted=float(np.average(contacts,weights=axial)),
                contacts_quantiles_nearest_rank=np.quantile(contacts,[.1,.25,.5,.75,.9,.99],method="inverted_cdf").astype(int).tolist(),
                count_vs_age_spearman=float(spearmanr(contacts,row["time"]-birth).statistic),
                count_vs_final_energy_spearman=float(spearmanr(contacts,energy).statistic),
                by_encounters=partition(subsets,contacts,energy,axial,birth),
                by_birth_quartile_oldest_first=partition(quartiles(birth),contacts,energy,axial,birth),
                by_final_energy_quartile_lowest_first=partition(quartiles(energy),contacts,energy,axial,birth))


def run(opts):
    with np.load(opts.history) as data:
        history=dict(data)
    if history["t"][-1]<7.05 or not np.allclose(history["K"],2.54,rtol=0.,atol=2.e-12):
        raise ValueError("Need a constant-K history covering the endpoints")
    orbits=runpy.run_path(str(opts.polarization/"orbits/pairs.py"))
    prior=runpy.run_path(str(opts.scattering/"live/evolve.py"))
    coarse=runpy.run_path(str(opts.polarization/"orbits/coarse.py"))
    coarse["thermal_cells"].__globals__["quad"]=prior["thermal_quadrature"]
    reference=[row for name in ("LIVE","CROSSINGS") for row in json.loads((HERE/(name+".json")).read_text())["runs"]
               if row["pole"] and row["mode"]=="particle" and row["time"] in (6.67,7.05)]
    with tempfile.TemporaryDirectory(prefix="tagged-encounters-") as folder:
        work=Path(folder)
        snapshots=work/"snapshots"; snapshots.mkdir()
        binary=wall["build"](opts.polarization,opts.scattering,work,record_contacts=True)
        context=history,orbits,coarse,prior,binary
        def one(ref):
            task=(ref["time"],ref["pairs"].bit_length()-1,ref["steps"],ref["seed"],"particle",True)
            row=wall["run_case"](context,task,work,snapshots,record_contacts=True)
            if row["events"]!=ref["events"] or not np.isclose(row["gas"],ref["gas"],rtol=1e-7):
                raise ValueError("Instrumented replay differs from recorded calculation")
            with np.load(snapshots/(row["label"]+".npz")) as sample:
                return summarize(row,sample)
        with ThreadPoolExecutor(opts.jobs) as pool:
            rows=list(pool.map(one,reference))
    opts.output.parent.mkdir(exist_ok=True,parents=True)
    opts.output.write_text(json.dumps({"runs":rows},indent=2,sort_keys=True,allow_nan=False)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ("polarization","scattering","history","output"):
        parser.add_argument("--"+key,type=Path,required=True)
    parser.add_argument("--jobs",type=int,default=2)
    run(parser.parse_args())
