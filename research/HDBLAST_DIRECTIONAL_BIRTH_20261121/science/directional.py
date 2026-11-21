"""Test angularly isotropic sharp-shell pair injection on pinned HDblast histories."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import runpy
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
COLUMNS=runpy.run_path(str(HERE/"history.py"))["STEP_COLUMNS"]
EDGES=[-1.,-.5,-.1,-.01,0.,.01,.1,.5,1.]


def analyze(prior, scratch, histories, jobs):
    exe=scratch/"angular"
    subprocess.run([os.environ.get("CXX","g++"),"-std=c++17","-O3",str(HERE/"angular.cpp"),"-o",str(exe)],check=True)
    previous=json.loads((prior/"kinetics/SWEEP.json").read_bytes())
    specs=[(name,grid,scale) for name,h in histories.items() for grid in h["runs"] for scale in ("action","floor")]
    def one(spec):
        name,grid,scale=spec
        h=histories[name]
        proc=subprocess.run([str(exe),h["runs"][grid]["path"],str(h["reference_rest"]),str(h["mass"]),str(h[scale])],
                            check=True,capture_output=True,text=True)
        windows={}
        for row in proc.stdout.splitlines():
            window,key,value=row.split()
            windows.setdefault(window,{})[key]=float(value)
        for window,d in windows.items():
            ref=previous["runs"][name]["runs"][grid]
            prefix=("extended_" if window=="extended" else "")+("main" if scale=="action" else "floor")+"_sweep_"
            for key,source in [("pooled_peak_demand_over_capacity","peak"),("pooled_denied_over_all_requested","deficit_fraction")]:
                if not math.isclose(d[key],ref[prefix+source],rel_tol=2e-5,abs_tol=1e-9):
                    raise ValueError(f"Prior aggregate disagreement in {spec}, {window}, {key}: {d[key]} vs {ref[prefix+source]}")
            fractions=[d.pop(f"bin_{i}_denied_over_all") for i in range(len(EDGES)-1)]
            if not math.isclose(sum(fractions),d["denied_over_all_requested"],abs_tol=2e-12,rel_tol=2e-9):
                raise ValueError("Angular bands do not sum to whole-sphere result")
            if d["denied_over_all_requested"]+1e-12<d["pooled_denied_over_all_requested"]:
                raise ValueError("Directional clipping should be at least as restrictive as pooling")
            d["angular_bins_denied_over_all"]=fractions
        return name,grid,scale,windows
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        rows=list(pool.map(one,specs))
    results={}
    for n,g,scale,w in rows:
        h=histories[n]
        case=results.setdefault(n,{k:h[k] for k in ("mass","q","reference_rest","action","floor")})
        gr=case.setdefault("runs",{}).setdefault(g,{k:h["runs"][g][k] for k in ("end_rest","hot_factor")})
        gr[scale]=w
    return dict(predecessor="dc99ff5c5894abcff3e11e02e4b97e19da31de00",
                source_sha256={k:hashlib.sha256((prior/k).read_bytes()).hexdigest() for k in (
                    "coupled/joined.py","coupled/joined_forward.cpp","envelope/build.py","envelope/ENVELOPE.json","kinetics/SWEEP.json")},
                particle_direction="positive daughter relative to its electric force",
                angular_cosine_edges=EDGES,delivery_windows={"central":[.01,.99],"extended":[.00001,.99999]},
                runs=results)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prior",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--native-output",type=Path,help="Keep intermediate native steps here for additional analyses.")
    parser.add_argument("--jobs",type=int,default=2)
    args=parser.parse_args()
    if args.jobs<1:parser.error("--jobs must be positive")
    build=runpy.run_path(str(HERE/"history.py"))["build_histories"]
    with tempfile.TemporaryDirectory(prefix="hdblast-directions-") as tmp:
        folder=args.native_output.resolve() if args.native_output else Path(tmp)
        folder.mkdir(parents=True,exist_ok=True)
        history=build(args.prior,folder,args.jobs)
        (folder/"traces.json").write_text(json.dumps(history,indent=2)+"\n")
        data=analyze(args.prior,Path(tmp),history,args.jobs)
    args.output.write_text(json.dumps(data,sort_keys=True,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    main()
