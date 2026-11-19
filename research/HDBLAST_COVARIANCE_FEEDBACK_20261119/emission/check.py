"""Check the archived pair experiment and optionally run both spectral solvers."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
STATES={"retain","erase","reflect","kick"}
ARMS={"plus":-1,"minus":1,"far":1}


def check(d):
    if d["proper_cut"]!=2.10 or d["switches"]!={"wide":[1.35,2.05],"narrow":[1.45,1.95]}:
        raise ValueError("Readout definition")
    if set(d["grids"])!={"coarse","fine"}:
        raise ValueError("Both native grids required")
    verified=0
    for grid, data in d["grids"].items():
        if (not 0 <= data["max_independent_difference_over_total"] < 2e-6 or
                max(data["original_vacuum_to_retained_fork_relative_discrepancy"].values()) >= 1e-5 or
                max(data["independent_contrast_error_pp"].values()) >= 2e-4 or
                data["independent_midpoint_error_fraction"] >= 1e-7):
            raise ValueError("Independent replay")
        walls=data["walls"]
        if set(walls)!=STATES:
            raise ValueError("Four preparations required")
        for state,runs in walls.items():
            if set(runs)!={"quiet",*ARMS}:
                raise ValueError("Native pulse arms")
            for arm,case in runs.items():
                for probe,windows in case["states"].items():
                    for w, bands in windows.items():
                        for band, observables in bands.items():
                            for key, values in observables.items():
                                incl=values["total"]
                                if (not incl>0 or values["vacuum"]<0 or
                                   not math.isclose(values["vacuum"]+values["excess"],incl,
                                                    rel_tol=1e-10,abs_tol=1e-13) or
                                   not math.isclose(case["a_final"]**3*values["total_density"],incl,
                                                    rel_tol=1e-10,abs_tol=1e-13)):
                                    raise ValueError(("Spectrum components",grid,state,arm,probe,w,band,key))
            for arm, sign in ARMS.items():
                for w in ("wide","narrow"):
                    for band in ("full","screen"):
                        row=data["contrasts"][state][arm][w][band]
                        for obs in ("number","energy"):
                            for key in ("total","excess","total_density"):
                                fixed=walls["retain"][arm]["states"][state][w][band][obs]
                                quiet=walls["retain"]["quiet"]["states"][state][w][band][obs]
                                active=runs[arm]["states"][state][w][band][obs]
                                active_quiet=runs["quiet"]["states"][state][w][band][obs]
                                norm=quiet["total_density" if key=="total_density" else "total"]
                                contrast=row[obs+"_"+key]
                                expected=sign*(-1 if key=="excess" else 1)
                                for name, val in (("fixed",(fixed[key]-quiet[key])*100/norm),
                                                  ("responding",(active[key]-active_quiet[key])*100/norm),
                                                  ("geometry_effect",contrast["responding"]-contrast["fixed"])):
                                    if (not math.isclose(contrast[name],val,rel_tol=2e-10,abs_tol=1e-12) or
                                            (name!="geometry_effect" and val*expected<=0)):
                                        raise ValueError(("Pulse contrast",grid,state,arm,w,band,obs,key,name))
                                verified+=1
        for arm,windows in data["midpoint_fraction"].items():
            for w,bands in windows.items():
                for band, checks in bands.items():
                    for obs in ("number","energy"):
                        v=checks[obs+"_total"]
                        if abs(v["fixed"])>1e-12 or not (-1e-5<v["responding"]<-2e-7):
                            raise ValueError(("Covariance-midpoint diagnostic",grid,arm,w,band,obs,v))
    if verified!=576:
        raise ValueError("Pulse coverage")
    return verified


def match(a,b,path=()):
    if isinstance(a,dict):
        if not isinstance(b,dict) or a.keys()!=b.keys():
            raise ValueError(("Keys",path))
        for k,v in a.items():
            match(v,b[k],(*path,k))
    elif isinstance(a,list):
        if not isinstance(b,list) or len(a)!=len(b):
            raise ValueError(("List",path))
        for i,v in enumerate(a):
            match(v,b[i],(*path,str(i)))
    elif isinstance(a,(int,float)):
        # The released data use finite mesh histories. Replaying the native
        # solver can move small differences of pulse counts by roundoff.
        absolute = 4e-7 if "contrasts" in path or "independent_contrast_error_pp" in path else 5e-9
        if not math.isclose(a,b,rel_tol=3e-6,abs_tol=absolute):
            raise ValueError(("Number",path,a,b))
    elif a!=b:
        raise ValueError(("Value",path,a,b))


def check_supplement(data,split,main_data):
    for grid in ("fine","coarse"):
        d=data[grid]
        tab=d["normalized_differences"]["wide"]["full"]["number_total"]
        if not (-.0015 < tab["state_difference_responding"] < -.0014
                and 1.5e-6 < tab["geometry_difference"] < 2e-6
                and abs(tab["midpoint_fixed"]) < 1e-12 and d["max_difference_between_methods"] < 1e-9):
            raise ValueError(("Source-derivative-matched readout",grid))
        for state in ("smoothplus","smoothminus","erase"):
            if d["values"][state][state]["wide"]["full"]["number"]["total"]<=0:
                raise ValueError("Supplement counts must be positive")
        original=main_data["grids"][grid]["walls"]["erase"]["quiet"]["states"]["erase"]
        erase=d["values"]["erase"]["erase"]["wide"]["full"]["number"]["total"]
        if not math.isclose(erase,original["wide"]["full"]["number"]["total"],rel_tol=1e-12):
            raise ValueError("Different erased-state readout in supplement")
        z=split[grid]
        control=z["symmetric_contributions"]["number"]["total"]
        if not (8e-7 < control["scale_factor"] < 1.1e-6
                and -2.8e-5 < control["mass"] < -2.5e-5 and z["max_method_difference"] < 1e-9):
            raise ValueError(("Controlled mass/expansion comparison",grid))
        for obs,row in z["symmetric_contributions"].items():
            for label,v in row.items():
                if not math.isclose(v["scale_factor"]+v["mass"],v["total"], rel_tol=1e-7,abs_tol=2e-14):
                    raise ValueError(("Four-corner identity",grid,obs,label))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("repository",type=Path)
    ap.add_argument("native",type=Path)
    ap.add_argument("--full",action="store_true")
    args=ap.parse_args()
    data=json.loads((HERE/"EMISSION.json").read_bytes())
    num=check(data)
    smooth=json.loads((HERE/"SMOOTH.json").read_bytes())
    split=json.loads((HERE/"BACKGROUND_SPLIT.json").read_bytes())
    check_supplement(smooth,split,data)
    if args.full:
        env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1")
        fresh=json.loads(subprocess.check_output([sys.executable,"-I","-B",str(HERE/"emission.py"),
                                                  str(args.repository),str(args.native)],env=env))
        check(fresh)
        match(data,fresh)
        rebuilt=[]
        for script,extra,saved in (("smooth_emission.py",[str(args.native.parent/"native_smooth")],smooth),
                                   ("background_split.py",[],split)):
            rows=json.loads(subprocess.check_output([sys.executable,"-I","-B",str(HERE/script),
                        str(args.repository),str(args.native),*extra],env=env))
            match(saved,rows)
            rebuilt.append(rows)
        check_supplement(*rebuilt,fresh)
    print(f"Verified {num} fixed and responding pulse comparisons; checked common-wall affine cancellation, matched derivatives and mass/scale-factor contributions. Full replay: {args.full}.")


if __name__=="__main__":
    main()
