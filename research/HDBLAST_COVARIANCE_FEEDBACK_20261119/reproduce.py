"""Verify the HDblast covariance feedback capsule and rebuild its calculations."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from zipfile import ZipFile


HERE = Path(__file__).resolve().parent
MANIFEST_DIGEST = "7e842a32937fa9890c4ee29bd6db7cfcd4d49d97ce00241b3d024cffc8ee8b6a"
BASE = "838834ff58748b559a14726a3a16b8716051be35"
GRIDS = ("coarse", "fine")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def close(value, target, label, relative=2e-7, absolute=1e-10):
    require(math.isfinite(value) and math.isfinite(target)
            and math.isclose(value, target, rel_tol=relative, abs_tol=absolute),
            f"{label}: {value!r} vs {target!r}")


def load(path):
    return json.loads((HERE/path).read_bytes())


def verify_sources(repository):
    data=(HERE/"SOURCES.json").read_bytes()
    require(hashlib.sha256(data).hexdigest()==MANIFEST_DIGEST, "Historical manifest SHA-256")
    manifest=json.loads(data)
    require(manifest["base"]==BASE and len(manifest["sources"])==13, "Historical revision and count")
    for row in manifest["sources"]:
        path=(repository/row["path"]).resolve()
        require(path.is_relative_to(repository), "Historical path")
        payload=path.read_bytes()
        require(len(payload)==row["bytes"] and hashlib.sha256(payload).hexdigest()==row["sha256"],
                f"Historical input changed: {row['path']}")


def run(script, args=(), *, numerical=False, capture=False):
    options=["-I","-B"]
    if not numerical:
        options.append("-S")
    if sys.flags.optimize:
        options.append("-O")
    command=[sys.executable,*options,str(HERE/script),*map(str,args)]
    env=dict(os.environ, OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1")
    print("Running",script,*args,file=sys.stderr,flush=True)
    result=subprocess.run(command,env=env,check=True,stdout=subprocess.PIPE if capture else None)
    return json.loads(result.stdout) if capture else None


def equal_report(actual, expected, key=""):
    if isinstance(expected, dict):
        require(isinstance(actual,dict) and actual.keys()==expected.keys(),f"Keys at {key}")
        for field in expected:
            equal_report(actual[field],expected[field],f"{key}/{field}")
    elif isinstance(expected, list):
        require(isinstance(actual,list) and len(actual)==len(expected),f"Length at {key}")
        for i,value in enumerate(expected):
            equal_report(actual[i],value,f"{key}/{i}")
    elif isinstance(expected, float):
        # One-sided second derivatives and differences cancelled to fifth order
        # amplify roundoff; these two field-specific tolerances cover that loss.
        if key.endswith("rho_double_prime") and "endpoint_differences" in key:
            close(actual,expected,key,relative=3e-5,absolute=.003)
        elif key.endswith("force_over_quintic"):
            close(actual,expected,key,relative=.0003,absolute=2e-5)
        else:
            close(actual,expected,key,relative=2e-6,absolute=1e-8)
    else:
        require(type(actual) is type(expected) and actual==expected,f"Value at {key}")


def verify_local_physics():
    native=load("native/FORK.json")
    fork=load("initial/FORK.json")
    corner=load("initial/CORNER.json")
    conditioning=load("moments/CONDITIONING.json")
    pulses=load("impulse/FINITE.json")
    with ZipFile(HERE/"native/FORK_INPUTS.zip") as inputs:
        require(set(inputs.namelist())=={g+".npz" for g in GRIDS}, "Native restart inventory")
        for grid in GRIDS:
            name=grid+".npz"
            digest=hashlib.sha256(inputs.read(name)).hexdigest()
            require(digest==fork[grid]["checkpoint_sha256"],"Native restart digest")
            source=native[grid]
            for name,field in (("a","a"),("phi","phi_b"),("s","H0tau"),("lapse","lapse")):
                close(fork[grid][name],source[field],f"{grid} fork {field}")
            states=source["states"]
            for i,key in enumerate(("X","P","C")):
                equal_report(fork[grid]["covariance"][i],states["retain"][key],"Native retained "+key)
            data=corner[grid]
            c=states["retain"]["C"]
            require(min(c)>0 and data["independent_257_min_C"]>.003, "Positive retained covariance")
            close(min(c),data["C_min"],"Lowest native C")
            close(max(c),data["C_max"],"Largest native C")
            require(data["max_relative_same_force_derivative_perturbation"] < 1e-6
                    and data["independent_max_C_error_native_nodes"]<1e-7, "Purity and independent evolution")
            for name, row in data["junction_derivatives"].items():
                state_name = "erase" if name=="dephase" else name
                before,after=states["retain"]["neumann_T"],states[state_name]["neumann_T"]
                close(after[1]-before[1],row["B_TZ_change"],f"{name} differentiated metric junction")
                close(after[2]-before[2],row["phi_TZ_change"],f"{name} differentiated scalar junction")
            smooth=data["smooth_control"]
            require(abs(smooth["delta_C_moment0"])<1e-13 and abs(smooth["delta_C_moment1"])<1e-12
                    and smooth["min_determinant"]>.25, "Matched mixed state")
            m2=smooth["delta_C_moment2"]
            close(smooth["pressure_third"],16*m2/(3*source["a"]**7),"Cubic pressure")
            close(smooth["force_fifth"],32e4*(source["phi_b"]-.5)*m2/source["a"]**7,"Quintic force")
            for key in ("pressure_over_cubic","force_over_quintic"):
                first,second=smooth["curve"][:2]
                require(abs(1-first[key])<abs(1-second[key])<.02,"Small-time response "+key)
            constructed=load("native_smooth/FORK.json")[grid]["states"]
            trajectory=load("native_smooth/SMOOTH.json")[grid]
            b0,b1,norm=smooth["polynomial_coefficients"]
            for q,old,plus,minus in zip(source["q"],c,constructed["smoothplus"]["C"],constructed["smoothminus"]["C"]):
                t=(q/8.)**2
                close(plus,.9*old*(t*t-b1*t-b0)/norm,"Constructed smooth state")
                close(minus,-plus,"Opposite smooth state")
            close(trajectory["moments"]["4"],m2,"Analytic and native fourth-momentum moment")
            for field in ("initial_delta_native_sources","initial_delta_native_T_derivatives",
                          "initial_delta_neumann","initial_delta_neumann_T"):
                require(max(map(abs,trajectory[field]))<1e-12,"Matched initial "+field)
            difference=trajectory["plus_minus_minus"]
            require(difference["H_over_H0"]<0 and difference["p_hat"]>0 and difference["j_hat"]<0,
                    "Source-derivative-matched quiet shell separation")
            for level in conditioning[grid]:
                ratios=level["first_twelve_ratios"]
                require(all(0<=y<=x for x,y in zip(ratios,ratios[1:])) and ratios[0]==1.,"Singular spectrum")
                for cutoff,count,field in ((.01,3,"above_one_percent"),(.001,4,"above_one_per_mille"),
                                           (1e-6,6,"above_one_part_per_million")):
                    require(sum(x>cutoff for x in ratios)==level[field]==count,"Force resolution")
            pulse=pulses[grid]
            require(pulse["max_determinant_drift"]<1e-13
                    and max(map(abs,pulse["ideal_source_changes"].values()))<1e-11,"Neutral canonical shear")
            close(pulse["integrated_C"],data["native_integral_C"],"Impulse input covariance")
            close(pulse["integrated_C_after"],-pulse["integrated_C"],"Impulse reflected mean covariance")
            previous=math.inf
            for test in pulse["finite_pulses"]:
                error=test["relative_cross_covariance_error_vs_ideal"]
                require(abs(test["external_work"])<1e-10 and test["delta_mass_squared"]>0
                        and 0<error<previous,"Finite-duration work and convergence")
                previous=error
            require(previous<.02,"Impulse approximation")


def full_replay(repository):
    cases=[
        ("audit/proper_time.py",[repository],"audit/PROPER.json"),
        ("audit/proper_time.py",[repository,"--nodes",48,"--samples",4097],"audit/REFINED.json"),
        ("initial/corner.py",[repository],"initial/CORNER.json"),
        ("moments/conditioning.py",[repository],"moments/CONDITIONING.json"),
    ]
    with tempfile.TemporaryDirectory() as folder:
        work=Path(folder)
        with ZipFile(HERE/"native/FORK_INPUTS.zip") as inputs:
            for grid in GRIDS:
                (work/(grid+".npz")).write_bytes(inputs.read(grid+".npz"))
        cases.append(("impulse/finite_pulse.py",[work/"coarse.npz",work/"fine.npz",repository],"impulse/FINITE.json"))
        for script,args,reference in cases:
            computed=run(script,args,numerical=True,capture=True)
            equal_report(computed,load(reference),reference)
    run("emission/check.py",[repository,HERE/"native","--full"],numerical=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository",type=Path,default=HERE.parent.parent)
    parser.add_argument("--full",action="store_true",help="Recompute the diagnostics and spectra on archived trajectories")
    parser.add_argument("--native",action="store_true",help="Also rebuild both preparation slices and rerun the native branches")
    parser.add_argument("--jobs",type=int,default=2,help="Simultaneous native jobs")
    args=parser.parse_args()
    require(args.jobs>=1,"Positive worker count required")
    repo=args.repository.resolve()
    verify_sources(repo)
    run("audit/check.py",[HERE/"audit/PROPER.json",HERE/"audit/REFINED.json"])
    equal_report(run("moments/force_jet.py",capture=True),load("moments/FORCE_JET.json"),"Rational force moments")
    verify_local_physics()
    run("native/reproduce.py",[repo])
    run("native_smooth/reproduce.py",[repo])
    run("emission/check.py",[repo,HERE/"native"])
    if args.full or args.native:
        full_replay(repo)
    if args.native:
        run("native/rebuild_forks.py",[repo,"--jobs",args.jobs],numerical=True)
        run("native/reproduce.py",[repo,"--full","--jobs",args.jobs],numerical=True)
        run("native_smooth/reproduce.py",[repo,"--full","--jobs",args.jobs],numerical=True)
    print(f"Verified 13 historical sources, native covariance and junction checks, force moments, "
          f"finite impulse and pair spectra. Stored-wall numerical replay: {args.full or args.native}. "
          f"Native recomputation: {args.native}.")


if __name__=="__main__":
    main()
