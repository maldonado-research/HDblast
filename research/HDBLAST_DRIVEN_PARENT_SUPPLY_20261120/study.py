"""Check the sourced neutral-parent study or independently reproduce its six numerical products."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parent
PRODUCTS=("coupled/PARENT","coupled/NARROW","analysis/COHERENT","replay/TRACE","replay/DRIVER","conformal/FORCED")
HISTORY="5dbed28e158584403a4cd238d392aab15c00480522cf489d1ee0a36f1d550f90"
REFERENCE="489e057448d5194959c84a08e317b3491d08c0b2e41c7eaa1d6a964540ac0901"

def require(test,message):
    if not test: raise ValueError(message)

def authenticate():
    manifest=ROOT/"SHA256SUMS"
    files={p.relative_to(ROOT).as_posix():p for p in ROOT.rglob("*")
           if p.is_file() and p!=manifest and "__pycache__" not in p.parts}
    entries={}
    for line in manifest.read_text().splitlines():
        checksum,name=line.split("  ",1)
        require(name in files and name not in entries and not files[name].is_symlink(),"Invalid manifest entry "+name)
        entries[name]=checksum
    require(entries.keys()==files.keys(),"Unexpected file inventory")
    for name,path in files.items():
        require(hashlib.sha256(path.read_bytes()).hexdigest()==entries[name],"Checksum failed "+name)
    require(entries["prior/CHARGED.json"]==HISTORY and entries["prior/parent_modes.py"]==REFERENCE,
            "The #188 history or #189 reference implementation changed")
    print(f"Verified {len(entries)} files including the two historical inputs.",flush=True)

def check(products):
    c=products["coupled/PARENT"]; runs=c["runs"]; r=runs["primary"]
    require(r["bounce_tau"]<r["drive_end_time"]<r["opening_time"] and
            1.513<r["drive_end_target_ratio"]<1.515 and 1.672<r["opening_target_ratio"]<1.675,
            "Empty-wall transfer or event ordering failed")
    require(r["min_bulk_mass"]>0 and r["min_lapse"]>0 and
            r["friedmann_error"]<1e-8 and r["driver_balance_error"]<1e-7 and
            r["oscillator_balance_error"]<1e-12 and r["tube_mass_error"]<1e-4,
            "Physical or numerical budget failed")
    require(r["opening_drive_work"]>7e1*r["opening_mass_work"] and
            abs(runs["flip_drive"]["opening_target_ratio"]-r["opening_target_ratio"])<1e-8 and
            abs(runs["fine"]["opening_target_ratio"]-r["opening_target_ratio"])<1e-4,
            "Source work, sign symmetry or mesh convergence failed")
    require(runs["no_drive"]["opening_target_ratio"]==0 and "opening_time" not in runs["no_pulse"]
            and "opening_time" not in runs["dose_high"],"A source/trajectory control changed")
    stable=runs["stable"]
    require(stable["end_H"]<0 and stable["constant_mass_half_cycles"]==11 and
            1.34<stable["half_cycles"][-1][3]/(.0008/1.9)<1.36 and
            stable["opening_time"]<stable["maximum_tau"],"Stable continuation failed")
    high_mass=products["coupled/NARROW"]["runs"]
    require(2.57<high_mass["fine"]["opening_target_ratio"]<2.59 and
            high_mass["fine"]["opening_mass_WKB"]<.035 and
            high_mass["wide"]["drive_end_target_ratio"]<1e-14 and
            "opening_time" not in high_mass["high"],"High-mass bandwidth control failed")
    audit=products["analysis/COHERENT"]; interval=audit["source_off_after_gate_expanding"]
    require(audit["opening"]["mass_sweep_over_M2"]>5 and interval["complete_periods"]>=5 and
            1.34<interval["action_over_target_mean"]<1.36 and interval["action_over_target_min"]>1.1,
            "Adiabaticity or pressure/action audit failed")
    trace=products["replay/TRACE"]
    require(trace["checks"]["max_absolute_ledger_residual"]<1e-8 and
            trace["pre_mass_to_opening_phase_independent_bounds"][0]>1.6,"Independent driven trajectory failed")
    drv=products["replay/DRIVER"]["runs"]
    require(1.44<drv["primary"]["parent_action_over_target"]<1.45 and
            drv["zero_displacement"]["parent_action_over_target"]==0,"Finite-driver benchmark failed")
    forced=products["conformal/FORCED"]
    n=forced["cases"]["nominal"]; h=forced["cases"]["half_amplitude"]
    require(forced["history_sha256"]==HISTORY and 1.416<n["ratio_to_fluid_target"]<1.417 and
            n["work_balance_error"]<1e-10 and abs(4*h["ratio_to_fluid_target"]-n["ratio_to_fluid_target"])<1e-8
            and n["source_amplitude_for_target"]>forced["nominal_necessary_amplitude_bound"]>0,
            "Conformal source, quadratic scaling or energy accounting failed")

def compare(actual,expected,path):
    if isinstance(expected,dict):
        require(isinstance(actual,dict) and actual.keys()==expected.keys(),"Changed keys "+path)
        for key in expected: compare(actual[key],expected[key],path+"/"+key)
    elif isinstance(expected,list):
        require(isinstance(actual,list) and len(actual)==len(expected),"Changed array length "+path)
        for i,(a,e) in enumerate(zip(actual,expected)): compare(a,e,path+"/"+str(i))
    elif type(expected) is float:
        key=path.rsplit("/",1)[-1]
        # Mass reconstructed across the bulk collar involves cancellation at the
        # large-radius censor boundary. Nominal conservation bounds are checked above.
        absolute={"tube_mass_error":2e-6,"end_tube_mass":2e-6,
                  "surface_balance_error":4e-8,"sheet_balance_error":4e-8}.get(key,2e-9)
        require(type(actual) in (float,int) and math.isfinite(actual) and
                math.isclose(actual,expected,rel_tol=1e-6,abs_tol=absolute),
                f"{path} has {actual}; expected {expected}")
    else:
        require(type(actual)==type(expected) and actual==expected,"Changed value "+path)

def run_all(out):
    out=out.resolve()
    require(not out.is_relative_to(ROOT) and (not out.exists() or out.is_dir() and not any(out.iterdir())),
            "Use a new or empty output directory outside the study")
    targets={n:out/(n+".json") for n in PRODUCTS}
    for p in targets.values():p.parent.mkdir(parents=True,exist_ok=True)
    def run(script,product,*extra):
        print("Replaying "+product,flush=True)
        subprocess.run([sys.executable,"-I","-B",str(ROOT/script),"--output",str(targets[product]),
                        *map(str,extra)],check=True)
    with ThreadPoolExecutor(max_workers=3) as pool:
        forced=pool.submit(run,"conformal/forced.py","conformal/FORCED")
        driver=pool.submit(run,"replay/driver.py","replay/DRIVER")
        run("coupled/production.py","coupled/PARENT","--jobs","2")
        narrow=pool.submit(run,"coupled/production.py","coupled/NARROW",
                           "--cases",ROOT/"coupled/narrow_cases.json","--jobs","2")
        subprocess.run([sys.executable,"-I","-B",str(ROOT/"coupled/verify.py"),
                        str(targets["coupled/PARENT"])],check=True)
        trace=pool.submit(run,"replay/trace.py","replay/TRACE","--source",targets["coupled/PARENT"])
        run("analysis/coherent.py","analysis/COHERENT","--source",targets["coupled/PARENT"])
        for future in (forced,driver,trace,narrow):future.result()
    p={k:json.loads(v.read_bytes()) for k,v in targets.items()}
    check(p)
    for name,result in p.items():compare(result,json.loads((ROOT/(name+".json")).read_bytes()),name)
    print("All six numerical products and physical tests passed.",flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify",action="store_true")
    parser.add_argument("--output-dir",type=Path)
    args=parser.parse_args()
    if args.output_dir and not args.verify:parser.error("--output-dir requires --verify")
    authenticate()
    check({k:json.loads((ROOT/(k+".json")).read_bytes()) for k in PRODUCTS})
    if args.verify:
        if args.output_dir:run_all(args.output_dir)
        else:
            with tempfile.TemporaryDirectory() as temp:run_all(Path(temp))
