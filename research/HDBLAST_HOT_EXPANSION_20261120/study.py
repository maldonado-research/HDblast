"""Authenticate the HDblast hot-expansion archive or recompute each numerical product."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
PRODUCTS = ("exit/EXIT","exit/PLATEAU","axion/GROWTH","axion/HANDOFF","radius/RADIUS","phase/PHASE",
            "criterion/CRITERION","helicity/HELICITY","helicity/EVEN","occupation/OCCUPATION")
PRIOR = {"DECAY.json":"39e3af48a019c3d0bc3d7d2f7818b0e8a68f0a260141b0b15fd56e518baf6347",
         "decay.py":"133ccf1f41360d4499cf9c7a98b38fcdca57af994c060602aeefec1acca526f4",
         "decay_forward.cpp":"e3e20e80d7558895ec189f46e0989c577acccdd81b46828b39aeb463577be611"}

def require(test, message):
    if not test:
        raise ValueError(message)

def authenticate():
    manifest=ROOT/"SHA256SUMS"
    files={p.relative_to(ROOT).as_posix():p for p in ROOT.rglob("*")
           if p.is_file() and p!=manifest and "__pycache__" not in p.parts}
    entries={}
    for line in manifest.read_text().splitlines():
        digest,name=line.split("  ",1)
        require(name in files and name not in entries and not files[name].is_symlink(),
                "Unexpected manifest path "+name)
        entries[name]=digest
    require(files.keys()==entries.keys(),"File inventory mismatch")
    for name,path in files.items():
        require(hashlib.sha256(path.read_bytes()).hexdigest()==entries[name],"Checksum mismatch "+name)
    require(all(entries["prior/"+n]==d for n,d in PRIOR.items()),"Predecessor provenance mismatch")
    print(f"Verified {len(entries)} files and all three predecessor objects.",flush=True)

def compare(actual, expected, path):
    if isinstance(expected,dict):
        require(isinstance(actual,dict) and actual.keys()==expected.keys(),"Changed keys "+path)
        for key in expected:
            compare(actual[key],expected[key],path+"/"+key)
    elif isinstance(expected,list):
        require(isinstance(actual,list) and len(actual)==len(expected),"Changed length "+path)
        for j,(a,b) in enumerate(zip(actual,expected)):
            compare(a,b,path+"/"+str(j))
    elif type(expected) is float:
        key=path.rsplit("/",1)[-1]
        # Subtracting two nearly equal enclosed masses on an extended large-radius
        # trajectory is sensitive to compiler arithmetic. Check the primary constraint
        # and conservation residuals separately below.
        absolute={"tube_mass_error":2e-5,"end_tube_mass":1.5e-5,
                  "surface_balance_error":8e-8,"sheet_balance_error":8e-8,
                  "peak_daughter_velocity":5e-9}.get(key,2e-9)
        require(isinstance(actual,(float,int)) and math.isfinite(actual) and
                math.isclose(actual,expected,abs_tol=absolute,rel_tol=1e-6),
                f"{path}: {actual} versus {expected}")
    else:
        require(type(actual) is type(expected) and actual==expected,"Changed value "+path)

def check(p):
    even=p["exit/EXIT"]["runs"]; odd=p["axion/GROWTH"]["runs"]
    e=even["primary"]; o=odd["primary"]
    for r in (e,o):
        require(0 < r["bounce_tau"] < r["exit_complete_tau"] < r["first_ninety_intrinsic_time"]
                < r["maximum_tau"] and r["end_H"]<0 and r["end_radiation"]==0,
                "Missing ordered events on a zero-residual-radiation run")
        require(min(r[k] for k in ("min_bulk_mass","min_lapse","min_radial_gradient_norm"))>0
                and r["friedmann_error"]<1e-8 and r["gas_balance_error"]<1e-8
                and r["tube_mass_error"]<1e-4 and r["sheet_balance_error"]<3e-6,
                "Primary constraint, admissibility or energy balance failed")
    require(e["ninety_intrinsic_intervals"]==1 and 3.33<e["last_ninety_intrinsic_radius"]/
            e["first_ninety_intrinsic_radius"]<3.35 and e["doubling_intrinsic_fraction"]>.985 and
            abs(e["doubling_pressure_over_gas_density"]-1/3)<1e-8, "Two-stage doubling failed")
    require(2.68<o["hot_growth"]<2.69 and o["hot_from_tau"]<o["first_double_tau"]<o["hot_to_tau"]
            and abs(o["hot_growth"]-odd["displacement"]["hot_growth"])<2e-5
            and abs(o["hot_growth"]-odd["late_refined"]["hot_growth"])<2e-4,
            "Parity-odd doubling or displacement/grid check failed")
    require(odd["even_energy_matched"]["maximum_intrinsic_fraction"]<.54 and
            odd["no_drag"]["maximum_intrinsic_fraction"]<.03 and
            "exit_complete_tau" not in odd["dose_plus"],"Physical controls changed")
    require(len(p["axion/HANDOFF"]["rows"])==1025,"Incomplete scalar transition")
    rad=p["radius/RADIUS"]["runs"]["larger_charge_finer_mesh"]
    require(1.037<math.exp(rad["ninety_intrinsic_expanding_efolds"])<1.039,
            "Calibrated larger-wall control changed")
    c=p["criterion/CRITERION"]
    require(c["previous"]["subsequent_radius_factor_upper_if_hot"]<1.01 and
            2.34<c["two_stage"]["stress_diagnostic"]["longest_hot_saved_span_by_total_pressure_tolerance"]
            ["0.05"]["radius_factor"]<2.37,"Curvature or radiation-pressure test failed")
    require(1.70<p["phase/PHASE"]["limiting_low_energy_radiation"]["max_growth_after_pure_magnetic_90"]<1.72
            and p["helicity/HELICITY"]["modes"]["8192"]["-1"]["max_energy_gain"]>1e22
            and 1e7<p["helicity/EVEN"]["modes"]["4"]["combined"]["max_energy_gain"]<1.2e7,
            "Passive phase or angular transfer test failed")
    occ=p["occupation/OCCUPATION"]
    require(2.8e29<occ["action_for_odd_selected_shells_equal_coherent"]<3e29
            and 9e7<occ["two_stage"]["action_for_selected_shells_equal_coherent"]<1.1e8,
            "Finite quantum shell sum changed")

def replay(out):
    out=out.resolve()
    require(not out.is_relative_to(ROOT) and
            (not out.exists() or (out.is_dir() and not any(out.iterdir()))),
            "Use an empty output directory outside the study")
    out.mkdir(parents=True,exist_ok=True)
    paths={n:out/(n+".json") for n in PRODUCTS}
    for target in paths.values(): target.parent.mkdir(exist_ok=True)
    def run(script,*args):
        print("Replaying "+script,flush=True)
        subprocess.run([sys.executable,"-I","-B",str(ROOT/script),*map(str,args)],check=True)
    with ThreadPoolExecutor(max_workers=3) as pool:
        a=pool.submit(run,"axion/growth.py","--output",paths["axion/GROWTH"],"--handoff",paths["axion/HANDOFF"],"--jobs","2")
        b=pool.submit(run,"exit/growth.py","--output",paths["exit/EXIT"],"--jobs","2")
        radius=pool.submit(run,"radius/radius.py","--output",paths["radius/RADIUS"],"--jobs","2")
        a.result(); b.result()
        even=pool.submit(run,"helicity/even.py","--output",paths["helicity/EVEN"])
        odd=pool.submit(run,"helicity/helicity.py","--trace",paths["axion/HANDOFF"],"--output",paths["helicity/HELICITY"])
        run("exit/growth.py","--cases",ROOT/"exit/plateau_cases.json","--output",paths["exit/PLATEAU"],"--jobs","2")
        run("exit/verify.py",paths["exit/EXIT"])
        run("criterion/criterion.py","--two-stage",paths["exit/EXIT"],"--dual",paths["axion/GROWTH"],
            "--handoff",paths["axion/HANDOFF"],"--output",paths["criterion/CRITERION"])
        run("phase/phase.py","--output",paths["phase/PHASE"],"--even",paths["exit/EXIT"],"--axion",paths["axion/GROWTH"])
        even.result(); odd.result()
        run("occupation/occupation.py","--helicity",paths["helicity/HELICITY"],"--even",paths["helicity/EVEN"],
            "--handoff",paths["axion/HANDOFF"],"--output",paths["occupation/OCCUPATION"])
        radius.result()
    products={n:json.loads(f.read_bytes()) for n,f in paths.items()}
    check(products)
    for n,obj in products.items(): compare(obj,json.loads((ROOT/(n+".json")).read_bytes()),n)
    print("Ten numerical products reproduced and physical checks passed.",flush=True)

if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify",action="store_true")
    parser.add_argument("--output-dir",type=Path)
    args=parser.parse_args()
    if args.output_dir and not args.verify: parser.error("--output-dir requires --verify")
    authenticate()
    check({n:json.loads((ROOT/(n+".json")).read_bytes()) for n in PRODUCTS})
    if args.verify:
        if args.output_dir: replay(args.output_dir)
        else:
            with tempfile.TemporaryDirectory(prefix="hdblast-hot-") as temp: replay(Path(temp))
