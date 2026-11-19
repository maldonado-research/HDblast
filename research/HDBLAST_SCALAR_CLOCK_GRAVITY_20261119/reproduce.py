"""Authenticate the offline capsule, replay released observations or regenerate all native runs."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

HERE = Path(__file__).resolve().parent
BASE = "ec29c348fdfd11248630b80d4f4cee06a6f0026d"
NATIVE = Path("research/HDBLAST_COVARIANCE_FEEDBACK_20261119/native")
GRIDS = ("coarse", "fine")
ARMS = ("quiet", "plus", "minus", "far")
ENV = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
           PYTHONDONTWRITEBYTECODE="1")


def require(condition, label):
    if not condition:
        raise ValueError(label)


def read(path):
    return json.loads(path.read_bytes())


def sha(path):
    with path.open("rb") as stream:
        digest = hashlib.sha256()
        for block in iter(lambda: stream.read(1024*1024), b""):
            digest.update(block)
        return digest.hexdigest()


def authenticate(repo):
    root_manifest = HERE/"SHA256SUMS"
    covered = set()
    for manifest in (root_manifest, *sorted(p for p in HERE.rglob("SHA256SUMS") if p != root_manifest)):
        listed = set()
        for line in manifest.read_text().splitlines():
            digest, name = line.split("  ", 1)
            path = (manifest.parent/name).resolve()
            require(path.is_relative_to(HERE) and path != manifest and path not in listed, name)
            require(sha(path) == digest, f"Hash mismatch in {manifest.name}: {name}")
            listed.add(path)
        covered.update(listed)
        covered.add(manifest)
    files = {p.resolve() for p in HERE.rglob("*") if p.is_file() and "__pycache__" not in p.parts}
    require(files == covered, f"Unmanifested or missing files: {sorted(str(p) for p in files ^ covered)[:10]}")
    declared = read(HERE/"SOURCES.json")
    require(declared["base"] == BASE, "Historical source revision")
    pinned = {}
    for item in declared["sources"]:
        path = (repo/item["path"]).resolve()
        require(path.is_relative_to(repo) and path.exists(), f"Missing source {path}")
        require(path.stat().st_size == item["bytes"] and sha(path) == item["sha256"], str(path))
        require(path.name not in pinned, f"Duplicate upstream source {path.name}")
        pinned[path.name] = path
    names = {"clock_test.py","budget.py","rolloff5d_v1.py","evolve_b1.py","shell_matter.py","static_w.py",
             "FORK.json","WALLS.zip","FORK_INPUTS.zip","rolling_data.py","quantum_shell.py"}
    require(set(pinned) == names, "Historical input inventory")
    vend = HERE/"bulk/vendor"
    for key in names-{"clock_test.py","budget.py","FORK.json","WALLS.zip"}:
        require((vend/key).read_bytes() == pinned[key].read_bytes(), f"Vendored upstream source {key}")
    fork = read(pinned["FORK.json"])
    refs = read(vend/"refs.json")
    require(refs == {g:fork[g]["reference_A"] for g in GRIDS}, "Reference scale factors")
    with zipfile.ZipFile(vend/"references.zip") as released, zipfile.ZipFile(pinned["WALLS.zip"]) as original:
        require(set(released.namelist()) == {f"{g}_{a}.csv" for g in GRIDS for a in ARMS}, "Retained wall inventory")
        for g in GRIDS:
            for a in ARMS:
                require(released.read(f"{g}_{a}.csv") == original.read(f"{g}/retain/{a}/wall.csv"),
                        f"Historical wall {g}/{a}")
    print(f"Authenticated {len(files)} capsule files and {len(pinned)} historical sources.", flush=True)


def run(script, *args, json_output=None):
    cmd = [sys.executable, "-IB", str(HERE/script), *map(str, args)]
    if json_output is None:
        subprocess.run(cmd, check=True, env=ENV)
    else:
        with json_output.open("w") as out:
            subprocess.run(cmd, check=True, env=ENV, stdout=out)


def numbers(actual, expected, path=(), *, relative=1e-7, absolute=5e-10, tolerance=None):
    if isinstance(expected, dict):
        require(isinstance(actual, dict) and actual.keys() == expected.keys(), f"Keys {path}")
        for k in expected:
            numbers(actual[k], expected[k], (*path,k), relative=relative, absolute=absolute, tolerance=tolerance)
    elif isinstance(expected, list):
        require(isinstance(actual, list) and len(actual) == len(expected), f"Length {path}")
        for i, (x,y) in enumerate(zip(actual,expected,strict=True)):
            numbers(x,y,(*path,str(i)),relative=relative,absolute=absolute,tolerance=tolerance)
    elif isinstance(expected, float):
        rel, gap = tolerance(path, relative, absolute) if tolerance else (relative, absolute)
        require(math.isclose(actual,expected,rel_tol=rel,abs_tol=gap),
                f"Numerical disagreement {'/'.join(path)}: obtained {actual}, expected {expected}, atol={gap}, rtol={rel}")
    else:
        require(actual == expected, f"Value {path}: {actual} / {expected}")


def authenticate_tables(parent, native, archived):
    manifest = read(parent/"SOURCES.json")
    expected = read(HERE/"source/SOURCES.json")
    if archived:
        require(manifest == expected, "Archived source manifest")
    require(manifest["native_inputs"] == {n:sha(native/n) for n in ("FORK.json","WALLS.zip")}, "Source ancestry")
    expect_cases = {f"{g}_{a}_{d}" for g in GRIDS for a in ARMS for d in (60,120)}
    require(set(manifest["runs"]) == expect_cases, "Source case inventory")
    seen = set()
    for case,row in manifest["runs"].items():
        require(set(row["files"]) == {"bare","static"}, f"Prescriptions {case}")
        for scheme,meta in row["files"].items():
            name = f"tables/{case}_{scheme}.npz"
            require(meta["path"] == name and sha(parent/name) == meta["sha256"], f"Source file {name}")
            seen.add((parent/name).resolve())
            numbers(meta["readout"],expected["runs"][case]["files"][scheme]["readout"],
                    (case,scheme),relative=2e-8,absolute=2e-8)
    require({f.resolve() for f in (parent/"tables").glob("*.npz")} == seen, "Source table coverage")


def check_csv(result_dir, regenerated):
    for filename in ("READOUTS","HISTORIES","DIAGNOSTICS"):
        def index(folder):
            with (folder/(filename+".csv")).open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            tags = ("grid","arm","cutoff","scheme","step") + (("s",) if filename == "HISTORIES" else ())
            data = {tuple(r[k] for k in tags):r for r in rows}
            require(len(data) == len(rows), f"Duplicate case in {filename}")
            return data,tags
        current,ids = index(result_dir)
        reference,_ = index(HERE/"bulk")
        require(current.keys() == reference.keys(), filename)
        for tag, row in current.items():
            match = reference[tag]
            require(row.keys() == match.keys(), f"CSV columns {filename}")
            for key in set(row)-set(ids):
                x,y = float(row[key]),float(match[key])
                if regenerated:
                    rel,ab = 2e-5,2e-8
                    if key in ("d_rho","d_p","d_j","d_E","d_work","d_pressure_work","d_Wy_over_H02"):
                        ab = 4e-5
                    elif key in ("physical_h_over_h","d_phi","d_A","d_a","d_a_over_a"):
                        ab = 3e-9 if key!="physical_h_over_h" else 4e-10
                    elif key.startswith("ledger_max"):
                        ab = 2e-7
                    elif key.startswith("d_h_minus_quiet") or key.startswith("d_h_over_h_minus_quiet"):
                        if abs(y)>2e-8:
                            require(x*y>0, f"Robust pulse contrast sign {tag} {key}")
                    elif key == "preturn_max_abs":
                        ab = 0.
                else:
                    rel,ab = 4e-7,2e-10
                    if key in ("d_rho","d_p","d_j","d_E","d_work","d_pressure_work"):
                        ab = 2e-8
                require(math.isclose(x,y,rel_tol=rel,abs_tol=ab), f"{filename} {tag} {key}: {x} / {y}")
    checks=read(result_dir/"BASELINE_CHECK.json")
    require(set(checks)=={g+"_"+a for g in GRIDS for a in ARMS}, "Eight baseline checks")
    if not regenerated:
        require(max(max(fields.values()) for fields in checks.values())<=5e-9, "Zero-source archival agreement")
    # For new integrations, bulk/summarize.py checks every field against a
    # geometry-, matter-, or constraint-specific tolerance before it writes.


def compare_histories(new, archived):
    import numpy as np
    from scipy.interpolate import CubicHermiteSpline, CubicSpline
    for file in sorted((archived/"tables").glob("*.npz")):
        with np.load(file,allow_pickle=False) as old, np.load(new/"tables"/file.name,allow_pickle=False) as current:
            sample = np.linspace(max(old["s"][0],current["s"][0]),2.10,1025)
            for k in ("rho","p","j","rad"):
                curves = [CubicHermiteSpline(d["s"],d[k],d["d"+k])(sample) for d in (old,current)]
                tolerance = 2e-8*max(1.,np.max(np.abs(curves[0])))
                np.testing.assert_allclose(curves[0],curves[1],atol=tolerance,rtol=2e-8,err_msg=file.name+"/"+k)
    old_runs = {p.name for p in (archived/"runs").glob("*.npz")}
    require(len(old_runs)==92 and old_runs == {p.name for p in (new/"runs").glob("*.npz")}, "Native run inventory")
    sample = np.linspace(1.31,2.10,257)
    for name in sorted(old_runs):
        with np.load(archived/"runs"/name,allow_pickle=False) as old, np.load(new/"runs"/name,allow_pickle=False) as now:
            for key,gap in (("h",3e-9),("A",5e-10),("B",5e-10),("phi",5e-10),("rho",2e-6),("p",4e-6),("j",1e-5)):
                pair=[CubicSpline(d["s"],d[key])(sample) for d in (old,now)]
                np.testing.assert_allclose(*pair,rtol=3e-9,atol=gap,err_msg=name+"/"+key)


def science(folder, native, jobs, cache, regenerated):
    tables,runs = folder/"tables",folder/"runs"
    print("Checking recomputed scientific observables." if regenerated else "Checking released histories and independent scientific diagnostics.", flush=True)
    run("source/verify.py","--native",native,"--tables",tables,"--output",folder/"verification.json")
    # verify.py checks Ward residuals both at and between knots, causal/spectral agreement and source component sums.
    actual=read(folder/"verification.json")
    numbers(actual,read(HERE/"source/VERIFICATION.json"),relative=1e-3,absolute=2e-8)
    run("theory/onset.py",native,json_output=folder/"onset.json")
    numbers(read(folder/"onset.json"),read(HERE/"theory/ONSET.json"),relative=3e-6,absolute=4e-9)
    run("theory/onset.py",native,"--tables",tables,json_output=folder/"injected.json")
    numbers(read(folder/"injected.json"),read(HERE/"theory/INJECTED.json"),relative=2e-4,absolute=2e-8)
    run("bulk/summarize.py","--runs",runs,"--output",folder,*(["--regenerated"] if regenerated else []))
    check_csv(folder,regenerated)
    # A new source file may have different bytes. The native runs must identify the tables they actually used.
    run("bulk/verify.py","--runs",runs,"--tables",tables,"--results",folder)
    for script,reference,options in (
        ("audit/audit.py","audit/AUDIT.json",[runs]),
        ("audit/stages.py","audit/STAGES.json",["--bulk",HERE/"bulk","--cache",cache]),
        ("radiation/radiation.py","radiation/RADIATION.json",[runs]),
        ("comparator/canonical.py","comparator/CANONICAL.json",["--runs",runs]),
        ("calibration/calibration.py","calibration/CALIBRATION.json",
         ["--bulk",HERE/"bulk","--tables",tables,"--runs",runs,"--cache",cache,"--jobs",jobs])):
        output=folder/Path(reference).name
        run(script,*options,"--output",output)
        def tolerances(path,rel,gap):
            if "zero_slope_kappa" in path:
                return 0.,.005
            if any(k in path for k in ("dh_dlambda","held_out_error","predicted_dh_kappa3")):
                return 0.,2e-8
            if "full_to_direct_h2_slope" in path:
                return 4e-4,3e-5
            if "h2_response_over_h02_per_lambda" in path:
                return 1e-4,4e-8
            if reference.endswith("/AUDIT.json") and "derivative_per_y_squared" in path:
                # Normalized signed Hamiltonian residuals can amplify roundoff at the fine wall.
                return 1e-3,3e-8 if regenerated else 2e-9
            return rel,gap
        numbers(read(output),read(HERE/reference),(reference,),relative=3e-6,
                absolute=2e-8 if regenerated else 2e-9,tolerance=tolerances)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument("--verify",action="store_true",help="Decode released arrays and rerun independent checks")
    modes.add_argument("--regenerate",action="store_true",help="Also regenerate 32 source tables and all 92 native evolutions")
    parser.add_argument("--jobs",type=int,default=2)
    parser.add_argument("--repo",type=Path,default=HERE.parents[1],help="Root of the offline HDblast repository")
    args=parser.parse_args()
    require(args.jobs>=1,"Positive worker count required")
    repo=args.repo.resolve()
    authenticate(repo)
    if not (args.verify or args.regenerate):
        return
    native=repo/NATIVE
    with tempfile.TemporaryDirectory(prefix="hdblast-live-") as temp:
        root=Path(temp)
        cache=root/"cache"
        old=root/"released"
        run("source/archive.py","extract","--archive",HERE/"source/TABLES.zip","--output",old)
        authenticate_tables(old,native,True)
        run("bulk/pack_runs.py","unpack","--archive",HERE/"bulk/RUNS.zip","--runs",old/"runs")
        require(len(list((old/"runs").glob("*.npz")))==92,"Expected 92 released native histories")
        science(old,native,args.jobs,cache,False)
        if args.regenerate:
            new=root/"regenerated"
            run("source/source.py","--native",native,"--output",new/"tables","--jobs",args.jobs)
            authenticate_tables(new,native,False)
            run("bulk/ensemble.py","--tables",new/"tables","--output",new/"runs","--cache",cache,"--jobs",args.jobs)
            compare_histories(new,old)
            science(new,native,args.jobs,cache,True)
    print("Regenerated and verified the complete live response." if args.regenerate else "Verified released arrays, scientific results and held-out integrations.",flush=True)


if __name__=="__main__":
    main()
