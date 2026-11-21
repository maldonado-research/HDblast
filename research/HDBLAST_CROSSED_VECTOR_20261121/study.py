"""Rebuild and verify the crossed-vector calculations from the pinned wall histories and elastic solver."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np


ROOT=Path(__file__).resolve().parent
RESEARCH=ROOT.parent
AXIAL=RESEARCH/"HDBLAST_AXIAL_CURRENT_20261121"
PHOTON=RESEARCH/"HDBLAST_COLLISION_PHOTONS_20261121"
POLARIZATION=RESEARCH/"HDBLAST_POLARIZATION_20261121"
SCATTERING=RESEARCH/"HDBLAST_SCATTERING_20261121"
PRODUCTS=("crossed-exploration/CROSSED.json","crossed-exploration/POPULATION.json",
          "crossed-stress/STRESS.json","unitarity/PARTIAL.json")
STRESS_INPUTS=("t6.67000-p14-n4096-r91283-particle-pole","t7.05000-p14-n8192-r91283-particle-pole",
               "t6.67000-p13-n4096-r70163-particle-pole","t7.05000-p13-n8192-r70163-particle-pole")


def authenticate():
    for base,manifest in ((RESEARCH.parent,ROOT/"parent-inputs.sha256"),(ROOT,ROOT/"SHA256SUMS")):
        names=set()
        for line in manifest.read_text().splitlines():
            digest,rel=line.split("  ",1)
            path=base/rel
            if rel in names or not path.resolve().is_relative_to(base) or path.is_symlink():
                raise ValueError(f"Duplicate or out-of-tree manifest entry {rel}")
            with path.open("rb") as stream:
                actual=hashlib.file_digest(stream,"sha256").hexdigest()
            if actual!=digest:
                raise ValueError(f"Checksum failure for {rel}")
            names.add(rel)
        if base==ROOT:
            files={p.relative_to(ROOT).as_posix() for p in ROOT.rglob("*") if p.is_file() and p!=manifest}
            if files!=names:
                raise ValueError(f"Scientific file inventory differs from manifest: {names^files}")
        print(f"Verified {len(names)} files in {manifest.name}",flush=True)


def compare(actual,expected,path=""):
    if isinstance(expected,dict):
        if isinstance(actual,dict) and set(actual)==set(expected):
            for key in expected:
                compare(actual[key],expected[key],f"{path}/{key}")
            return
    elif isinstance(expected,list):
        if isinstance(actual,list) and len(actual)==len(expected):
            for i,(v,ref) in enumerate(zip(actual,expected)):
                compare(v,ref,f"{path}/{i}")
            return
    elif type(expected) is float and type(actual) in (int,float):
        # Near-zero residuals describe integration or conservation checks. Physical rates
        # and signed derivatives, including very small ones, must match in relative units.
        residual=any(token in path for token in ("_error","_residual","midpoint_checks"))
        if math.isfinite(actual) and math.isclose(actual,expected,rel_tol=2.e-5,abs_tol=2.e-10 if residual else 0.):
            return
    elif type(actual) is type(expected) and actual==expected:
        return
    raise ValueError(f"{path}: regenerated {actual!r}, reference {expected!r}")


def validate(dir):
    for product in PRODUCTS:
        compare(json.loads((dir/product).read_bytes()),json.loads((ROOT/product).read_bytes()),product)
        print(f"Reproduced {product}",flush=True)


def run(script,*args):
    env=dict(os.environ,OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1")
    print(f"Running {script.parent.name}/{script.name}",flush=True)
    subprocess.run([sys.executable,"-I","-B",str(script),*map(str,args)],env=env,check=True)


def endpoints(events,temp):
    # The instrumented parent writes the final matrix in the same stable birth
    # order as the wall-current snapshots, with the positive track first in each pair.
    expected={r["label"]:r for r in json.loads((PHOTON/"events/EVENTS.json").read_bytes())["runs"]}
    if {p.stem for p in events.glob("*.npz")}!=set(expected):
        raise ValueError("Expected the eight archived dedicated-interaction collision histories")
    for label,ref in expected.items():
        with np.load(events/(label+".npz"),allow_pickle=False) as z:
            p=z["p_final"]
            n=len(p)
            e=np.hypot(1.,np.linalg.norm(p,axis=1))
            if (str(z["label"])!=label or n!=2*ref["pairs"] or len(z["s"])!=ref["events"] or
                    p.shape!=(n,3) or not np.isclose(float(z["gas"]),ref["gas"],rtol=2.e-8) or
                    not np.isclose(float(z["weight"])*float(z["endpoint_a"])*e.sum(),ref["gas"],rtol=2.e-8)):
                raise ValueError(f"Endpoint data disagree with the checked parent run {label}")
            if label in STRESS_INPUTS:
                signs=np.tile([1.,-1.],n//2)
                np.savez(temp/(label+".npz"),p=p,sign=signs)
    return [temp/(label+".npz") for label in STRESS_INPUTS]


def replay(output,jobs,existing_events,existing_history):
    if any(output.is_relative_to(study) for study in (ROOT,AXIAL,SCATTERING,POLARIZATION,PHOTON)):
        raise ValueError("Write regenerated outputs outside the input studies")
    if output.exists() and (not output.is_dir() or any(output.iterdir())):
        raise ValueError("Choose an empty or new output directory")
    for product in PRODUCTS:
        (output/product).parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="hdblast-crossed-vector-") as folder:
        temp=Path(folder)
        events=existing_events if existing_events is not None else temp/"events"
        history=existing_history if existing_history is not None else temp/"dedicated_half.npz"
        if existing_events is None:
            run(AXIAL/"fields/extract.py","--study",POLARIZATION,"--cases","dedicated_half","--end","7.0505",
                "--output-dir",temp,"--output",temp/"fields.json","--jobs",jobs)
            run(PHOTON/"events/replay.py","--axial",AXIAL,"--scattering",SCATTERING,"--polarization",POLARIZATION,
                "--history",history,"--artifacts",events,"--output",temp/"events.json","--jobs",jobs)
            compare(json.loads((temp/"events.json").read_bytes()),json.loads((PHOTON/"events/EVENTS.json").read_bytes()),"events")
        snapshots=endpoints(events,temp)
        fixed=ROOT/"crossed-exploration/check.py"
        run(fixed,"--output",output/PRODUCTS[0])
        run(ROOT/"unitarity/partial_wave.py","--events",events,"--output",output/PRODUCTS[3])
        with ThreadPoolExecutor(max_workers=min(2,jobs)) as pool:
            census=pool.submit(run,ROOT/"crossed-exploration/population.py","--artifacts",events,"--jobs",max(1,jobs-1),
                              "--output",output/PRODUCTS[1])
            stress=pool.submit(run,ROOT/"crossed-stress/generator.py","--snapshots",*snapshots,"--history",history,
                               "--kernel",fixed,"--power",14,"--batches",8,"--polar",10,"--azimuths",8,
                               "--output",output/PRODUCTS[2])
            census.result()
            stress.result()
    validate(output)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group()
    mode.add_argument("--verify",action="store_true",help="Regenerate results from the pinned parents and compare every output")
    mode.add_argument("--validate",metavar="DIRECTORY",type=Path,help="Compare four previously regenerated JSON products")
    mode.add_argument("--authenticate-only",action="store_true",help="Check the study and dependency file manifests (default)")
    parser.add_argument("--output-dir",type=Path,help="Keep regenerated JSON products in a new or empty directory")
    parser.add_argument("--jobs",type=int,default=2)
    parser.add_argument("--events",type=Path,help="Reuse eight instrumented parent collision files")
    parser.add_argument("--history",type=Path,help="Dedicated extended wall history regenerated with the same parent")
    opt=parser.parse_args()
    if opt.jobs<1 or bool(opt.events)!=bool(opt.history) or (not opt.verify and (opt.events or opt.output_dir)):
        parser.error("Use --verify for reconstruction options and supply --events and --history together")
    authenticate()
    if opt.validate:
        validate(opt.validate)
    if opt.verify:
        if opt.output_dir:
            replay(opt.output_dir.resolve(),opt.jobs,opt.events,opt.history)
        else:
            with tempfile.TemporaryDirectory(prefix="checked-crossed-vector-") as folder:
                replay(Path(folder),opt.jobs,opt.events,opt.history)


if __name__=="__main__":
    main()
