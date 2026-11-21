"""Tabulate low-drift, hot gas temperatures and maximal required elastic rate by temperature bin."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import runpy
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent


def replay(base, names):
    envelope=json.loads((base/"envelope/ENVELOPE.json").read_text())
    transform=runpy.run_path(str(base/"envelope/build.py"))["transform"]
    runner=runpy.run_path(str(base/"coupled/joined.py"))["run"]
    source=transform((base/"coupled/joined_forward.cpp").read_text(), centered=True)
    signature=hashlib.sha256(source.encode()).hexdigest()
    if signature!=envelope["centered_solver_sha256"]:
        raise ValueError("Wrong parent solver")
    anchor="        const double margin=gas.intrinsic_random-.9*wnew.e;"
    if source.count(anchor)!=1:
        raise ValueError("Incompatible solver")
    with tempfile.TemporaryDirectory(prefix="hdblast-resonance-") as scratch:
        cpp=Path(scratch)/"measure.cpp"; exe=cpp.with_suffix("")
        cpp.write_text(source.replace(anchor,(HERE/"thermal.inc").read_text()+anchor))
        subprocess.run(["g++","-O3","-std=c++17",str(cpp),"-o",str(exe)],check=True)
        def measure(name):
            params=envelope["parameters"][name]
            result=runner(exe,dict(params,sample=0))
            reference=envelope["runs"][name]
            for key in ("end_rest","hot_efolds","low_drift_drag_floor"):
                if not math.isclose(result[key],reference[key],rel_tol=3e-8):
                    raise ValueError("Readout changed "+key)
            bins={}
            for k,val in result.items():
                if k.startswith("resonance_bin_"):
                    index,metric=k[len("resonance_bin_"):].split("_",1)
                    bins.setdefault(index,{})[metric]=val
            rows=[bins[i] for i in sorted(bins,key=int)]
            return dict(mu=params['daughter_mass'],desired_drag=params['drag'],bins=rows,
                        maximum_action_floor=max(r["floor"] for r in rows),
                        duration=sum(r["duration"] for r in rows))
        with ThreadPoolExecutor(len(names)) as pool:
            measured=dict(zip(names,pool.map(measure,names)))
    return dict(parent_sha256=signature,runs=measured)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--base",type=Path,required=True)
    cli.add_argument("--output",type=Path,required=True)
    args=cli.parse_args()
    obj=replay(args.base,["mass_4e12","mass_4e12_half_step"])
    args.output.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n")
