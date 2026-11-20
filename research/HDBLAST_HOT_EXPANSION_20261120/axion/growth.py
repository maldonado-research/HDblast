"""Replay the wall with a transient theta F dual(F) interaction and a charged gas."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import runpy
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
prior = runpy.run_path(str(HERE.parent / "prior/decay.py"))
COLUMNS = prior["COLUMNS"] + ["comoving_odd_to_EM"]
HANDOFF_COLUMNS = "t a psi pi K B X theta".split()
NORM, STOP = prior["waveform_constants"]()
BASE = dict(cells=64, late=6, cfl=.08, end=65., K=6000., theta=28000.,
            exit=.369, span=.0006, gate=.371, rate=100., nu=6., Q=1e8,
            D=.0001, dose=.123, canonical=False)
VARIANTS = {
    "primary": {},
    "coarse": {"cells":24,"late":3},
    "medium": {"cells":48,"late":4},
    "late_refined": {"cells":64,"late":10,"cfl":.07},
    "displacement": {"canonical":True},
    "theta_off": {"theta":0.},
    "even_energy_matched": {"theta":0., "K":6000. + 28000.**2/6000.04},
    "reversed": {"theta":-28000.},
    "no_charge": {"Q":0.}, "no_decay": {"rate":0.}, "no_drag": {"nu":0.},
    "broad_exit": {"span":.002}, "earlier_exit":{"exit":.36},
    "late_decay_gate": {"gate":.77}, "fast_decay_drag":{"rate":1500.,"nu":14.},
    "dose_minus": {"dose":.12295}, "dose_plus":{"dose":.12305},
}

def calculate(exe, changes, sample=False):
    p = BASE | changes
    args = [p["cells"],.002,p["dose"],p["end"],NORM,STOP,p["cfl"],p["late"],.08,
            p["K"],int(sample),p["D"],p["Q"],p["nu"],p["gate"],.03,p["rate"],
            .95,1.05,p["exit"],p["theta"],p["span"],int(p["canonical"])]
    raw = subprocess.check_output([str(exe),*map(str,args)],text=True)
    metrics, trace, handoff = {}, [], []
    for line in raw.splitlines():
        name,*values = line.split()
        numbers = list(map(float,values))
        if not all(map(math.isfinite,numbers)):
            raise ValueError(line)
        if name in ("@","#"):
            if len(numbers) != len(COLUMNS if name == "@" else HANDOFF_COLUMNS):
                raise ValueError(line)
            (trace if name == "@" else handoff).append(numbers)
        else:
            metrics[name], = numbers
    # Keep the run through its first expansion maximum with a short contracting tail.
    if sample:
        if len(handoff)!=1025:
            raise ValueError("Incomplete scalar-exit trajectory")
        # Keep full precision for field-crossing rows and regularly spaced times.
        # The source phase, scalar exit and heating interval have different time scales.
        metrics["trace"] = trace
    return metrics, handoff

def dump(path,data):
    path.write_text(json.dumps(data,sort_keys=True,separators=(",",":"),allow_nan=False)+"\n")

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--handoff",type=Path,required=True)
    parser.add_argument("--jobs",type=int,default=2)
    args=parser.parse_args()
    with tempfile.TemporaryDirectory() as folder:
        exe=Path(folder)/"forward"
        subprocess.run(["g++","-O3","-std=c++17",str(HERE/"forward.cpp"),"-o",str(exe)],check=True)
        def evaluate(name):
            # Controls resolve the source with 48 elements and the later stage with 4.
            params=VARIANTS[name]
            if name not in {"primary","coarse","medium","late_refined","displacement"}:
                params={"cells":48,"late":4}|params
            return BASE|params, calculate(exe,params,sample=name=="primary")
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            products=dict(zip(VARIANTS,pool.map(evaluate,VARIANTS)))
        dump(args.output,{"parameters":{k:p for k,(p,_) in products.items()},
             "trace_columns":COLUMNS,"runs":{k:r for k,(_,(r,_)) in products.items()}})
        dump(args.handoff,{"columns":HANDOFF_COLUMNS,"rows":products["primary"][1][1]})

if __name__=="__main__":
    main()
