"""Test cumulative occupation on the paired canonical labels swept by #192's birth spheres."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
import os
from pathlib import Path
import runpy
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
FRACTIONS = (.00001,.99999)
# A time quantile offset, Gaussian order, and native step size test separate errors.
CHECKS = {"native":("native",10,16384,.5,8),
          "new_time_offset":("native",10,16384,.173,8),
          "more_samples":("native",10,32768,.5,8),
          "coarse_integral":("native",7,8192,.5,4),
          "half_step":("half_step",10,16384,.5,8)}
FIELDS = ("local_hard", "total_weight", "sampled_local_hard", "history_extra_hard",
          "multiple_visit_weight", "global_hard_direct", "one_way_correlated",
          "one_way_uncorrelated","history_extra_at_reset_rate_18")


def evaluate(exe, traces, name, check):
    grid, depth, points, phase_t, order = CHECKS[check]
    record=traces[name]
    mu=record["mass"]
    load=math.pi**2*record["action"]/(6*mu**4)
    lo,hi=(record["reference_rest"]*v for v in FRACTIONS)
    args=[record["runs"][grid]["path"],mu,load,lo,hi,depth,points,phase_t,order]
    result=subprocess.check_output([str(exe),*map(str,args)],text=True)
    row=dict(zip(FIELDS,map(float,result.split()),strict=True))
    row["fraction_sampled_with_repeated_label"]=row.pop("multiple_visit_weight")/row["total_weight"]
    row["global_hard_decomposed"]=row["local_hard"]+row["history_extra_hard"]
    row["hard_with_reset_18"]=row["local_hard"]+row["history_extra_at_reset_rate_18"]
    if (not 0<=row["history_extra_at_reset_rate_18"]<=row["history_extra_hard"]*(1+1e-10) or
        abs(row["total_weight"]-1)>1.e-3 or
        abs(row["global_hard_direct"]/row["global_hard_decomposed"]-1)>.025 or
        not row["global_hard_decomposed"]<row["one_way_correlated"]<row["one_way_uncorrelated"]):
        raise ValueError(f"Label integration accuracy or occupancy ordering failed on {name}/{check}")
    return row


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    inputs=parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--prior",type=Path,help="Path to the merged #192 study")
    inputs.add_argument("--traces",type=Path,help="Readout manifest from science/history.py")
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--jobs",type=int,default=2)
    args=parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="hdblast-labels-") as folder:
        work=Path(folder)
        if args.traces:
            traces=json.loads(args.traces.read_bytes())
        else:
            build=runpy.run_path(str(HERE.parent/"science/history.py"))["build_histories"]
            traces=build(args.prior,work,args.jobs)
        exe=work/"labels"
        subprocess.run([os.environ.get("CXX","g++"),"-O3","-std=c++17",str(HERE/"labels.cpp"),"-o",str(exe)],check=True)
        tasks=[(n,c) for n in ("mass_4e10","mass_4e12") for c in CHECKS]
        records={n:{} for n in traces}
        with ThreadPoolExecutor(args.jobs) as pool:
            for (name,check),row in zip(tasks,pool.map(lambda pair:evaluate(exe,traces,*pair),tasks)):
                records[name][check]=row
    doc={"window_fraction":FRACTIONS,"fraction_denominator":"requested births in the measured window",
         "numerical_controls":CHECKS,"runs":records}
    args.output.write_text(json.dumps(doc,indent=2,sort_keys=True,allow_nan=False)+"\n")

if __name__=="__main__":
    main()
