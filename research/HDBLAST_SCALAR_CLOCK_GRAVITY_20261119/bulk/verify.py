"""Check the portable signed-run archive, source hashes, matched-time slopes and native controls."""
import argparse
import csv
import hashlib
from pathlib import Path
import sys

import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import ensemble
import summarize


def verify(runs,tables,results):
    slopes={}
    max_diff=0.
    for g,a,l,d,k in ensemble.cases():
        stem=f"{g}_{a}_zero" if l==0 else f"{g}_{a}_{d}_{k}_{'p' if l>0 else 'm'}{abs(l):g}"
        with np.load(runs/(stem+".npz"),allow_pickle=False) as data:
            if l:
                digest=hashlib.sha256((tables/f"{g}_{a}_{d}_{k}.npz").read_bytes()).hexdigest()
                if digest!=str(data["source_sha256"]):
                    raise ValueError("Table hash does not match "+stem)
            if l<=0:
                continue
            h=summarize.CubicSpline(data["s"],data["h"])(2.10)
        nm=f"{g}_{a}_{d}_{k}_m{abs(l):g}.npz"
        with np.load(runs/nm,allow_pickle=False) as neg:
            hm=summarize.CubicSpline(neg["s"],neg["h"])(2.10)
        with np.load(runs/f"{g}_{a}_zero.npz",allow_pickle=False) as z:
            hz=summarize.CubicSpline(z["s"],z["h"])(2.10)
        slopes[g,a,d,k,l]=float((h-hm)/(2*l*hz))
    with (results/"READOUTS.csv").open() as f:
        rows=list(csv.DictReader(f))
    if len(rows)!=len(slopes):
        raise ValueError("Case count")
    for row in rows:
        key=(row["grid"],row["arm"],int(row["cutoff"]),row["scheme"],float(row["step"]))
        err=abs(float(row["d_h_over_h"])-slopes[key])
        max_diff=max(max_diff,err)
    if max_diff>1e-13:
        raise ValueError("Archive disagrees with published response")
    print("Verified",len(slopes),"signed pairs and their source hashes. Max fractional-slope discrepancy",max_diff)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs",type=Path,default=HERE/"runs")
    ap.add_argument("--tables",type=Path,default=HERE.parent/"source/tables")
    ap.add_argument("--results",type=Path,default=HERE)
    args=ap.parse_args()
    verify(args.runs,args.tables,args.results)
