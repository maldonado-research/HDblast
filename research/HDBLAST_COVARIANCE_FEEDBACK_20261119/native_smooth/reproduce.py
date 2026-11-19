"""Check or rebuild the four first-derivative-matched quiet histories."""
import argparse
import csv
import io
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

HERE=Path(__file__).resolve().parent


def verify(folder):
    fork=json.loads((folder/"FORK.json").read_bytes())
    report=json.loads((folder/"SMOOTH.json").read_bytes())
    with zipfile.ZipFile(folder/"WALLS.zip") as z:
        for grid in ("coarse","fine"):
            f=fork[grid]
            d=report[grid]
            signed=[]
            for s in ("smoothplus","smoothminus"):
                state=f["states"][s]
                if state["min_determinant"]<.25:raise ValueError("Positive covariance")
                series=list(csv.DictReader(io.StringIO(z.read(f"{grid}/{s}/quiet/wall.csv").decode())))
                i=next(i for i,r in enumerate(series) if float(r["H0tau"])>=2.10)
                a,b=series[i-1],series[i]
                t=(2.10-float(a["H0tau"]))/(float(b["H0tau"])-float(a["H0tau"]))
                h=float(a["H_over_H0"])*(1-t)+t*float(b["H_over_H0"])
                signed.append(h)
                if abs(h-d[s]["H_over_H0"])>1e-13:raise ValueError("Proper readout")
            if abs(signed[0]-signed[1]-d["plus_minus_minus"]["H_over_H0"])>1e-13:raise ValueError("Delta H")
            if not -4.7e-7 < signed[0]-signed[1] < -4.3e-7:raise ValueError("Expansion difference")
            if abs(d["moments"]["0"])>1e-13 or abs(d["moments"]["2"])>1e-12:raise ValueError("Initial covariance moments")
            for key in ("native_sources","native_T_derivatives","neumann","neumann_T"):
                if max(abs(x) for x in d["initial_delta_"+key])>1e-14:raise ValueError("Initial matching")
    return report


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("repo",type=Path)
    ap.add_argument("--full",action="store_true")
    ap.add_argument("--jobs",type=int,default=2)
    args=ap.parse_args()
    gold=verify(HERE)
    if args.full:
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            subprocess.run([sys.executable,"-I","-B",str(HERE/"smooth.py"),str(args.repo.resolve()),str(root/"work"),
                            str(root/"out"),"--jobs",str(args.jobs)],check=True)
            actual=verify(root/"out")
            for g in ("coarse","fine"):
                for key in ("H_over_H0","phi_b","rho_hat","p_hat","j_hat"):
                    if not math.isclose(gold[g]["plus_minus_minus"][key],actual[g]["plus_minus_minus"][key],
                                        rel_tol=.001,abs_tol=3e-12):raise ValueError("Cross-platform smooth response")
    print(f"Verified four quiet histories with matched initial sources and derivatives. Full numerical run: {args.full}")


if __name__=="__main__":main()
