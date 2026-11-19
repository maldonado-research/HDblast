"""Check the native results with Python, or run the complete matter/bulk evolution."""
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import subprocess
import sys
import zipfile

HERE=Path(__file__).resolve().parent
GRIDS=("coarse","fine")
STATES=("retain","erase","reflect","kick")
ARMS=("quiet","plus","minus","far")


def require(ok,label):
    if not ok:raise ValueError(label)


def quick():
    forks=json.loads((HERE/"FORK.json").read_bytes())
    data=json.loads((HERE/"NATIVE.json").read_bytes())
    with zipfile.ZipFile(HERE/"WALLS.zip") as z:
        for g in GRIDS:
            f=forks[g];w=f["weight"]
            old=f["states"]["retain"]
            kappa=-2*sum(b*c for b,c in zip(w,old["C"]))/sum(b*x for b,x in zip(w,old["X"]))
            kick=f["states"]["kick"]
            require(math.isclose(kappa,kick["kappa"],rel_tol=1.e-13),"Canonical impulse")
            for state in STATES:
                row=f["states"][state]
                for lhs,rhs in zip(row["native_sources"],old["native_sources"]):
                    require(abs(lhs-rhs)<1e-17,"Initial matter sources")
                require(row["min_determinant"]>=.25-2e-10,"Initial positivity")
                for arm in ARMS:
                    name=f"{g}/{state}/{arm}"
                    rows=list(csv.DictReader(io.StringIO(z.read(name+"/wall.csv").decode())))
                    lo=next(i for i,r in enumerate(rows) if float(r["H0tau"])>=2.1)-1
                    t0,t1=map(float,(rows[lo]["H0tau"],rows[lo+1]["H0tau"]))
                    u=(2.1-t0)/(t1-t0)
                    for key in ("H_over_H0","phi_b","rho_hat","p_hat","j_hat","comoving_energy"):
                        a,b=map(float,(rows[lo][key],rows[lo+1][key]))
                        require(math.isclose((1-u)*a+u*b,data["rows"][name][key],rel_tol=1e-10,abs_tol=1e-10),
                                name+"/"+key)
    require(len(data["rows"])==32,"Native experiment coverage")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("repository",type=Path)
    ap.add_argument("--jobs",type=int,default=2)
    ap.add_argument("--full",action="store_true")
    ap.add_argument("--prefork",action="store_true",help="Also reproduce preparation checkpoints from T=0")
    args=ap.parse_args()
    quick()
    if args.full:
        subprocess.run([sys.executable,"-I","-B",str(HERE/"verify_native.py"),str(args.repository.resolve()),
                        "--full","--jobs",str(args.jobs)],check=True)
    if args.prefork:
        subprocess.run([sys.executable,"-I","-B",str(HERE/"rebuild_forks.py"),str(args.repository.resolve()),"--jobs",str(args.jobs)],check=True)
    print("Verified four preparations and 32 native wall histories.")

if __name__=="__main__":main()
