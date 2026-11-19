"""Generate the archived controls and the full signed-strength response matrix."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
GRIDS = ("fine", "coarse")
ARMS = ("quiet", "plus", "minus", "far")


def cases():
    for g in GRIDS:
        for a in ARMS:
            yield g,a,0,120,"static"
    for d in (120,60):
        for k in ("static","bare"):
            for a in ARMS:
                for g in GRIDS:
                    for strength in (.01,-.01):
                        yield g,a,strength,d,k
    for d in (120,60):
        for k in ("static","bare"):
            for g in GRIDS:
                for arm in (("quiet","far") if (d,k)==(120,"static") else ("quiet",)):
                    for strength in (.005,-.005):
                        yield g,arm,strength,d,k


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--batch",choices=("zero","nonzero","all"),default="all")
    ap.add_argument("--output",type=Path,default=HERE/"runs")
    ap.add_argument("--tables",type=Path,default=HERE.parent/"source/tables")
    ap.add_argument("--jobs",type=int,default=2)
    ap.add_argument("--resume",action="store_true")
    ap.add_argument("--cache",type=Path,default=Path(tempfile.gettempdir())/"hdblast-live-static-cache")
    args=ap.parse_args()
    env=dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    def execute(case):
        g,a,l,d,k=case
        name=f"{g}_{a}_zero" if l==0 else f"{g}_{a}_{d}_{k}_{'p' if l>0 else 'm'}{abs(l):g}"
        output_path=args.output/(name+".npz")
        if args.resume and output_path.exists():
            if not l:
                return
            table=args.tables/f"{g}_{a}_{d}_{k}.npz"
            digest=hashlib.sha256(table.read_bytes()).hexdigest()
            with np.load(output_path,allow_pickle=False) as saved:
                if "source_sha256" in saved and str(saved["source_sha256"])==digest:
                    return
        subprocess.run([sys.executable,"-IB",str(HERE/"live.py"),*map(str,case),"--output",str(args.output),
                        "--tables",str(args.tables),"--cache",str(args.cache)],env=env,check=True)
    # Warm the static geometry cache before any child processes open it.
    selected=[c for c in cases() if args.batch=="all" or (c[2]==0)==(args.batch=="zero")]
    if not list((args.cache/"runs/shell_cache").glob("*.pkl")):
        execute(selected.pop(0))
    with ThreadPoolExecutor(args.jobs) as pool:
        list(pool.map(execute,selected))

if __name__=="__main__":
    main()
