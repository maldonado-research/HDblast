"""Evolve two quiet preparations sharing initial sources and their first derivatives."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import zipfile
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/"native"))
import native
import summarize

STATES=("smoothplus","smoothminus")


def transform(cov,state,weight):
    result=cov.copy()
    x,p,c=result[3:].reshape(3,-1)
    t=((1+np.polynomial.legendre.leggauss(len(weight))[0])/2)**2
    basis=np.column_stack((np.ones_like(t),t))
    measure=weight*c
    b0,b1=np.linalg.solve(basis.T@(measure[:,None]*basis),basis.T@(measure*t*t))
    extrema=np.array([0.,1.,np.clip(b1/2,0,1)])
    norm=np.max(abs(extrema*extrema-b1*extrema-b0))
    sign=1 if state=="smoothplus" else -1
    c*=sign*.9*(t*t-b1*t-b0)/norm
    return result,0.


def package(work,out):
    forks={g:json.loads((work/g/"FORK.json").read_bytes()) for g in native.GRIDS}
    (out/"FORK.json").write_text(summarize.json_text(forks))
    with zipfile.ZipFile(out/"WALLS.zip","w") as z:
        for g in native.GRIDS:
            for state in STATES:
                name=f"{g}/{state}/quiet/wall.csv"
                summarize.store(z,name,(work/name).read_bytes())


def analyze(folder):
    forks=json.loads((folder/"FORK.json").read_bytes())
    settings={}
    keys=summarize.READOUT
    with zipfile.ZipFile(folder/"WALLS.zip") as walls,zipfile.ZipFile(HERE.parent/"native"/"WALLS.zip") as comparison:
        for g in native.GRIDS:
            q=np.array(forks[g]["q"])
            w=np.array(forks[g]["weight"])
            states=forks[g]["states"]
            delta=np.array(states[STATES[0]]["C"])-np.array(states[STATES[1]]["C"])
            report={}
            for s in (*STATES,"erase"):
                bundle=comparison if s=="erase" else walls
                tab=np.genfromtxt(io.BytesIO(bundle.read(f"{g}/{s}/quiet/wall.csv")),names=True,delimiter=",")
                report[s]={k:float(np.interp(summarize.CUT,tab["H0tau"],tab[k])) for k in keys}
            for s in STATES:
                report[s+"_minus_erase"]={k:report[s][k]-report["erase"][k] for k in keys}
            report["plus_minus_minus"]={k:report[STATES[0]][k]-report[STATES[1]][k] for k in keys}
            report["moments"]={str(k):float(w@(q**k*delta)) for k in (0,2,4)}
            for key in ("native_sources","native_T_derivatives","neumann","neumann_T"):
                diff=np.array(states[STATES[0]][key])-states[STATES[1]][key]
                if np.max(abs(diff))>5.e-15:raise ValueError("First-order junction matching")
                report["initial_delta_"+key]=diff.tolist()
            if abs(report["moments"]["0"])>1e-14 or abs(report["moments"]["2"])>1e-12:
                raise ValueError("Initial derivative moments")
            settings[g]=report
    return settings


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("repo",type=Path)
    ap.add_argument("work",type=Path)
    ap.add_argument("output",type=Path)
    ap.add_argument("--jobs",type=int,default=2)
    ap.add_argument("--one",nargs=2,metavar=("GRID","STATE"))
    args=ap.parse_args()
    repo,work,out=args.repo.resolve(),args.work.resolve(),args.output.resolve()
    if args.one:
        native.worker(repo,work,*args.one,"quiet")
        return
    native.install(repo,work)
    out.mkdir(exist_ok=True,parents=True)
    for g in native.GRIDS:native.fork(repo,work,g,STATES,transform)
    env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1")
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        tasks=[pool.submit(subprocess.run,[sys.executable,"-I","-B",str(HERE/"smooth.py"),str(repo),str(work),str(out),
                                           "--one",g,s],env=env,check=True) for g in native.GRIDS for s in STATES]
        for task in tasks:task.result()
    package(work,out)
    (out/"SMOOTH.json").write_text(summarize.json_text(analyze(out)))

if __name__=="__main__":main()
