"""Audit source-matched native histories and optionally rerun the 32 configurations."""
import argparse
import hashlib
import io
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import native
import summarize

EARLIER="7f3289631a80291f6180e51f5b648a135c6e33c6366d4f97c32a4b90630cab15"


def require(ok,label):
    if not ok:raise ValueError(label)


def close(a,b,label,rtol=2e-7,atol=5e-10):
    require(np.allclose(a,b,rtol=rtol,atol=atol,equal_nan=True),label)


def nested(a,b,path=""):
    if isinstance(b,dict):
        require(set(a)==set(b),path+" keys")
        for k in b:nested(a[k],b[k],path+"/"+k)
    else:
        close(a,b,path)


def check(folder,repo):
    doc=json.loads((folder/"NATIVE.json").read_bytes())
    forks=json.loads((folder/"FORK.json").read_bytes())
    require(len(doc["rows"])==32 and len(doc["pulse_minus_quiet"])==24,"Run coverage")
    with zipfile.ZipFile(folder/"FORK_INPUTS.zip") as z:
        for g in native.GRIDS:native.checked(z.read(g+".npz"),native.FORK_PIN[g])
    earlier=native.checked((repo/native.OLD/"WALLS.zip").read_bytes(),EARLIER)
    with zipfile.ZipFile(io.BytesIO(earlier)) as archive,zipfile.ZipFile(folder/"WALLS.zip") as walls:
        for g in native.GRIDS:
            f=forks[g]
            wf=np.array(f["weight"])
            q=np.array(f["q"])
            gauss,w=np.polynomial.legendre.leggauss(24)
            close(q,4*(1+gauss),"Band quadrature")
            close(wf,4*w*q*q/(2*np.pi**2),"Mode density")
            states=f["states"]
            r=states["retain"]
            x,p,c=[np.array(r[k]) for k in ("X","P","C")]
            kick=-2*wf@c/(wf@x)
            for s in native.STATES:
                row=states[s]
                close(row["X"],x,"Initial field variance",atol=1e-13)
                close(row["native_sources"],r["native_sources"],"Initial sources",atol=1.e-17)
                close(row["neumann"],r["neumann"],"Initial junctions",atol=1.e-15,rtol=1e-15)
                close(row["neumann_T"][0],r["neumann_T"][0],"Initial metric derivative")
                close(row["shell_constraint"],r["shell_constraint"],"Initial H and M",atol=1.e-12)
                if s=="kick":
                    close(row["kappa"],kick,"Kick")
                    close(row["C"],c+kick*x,"Sheared covariance")
                    close(row["P"],p+2*kick*c+kick*kick*x,"Sheared momentum")
                else:
                    close(row["P"],p,"Fixed momentum")
                    close(row["C"],{"reflect":-c,"retain":c,"erase":0.}[s],"Prescribed covariance")
                for a in native.ARMS:
                    key=f"{g}/{s}/{a}"
                    result=doc["rows"][key]
                    require(result["min_determinant"]>.24999999 and result["max_energy_ledger_residual"]<3e-6,
                            "Positivity and conservation "+key)
                    require(result["max_shell_constraint"]["Hshell"]<1e-4 and result["max_shell_constraint"]["Mshell"]<2e-4,
                            "Shell constraint residual "+key)
                    if a!="quiet":
                        d=doc["pulse_minus_quiet"][key]
                        require((d["H_over_H0"]>0)==(a=="minus"),"Expansion pulse "+key)
                        require((d["rho_hat"]>0)==(a!="plus"),"Parent density pulse "+key)
                    if s!="retain":
                        require(doc["phase_minus_retained"][key]["H_over_H0"]>0,"Expansion phase "+key)
                    if s=="retain":
                        v=np.genfromtxt(io.BytesIO(walls.read(key+"/wall.csv")),names=True,delimiter=",")
                        old=np.genfromtxt(io.BytesIO(archive.read(f"{g}/y0_b1/{a}.csv")),names=True,delimiter=",")
                        for name in set(v.dtype.names)&set(old.dtype.names):
                            close(v[name],old[name],"Historical retained "+name,rtol=1e-8,atol=1e-9)
    # Recreate every interpolated contrast and diagnostic from the public time series.
    with tempfile.TemporaryDirectory() as scratch:
        path=Path(scratch)
        for name in ("WALLS.zip","FORK.json"):
            (path/name).symlink_to((folder/name).resolve())
        nested(summarize.analyze(path),doc)
    return doc,forks


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("repository",type=Path)
    ap.add_argument("--full",action="store_true",help="Regenerate all 32 histories from bundled restart states")
    ap.add_argument("--jobs",type=int,default=2)
    args=ap.parse_args()
    repo=args.repository.resolve()
    original,forks=check(HERE,repo)
    if args.full:
        with tempfile.TemporaryDirectory() as scratch:
            work=Path(scratch)/"native"
            subprocess.run([sys.executable,"-I","-B",str(HERE/"native.py"),str(repo),str(work),
                            "--jobs",str(args.jobs)],check=True)
            out=Path(scratch)/"report"
            out.mkdir()
            summarize.package(work,out)
            nested(summarize.analyze(out),original)
            nested(json.loads((out/"FORK.json").read_bytes()),forks)
    print(f"Verified 32 native histories and 24 pulse contrasts; full evolution: {args.full}")

if __name__=="__main__":main()
