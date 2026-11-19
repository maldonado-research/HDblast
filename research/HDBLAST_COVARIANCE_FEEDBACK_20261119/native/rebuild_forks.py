"""Rebuild the two T=1.11 native restart states from the earlier public 5D solver."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import io
import json
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


def rebuild(repo,tmp,grid,zsource):
    # The old launcher verifies the bulk and pulse scripts before staging them.
    # Verify its dynamics and normalization data as well before running.
    driver=repo/native.OLD/"run_coupled.py"
    native.checked(driver.read_bytes(), native.PINS["run_coupled.py"])
    native.checked((repo/native.OLD/"quantum_shell.py").read_bytes(), native.PINS["quantum_shell.py"])
    native.checked((repo/native.CLASSICAL/"runs"/grid/"quiet"/"rolling_summary.json").read_bytes(),native.REF_PIN[grid])
    proc_env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1")
    for stage in ("initial","fork"):
        subprocess.run([sys.executable,"-I","-B",str(driver),str(repo),str(tmp),"--grid",grid,"--branch",stage,
                        "--nodes","24","--decay","0","--prep-T","0.70"],env=proc_env,check=True)
    new=tmp/grid/"k24_y0_b1_t0.700"/"fork.npz"
    delta={}
    with np.load(new,allow_pickle=False) as actual,np.load(io.BytesIO(zsource),allow_pickle=False) as reference:
        if set(actual.files)!=set(reference.files):raise ValueError("Checkpoint fields")
        for key in reference:
            a,b=actual[key],reference[key]
            if a.shape!=b.shape:raise ValueError((grid,key,"shape"))
            if a.dtype.kind in "fci":
                diff=float(np.max(abs(a-b)))
                scale=float(np.max(abs(b)))
                if not np.allclose(a,b,rtol=1.e-8,atol=4.e-9):
                    raise ValueError((grid,key,diff,scale))
                delta[key]=diff
            elif not np.array_equal(a,b):raise ValueError((grid,key))
    return delta


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("repository",type=Path)
    ap.add_argument("--jobs",type=int,default=2)
    args=ap.parse_args()
    with tempfile.TemporaryDirectory() as td,zipfile.ZipFile(HERE/"FORK_INPUTS.zip") as bundle:
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            pending={g:pool.submit(rebuild,args.repository.resolve(),Path(td),g,
                                   native.checked(bundle.read(g+".npz"),native.FORK_PIN[g])) for g in native.GRIDS}
            print(json.dumps({g:f.result() for g,f in pending.items()},indent=2))

if __name__=="__main__":
    main()
