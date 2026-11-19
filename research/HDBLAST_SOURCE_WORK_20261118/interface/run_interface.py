"""Regenerate the two-sided bulk-flux diagnostic using strictly interior scalar samples."""
import argparse
import concurrent.futures
import csv
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import zipfile
import numpy as np


ARMS = ("quiet","plus","minus","far")
GRIDS = ("coarse","fine")
SET = "k24_y1_b1_t0.700"
COLS = ("T","H0tau","A_b","B_b","phi_b","v_over_H0","H_over_H0","j_hat","scalar_work",
        "J0","sigma1","rho_hat","rad_hat","comoving_energy","pressure_work","bulk_phi_T",
        "gradient_wall_z", *(f"inside_phi_{i}" for i in range(1,6)))
PATCH = """
            row['bulk_phi_T'] = float(F['phi_T'][-1])
            row['gradient_wall_z'] = float(self.z[-1])
            if self.mode != 'full' or abs(self.z[-1]) > 1e-12:
                raise ValueError('The postfork observable requires the full conformal chart')
            for j in range(1,6):
                row['inside_phi_'+str(j)] = float(F['phi'][-1-j])
                row['inside_z_'+str(j)] = float(self.z[-1-j])
"""


def compact(root, parameters, output):
    with zipfile.ZipFile(output,"w",compression=zipfile.ZIP_DEFLATED) as archive:
        def write(name,content):
            item=zipfile.ZipInfo(name,(2026,11,18,0,0,0))
            item.compress_type=zipfile.ZIP_DEFLATED
            item.external_attr=0o100644 << 16
            archive.writestr(item,content)
        for grid in GRIDS:
            for arm in ARMS:
                f=root/grid/SET/arm/"wall_timeseries.npz"
                with np.load(f,allow_pickle=False) as z:
                    vals=dict(zip(z["cols"],z["rec"].T,strict=True))
                    tab=io.StringIO()
                    writer=csv.writer(tab,lineterminator="\n")
                    writer.writerow(COLS)
                    for row in zip(*(vals[k] for k in COLS),strict=True):
                        writer.writerow([format(x,".17g") for x in row])
                    write(f"{grid}/{arm}.csv",tab.getvalue())
                    coords=[float(vals[f"inside_z_{i}"][0]) for i in range(1,6)]
                    for i,x in enumerate(coords,1):
                        if np.any(vals[f"inside_z_{i}"]!=x):raise ValueError("Moving spatial grid")
                    if grid in parameters["interior_z"] and parameters["interior_z"][grid]!=coords:
                        raise ValueError("Mismatched arm coordinates")
                    parameters["interior_z"][grid]=coords
        write("parameters.json",json.dumps(parameters,indent=2,sort_keys=True)+"\n")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repo",type=Path,help="Input repository root containing previous B1 sources")
    parser.add_argument("output",type=Path,help="ZIP file for measured interface histories")
    parser.add_argument("--coupled",type=Path,default=Path(__file__).resolve().parents[1]/"coupled")
    parser.add_argument("--forks",type=Path,help="Optionally reuse the preceding coupled runner's local fork checkpoints")
    parser.add_argument("--jobs",type=int,default=2)
    args=parser.parse_args()
    repo=args.repo.resolve()
    coupled=args.coupled.resolve()
    spec=importlib.util.spec_from_file_location("coupled_driver",coupled/"run_coupled.py")
    driver=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(driver)
    with tempfile.TemporaryDirectory() as tmp:
        work=Path(tmp)
        env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1")
        def run(grid,arm):
            subprocess.run([sys.executable,"-I","-B",coupled/"run_coupled.py",repo,work,
                            "--grid",grid,"--branch",arm,"--decay","1","--nodes","24"],env=env,check=True)
        params={"b":6.49e-6,"interior_z":{},"rb":{},"A_ref":{},"fork_digests":{}}
        def setup(grid):
            href=repo/driver.OLD/"runs"/grid/"quiet/rolling_summary.json"
            hist=json.loads(href.read_bytes())
            params["rb"][grid]=hist["rho_b"]
            params["A_ref"][grid]=hist["matter"]["events"][0]["A_star"]
            target=work/grid/SET/"fork.npz"
            if args.forks:
                src=args.forks.resolve()/grid
                # Reconstruct the coupled solver from authenticated historical inputs.
                solver=work/grid/"solver"
                solver.mkdir(parents=True)
                for path,digest in driver.PINS.items():
                    data=(repo/path).read_bytes()
                    if hashlib.sha256(data).hexdigest()!=digest:raise ValueError(path)
                    if path.name=="evolve_b1.py":data=driver.patched(data.decode()).encode()
                    (solver/path.name).write_bytes(data)
                b=(repo/driver.OLD/"rolling_data.py").read_bytes()
                if hashlib.sha256(b).hexdigest()!=driver.MATCH_PIN:raise ValueError("Pulse data builder")
                (solver/"rolling_data.py").write_bytes(b)
                target.parent.mkdir()
                shutil.copyfile(src/SET/"fork.npz",target)
            else:
                run(grid,"initial")
                run(grid,"fork")
            params["fork_digests"][grid]=hashlib.sha256(target.read_bytes()).hexdigest()
            p=work/grid/"solver/evolve_b1.py"
            s=p.read_text()
            anchor="            return row, F"
            if s.count(anchor)!=1:raise ValueError("Unknown record function")
            # Record postfork every native RK4 step to resolve oscillations in the
            # strictly interior derivative. This changes the readout cadence only.
            cadence = "every = max(1, int(round(rec_dT/dt)))"
            if s.count(cadence)!=1:raise ValueError("Unknown snapshot cadence")
            p.write_text(s.replace(anchor,PATCH+anchor).replace(cadence,"every = 1"))
            run(grid,"quiet")
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
            list(pool.map(setup,GRIDS))
            futures=[pool.submit(run,grid,arm) for grid in GRIDS for arm in ARMS[1:]]
            for fut in futures:fut.result()
        compact(work,params,args.output)


if __name__=="__main__":
    main()
