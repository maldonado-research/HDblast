"""Continue the closed finite-band Gaussian B1 solution from its native T=1.11 fork."""
import argparse
import csv
import hashlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import zipfile

import numpy as np

HERE = Path(__file__).resolve().parent
OLD = Path("research/HDBLAST_SOURCE_WORK_20261118/coupled")
CLASSICAL = Path("research/HDBLAST_PULSE_PARTICLE_RESPONSE_20261118/classical")
GRIDS = ("coarse","fine")
STATES = ("retain","erase","reflect","kick")
ARMS = ("quiet","plus","minus","far")
PINS = {"run_coupled.py":"05df56cef9d16c2e7c71edbc57f6c6a0ac86043cea79d0c0d9296292bce872fb",
        "quantum_shell.py":"b056bd57908834ddafc2f0b44531d20619c375e3be1eafbaabc4456dd1536d97"}
REF_PIN = {"coarse":"91601211fce4a617838d3a709849a138ec29a23998cc7c3ba2416ec0b0b8aa8c",
           "fine":"abec0887ee4de00ba0757e7d96987610f36b7157a85bc4f00696fee8064a6810"}
FORK_PIN = {"fine":"d88e28d8d60eaa3803c1db54b8ac6c183a49aba9a21f71256ad0ae7e1f3613c5",
            "coarse":"02f3937aefa031ba24ba29bb79f16aa291512f708b9a5885dc908c338849a5b8"}
COLS = ("T","H0tau","A_b","B_b","phi_b","H_over_H0","v_over_H0","rho_hat","p_hat","j_hat",
        "comoving_number","scalar_work","pressure_work","comoving_energy","min_uncertainty","min_occupation",
        "Hshell","Mshell","Hmax_near","Mmax_near")


def checked(data, digest):
    if hashlib.sha256(data).hexdigest() != digest:
        raise ValueError("Unexpected input SHA-256: "+digest)
    return data


def dump(path,data):
    path.parent.mkdir(exist_ok=True,parents=True)
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")


def install(repo,work):
    solver=work/"solver"
    solver.mkdir(exist_ok=True,parents=True)
    for name,pin in PINS.items():
        (solver/name).write_bytes(checked((repo/OLD/name).read_bytes(),pin))
    sys.path.insert(0,str(solver))
    import run_coupled as driver
    for path,pin in driver.PINS.items():
        content=checked((repo/path).read_bytes(),pin)
        if path.name=="evolve_b1.py":
            content=driver.patched(content.decode()).encode()
        driver.install(solver/path.name,content)
    driver.install(solver/"rolling_data.py",checked((repo/CLASSICAL/"rolling_data.py").read_bytes(),driver.MATCH_PIN))
    return solver


def run_instance(repo,work,grid):
    sys.path.insert(0,str(work/"solver"))
    from quantum_shell import QuantumShell
    from evolve_b1 import Run
    refpath=repo/CLASSICAL/"runs"/grid/"quiet"/"rolling_summary.json"
    ref=json.loads(checked(refpath.read_bytes(),REF_PIN[grid]))["matter"]["events"][0]["A_star"]
    fine=grid=="fine"
    return Run(delta=.1,c=.5975949350280132,d=-3.106933495673783,dc=.01,dzf=.0005 if fine else .001,
               dzc=.002 if fine else .004,L=16.,kappa=10.,project=None,cfl=.5,order=4,xc=1.9,
               matter=lambda rb: QuantumShell(rb,ref,nodes=24,feedback=True,decay=0.),
               log=lambda msg:None)


def transform(cov,state,weight):
    out=cov.copy()
    x,p,c=out[3:].reshape(3,-1)
    kappa=0.
    if state=="kick":
        kappa=-2*np.dot(weight,c)/np.dot(weight,x)
        p+=2*kappa*c+kappa*kappa*x
        c+=kappa*x
    elif state=="erase":
        c[:]=0.
    elif state=="reflect":
        c*=-1
    return out,float(kappa)


def fork(repo,work,grid,labels=STATES,transform_cov=transform):
    run=run_instance(repo,work,grid)
    src=HERE/"FORK_INPUTS.zip"
    with zipfile.ZipFile(src) as z:
        data=checked(z.read(grid+".npz"),FORK_PIN[grid])
    with np.load(io.BytesIO(data),allow_pickle=False) as saved:
        cp=dict(saved)
    t=float(cp['T'])
    y=cp['bulk']
    run.mode=str(cp['mode'].item())
    a_log,b_log,phi,at,vt=run.shell_vals(t,y)
    matter=run.matter
    a,_,_,_,energy=matter.background(a_log,phi)
    states={}
    for state in labels:
        m,kappa=transform_cov(cp['matter'],state,matter.weight)
        x,p,c=matter.view(m)
        dm=matter.rhs(m,a_log,phi,at,vt,math.exp(b_log))
        rho,pr,j=matter.totals(m,a_log,phi)
        dt_sources=matter.totals_T(m,dm,a_log,phi,at,vt,math.exp(b_log))
        f=run.fields(t,y,m)
        phitt=run.accel(f)[2][-1]
        h,mo,sh,sm=run.constraints(f,y[3],m,dm,phitt,t)
        occupation=.5*(a**3*energy*x+p/(a**3*energy)-1)
        real=.5*(a**3*energy*x-p/(a**3*energy))
        states[state]=dict(X=x.tolist(),P=p.tolist(),C=c.tolist(),n=occupation.tolist(),Re_c=real.tolist(),kappa=kappa,
                           native_sources=[rho,pr,j],native_T_derivatives=list(dt_sources),
                           neumann=list(run.neumann(a_log,b_log,phi,vt,m,t)),
                           neumann_T=list(run.neumann_t(a_log,b_log,float(f["B_T"][-1]),at,phi,vt,phitt,m,dm,t)),
                           shell_constraint=[float(np.max(abs(h[-6:])/sh[-6:])),float(np.max(abs(mo[-6:])/sm[-6:]))],
                           min_determinant=float(np.min(x*p-c*c)))
        target=work/grid/state
        target.mkdir(parents=True,exist_ok=True)
        np.savez(target/"fork.npz",**{**cp,"matter":m})
    doc=dict(T=t,H0tau=float(cp['tau']),A_b=a_log,B_b=b_log,phi_b=phi,a=a,H_over_H0=at*run.rb/math.exp(b_log),
             v_over_H0=vt*run.rb/math.exp(b_log),reference_A= matter.ref,lapse=math.exp(b_log)/run.rb,
             q=matter.k.tolist(),weight=matter.weight.tolist(),frequency=energy.tolist(),states=states)
    dump(work/grid/"FORK.json",doc)


def write_wall(path):
    with np.load(path/"wall_timeseries.npz",allow_pickle=False) as d:
        names=d["cols"].tolist()
        table=d["rec"][:,[names.index(k) for k in COLS]]
    with (path/"wall.csv").open("w") as f:
        np.savetxt(f,table,header=",".join(COLS),comments="",delimiter=",",fmt="%.16g")


def worker(repo,work,grid,state,arm):
    run=run_instance(repo,work,grid)
    run.branch=arm
    run.fork=work/grid/state/"fork.npz"
    root=work/grid/state/arm
    run.evolve(T_final=1.62,rec_dT=.001,tag="wall",out_dir=str(root),wall_limit=2000)
    d=json.loads((root/"wall_summary.json").read_bytes())
    if d["stop_reason"]!="T_final":
        raise ValueError(d["stop_reason"])
    if hasattr(run,"match"):
        dump(root/"match.json",run.match)
    write_wall(root)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("repo",type=Path)
    ap.add_argument("work",type=Path)
    ap.add_argument("--prepare",action="store_true")
    ap.add_argument("--one",nargs=3,metavar=("GRID","STATE","ARM"))
    ap.add_argument("--jobs",type=int,default=2)
    args=ap.parse_args()
    repo,work=args.repo.resolve(),args.work.resolve()
    if args.one:
        worker(repo,work,*args.one)
        return
    install(repo,work)
    for grid in GRIDS:
        fork(repo,work,grid)
    if args.prepare:
        return
    from concurrent.futures import ThreadPoolExecutor
    env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1")
    cases=[(g,s,a) for s in STATES for a in ARMS for g in GRIDS]
    with ThreadPoolExecutor(args.jobs) as pool:
        jobs=[pool.submit(subprocess.run,[sys.executable,"-I","-B",str(HERE/"native.py"),str(repo),str(work),"--one",*case],
                          check=True,env=env) for case in cases]
        for case,job in zip(cases,jobs):
            job.result()
            print("/".join(case),flush=True)


if __name__=="__main__":
    main()
