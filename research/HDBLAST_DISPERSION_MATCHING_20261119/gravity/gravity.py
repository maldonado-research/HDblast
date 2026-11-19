"""Replay the matched-static, two-cutoff gravitational prediction from pinned upstream data."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

HERE=Path(__file__).resolve().parent
REV="2cbafaefde87652135d8222e50a53da4eb704c41"
OLD=Path("research/HDBLAST_SCALAR_CLOCK_GRAVITY_20261119")
NATIVE=Path("research/HDBLAST_COVARIANCE_FEEDBACK_20261119/native")
PINS={
    OLD/"bulk/READOUTS.csv":"acd967b4d2ed780a13acb2409a36a6e663561428886106ff133cf75fd8716fce",
    OLD/"bulk/live.py":"96b2613f413645136100b80a8ecb1e35f73bc164ca583dadbf09ec3fbb1b4409",
    OLD/"bulk/vendor/SHA256SUMS":"20119edbbe6be623a4b5bfec2bd3a469a3f078417910f81eb67a4877b01ce8ff",
    OLD/"source/source.py":"315f602a242c069d06d641b70a00db651d50bd4164436828503f078427a06aad",
    OLD/"source/SOURCES.json":"975910253a8de127f51862bbd6f510de9f38a1bb4d36f13b4be26696b49aa729",
    NATIVE/"FORK.json":"d1028e10aa6440c98d71e5c99cf2c18fd336606d3068816f02d01db499918038",
    NATIVE/"WALLS.zip":"a18041b6ab8591cc88e7613e0759c812d00da5655953b5cb3a25fbb36bc9d05d",
}
ALPHA=1.76
CUTS=(60,120)
GRIDS=("coarse","fine")
EPS=(.01,.005)
SOURCES=("rho","p","j","drho","dp","dj")


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def authenticate(repo):
    for p,digest in PINS.items():
        if hashlib.sha256((repo/p).read_bytes()).hexdigest()!=digest:
            raise ValueError(f"Historical input differs at {p} for revision {REV}")
    for line in (repo/OLD/"bulk/vendor/SHA256SUMS").read_text().splitlines():
        digest,name=line.split()
        if hashlib.sha256((repo/OLD/"bulk/vendor"/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(name)


def moment(q, cutoff):
    # The same 64-node pair rule that defines the finite Hamiltonians in the released study.
    x,w=np.polynomial.legendre.leggauss(64)
    om=q[:,None]+(cutoff-q[:,None])*(x+1)/2
    weight=(cutoff-q[:,None])*w/2
    return np.sum((om**2-q[:,None]**2)*weight/(4*np.pi**2*om),axis=1)


def local_response(q,weights,initial,background,clock,t):
    """Independent DOP853 solve of the effective local quadratic perturbation for every parent mode."""
    _,loga,mass=background
    K=ALPHA*moment(q,120)
    def rhs(s,state):
        X,C,P,dX,dC,dP=state.reshape(6,-1)
        a=np.exp(loga(s)); m=mass(s); f=clock(m/100+.5)[0]
        inv=a**-3
        b=a*q*q+a**3*m*m
        contact=a*f*f*K
        return np.array([2*inv*C,inv*P-b*X,-2*b*C,2*inv*dC,
                         inv*dP-b*dX-contact*X,-2*b*dC-2*contact*C]).ravel()
    free=np.array(initial)
    sol=solve_ivp(rhs,(t[0],t[-1]),np.vstack([free,np.zeros_like(free)]).ravel(),t_eval=t,
                  method="DOP853",rtol=3e-12,atol=1e-14,max_step=.0015)
    if not sol.success:
        raise RuntimeError(sol.message)
    X,C,P,dX,dC,dP=sol.y.reshape(6,len(q),len(t))
    a=np.exp(loga(t)); h=loga(t,1); m=mass(t); md=mass(t,1); v=md/100
    f,fp,fpp=clock(m/100+.5); fd=fp*v
    k2=(q[:,None]/a)**2
    stiffness=a*q[:,None]**2+a**3*m*m
    K=K[:,None]
    dXdt=2*dC/a**3
    dPdt=-2*stiffness*dC-2*a*f*f*K*C
    Xdt=2*C/a**3
    density=(dP/a**6+(m*m+k2)*dX)/2
    pressure=(dP/a**6-(m*m+k2/3)*dX)/2
    current=100*m*dX
    dr=((dPdt-6*h*dP)/a**6+(2*m*md-2*h*k2)*dX+(m*m+k2)*dXdt)/2
    dp=((dPdt-6*h*dP)/a**6-(2*m*md-2*h*k2/3)*dX-(m*m+k2/3)*dXdt)/2
    dj=100*(md*dX+m*dXdt)
    rc=K*f*f*X/(2*a*a)
    rcprime=K/a**2*(f*fd*X+.5*f*f*Xdt-h*f*f*X)
    jc=K*f*fp*X/a**2
    jcprime=K/a**2*(v*(fp*fp+f*fpp)*X+f*fp*Xdt-2*h*f*fp*X)
    all_fields=(density+rc,pressure-rc/3,current+jc,dr+rcprime,dp-rcprime/3,dj+jcprime)
    return {k:weights@field for k,field in zip(SOURCES,all_fields,strict=True)}


def old_slopes(repo):
    with (repo/OLD/"bulk/READOUTS.csv").open() as f:
        return {(r["grid"],float(r["step"]),int(r["cutoff"]),r["scheme"]):float(r["d_h"])
                for r in csv.DictReader(f) if r["arm"]=="quiet"}


def prediction(slopes,grid,step):
    gain={n:slopes[grid,step,n,"static"]-slopes[grid,step,n,"bare"] for n in CUTS}
    static={n:slopes[grid,step,n,"static"] for n in CUTS}
    base={n:slopes[grid,step,n,"bare"] for n in CUTS}
    return dict(h_slope_static=static,h_slope_bare=base,contact_direction_gain=gain,
                critical_reference_kappa={n:1-static[n]/gain[120] for n in CUTS},
                matched_h_slope={n:static[n]+ALPHA*gain[120] for n in CUTS},
                predicted_matched_60_minus_120=static[60]-static[120])


def prepare(repo,work,grid):
    src=module("prior_source",repo/OLD/"source/source.py")
    native=repo/NATIVE
    q,w,states,bg=src.read_native(native,grid)
    clock=src.ScalarClock(src.read_native(native,"fine")[-1][-1])
    saved=json.loads((repo/OLD/"source/SOURCES.json").read_bytes())
    out=work/"tables"; out.mkdir(exist_ok=True,parents=True)
    for cutoff in CUTS:
        row=src.integrate(native,grid,"quiet",cutoff,output=out)
        for k in ("bare","static"):
            expected=saved["runs"][f"{grid}_quiet_{cutoff}"]["files"][k]["readout"]
            for field,x in row["files"][k]["readout"].items():
                np.testing.assert_allclose(x,expected[field],rtol=3e-8,atol=2e-8)
    def read(cut,k):
        with np.load(out/f"{grid}_quiet_{cut}_{k}.npz",allow_pickle=False) as d:
            return dict(d)
    ref=read(120,"static")
    bare=read(120,"bare")
    extra=local_response(q,w,states["retain"],bg,clock,ref["s"])
    discrepancy={}
    for k in SOURCES:
        check=ALPHA*(ref[k]-bare[k])
        discrepancy[k]=float(np.max(abs(extra[k]-check))/max(np.max(abs(check)),1e-25))
        if discrepancy[k]>2e-5:
            raise ValueError(f"Independent local contact mismatch {grid} {k} {discrepancy[k]}")
    ward={}
    for cut in CUTS:
        data=read(cut,"static")
        np.testing.assert_array_equal(ref["s"],data["s"])
        new={k:data[k]+extra[k] for k in SOURCES}
        ward[cut]=float(np.max(abs(new["drho"]+3*data["h"]*(new["rho"]+new["p"])-new["j"]*data["v"])))
        if ward[cut]>5e-8: raise ValueError(f"Local Ward mismatch {ward[cut]}")
        np.savez_compressed(out/f"{grid}_quiet_{cut}_matched.npz",**new,**{k:data[k] for k in ("s","rad","drad")})
    ratio=moment(q,120)/moment(q,60)
    center=len(q)//2
    scalar_fit=1+ALPHA*ratio[center]
    mismatch=(scalar_fit-1)*moment(q,60)/(ALPHA*moment(q,120))-1
    return dict(q=q.tolist(),effective_kappa_60=(1+ALPHA*ratio).tolist(),
                constant_kappa_pivot_test=dict(q=float(q[center]),kappa=scalar_fit,
                    residual_relative_range=[float(mismatch.min()),float(mismatch.max())]),
                direct_contact_vs_archived_peak_relative_error=discrepancy,
                max_abs_ward=ward,table_pair_energy={n:float(read(n,"static")["daughter_comoving"][
                    np.argmin(abs(ref["s"]-2.10))]) for n in CUTS})


def worker(args):
    grid,cut,strength=args.worker
    solver=module("prior_live",args.repo/OLD/"bulk/live.py")
    solver.run_case(grid,"quiet",float(strength),int(cut),"matched",args.work/"tables",args.work/"runs",args.cache)


def readout(work,grid,cut,eps):
    signed=[]
    for tag in ("p","m"):
        p=work/"runs"/f"{grid}_quiet_{cut}_matched_{tag}{eps:g}.npz"
        with np.load(p,allow_pickle=False) as d:
            s=d["s"];vol=np.exp(3*(d["A"]-d["Aref"]))
            columns={"h":d["h"],"A":d["A"],"phi":d["phi"],"pair":vol*d["rad"]}
            vals={k:float(CubicSpline(s,y)(2.1)) for k,y in columns.items()}
            vals["h_linear"]=float(np.interp(2.10,s,d["h"]))
            signed.append(vals)
    return {k:(signed[0][k]-signed[1][k])/(2*eps) for k in signed[0]}


def run(args):
    authenticate(args.repo)
    slopes=old_slopes(args.repo)
    report=dict(source_revision=REV,alpha=ALPHA,reference_cutoff=120,readout_proper_time=2.10,
                arm="quiet",grids={g:{"sources":prepare(args.repo,args.work,g),"steps":{}} for g in GRIDS})
    jobs=[(g,c,l) for g in GRIDS for c in CUTS for e in EPS for l in (e,-e)]
    def launch(case):
        cmd=[sys.executable,"-IB",str(HERE/"gravity.py"),"--repo",str(args.repo),"--work",str(args.work),
             "--cache",str(args.cache),"--worker",*map(str,case)]
        subprocess.run(cmd,check=True,env={**os.environ,"OPENBLAS_NUM_THREADS":"1","OMP_NUM_THREADS":"1","MKL_NUM_THREADS":"1"})
    if not list((args.cache/"runs/shell_cache").glob("*.pkl")): launch(jobs.pop(0))
    with ThreadPoolExecutor(args.jobs) as pool:
        list(pool.map(launch,jobs))
    for g in GRIDS:
        for eps in EPS:
            row=prediction(slopes,g,eps)
            measured={c:readout(args.work,g,c,eps) for c in CUTS}
            row["observed"]=measured
            row["heldout_error"]={c:measured[c]["h"]-row["matched_h_slope"][c] for c in CUTS}
            row["observed_60_minus_120"]=measured[60]["h"]-measured[120]["h"]
            if not measured[60]["h"]<0<measured[120]["h"]:
                raise ValueError(f"Held-out opposite signs not recovered: {g} {eps} {measured}")
            if max(abs(e) for e in row["heldout_error"].values())>2e-8:
                raise ValueError(f"Held-out native slopes disagree: {row['heldout_error']}")
            for c in CUTS:
                np.testing.assert_allclose(measured[c]["pair"],report["grids"][g]["sources"]["table_pair_energy"][c],rtol=1e-7)
            report["grids"][g]["steps"][str(eps)]=row
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo",type=Path,required=True)
    ap.add_argument("--jobs",type=int,default=2)
    ap.add_argument("--work",type=Path,help="Preserve regenerated tables and sixteen native wall series here")
    ap.add_argument("--cache",type=Path,help="Optional trusted native static-shell cache")
    ap.add_argument("--output",type=Path,default=HERE/"GRAVITY.json")
    ap.add_argument("--worker",nargs=3,help=argparse.SUPPRESS)
    args=ap.parse_args()
    if args.worker: worker(args)
    else:
        with tempfile.TemporaryDirectory() as tmp:
            args.work=(args.work or Path(tmp)).resolve()
            args.cache=(args.cache or Path(tmp)/"cache").resolve()
            args.repo=args.repo.resolve()
            run(args)
