"""Compare the native gravitational response to the full and compressed quantum bath."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from scipy.interpolate import CubicSpline

HERE=Path(__file__).resolve().parent
PRIOR=Path("research/HDBLAST_SCALAR_CLOCK_GRAVITY_20261119/bulk")
BASE="67b2e6ab85bb14d09f248415fdf31b7b66326c31"
PINS={"live.py":"96b2613f413645136100b80a8ecb1e35f73bc164ca583dadbf09ec3fbb1b4409",
      "vendor/SHA256SUMS":"20119edbbe6be623a4b5bfec2bd3a469a3f078417910f81eb67a4877b01ce8ff",
      "READOUTS.csv":"acd967b4d2ed780a13acb2409a36a6e663561428886106ff133cf75fd8716fce"}
CASES=(("fine","quiet"),("coarse","quiet"),("fine","minus"))
MODELS=("low","N4","N64")
FIELDS=("h","A","phi","v","rho","p","j","pair","E")
TIMES=np.linspace(1.35,2.10,1201)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def authenticate(repo):
    prior=repo/PRIOR
    for file,expected in PINS.items():
        if digest(prior/file)!=expected:
            raise ValueError(f"Unexpected prior input {file}")
    for line in (prior/"vendor/SHA256SUMS").read_text().splitlines():
        expected,name=line.split()
        if digest(prior/"vendor"/name)!=expected: raise ValueError(name)
    return prior


def steps(grid,arm):
    return (.01,.005) if (grid,arm)==("fine","quiet") else (.01,)


def signed_name(grid,arm,model,strength):
    return f"{grid}_{arm}_120_{model}_{'p' if strength>0 else 'm'}{abs(strength):g}.npz"


def worker(prior,work,cache,case):
    grid,arm,model,strength=case
    spec=importlib.util.spec_from_file_location("native_wall",prior/"live.py")
    native=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(native)
    native.run_case(grid,arm,float(strength),120,model,work/"sources",work/"runs",cache)


def prepare(tables,work):
    # Keep the exact input bytes. Renaming only supplies the native loader's
    # established grid/arm/cutoff/scheme convention. Source files retain provenance.
    dest=work/"sources"; dest.mkdir(exist_ok=True,parents=True)
    hashes={}
    for grid,arm in CASES:
        for model in MODELS:
            input_name=f"{grid}_{arm}_{model}_static.npz"
            src=tables/input_name
            target=dest/f"{grid}_{arm}_120_{model}.npz"
            shutil.copyfile(src,target)
            hashes[input_name]=digest(target)
            with np.load(target,allow_pickle=False) as d:
                if any(np.any(d[k][d["s"]<1.34]!=0) for k in ("rho","p","j","rad","drho","dp","dj","drad")):
                    raise ValueError(f"New matter before gate {input_name}")
    local=dest/"fine_quiet_120_difference.npz"
    with np.load(dest/"fine_quiet_120_N4.npz",allow_pickle=False) as four, np.load(dest/"fine_quiet_120_N64.npz",allow_pickle=False) as full:
        np.testing.assert_array_equal(four["s"],full["s"])
        fields=("rho","p","j","rad","drho","dp","dj","drad")
        np.savez_compressed(local,s=four["s"],**{k:four[k]-full[k] for k in fields})
    return hashes


def expansion_from_trace(d,a,cut):
    """Integrate d(a²h)/ds using junction sources and wall fields (previous native trace diagnostic)."""
    rb=float(d["rb"])
    convert=float(d["b"])/rb
    rho,p,j=(convert*d[k] for k in ("rho","p","j"))
    phi=d["phi"]; sigma=d["sigma"]; grad=d["sigma1"]
    potential=.5*(phi**2-1)**2-2*(1-phi+phi**3/3)**2/3
    vacuum=sigma*sigma/18-grad*grad/24+potential/3
    stress=sigma*(rho-3*p)/36-rho*(rho+3*p)/36-j*(2*grad+j)/24
    rate=a**2*(rb**2*(vacuum+stress)-d["v"]**2/6)
    history=CubicSpline(d["s"],rate)
    boundary=float(CubicSpline(d["s"],a)(cut))
    return (a[0]**2*d["h"][0]+history.integrate(d["s"][0],cut))/boundary**2


def load_case(work,grid,arm,model,strength,pre):
    path=work/"runs"/signed_name(grid,arm,model,strength)
    with np.load(path,allow_pickle=False) as d:
        if str(d["source_sha256"])!=digest(work/"sources"/f"{grid}_{arm}_120_{model}.npz"):
            raise ValueError(f"Mismatched source {path.name}")
        if float(d["lam"])!=strength or not np.all(np.diff(d["s"])>0) or d["s"][-1]<TIMES[-1]:
            raise ValueError(f"Native time or strength {path.name}")
        subset=d["s"]<1.34
        control=np.stack([d[k][subset] for k in ("T","s","A","B","h","v","phi","rho","p","j")])
        if pre: np.testing.assert_array_equal(control,pre[0])
        else: pre.append(control)
        a=np.exp(d["A"]-float(d["Aref"]))
        data={k:(a**3*d["rad"] if k=="pair" else d[k]) for k in FIELDS}
        series={k:CubicSpline(d["s"],val)(TIMES) for k,val in data.items()}
        series["h_trace_cut"]=float(expansion_from_trace(d,a,TIMES[-1]))
        series["h_linear_cut"]=float(np.interp(2.10,d["s"],d["h"]))
        norms={k:float(np.max(d[k][d["s"]<=2.1])) for k in ("Hshell","Mshell","Hnear","Mnear")}
    return series,norms


def norm(y,t=TIMES):
    return float(np.sqrt(np.trapezoid(y*y,x=t)))


def relative(num,den):
    return None if den==0 else float(num/den)


def analyze(work,prior):
    with (prior/"READOUTS.csv").open() as f:
        old={(r["grid"],r["arm"],int(r["cutoff"]),float(r["step"])):r for r in csv.DictReader(f)
             if r["scheme"]=="static"}
    results={}
    for grid,arm in CASES:
        pre=[]
        group={}
        for eps in steps(grid,arm):
            models={}; constraints={}
            for model in MODELS:
                positive,cp=load_case(work,grid,arm,model,eps,pre)
                negative,cm=load_case(work,grid,arm,model,-eps,pre)
                models[model]={key:(positive[key]-negative[key])/(2*eps) for key in positive}
                constraints[model]={key:max(cp[key],cm[key]) for key in cp}
            low,four,exact=(models[k] for k in MODELS)
            endpoint={}; curves={}
            for field in FIELDS:
                hi=exact[field]-low[field]
                reproduced=four[field]-low[field]
                error=four[field]-exact[field]
                endpoint[field]={"low":float(low[field][-1]),"full":float(exact[field][-1]),
                    "four":float(four[field][-1]),"high_increment":float(hi[-1]),
                    "four_increment":float(reproduced[-1]),"four_minus_full":float(error[-1]),
                    "absolute_error_over_high":relative(abs(error[-1]),abs(hi[-1]))}
                n=norm(hi)
                curves[field]={"high_norm":n,"error_norm":norm(error),"relative_to_high":relative(norm(error),n),
                            "relative_to_full":relative(norm(error),norm(exact[field])),
                            "relative_to_high_stride_two":relative(norm(error[::2],TIMES[::2]),norm(hi[::2],TIMES[::2]))}
            high=exact["h"]-low["h"]; error=four["h"]-exact["h"]
            linhi=exact["h_linear_cut"]-low["h_linear_cut"]
            liner=four["h_linear_cut"]-exact["h_linear_cut"]
            ref=old[grid,arm,120,eps]; low_ref=old[grid,arm,60,eps]
            grab=lambda a:[float(v) for v in a[::8]]
            group[str(eps)]=dict(endpoint=endpoint,history_error=curves,constraint_maxima=constraints,
                interpolation=dict(linear_endpoint_high=float(linhi),linear_endpoint_error=float(liner),
                    cubic_minus_linear_high=float(high[-1]-linhi),cubic_minus_linear_error=float(error[-1]-liner)),
                trace_check=dict(high_from_trace=float(exact["h_trace_cut"]-low["h_trace_cut"]),
                    error_from_trace=float(four["h_trace_cut"]-exact["h_trace_cut"]),
                    high_metric_minus_trace=float(high[-1]-(exact["h_trace_cut"]-low["h_trace_cut"])),
                    error_metric_minus_trace=float(error[-1]-(four["h_trace_cut"]-exact["h_trace_cut"]))),
                archive_control=dict(zero_expansion_at_cut=float(ref["zero_h"]),
                    low_minus_previous_low=float(low["h"][-1]-float(low_ref["d_h"])),
                    full_minus_previous_unsplit_band=float(exact["h"][-1]-float(ref["d_h"]))),
                history=dict(s=grab(TIMES),high=grab(high),four_minus_low=grab(four["h"]-low["h"])))
        results[f"{grid}_{arm}"]=group
    residual={}
    pre=[]
    for eps in (1.,.5):
        pos,_=load_case(work,"fine","quiet","difference",eps,pre)
        neg,_=load_case(work,"fine","quiet","difference",-eps,pre)
        dh=(pos["h"]-neg["h"])/(2*eps)
        trace=(pos["h_trace_cut"]-neg["h_trace_cut"])/(2*eps)
        linear=(pos["h_linear_cut"]-neg["h_linear_cut"])/(2*eps)
        ref=results["fine_quiet"]["0.01"]
        residual[str(eps)]={"h":float(dh[-1]),"h_from_trace":float(trace),
            "h_via_linear_interpolation":float(linear),
            "relative_endpoint_to_high":relative(abs(dh[-1]),abs(ref["endpoint"]["h"]["high_increment"])),
            "relative_history_to_high":relative(norm(dh),ref["history_error"]["h"]["high_norm"]),
            "pair":float((pos["pair"][-1]-neg["pair"][-1])/(2*eps)),
            "delta_h_history":[float(x) for x in dh[::8]]}
    fine=results["fine_quiet"]
    one,two=[fine[str(e)]["endpoint"]["h"] for e in (.01,.005)]
    return dict(readout=2.10,native_runs=24,exploratory_native_runs=4,
                residual_control=dict(auxiliary_source_sha256=digest(work/"sources"/"fine_quiet_120_difference.npz"),
                    runs=residual,half_minus_full_h=residual["0.5"]["h"]-residual["1.0"]["h"]),relative_error_norm="L2 in native proper time, 1.35 to 2.10",
                matched_static_residual=0,results=results,
                step_half_minus_standard={k:two[k]-one[k] for k in ("low","full","four","high_increment","four_minus_full")})


def run(args):
    prior=authenticate(args.repo)
    if args.worker:
        worker(prior,args.work,args.cache,args.worker)
        return
    hashes=prepare(args.tables,args.work)
    cases=[(g,a,k,x) for g,a in CASES for e in steps(g,a) for k in MODELS for x in (e,-e)]
    cases.extend(("fine","quiet","difference",l) for l in (1.,-1.,.5,-.5))
    def launch(c):
        out=args.work/"runs"/signed_name(*c)
        if args.resume and out.exists():
            g,a,k,_=c
            with np.load(out,allow_pickle=False) as d:
                if str(d["source_sha256"])==digest(args.work/"sources"/f"{g}_{a}_120_{k}.npz"):return
        cmd=[sys.executable,"-IB",str(HERE/"gravity.py"),"--repo",str(args.repo),
             "--work",str(args.work),"--cache",str(args.cache),"--worker",*map(str,c)]
        subprocess.run(cmd,check=True,env={**os.environ,"OPENBLAS_NUM_THREADS":"1","OMP_NUM_THREADS":"1","MKL_NUM_THREADS":"1"})
    if not list((args.cache/"runs/shell_cache").glob("*.pkl")):
        launch(cases.pop(0))
    with ThreadPoolExecutor(args.jobs) as pool:
        list(pool.map(launch,cases))
    result={"upstream_revision":BASE,"source_sha256":hashes,**analyze(args.work,prior)}
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--repo",type=Path,required=True)
    ap.add_argument("--tables",type=Path,default=HERE.parent/"source/tables")
    ap.add_argument("--output",type=Path,default=HERE/"GRAVITY.json")
    ap.add_argument("--jobs",type=int,default=2)
    ap.add_argument("--work",type=Path,help="Keep source aliases and full native output here")
    ap.add_argument("--cache",type=Path,help="Optional trusted static geometry cache")
    ap.add_argument("--resume",action="store_true")
    ap.add_argument("--worker",nargs=4,help=argparse.SUPPRESS)
    args=ap.parse_args()
    with tempfile.TemporaryDirectory() as t:
        args.repo=args.repo.resolve()
        args.work=(args.work or Path(t)/"work").resolve()
        args.cache=(args.cache or Path(t)/"cache").resolve()
        run(args)
