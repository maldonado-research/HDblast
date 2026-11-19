"""Reconstruct equal-proper-time responses from independently evolved signed-strength runs."""
import argparse
import csv
import io
import json
from pathlib import Path
import zipfile

import numpy as np
from scipy.interpolate import CubicSpline

HERE = Path(__file__).resolve().parent
ARMS = ("quiet","plus","minus","far")
GRIDS = ("coarse","fine")
KEYS = ("T","A","a","phi","h","v","rho","p","j","E","rad","work","pressure_work","Wy_over_H02")
CLOCKS = [1.34,1.45,1.55,1.65,1.75,1.85,1.95,2.05,2.10]
CONTROL = {"T":"T","s":"H0tau","A":"A_b","B":"B_b","phi":"phi_b","h":"H_over_H0","v":"v_over_H0",
           "rho":"rho_hat","p":"p_hat","j":"j_hat","work":"scalar_work","pressure_work":"pressure_work",
           "E":"comoving_energy","Hnear":"Hmax_near","Mnear":"Mmax_near","Hshell":"Hshell","Mshell":"Mshell"}


def load(path):
    with np.load(path,allow_pickle=False) as d:
        data={k:d[k] for k in d}
    data["a"] = np.exp(data["A"]-data["Aref"])
    return data


def readouts(data, clocks):
    if min(clocks)<min(data["s"]) or max(clocks)>max(data["s"]):
        raise ValueError("Missing readout coverage")
    return {k:CubicSpline(data["s"],data[k])(clocks) for k in KEYS}


def check_baseline(data,grid,arm,regenerated=False):
    with zipfile.ZipFile(HERE/"vendor/references.zip") as z:
        ref=np.genfromtxt(io.BytesIO(z.read(f"{grid}_{arm}.csv")),delimiter=",",names=True)
    indices=np.rint((ref["T"]-data["T"][0])/float(data["dt"])).astype(int)
    valid=(indices>=0)&(indices<len(data["T"]))
    ix=indices[valid]
    errors={}
    for new,old in CONTROL.items():
        test=ref[old][valid]
        mask=np.isfinite(test)
        errors[new]=float(np.max(abs(data[new][ix[mask]]-test[mask])))
    # Regenerated PDE trajectories inherit platform-dependent cancellation in
    # momentum stencils and in large oscillatory currents. Raw archived arrays
    # have only the original CSV writing error.
    limits={k:5e-9 for k in errors}
    if regenerated:
        limits.update(A=5e-10,B=5e-10,phi=5e-10,h=3e-9,v=3e-9,rho=2e-6,p=4e-6,j=1e-5,
                      E=2e-6,work=2e-6,pressure_work=5e-7,Hnear=2e-7,Mnear=2e-7,
                      Hshell=2e-7,Mshell=2e-7)
    if any(error>limits[k] for k,error in errors.items()):
        raise ValueError(f"Zero source fails archive check {grid} {arm}: {errors}")
    return errors


def summarize(run_dir, out_dir, regenerated=False):
    out_dir.mkdir(exist_ok=True,parents=True)
    baselines={}
    rows=[]
    diagnostic_rows=[]
    histories=[]
    for g in GRIDS:
        for arm in ARMS:
            zero=load(run_dir/f"{g}_{arm}_zero.npz")
            baselines[f"{g}_{arm}"]=check_baseline(zero,g,arm,regenerated)
            z=readouts(zero,CLOCKS)
            for cutoff in (60,120):
                for scheme in ("bare","static"):
                    for strength in ((.01,.005) if arm=="quiet" or (arm,cutoff,scheme)==("far",120,"static") else (.01,)):
                        prefix=f"{g}_{arm}_{cutoff}_{scheme}_"
                        plus=load(run_dir/(prefix+f"p{strength:g}.npz"))
                        minus=load(run_dir/(prefix+f"m{strength:g}.npz"))
                        pos=readouts(plus,CLOCKS)
                        neg=readouts(minus,CLOCKS)
                        D={k:(pos[k]-neg[k])/(2*strength) for k in KEYS}
                        info=dict(grid=g,arm=arm,cutoff=cutoff,scheme=scheme,step=strength)
                        for i,s in enumerate(CLOCKS):
                            row=dict(**info,s=s,**{f"d_{k}":float(D[k][i]) for k in KEYS},
                                     d_h_over_h=float(D["h"][i]/z["h"][i]),
                                     d_a_over_a=float(D["a"][i]/z["a"][i]),
                                     d_planar_area_over_area=float(2*D["A"][i]),
                                     d_volume_over_volume=float(3*D["A"][i]),
                                     physical_h_over_h=float((pos["h"][i]-z["h"][i])/z["h"][i]),
                                     zero_h=float(z["h"][i]),zero_phi=float(z["phi"][i]))
                            histories.append(row)
                            if s==2.10:
                                rows.append(row.copy())
                        # Before any gate opens, all six equations and native matter coincide pointwise.
                        n=np.count_nonzero(zero["s"]<1.345)
                        err=max(float(np.max(np.abs(zero[k][:n]-signed[k][:n]))) for k in
                                ("A","B","phi","h","v","j","rho","p","junction","junction_T","inside_phi")
                                for signed in (plus,minus))
                        if err!=0.0:
                            raise ValueError(f"Pre-gate control changed: {prefix} {err}")
                        diag=dict(**info,preturn_max_abs=err,
                                  h_step_slope_cubic=float(D["h"][-1]),
                                  h_step_slope_linear=float((np.interp(2.1,plus["s"],plus["h"])-
                                                             np.interp(2.1,minus["s"],minus["h"]))/(2*strength)))
                        # Compute absolute maxima before the common cut for each signed solution separately.
                        for nm in ("Hnear","Mnear","Hshell","Mshell"):
                            for label,dat in (("zero",zero),("plus",plus),("minus",minus)):
                                diag[nm+"_"+label]=float(np.max(dat[nm][dat["s"]<=2.10]))
                        # Signed closure as energy equation in the integrated live system.
                        for label,dat in (("plus",plus),("minus",minus)):
                            ledger=(dat["E"]-dat["E"][0])+(dat["pressure_work"]-dat["pressure_work"][0])-(dat["work"]-dat["work"][0])
                            mask=dat["s"]<=2.10
                            diag["ledger_max_"+label]=float(np.max(abs(ledger[mask])))
                        diagnostic_rows.append(diag)
    for group in ((g,c,k,e) for g in GRIDS for c in (60,120) for k in ("bare","static") for e in (.01,.005)):
        selected=[r for r in rows if (r["grid"],r["cutoff"],r["scheme"],r["step"])==group]
        q=next(r for r in selected if r["arm"]=="quiet")
        for r in selected:
            r["d_h_over_h_minus_quiet"]=r["d_h_over_h"]-q["d_h_over_h"]
            r["d_h_minus_quiet"]=r["d_h"]-q["d_h"]
    for name,tab in (("READOUTS",rows),("HISTORIES",histories),("DIAGNOSTICS",diagnostic_rows)):
        with (out_dir/(name+".csv")).open("w") as f:
            writer=csv.DictWriter(f,fieldnames=list(tab[0]),lineterminator="\n")
            writer.writeheader()
            writer.writerows(tab)
    (out_dir/"BASELINE_CHECK.json").write_text(json.dumps(baselines,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--runs",type=Path,default=HERE/"runs")
    ap.add_argument("--output",type=Path,default=HERE)
    ap.add_argument("--regenerated",action="store_true",help="Use architecture-aware field tolerances on new native integrations")
    a=ap.parse_args()
    summarize(a.runs,a.output,a.regenerated)
