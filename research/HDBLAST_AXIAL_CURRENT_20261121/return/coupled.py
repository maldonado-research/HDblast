"""Iterate free Maxwell propagation in the curved collar with the three-current wall."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import tempfile
import shutil

import numpy as np


HERE=Path(__file__).resolve().parent
response=runpy.run_path(str(HERE/"response.py"))
wall=runpy.run_path(str(HERE.parent/"live/evolve.py"))


def readout(row):
    return {key:row[key] for key in ("events","hit","mean_encounters","tv20","end_E","first_field_reversals",
        "rest_pressures_over_energy","gas","max_field_residual_over_scale","matter_residual_over_scale")}


def terminal_rms(t,e):
    mask=t>=t[-1]-.01
    return float(np.sqrt(np.trapezoid(e[mask]**2,t[mask])/(t[-1]-t[mask][0])))


def check_resolution(g, temp, t, delta, standard):
    controls=[]
    for refinement,cfl in ((1,.03),(2,.06)):
        field,ledger=response["calculate"](g,temp,t,delta,refinement,cfl)
        controls.append(dict(cells=ledger["cells"],cfl=cfl,
                             peak_absolute_incoming_difference=np.max(np.abs(field-standard),axis=0).tolist(),
                             peak_bulk_balance_residual_over_scale=ledger["peak_bulk_balance_residual_over_scale"]))
    return controls


def check_coupled_mesh(ctx,g,work,task,t,delta,initial_inflow):
    """Increase the Maxwell grid to 12 cells and repeat the charged-particle simulation."""
    old=initial_inflow
    extra,ledger=response["calculate"](g,work,t,delta,refine=2)
    sweeps=[]
    for _ in range(8):
        row=wall["run_case"](ctx,task,work,work,{"t":t,"delta_incoming":extra})
        with np.load(work/(row["label"]+".npz")) as sample:
            new,_=response["calculate"](g,work,t,sample["delta_E"],refine=2)
            terminal=terminal_rms(t,sample["E"][:,2])
        residual=np.max(abs(new-extra),axis=0)
        sweeps.append(dict(cells=ledger["cells"],inflow_change_from_previous=np.max(abs(extra-old),axis=0).tolist(),
                           inflow_residual=residual.tolist(),wall=readout(row),terminal_window_axial_E_rms=terminal))
        if max(residual)<1e-14:
            break
        old=extra
        extra=new
    return sweeps


def sweep(polarization,scattering,history,live_output,geometry,out,jobs,extended):
    reference_file, snapshots = ("CROSSINGS.json", "cross-snapshots") if extended else ("LIVE.json", "snapshots")
    refs=json.loads((live_output/reference_file).read_text())["runs"]
    with np.load(history) as z:h=dict(z)
    orbits=runpy.run_path(str(polarization/"orbits/pairs.py"))
    coarse=runpy.run_path(str(polarization/"orbits/coarse.py"))
    live=runpy.run_path(str(scattering/"live/evolve.py"))
    coarse["thermal_cells"].__globals__["quad"]=live["thermal_quadrature"]
    with tempfile.TemporaryDirectory(prefix="self-consistent-collar-") as tmp:
        root=Path(tmp)
        g=response["prepare"](geometry,root)
        ctx=(h,orbits,coarse,live,wall["build"](polarization,scattering,root))
        test_t=g["t"]
        test_inputs=np.zeros((len(test_t),3))
        test_inputs[:,2]=g["drive"]
        controls=[]
        for refine,cfl in ((1,.06),(1,.03),(2,.06)):
            ref_i,check=response["calculate"](g,root,test_t,test_inputs,refine,cfl)
            check["peak_error_against_independent_solver"]=float(max(abs(ref_i[:,2]-g["reference"])))
            controls.append(check)
        # Independent zero input is also propagated in both inactive components of these controls.
        labels=[(6.45,12,1024,18431,True), (6.67,12,2048,18431,True),
                (6.67,12,2048,18431,False), (6.67,13,2048,70163,True),
                (6.67,14,4096,91283,True)]
        if extended:
            labels=[(7.05,13,4096,70163,True), (7.05,13,8192,70163,True), (7.05,14,8192,91283,True)]
        def one(spec):
            tm,power,steps,seed,pole=spec
            task=tm,power,steps,seed,"particle",pole
            source=next(r for r in refs if abs(r["time"]-tm)<1e-9 and r["pairs"]==2**power
                and r["steps"]==steps and r["seed"]==seed and r["mode"]=="particle" and r["pole"]==pole)
            label=source["label"]
            with tempfile.TemporaryDirectory(dir=root) as temp:
                temp=Path(temp)
                for file in ("geometry","collar"):
                    (temp/file).symlink_to(root/file)
                archived=live_output/snapshots/(label+".npz")
                snapshot=temp/(label+".npz")
                shutil.copyfile(archived,snapshot)
                iterations=[]
                with np.load(snapshot) as initial:
                    fixed=initial["E"]
                old_i=np.zeros_like(fixed)
                previous=fixed
                for iteration in range(8 if extended else 3):
                    with np.load(snapshot) as data:
                        t=data["t"]; delta=data["delta_E"]
                    reflected,ledger=response["calculate"](g,temp,t,delta)
                    # Force the numerical wall with the last iteration's physical incident characteristic.
                    data=wall["run_case"](ctx,task,temp,temp,{"t":t,"delta_incoming":reflected})
                    with np.load(snapshot) as fresh:
                        E=fresh["E"]
                    iterations.append(dict(number=iteration+1,peak_inflow_update=np.max(abs(reflected-old_i),axis=0).tolist(),
                        peak_field_update=np.max(abs(E-previous),axis=0).tolist(),
                        peak_field_change_from_open=np.max(abs(E-fixed),axis=0).tolist(),
                        collar_ledger=ledger,wall=readout(data)))
                    done = iteration >= 2 and np.max(abs(reflected-old_i)) < 1e-14
                    previous=E
                    old_i=reflected
                    if done:
                        break
                with np.load(snapshot) as current:
                    remaining,_=response["calculate"](g,temp,current["t"],current["delta_E"])
                    residual=np.max(abs(remaining-old_i),axis=0)
                    resolution=check_resolution(g,temp,current["t"],current["delta_E"],remaining) if extended else []
                    if out is not None:
                        out.mkdir(parents=True,exist_ok=True)
                        np.savez_compressed(out/(label+".npz"),t=current["t"],E=current["E"],B=current["B"],
                                            p=current["p"],sign=current["sign"],delta_incoming=old_i)
                if extended:
                    wall["add_common_target"]([data],temp,scattering)
                    iterations[-1]["wall"]["common_gas_tv20"]=data["common_gas_tv20"]
                    iterations[-1]["wall"]["terminal_window_axial_E_rms"]=terminal_rms(t,previous[:,2])
                initial_readout=readout(source)
                if extended:
                    initial_readout["common_gas_tv20"]=source["common_gas_tv20"]
                    initial_readout["terminal_window_axial_E_rms"]=terminal_rms(t,fixed[:,2])
                record=dict(label=label,open=initial_readout,iterations=iterations,final_inflow_residual=residual.tolist())
                if extended:
                    record["final_refinement_controls"]=resolution
                    if power==14:
                        with np.load(snapshot) as replay:
                            final_delta=replay["delta_E"]
                        record["coupled_mesh_control"]=check_coupled_mesh(ctx,g,temp,task,t,final_delta,old_i)
                return record
        with ThreadPoolExecutor(jobs) as pool:
            cases=list(pool.map(one,labels))
        return dict(cases=cases,manufactured=controls,
                    control_peak_ratio=float(np.max(abs(g["reference"]))/np.max(abs(g["drive"]))),
                    early_characteristic_bound=response["continuum_cap"](
                        {**g,"geo":g["geo"][g["t"]<=6.4505]},6.45-6.248),
                    characteristic_bound=response["continuum_cap"](g,max(g["t"])-min(g["t"])))


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for key in ("polarization","scattering","history","live-output","geometry","output"):
        p.add_argument("--"+key,type=Path,required=True)
    p.add_argument("--extended",action="store_true")
    p.add_argument("--snapshots",type=Path)
    p.add_argument("--jobs",type=int,default=2)
    a=p.parse_args()
    results=sweep(a.polarization.resolve(),a.scattering.resolve(),a.history,a.live_output,a.geometry,a.snapshots,a.jobs,a.extended)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(results,indent=2,sort_keys=True,allow_nan=False)+"\n")
