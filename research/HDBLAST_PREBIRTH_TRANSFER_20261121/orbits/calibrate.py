"""Refine source amplitudes for all broadly responsive wavelengths and check an inherited control."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import subprocess
import tempfile
import numpy as np

import sys
# Also work with isolated Python (-I).
sys.path.insert(0,str(Path(__file__).resolve().parent))
from run import GRIDS,HERE,make_trace,simulate

def standing_control(exe,prior,tmp,write,history,meta,label):
    # Equal starting Ey cosine, no incoming drive or feedback, source-free wave in a varying scale factor.
    grid=GRIDS["base"]; k=50.; amp=.005
    base=tmp/"prior_standing.txt"; write(history,meta,label,"0.01",grid,base)
    old=tmp/"original"
    subprocess.run(["g++","-std=c++17","-O3",str(prior/"particles/patch.cpp"),"-o",str(old)],check=True)
    archived=runpy.run_path(str(prior/"particles/run.py"))["run_case"](old,base,grid,k,amp,False,1)
    with base.open() as f:
        header=f.readline(); a0=float(header.split()[4]); nb=int(header.split()[1])
        births=[f.readline() for _ in range(nb)]
        states=np.loadtxt(f)
    eta=np.cumsum(states[:,3])-.5*states[:,3]
    factor=(a0/states[:,1])**2*amp
    ext=(factor*np.cos(k*a0*eta),factor*np.sin(k*a0*eta))
    out=tmp/"new_standing.txt"
    with out.open("w") as f:
        f.write(header); f.writelines(births)
        np.savetxt(f,np.column_stack((states,*ext)),fmt=["%.17g"]*6+["%d"]+["%.17g"]*2)
    new=simulate(exe,out,grid,k,1.,False)
    names={"cdf_gap":"transverse_cdf_distance","gas_conformal":"gas_conformal","stress_zz":"wall_axial_stress_share",
           "transverse_rms":"transverse_rms","supplied_wave_work":"spatial_field_work"}
    errors={key:abs(new[key]-archived[oldkey])/max(abs(archived[oldkey]),1e-10) for key,oldkey in names.items()}
    if max(errors.values())>5.e-10: raise ArithmeticError(errors)
    return {"previous":{v:archived[v] for v in names.values()},"new":{k:new[k] for k in names},"relative_differences":errors}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prior",type=Path,required=True); p.add_argument("--transfer",type=Path,required=True)
    p.add_argument("--orbit-product",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
    p.add_argument("--jobs",type=int,default=2)
    args=p.parse_args()
    prior=args.prior.resolve();meta=json.loads((prior/"history/markers.json").read_text())
    grid=json.loads(args.orbit_product.read_text())["runs"]
    write=runpy.run_path(str(prior/"particles/run.py"))["write_trace"]
    result={}
    with tempfile.TemporaryDirectory(prefix="hd-orbit-scale-") as folder:
        tmp=Path(folder);exe=tmp/"trajectory"
        subprocess.run(["g++","-std=c++17","-O3",str(HERE/"trajectory.cpp"),"-o",str(exe)],check=True)
        for label in ("dedicated_half","shared_half"):
            with np.load(prior/"history"/(label+".npz")) as f:h={k:f[k] for k in f.files}
            with np.load(args.transfer/(label+".npz")) as f:wave={k:f[k] for k in f.files}
            tested=[r for r in grid if r["history"]==label and r["resolution"]=="base" and r["feedback"]]
            order=sorted(tested,key=lambda x:(round(x["result"]["cdf_gap"],9),x["drive_amplitude"],x["mode_index"]))
            best={}
            for row in order:
                if row["drive_amplitude"]>0 and row["result"]["cdf_gap"]<.95:
                    best.setdefault(row["mode_index"],row["drive_amplitude"])
            choices=list(best)
            def batch(params):
                traces={}; spec=list(dict.fromkeys(params))
                for ix,amp,res in spec:
                    if (ix,res) not in traces:
                        traces[ix,res]=make_trace(tmp,write,h,wave,meta,label,"0.01",ix,res)
                def simulate_one(params):
                    ix,amp,res=params;k=wave["k"][ix]
                    stats=simulate(exe,traces[ix,res],GRIDS[res],k,amp,True)
                    return dict(mode_index=ix,k_physical=float(k),amplitude=float(amp),resolution=res,
                                result={key:stats[key] for key in ("cdf_gap","cdf_momentum_cut","mean_gamma","energy_over_thermal",
                                                                   "stress_eigen_min","stress_eigen_mid","stress_eigen_max",
                                                                   "mean_x_momentum_over_energy")})
                with ThreadPoolExecutor(args.jobs) as pool:return list(pool.map(simulate_one,spec))
            first=batch((ix,min(1.,best[ix]*2**p), "base") for ix in choices for p in np.linspace(-1,1,13))
            by_mode=[min((row for row in first if row['mode_index']==ix),key=lambda row:row['result']['cdf_gap']) for ix in choices]
            final=batch((r["mode_index"],min(1.,r["amplitude"]*2**p),"fine") for r in by_mode for p in (-1/12,0,1/12))
            preferred=min(final,key=lambda row:row['result']['cdf_gap'])
            check=batch((preferred['mode_index'],preferred['amplitude'],res) for res in ('spectral','angle'))
            result[label]=dict(responsive_broad_mode_indices=choices,coarse_amplitude_refinement=first,
                               fine_amplitude_refinement=final,best_fine=preferred,independent_resolution=check)
            if label=="dedicated_half":
                result["source_free_standing_wave_control"]=standing_control(exe,prior,tmp,write,h,meta,label)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
if __name__=="__main__":main()
