"""Continue the same pulse until the kinetic axial field reverses again and most reference work is done."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import tempfile

import numpy as np


HERE=Path(__file__).resolve().parent
base=runpy.run_path(str(HERE/"evolve.py"))


def follow(args):
    with np.load(args.history) as data:
        h=dict(data)
    if h["t"][-1]<7.05 or not np.allclose(h["K"],2.54,rtol=0.,atol=2.e-12):
        raise ValueError("Need an actual history at fixed K through wall time 7.05")
    scatter=runpy.run_path(str(args.scattering/"live/evolve.py"))
    orbits=runpy.run_path(str(args.polarization/"orbits/pairs.py"))
    coarse=runpy.run_path(str(args.polarization/"orbits/coarse.py"))
    coarse["thermal_cells"].__globals__["quad"]=scatter["thermal_quadrature"]
    snapshots=args.output.parent/"cross-snapshots"
    snapshots.mkdir(exist_ok=True,parents=True)
    tasks=[]
    for tm in (6.73,6.90,7.05):
        for power,steps,seed in ((12,4096,18431),(13,4096,70163),(13,8192,70163),(14,8192,91283)):
            for mode in ("fluid","particle"):
                tasks.append((tm,power,steps,seed,mode,True))
            if power==13 and steps==4096:
                tasks.append((tm,power,steps,seed,"particle",False))
    with tempfile.TemporaryDirectory(prefix="second-reversal-") as wd:
        work=Path(wd)
        context=h,orbits,coarse,scatter,base["build"](args.polarization,args.scattering,work)
        with ThreadPoolExecutor(args.jobs) as pool:
            rows=list(pool.map(lambda task:base["run_case"](context,task,work,snapshots),tasks))
    base["add_common_target"](rows,snapshots,args.scattering)
    for r in rows:
        with np.load(snapshots/(r["label"]+".npz")) as sample:
            t,e=sample["t"],sample["E"][:,2]
            mask=t>=t[-1]-.01
            r["terminal_window_axial_E_rms"]=float(np.sqrt(np.trapezoid(e[mask]**2,t[mask])/(t[-1]-t[mask][0])))
            r["terminal_window_axial_E_rms_over_first_field"]=r["terminal_window_axial_E_rms"]/abs(e[0])
        # The individually timed late crossings of a tiny fluctuating field vary
        # across particle samples. Record the first two and the terminal RMS.
        for key in ("first_field_reversals","first_reference_field_reversals","first_current_reversals","first_fluid_current_reversals"):
            r[key]=r[key][:2]
    args.output.write_text(json.dumps(dict(runs=rows),indent=2,sort_keys=True,allow_nan=False)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for k in ("scattering","polarization","history","output"):
        parser.add_argument("--"+k,type=Path,required=True)
    parser.add_argument("--jobs",type=int,default=2)
    follow(parser.parse_args())
