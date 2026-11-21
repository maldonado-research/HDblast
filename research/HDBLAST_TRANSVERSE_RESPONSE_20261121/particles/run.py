"""Deterministic nonlinear particle and Fourier-Maxwell response on an expanding local patch."""
import argparse
import json
import math
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid

HERE = Path(__file__).resolve().parent
FIELDS = ("representatives gas_conformal field_conformal initial_seed_energy mean_field_work "
          "spatial_field_work particle_energy_residual_relative max_single_step_field_residual "
          "peak_spatial_field_rms peak_induced_field_rms peak_density_mode transverse_rms transverse_median "
          "transverse_p90 transverse_force_integral_rms birth_tube_share wall_axial_stress_share "
          "expansion_mass_term injected_conformal transverse_cdf_distance").split()


def birth_directions(c, d, target, n):
    # Piecewise-linear rate with exact boundaries. Its integrated CDF is quadratic on each piece.
    if target <= 0 or d == 0:
        return np.linspace(-1 + 1/n, 1 - 1/n, n)
    cuts = np.unique(np.clip([-1, (c-target)/d, c/d, (c+target)/d, 1], -1, 1))
    w = np.minimum(1, np.abs(c - d*cuts) / target)
    areas = .5 * (w[1:] + w[:-1]) * np.diff(cuts)
    cumul = np.r_[0., np.cumsum(areas)]
    quantiles = (np.arange(n)+.5)/n*cumul[-1]
    cells = np.searchsorted(cumul[1:], quantiles, side="right")
    widths = np.diff(cuts)[cells]
    area = (quantiles-cumul[cells])/widths
    initial, final = w[cells], w[cells+1]
    coordinate = 2*area/(initial + np.sqrt(np.maximum(0., initial**2+2*(final-initial)*area)))
    return cuts[cells] + widths*coordinate


def write_trace(history, meta, label, marker, grid, target):
    cohorts, nangle, naz, nx, harmonics, nsteps = grid
    j = meta["runs"][label]["markers"][marker]
    a = history
    interp = lambda k, t: np.interp(t, a["t"], a[k])
    weights = j["D"]/cohorts
    times = np.interp(weights*(np.arange(cohorts)+.5), a["D"], a["t"])
    edges = np.unique(np.r_[np.linspace(a["t"][0], j["time"], nsteps+1), times])
    centers = .5*(edges[1:]+edges[:-1])
    mids = interp("a", centers)
    rows = []
    for t in times:
        scale = interp("a",t)
        pstar = math.sqrt(max(0.,interp("R",t)**2-1.))
        for z in birth_directions(interp("c",t),interp("d",t),interp("target",t),nangle):
            for ai in range(naz):
                rows.append([t,scale,pstar*z,pstar*math.sqrt(max(0.,1-z*z)),2*math.pi*(ai+.5)/naz])
    states = np.column_stack((interp("a",edges[:-1]),mids,interp("a",edges[1:]),
                              np.diff(edges)/mids,interp("X",centers),interp("B",centers),
                              nangle*naz*np.searchsorted(times,edges[:-1],side='right')))
    with target.open("w") as f:
        print(len(states),len(rows),nangle*naz,weights/(nangle*naz),a['a'][0],j['theta'], file=f)
        np.savetxt(f, rows, fmt="%.17g")
        np.savetxt(f,states,fmt=["%.17g"]*6+["%d"])


def run_case(exe, path, grid, k, seed, coupled, polarization):
    cmd = [str(exe),str(path),str(grid[4]),str(grid[3]),str(k),str(seed),str(int(coupled)),str(polarization)]
    values = [float(x) for x in subprocess.check_output(cmd, text=True).split()]
    if len(values)!=len(FIELDS) or not np.isfinite(values).all():
        raise ValueError("Invalid nonlinear result")
    r = dict(zip(FIELDS,values))
    if abs(r["particle_energy_residual_relative"])>5e-10 or r["max_single_step_field_residual"]>1e-9:
        raise ValueError("Nonlinear step violates matter or field energy ledger")
    r["integrated_spatial_balance_relative"]=(
        (r["field_conformal"]-r["initial_seed_energy"]+r["spatial_field_work"])/r["gas_conformal"]
        if coupled else None)
    if coupled and abs(r["integrated_spatial_balance_relative"])>5e-10:
        raise ValueError("Integrated wave-particle balance")
    return r


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--history",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--pilot",action="store_true")
    p.add_argument("--jobs",type=int,default=2)
    args=p.parse_args()
    meta=json.loads((args.history/"markers.json").read_bytes())
    grids={"base":(36,2,4,16,3,160),"fine":(72,2,4,24,3,320),
           "spectral":(72,2,4,32,7,320),"time":(36,2,4,16,3,320),
           "angular":(36,4,8,16,3,160)}
    variants=[]
    def add(label,k,seed,coupled=True,grid="base",pol=0):
        variants.append((label,k,seed,coupled,grid,pol))
    for channel in ("dedicated","shared"):
        name=channel+"_native"
        for k in (5.,50.,250.,1000.):
            for field in (1e-6,.005,.04):
                add(name,k,field)
        for field in (0.,.005,.04):
            add(name,50.,field,False)
        for field in (.005,.04):
            for grid in ("fine","spectral","time","angular"):
                add(channel+"_half",50.,field,True,grid)
            add(name,50.,field,True,"base",1)
        add(name,50.,.04,True,"base",2)
        add(name,50.,0.)
    if args.pilot:
        variants=[("dedicated_native",50.,v,True,"base",0) for v in (0,1e-6,.005,.04)]
    with tempfile.TemporaryDirectory(prefix="hd-transverse-particles-") as folder:
        tmp=Path(folder);exe=tmp/"patch"
        subprocess.run(["g++","-std=c++17","-O3",str(HERE/"patch.cpp"),"-o",str(exe)],check=True)
        for label,g in sorted({(v[0],v[4]) for v in variants}):
            with np.load(args.history/(label+".npz")) as data:
                write_trace(data,meta,label,"0.01",grids[g],tmp/(label+g))
        def one(v):
            label,k,seed,coupled,g,pol=v
            result=run_case(exe,tmp/(label+g),grids[g],k,seed,coupled,pol)
            return dict(history=label,k_physical_opening=k,initial_peak_field=seed,feedback=coupled,
                        resolution=g,polarization=["By","Ey","Ez"][pol],result=result)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            runs=list(pool.map(one,variants))
    args.output.write_text(json.dumps(dict(grids=grids,runs=runs),indent=2,sort_keys=True,allow_nan=False)+'\n')

if __name__ == "__main__":
    main()
