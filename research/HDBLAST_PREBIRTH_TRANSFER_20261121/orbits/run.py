"""Transport the charged daughters through measured prebirth wave-transfer histories."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import subprocess
import tempfile

import numpy as np
from scipy.special import kv

HERE = Path(__file__).resolve().parent
FIELDS = (
    "representatives cdf_gap cdf_momentum_cut transverse_rms transverse_median "
    "stress_zz stress_xx stress_yy stress_yz gas_conformal injected expansion_work "
    "axial_work supplied_wave_work induced_wave_work induced_field_energy induced_leakage "
    "peak_induced_field_rms peak_first_density_mode particle_energy_residual_rel induced_energy_residual_rel "
    "stress_eigen_min stress_eigen_mid stress_eigen_max mean_x_momentum_over_energy mean_gamma energy_over_thermal"
).split()
GRIDS = {
    "screen": (18, 2, 4, 8, 1, 96),
    "base": (36, 2, 4, 16, 3, 160),
    "fine": (72, 2, 4, 24, 3, 320),
    "spectral": (72, 2, 4, 32, 7, 320),
    "angle": (48, 4, 8, 24, 3, 240),
}
AMPLITUDES = (0., .005, .01, .02, .04, .08, .16, .32, .64, 1.)

def simulate(exe, trace, grid, k, amplitude, feedback):
    cmd = list(map(str, [exe, trace, grid[4], grid[3], float(k), float(amplitude), int(feedback)]))
    values = [float(x) for x in subprocess.check_output(cmd, text=True).split()]
    if len(values) != len(FIELDS) or not np.isfinite(values).all():
        raise ValueError("Incomplete or nonfinite trajectory output")
    row = dict(zip(FIELDS, values))
    if max(abs(row[s]) for s in ("particle_energy_residual_rel", "induced_energy_residual_rel")) > 5e-9:
        raise ArithmeticError("Orbit energy ledger failed")
    eig = [row["stress_eigen_"+s] for s in ("min", "mid", "max")]
    if abs(sum(eig)-1)>1e-10 or eig[0] < -1e-12 or eig[2] > 1+1e-12:
        raise ArithmeticError("Charge-summed stress failed")
    return row

def make_trace(tmp, write_trace, reference, waves, meta, label, marker, index, grid_name):
    dest = tmp / f"{label}-{marker}-{index}-{grid_name}.txt"
    if dest.exists():
        return dest
    base = tmp / f"{label}-{marker}-{grid_name}-base.txt"
    if not base.exists():
        write_trace(reference, meta, label, marker, GRIDS[grid_name], base)
    with base.open() as f:
        header = f.readline()
        steps, n_birth = map(int, header.split()[:2])
        birth_lines = [f.readline() for _ in range(n_birth)]
        data = np.loadtxt(f, ndmin=2)
    end = meta["runs"][label]["markers"][marker]["time"]
    # Input stores Δη=Δt/a_mid and splits the grid at each birth.
    dt = data[:, 3]*data[:, 1]
    mids = reference["t"][0] + np.cumsum(dt) - .5*dt
    if (len(mids) != steps or abs(reference["t"][0]+dt.sum()-end)>1e-10 or
            waves["t"][0]>mids[0]+1e-10 or waves["t"][-1]<mids[-1]-1e-10):
        raise ValueError("Transferred waveform does not cover the requested particle times")
    ext = [np.interp(mids,waves["t"],waves[key][index]) for key in ("E_cos","B_sin")]
    with dest.open("w") as f:
        f.write(header)
        f.writelines(birth_lines)
        np.savetxt(f, np.column_stack((data,*ext)), fmt=["%.17g"]*6+["%d"]+["%.17g"]*2)
    return dest

def check_uniform_drive(exe,tmp):
    # Polarized uniform electric pulse. Zero background and B produce a known impulse for every initial phase.
    nx=32
    h=.002
    target=tmp/"electric.txt"
    target.write_text(f"1 1 1 1 1 1\n0 1 0 0 0\n1 1 1 {h} 0 0 1 .0001 0\n")
    row=simulate(exe,target,(1,1,1,nx,1,1),0,1.,False)
    phi=2*np.pi*(np.arange(nx)+.61803398874989485)/nx
    u=1.e7*h*.0001*np.cos(phi)
    g=np.hypot(1.,u)
    r=np.sort(abs(u))
    F=1-(1+r*r)*kv(2,np.hypot(1.,r))/kv(2,1.)
    ks=max(np.max(abs(F-np.arange(nx)/nx)),np.max(abs(F-(1+np.arange(nx))/nx)))
    exact={"transverse_rms":np.sqrt(np.mean(u*u)),"mean_gamma":g.mean(),"stress_eigen_max":1.,
           "energy_over_thermal":g.mean()/(3+kv(1,1.)/kv(2,1.)),"cdf_gap":ks}
    err=max(abs(row[key]-value) for key,value in exact.items())
    if err>2e-12:
        raise ArithmeticError("Analytic transverse electric impulse failed")
    return float(err)

def analyze(prior, profile_dir, output, jobs):
    history=prior/"history"
    meta=json.loads((history/"markers.json").read_text())
    reference={}; waves={}
    for label in sorted(meta["runs"]):
        with np.load(history/(label+".npz")) as tab:
            reference[label] = {key: tab[key] for key in tab.files}
        with np.load(profile_dir/(label+".npz")) as tab:
            waves[label]={key: tab[key] for key in ("t","k","E_cos","B_sin")}
        wave=waves[label]
        if (len(wave["t"])<2 or not (np.diff(wave["t"])>0).all() or len(np.unique(wave["k"]))!=len(wave["k"])
                or not np.isfinite(wave["k"]).all() or np.any(wave["k"]<=0)
                or any(wave[key].shape!=(len(wave["k"]),len(wave["t"])) or not np.isfinite(wave[key]).all()
                       for key in ("E_cos","B_sin"))):
            raise ValueError("Invalid wave profile for "+label)
    helper=runpy.run_path(str(prior/"particles/run.py"))["write_trace"]
    rows=[]
    with tempfile.TemporaryDirectory(prefix="hd-prebirth-orbits-") as folder:
        tmp=Path(folder);exe=tmp/"trajectory"
        subprocess.run(["g++","-std=c++17","-O3",str(HERE/"trajectory.cpp"),"-o",str(exe)],check=True)
        analytic=check_uniform_drive(exe,tmp)
        def batch(variants):
            variants=list(dict.fromkeys(variants))
            traces={}
            for label,mark,ix,amp,fb,res in variants:
                key=(label,mark,ix,res)
                if key not in traces:
                    traces[key]=make_trace(tmp,helper,reference[label],waves[label],meta,label,mark,ix,res)
            def evaluate(v):
                label,mark,ix,amp,fb,res=v
                k=float(waves[label]["k"][ix])
                stats=simulate(exe,traces[label,mark,ix,res],GRIDS[res],k,amp,fb)
                return dict(history=label,marker=mark,mode_index=ix,k_physical_birth=k,
                            drive_amplitude=amp,feedback=fb,resolution=res,result=stats)
            with ThreadPoolExecutor(max_workers=jobs) as pool:
                result=list(pool.map(evaluate,variants))
            rows.extend(result)
            return result
        # Screen a predeclared amplitude grid and every supplied mode on the finer background.
        tests=[]
        for label in ("dedicated_half","shared_half"):
            for ix in range(len(waves[label]["k"])):
                for amp in AMPLITUDES:
                    if amp==0 and ix!=0:continue
                    tests.append((label,"0.01",ix,amp,False,"screen"))
        screen=batch(tests)
        base=[]
        for label in ("dedicated_half","shared_half"):
            subset=[r for r in screen if r["history"]==label]
            # Recalculate each mode's best screen amplitude and the three strongest overall
            # candidates. Runners-up at neighboring amplitudes reduce grid-selection artifacts.
            choose=lambda r:(round(r["result"]["cdf_gap"],9),r["drive_amplitude"],r["mode_index"])
            winners=[min((r for r in subset if r["mode_index"]==ix),key=choose)
                     for ix in range(len(waves[label]["k"]))]
            for r in winners+sorted(subset,key=choose)[:3]:
                ix=r["mode_index"]; j=AMPLITUDES.index(r["drive_amplitude"])
                for amp in AMPLITUDES[max(0,j-1):min(len(AMPLITUDES),j+2)]:
                    base.append((label,"0.01",ix,amp,True,"base"))
            base.append((label,"0.01",0,0.,True,"base"))
        base_results=batch(base)
        details=[]; bests={}
        for label in ("dedicated_half","shared_half"):
            subset=[r for r in base_results if r["history"]==label]
            best=min(subset,key=lambda r:(round(r["result"]["cdf_gap"],9),r["drive_amplitude"],r["mode_index"]))
            ix=best["mode_index"]; amp=best["drive_amplitude"]
            j=AMPLITUDES.index(amp);bests[label]=dict(mode_index=ix,amplitude=amp)
            for grid in ("fine","spectral","angle"):
                details.append((label,"0.01",ix,amp,True,grid))
            # Test convergence and amplitude near the preferred tabulated setting.
            native=label.replace("_half","_native")
            if len(waves[native]["k"])!=len(waves[label]["k"]):
                raise ValueError("Incompatible native/half modes")
            details.extend([(native,"0.01",ix,amp,True,"base"),(label,"0.01",ix,amp,False,"fine"),
                            (label,"0.01",ix,0.,True,"fine")])
            for neighbor in AMPLITUDES[max(0,j-1):min(len(AMPLITUDES),j+2)]:
                details.append((label,"0.01",ix,neighbor,True,"fine"))
            for marker in ("0.001","0.05"):
                details.append((label,marker,ix,amp,True,"fine"))
        batch(details)
    output.write_text(json.dumps(dict(grids=GRIDS,amplitude_grid=AMPLITUDES,
        analytic_uniform_drive_max_abs_error=analytic,best_base=bests,runs=rows),
        indent=2,sort_keys=True,allow_nan=False)+'\n')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prior",type=Path,required=True)
    p.add_argument("--transfer",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--jobs",type=int,default=2)
    args=p.parse_args();analyze(args.prior.resolve(),args.transfer,args.output,args.jobs)
if __name__=="__main__":
    main()
