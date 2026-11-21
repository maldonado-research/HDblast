"""Extract the inward wall characteristics and reconstruct the three fields of the parent solution."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import shlex
import subprocess
import tempfile

import numpy as np
from scipy.integrate import cumulative_trapezoid

CHANNELS = ((2, "early"), (2, "late_birth"))
BASE = "t a K X B incoming D phi tube_time emitted".split()
EXTRA = "H theta u R tried_rest e intrinsic c d target current net_native pump_native work_native".split()
WAVES = "E B I net pump emitted square".split()
END = 6.6805

# These expressions are the directional birth weights exported by #198's pinned ancestor.
BIRTH_RECORD = """const auto ch = oscillator(y);
            const double R = ch.mass_ratio;
            const double z = R>1. ? std::sqrt((R-1.)*(R+1.)) : 0.;
            const double c = z>0. ? h*z+R*ch.mass_slope*y[PI]/(2.*daughter_mass*z) : 0.;
            const double d = charge_mass*wnew.x;
            const double target = envelope_active && ch.power>0. && z>0. ?
                std::pow(std::acos(-1.),2)/6.*feedback_action/std::pow(daughter_mass,4)*ch.power/(R*z*z) : 0.;
            const double extra[] = {h,gas.theta,gas.u,R,y[TRIED_REST],wnew.e,gas.intrinsic_random,
                c,d,target,gas.current,y[WAVE_TO_SHEET],y[SCALAR_TO_SHEET],y[ELECTRIC_WORK_A4]};
            trace.write(reinterpret_cast<const char*>(extra),sizeof(extra));"""


def build(study, folder):
    launcher = runpy.run_path(str(study / "fields/launch.py"))
    prior, exe = launcher["build"](study / "parent", folder, CHANNELS)
    source = folder / "probe.cpp"
    code = source.read_text()
    anchor = "trace.write(reinterpret_cast<const char*>(record),sizeof(record));"
    if code.count(anchor) != 1:
        raise ValueError("Parent radial recorder differs")
    source.write_text(code.replace(anchor, anchor+"\n            "+BIRTH_RECORD))
    subprocess.run(["g++", "-std=c++17", "-O3", "-ffp-contract=off", "-I", str(prior/"boost"), str(source), "-o", str(exe)], check=True)
    return prior, exe


def reconstruct(t, a, k, x, b, inward, current):
    """Integrate conformal fields in proper time, treating each recorded driver independently."""
    electric = a*a*x
    magnetic = a*a*b
    f, g = electric[0], magnetic[0]
    peak_f = peak_g = 0.
    dt = np.diff(t)
    mid = lambda v: .5*(v[1:]+v[:-1])
    radius = mid(a)
    kin = mid(k)
    incoming = mid(inward)
    source = (4.*radius**2*incoming-mid(current)/radius)/kin
    # Implicit midpoint for F_dot=-2F/K-2G/a+source and G_dot=2F/a.
    for step, curl, damp, drive, ref_f, ref_g in zip(dt, dt/radius, dt/kin, source, electric[1:], magnetic[1:]):
        mid_f = (f-curl*g+.5*step*drive)/(1.+damp+curl**2)
        f, g = 2.*mid_f-f, g+2.*curl*mid_f
        peak_f = max(peak_f, abs(f-ref_f))
        peak_g = max(peak_g, abs(g-ref_g))
    scale = max(np.hypot(electric,magnetic))
    stored = .5*k*(electric**2+magnetic**2)
    incident = float(np.trapezoid(2.*a**4*inward**2, t))
    returned = float(np.trapezoid(2.*a**4*(x-inward)**2, t))
    work = float(np.trapezoid(a*x*current, t))
    ledger = stored[-1]-stored[0]-incident+returned+work
    report = dict(max_F_error_over_peak=float(peak_f/scale), max_G_error_over_peak=float(peak_g/scale),
                  stored_initial=float(stored[0]), stored_final=float(stored[-1]),
                  inward_energy=incident, outward_energy=returned, work_to_current=work,
                  ledger_residual_over_peak=float(abs(ledger)/max(stored)))
    return report


def energy_windows(h):
    """Signed exchange for each interval bounded by a fluid-current or electric-field reversal."""
    t, a, x, j, k, inc = (h[key] for key in ("t", "a", "X", "J", "K", "incoming"))
    stored = .5*k*a**4*(x*x+h["B"]**2)
    clocks = [(float(t[0]), "birth")]
    for field, name in ((x, "field_zero"), (j, "current_zero")):
        ix = np.flatnonzero(field[:-1]*field[1:] < 0)
        for i in ix:
            frac = -field[i]/(field[i+1]-field[i])
            clocks.append((float(t[i]+frac*(t[i+1]-t[i])), name))
    clocks.append((float(t[-1]), "end"))
    clocks.sort()
    powers = {"work_to_current": a*x*j, "inward_energy": 2.*a**4*inc**2,
              "outward_energy": 2.*a**4*(x-inc)**2}
    rows = []
    for (left, before), (right, after) in zip(clocks, clocks[1:]):
        keep = (t>left)&(t<right)
        knots = np.r_[left, t[keep], right]
        energy = {key:float(np.trapezoid(np.interp(knots,t,value),knots)) for key,value in powers.items()}
        change = float(np.interp(right,t,stored)-np.interp(left,t,stored))
        residual = change-energy["inward_energy"]+energy["outward_energy"]+energy["work_to_current"]
        rows.append(dict(start=left,end=right,first_event=before,last_event=after,stored_change=change,
                         residual_over_global_peak=float(abs(residual)/max(stored)),**energy))
    work = [row["work_to_current"] for row in rows]
    return dict(windows=rows,gross_work_to_particles=sum(max(w,0.) for w in work),
                gross_work_returned_by_particles=-sum(min(w,0.) for w in work))


def generate(label, study, prior, exe, metadata, output_dir, tmp, end):
    cfg = json.loads((prior/"prior/envelope/ENVELOPE.json").read_bytes())["parameters"]["mass_4e12"]
    spec = metadata["runs"][label]
    trace = tmp/(label+".bin")
    command = tmp/label
    tail = [2,metadata["action"],prior/"boost/BOOST_TABLE.dat",1,spec["channel_fraction"],trace]
    command.write_text("#!/bin/sh\nexec "+shlex.quote(str(exe))+' "$@" '+shlex.join(list(map(str,tail)))+"\n")
    command.chmod(0o755)
    diag = runpy.run_path(str(prior/"prior/coupled/joined.py"))["run"](
        command,dict(cfg,cfl=cfg["cfl"]*0.5,sample=0,end=end))
    nbase = len(BASE)+len(EXTRA)
    raw = np.fromfile(trace, dtype="=f8").reshape(-1,nbase+3*len(WAVES))
    if abs(diag["end_tau"]-end)>2e-12 or not np.isfinite(raw).all():
        raise ValueError("Invalid radial trace")
    first = np.searchsorted(raw[:,6],0.,side="right")-1
    if first < 0 or raw[first,0] < 6.18 or raw[first,6] != 0:
        raise ValueError("Missing recorded birth onset")
    all_waves = raw[:,nbase:].reshape(-1,3,len(WAVES))
    if np.any(all_waves[:,2,:] != 0):
        raise ValueError("Unforced radial characteristic changed")
    scale = np.sqrt(raw[-1,BASE.index("emitted")]/all_waves[-1,:2,WAVES.index("emitted")])
    rows = raw[first:]
    h = {k:rows[:,i] for i,k in enumerate(BASE+EXTRA)}
    t, a, d, k = (h[key] for key in ("t","a","D","K"))
    if np.any(np.diff(t)<=0) or np.any(np.diff(d)<0) or np.any(k!=2.54):
        raise ValueError(("Invalid gauge plateau or cumulative birth history",label,float(k.min()),float(k.max()),float(t[np.argmin(k)])))
    current = a**3*h["current"]
    q = metadata["charge_to_mass"]
    fluid = q*d*h["u"]/np.hypot(1.,h["u"])
    current_error = float(max(abs(current-fluid))/max(abs(current)))
    if current_error>2e-12:
        raise ValueError("Inherited fluid current normalization mismatch")
    h["J"] = current
    h["I"] = q*cumulative_trapezoid(a*h["X"],t,initial=0.)
    h["work_fraction"] = cumulative_trapezoid(abs(current*h["X"]),t,initial=0.)/spec["absolute_work"]
    for index in range(2):
        for var in ("E","B","I"):
            h[var+str(index+1)] = scale[index]*all_waves[first:,index,WAVES.index(var)]
    archive = json.loads((study/"fields/FIELDS.json").read_bytes())[label]["channels"]
    comparisons={}
    with np.load(study/"parent/parent/history"/(label+".npz")) as history:
        days = history["t"]
        for key in history.files:
            delta = np.max(np.abs(np.interp(days,t,h[key])-history[key]))
            # #198 independently assembled B via the full-history aX quadrature.
            tolerance = 1e-8 + 2e-6*max(abs(history[key]))
            if delta>tolerance:
                raise ValueError(("Archived particle input differs",label,key,delta,tolerance))
            comparisons[key]=float(delta/max(1.,max(abs(history[key]))))
    for axis, (_,profile) in enumerate(CHANNELS,1):
        ref = archive["2_"+profile]
        if not np.isclose(scale[axis-1],ref["flux_normalizing_multiplier"],rtol=1e-11):
            raise ValueError("Source normalization changed")
        for key, marker in spec["markers"].items():
            for var in ("E","B"):
                value = np.interp(marker["time"],t,h[var+str(axis)])
                if not np.isclose(value,ref["markers"][key][var],rtol=2e-6,atol=1e-8):
                    raise ValueError(("Early fields differ",label,axis,key,var))
    for key, value in (("action",metadata["action"]),("charge_to_mass",q),
                       ("daughter_mass_L",metadata["daughter_mass_L"]),("channel_fraction",spec["channel_fraction"])):
        h[key] = value
    np.savez_compressed(output_dir/(label+".npz"),**h)
    zeros = np.zeros_like(current)
    replay = {"axial":reconstruct(t,a,k,h["X"],h["B"],h["incoming"],current)}
    for axis in (1,2):
        replay["polarization_"+str(axis)] = reconstruct(t,a,k,h["E"+str(axis)],h["B"+str(axis)],h["I"+str(axis)],zeros)
    for answer in replay.values():
        if max(answer["max_F_error_over_peak"],answer["max_G_error_over_peak"],answer["ledger_residual_over_peak"])>1e-5:
            raise ValueError(("Maxwell recovery or energy check",label,answer))
    stored = .5*k*a**4*(h["X"]**2+h["B"]**2)
    budget = stored-h["net_native"]-h["pump_native"]+h["work_native"]
    return dict(start=float(t[0]),end=float(t[-1]),rows=len(t),channel_fraction=spec["channel_fraction"],
        gauge=float(k[0]),fraction_of_reference_work=float(h["work_fraction"][-1]),
        relative_fluid_current_error=current_error,max_archive_difference_over_scale=comparisons,
        native_axial_balance_error_over_peak=float(max(abs(budget))/max(stored)),
        launch_normalizations=scale.tolist(),maxwell_recovery=replay,axial_energy_windows=energy_windows(h))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study",required=True,type=Path)
    parser.add_argument("--output-dir",required=True,type=Path)
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--jobs",default=2,type=int)
    parser.add_argument("--end",type=float,default=END)
    parser.add_argument("--cases",nargs="+",choices=("dedicated_half","shared_half"),default=("dedicated_half","shared_half"))
    opt = parser.parse_args()
    study = opt.study.resolve()
    opt.output_dir.mkdir(parents=True,exist_ok=True)
    metadata = json.loads((study/"parent/parent/history/markers.json").read_bytes())
    labels = opt.cases
    with tempfile.TemporaryDirectory(prefix="hd-axial-drivers-") as tmp:
        folder = Path(tmp)
        prior, exe = build(study, folder)
        def one(label):
            return generate(label,study,prior,exe,metadata,opt.output_dir,folder,opt.end)
        with ThreadPoolExecutor(opt.jobs) as pool:
            runs = dict(zip(labels,pool.map(one,labels)))
    opt.output.parent.mkdir(parents=True,exist_ok=True)
    opt.output.write_text(json.dumps(runs,indent=2,sort_keys=True,allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
