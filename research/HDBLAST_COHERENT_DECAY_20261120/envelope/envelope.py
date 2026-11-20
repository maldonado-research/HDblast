"""Compare a cycle-averaged parent with the resolved run, then scan finite daughter masses."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import math
from pathlib import Path
import runpy
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
STUDY = HERE.parent
COUPLED = STUDY / "coupled"


def screening(p, r):
    """Necessary coupling and thermal-rate conditions at one action normalization."""
    mu, charge, kappa = (p[k] for k in ("daughter_mass", "charge_to_mass", "kappa"))
    y_squared = 4 * math.pi * kappa / (3 * mu)
    # The base K is bounded below by .001/40. Both added plateaux are positive.
    gauge_K_floor = .001 / 40.
    alpha_limit = .1
    mass_fraction_limit = .005
    gauge_floor = 3 * (charge*mu)**2 / (4 * math.pi * gauge_K_floor * alpha_limit)
    mass_ceiling = (mass_fraction_limit*mu/r["whole_chi_bound"])**2 / y_squared
    drag_floor = r.get("low_drift_drag_floor")
    result = {
        "alpha_gauge_limit": alpha_limit, "relative_mass_shift_limit": mass_fraction_limit,
        "gauge_K_floor": gauge_K_floor, "gauge_action_floor": gauge_floor,
        "yukawa_alpha_each": y_squared/(4*math.pi), "mass_shift_action_ceiling": mass_ceiling,
        "max_electric_force_over_m2": charge*r["whole_abs_field_max"]/mu,
        "gauge_to_mass_window_ratio": mass_ceiling/gauge_floor,
        "birth_rms_shift_at_local_unit_gauge_alpha": charge*math.sqrt(
            3*kappa/(3*mu)*r["end_birth_chi2_over_K"]/r["end_rest"]),
    }
    if drag_floor is not None:
        lower = max(gauge_floor, drag_floor)
        result["elastic_drag_action_floor"] = drag_floor
        result["three_condition_window_ratio"] = mass_ceiling/lower
        if lower <= mass_ceiling:
            witness = math.sqrt(lower*mass_ceiling)
            result["witness_action"] = witness
            result["witness_max_gauge_alpha_bound"] = alpha_limit*gauge_floor/witness
            result["witness_max_mass_shift_bound"] = mass_fraction_limit*math.sqrt(witness/mass_ceiling)
            result["witness_drag_capacity_over_needed"] = witness/drag_floor
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",required=True,type=Path)
    parser.add_argument("--jobs",type=int,default=2)
    parser.add_argument("--joined",type=Path,default=COUPLED/"JOINED.json",
                        help="Resolved calculation used for comparison and for the incoming velocity.")
    args=parser.parse_args()
    solver=runpy.run_path(str(COUPLED/"joined.py"))
    transform=runpy.run_path(str(HERE/"build.py"))["transform"]
    native_src=(COUPLED/"joined_forward.cpp").read_text()
    src=transform(native_src)
    centered_src=transform(native_src, centered=True)
    base=json.loads((COUPLED/"cases.json").read_text())["defaults"]
    # mu*v0 and mu*drive width fixed; geometry is allowed to respond in every case.
    configs={"resolved_comparison":(400,2e8,24,6,.08),
             "reduced_charge":(400,1e7,24,6,.08),
             "mass_4e3":(4e3,1e7,24,6,.08), "mass_4e4":(4e4,1e7,24,6,.08),
             "mass_4e6":(4e6,1e7,24,6,.08), "mass_4e8":(4e8,1e7,24,6,.08),
             "mass_4e10":(4e10,1e7,24,6,.08),
             "mass_4e10_centered":(4e10,1e7,24,6,.08),
             "mass_4e12":(4e12,1e7,24,6,.08),
             "mass_4e12_fine":(4e12,1e7,48,8,.08),
             "mass_4e12_half_step":(4e12,1e7,24,6,.04),
             "mass_4e10_fine":(4e10,1e7,48,8,.08), "mass_4e10_half_step":(4e10,1e7,24,6,.04),
             "mass_4e8_fine":(4e8,1e7,48,6,.08), "mass_4e8_half_step":(4e8,1e7,24,6,.04),
             "mass_4e8_lower_charge":(4e8,3e6,24,6,.08)}
    params={}
    for name,(mu,charge,early,late,cfl) in configs.items():
        params[name]=base|dict(daughter_mass=mu,charge_to_mass=charge,drive_amplitude=.008/mu,
            drive_width=.24/mu,early_cells=early,late_cells=late,cfl=cfl,final_gauge_increment=2.5)
    with tempfile.TemporaryDirectory(prefix="hdblast-averaged-") as folder:
        program=Path(folder)/"avg.cpp"
        exe=Path(folder)/"avg"
        program.write_text(src)
        subprocess.run(["g++","-std=c++17","-O3",str(program),"-o",str(exe)],check=True)
        centered_exe=Path(folder)/"centered"
        program.write_text(centered_src)
        subprocess.run(["g++","-std=c++17","-O3",str(program),"-o",str(centered_exe)],check=True)
        def evolve(name):
            use_centered=name.endswith("_centered") or params[name]["daughter_mass"]>=1e11
            return solver["run"](centered_exe if use_centered else exe,params[name])
        with ThreadPoolExecutor(args.jobs) as pool:
            runs=dict(zip(params,pool.map(evolve,params)))
    checks={}
    for name,r in runs.items():
        required=["handover_action","birth_min_K","whole_chi_bound","first_hot_time","stop_turnaround"]
        if not all(r.get(k,0)>0 for k in required): raise ValueError("Missing model milestone: "+name)
        if max(r["gas_balance_error"],r["oscillator_balance_error"],r["surface_balance_error"])>1e-5:
            raise ValueError("Energy ledger failed: "+name)
        checks[name]=screening(params[name],r)
    joined=json.loads(args.joined.read_text())
    inner=runpy.run_path(str(HERE/"inner.py"))["audit"](joined,runs)
    native=joined["runs"]["coarse"]
    kinetics=runs["resolved_comparison"]
    comparison={k:kinetics[k]/native[k]-1 for k in [
        "opening_time","hot_efolds","first_hot_radius","end_rest","end_decay_energy",
        "end_electric_work","end_birth_chi2_over_K","end_birth_amplitude2_over_K","maximum_radius"]}
    result=dict(transformed_solver_sha256=hashlib.sha256(src.encode()).hexdigest(),
                centered_solver_sha256=hashlib.sha256(centered_src.encode()).hexdigest(),
                parameters=params,runs=runs,screen=checks,inner_source=inner,averaged_minus_resolved_relative=comparison)
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+"\n")

if __name__=="__main__":
    main()
