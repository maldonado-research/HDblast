"""Reproduce the finite-collar coherent decay and charged-gas heating runs (Python 3.9+ and C++17)."""
import argparse
import hashlib
import json
import math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
COLUMNS = "t a H psi pi rad sheet e V source_time emitted returned source_radius K B X mass p_total comoving_wave_to_EM comoving_scalar_to_EM chi chi_dot v M rho_chi p_chi drive_work mass_work chi_pressure_work dU_dpsi_drive dU_dpsi_mass psi_em_work psi_pressure_work Gamma decay_power D S U gas pgas intrinsic theta velocity current decay_energy pairs birth_kinetic electric_work electric_work_a4 gas_pressure_work drag_integral mixing_integral field_exposure formation_phase optical_depth birth_chi2_over_K birth_amplitude2_over_K".split()
def waveform_constants():
    stop = 3 + math.log(39) / 4 + .0625
    def logistic(z):
        val = math.exp(-abs(z))
        return 1/(1+val) if z >= 0 else val/(1+val)
    def step(v):
        if v <= 0: return 0.
        if v >= 1: return 1.
        return logistic(1/(1-v)-1/v)
    def square(t):
        z=4*(t-3)
        return logistic(z)*logistic(-z)**2*step(t/.35)**2*step((stop-t)/.125)**2
    n=32768
    total=math.fsum((2 if j%2==0 else 4)*square(j*stop/n) for j in range(1,n))
    return math.sqrt(2*stop*total/(3*n)),stop


def run(executable, p):
    norm, stop = waveform_constants()
    args = [p["early_cells"], p["depth"], p["dose"], p["end"], norm, stop,
            p["cfl"], p["late_cells"], p["remesh_time"], p["final_gauge_increment"], p["sample"],
            p["decay_threshold"],
            p["decay_width"], p["initial_mass_ratio"], p["final_mass_ratio"],
            p["exit_field"], p["exit_width"], p["preparation_field"], p["preparation_width"],
            p["preparation_gauge_increment"], p["preparation_energy"], p["radius_limit"], p["daughter_mass"], p["drive_amplitude"],
            p["drive_start"], p["drive_width"], p["kappa"], p["charge_to_mass"], p["drag"]]
    proc = subprocess.run([str(executable), *map(str, args)], capture_output=True, text=True, check=True)
    result = {}
    trace = []
    for line in proc.stdout.splitlines():
        key, *values = line.split()
        values = [float(item) for item in values]
        if any(not math.isfinite(x) for x in values):
            raise ValueError(f"Nonfinite value in {key}")
        if key == "@":
            if len(values) != len(COLUMNS):
                raise ValueError("Inconsistent trace length")
            trace.append(values)
        else:
            result[key] = values[0]
    if trace:
        result["trace"] = trace
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cases", type=Path, default=HERE/"cases.json")
    parser.add_argument("--only", nargs="+", help="Run only the named cases")
    parser.add_argument("--jobs", type=int, default=2)
    args = parser.parse_args()
    config = json.loads(args.cases.read_text())
    names = args.only or list(config["cases"])
    parameters = {key: config["defaults"] | config["cases"][key] for key in names}
    source = HERE/"joined_forward.cpp"
    output = {"solver_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
              "trace_columns": COLUMNS, "parameters": parameters}
    with tempfile.TemporaryDirectory(prefix="hdblast-growth-") as temp:
        executable = Path(temp)/"growth"
        subprocess.run(["g++", "-std=c++17", "-O3", str(source), "-o", str(executable)], check=True)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            output["runs"] = dict(zip(names, pool.map(lambda n: run(executable, parameters[n]), names)))
    args.output.write_text(json.dumps(output, separators=(",", ":"), sort_keys=True, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
