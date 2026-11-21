"""Add direct vector-pole collisions to the early homogeneous #198 daughter calculation."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import shlex
import subprocess
import tempfile

import numpy as np
from scipy.integrate import quad


HERE = Path(__file__).resolve().parent


def build(study, work):
    code = (study / "feedback/homogeneous.cpp").read_text()
    changes = [
        ("int main(int argc", '#include "collide.hpp"\nint main(int argc'),
        ("if(argc!=7)", "if(argc!=11)"),
        ("Wave wave;size_t n=0;", """Wave wave;size_t n=0;
    PoleCollisions pole(p.size(),std::stod(argv[7]),std::stoull(argv[8]),std::stoi(argv[9])!=0);"""),
        ("const std::array<double,10> trace=", "pole.step(p,n,weight,l[0],r[0],r[1]);\n        const std::array<double,10> trace="),
        ('if(!raw||!hist)throw std::runtime_error("output write");', 'pole.save(argv[10]);\n    if(!raw||!hist)throw std::runtime_error("output write");\n    return 0;'),
    ]
    for old, new in changes:
        if code.count(old) != 1:
            raise ValueError("Unexpected parent homogeneous source at: " + old)
        code = code.replace(old, new)
    src = work / "evolve.cpp"
    src.write_text(code)
    binary = work / "evolve"
    subprocess.run(["g++", "-std=c++17", "-O3", "-I", str(HERE), str(src), "-o", str(binary)], check=True)
    return binary


def combined_stress(p):
    """Principal pressures of both signs in the rest frame of their combined stress tensor."""
    gamma = np.hypot(1., np.linalg.norm(p, axis=1))
    momentum = np.column_stack((gamma, p))
    stress = (momentum.T / gamma) @ momentum / len(p)
    mixed = stress * np.array([-1., 1., 1., 1.])
    values = np.linalg.eigvals(mixed)
    if np.max(np.abs(values.imag)) > 1e-10 or np.sum(values.real < 0) != 1:
        raise ValueError("Combined stress has no unique timelike energy eigenvalue")
    values = np.sort(values.real)
    rho = -values[0]
    if abs((rho - sum(values[1:])) - np.mean(1/gamma)) > 5e-11 * rho:
        raise ValueError("Stress trace identity failed")
    return (values[1:] / rho).tolist()


def absolute_clock_control(work):
    control = work/"clock_control"
    subprocess.run(["g++", "-std=c++17", "-O3", "-I", str(work), "-I", str(HERE),
                    str(HERE/"control.cpp"), "-o", str(control)], check=True)
    rows = []
    for line in subprocess.check_output([str(control)], text=True).splitlines():
        mode, n, events, hits = map(int, line.split())
        mean = (0.8, 0.4, 0.)[mode]
        prob = -np.expm1(-mean)
        if abs(events/n - mean) > 5*np.sqrt(mean/n) or abs(hits/n - prob) > 5*np.sqrt(prob*(1-prob)/n):
            raise ValueError("Constant-rate resonant-beam control failed")
        rows.append(dict(mode=mode, trials=n, expected_mean=mean, measured_mean=events/n,
                         expected_hit_fraction=float(prob), measured_hit_fraction=hits/n))
    return rows


def thermal_quadrature(integrand, lo, hi, **options):
    """The substitution x=u² smooths the thermal density at zero kinetic energy."""
    return quad(lambda u: 2*u*integrand(u*u), np.sqrt(lo), np.sqrt(hi),
                epsabs=options.get("epsabs", 1e-12), epsrel=2e-13)


def simulate(study, fields, jobs):
    feedback = runpy.run_path(str(study / "feedback/currents.py"))
    # The inherited target routine is loaded into a private namespace for this run.
    feedback["coarse"]["thermal_cells"].__globals__["quad"] = thermal_quadrature
    orbits = feedback["pairs"]
    history_dir = study / "parent/parent/history"
    histories = orbits["read_inputs"](history_dir, fields)
    metadata = json.loads((history_dir / "markers.json").read_bytes())
    density = metadata["action"] / (6 * metadata["daughter_mass_L"]**3)
    label = "dedicated_half"
    amplitude = 0.04/np.sqrt(2)
    sources = {"pressure_fit": orbits["uniform_pair"](orbits["FIT"][label][1]),
               "small_pulse": (2, 2, "early", "late_birth", amplitude, amplitude, np.pi/2)}
    snapshots = {}
    for fraction in ("0.001", "0.01", "0.05"):
        h, f, _ = histories[label]
        snapshots[fraction] = {label: (h, f, metadata["runs"][label]["markers"][fraction])}
    jobspecs = []
    for name in sources:
        for power, steps, seed in ((12, 512, 1), (13, 512, 71), (14, 512, 121), (13, 1024, 71)):
            jobspecs.append((name, power, steps, seed, "off", 0))
            for scatter_seed in (18431, 70163, 91283):
                jobspecs.append((name, power, steps, seed, "vector", scatter_seed))
            if power == 14:
                jobspecs.append((name, power, steps, seed, "isotropic", 18431))
    jobspecs = [(spec, "0.01") for spec in jobspecs]
    # Reuse exactly the same emitted waveforms at the other recorded work markers.
    for mark in ("0.001", "0.05"):
        for name in sources:
            for power, steps, birth_seed in ((13, 512, 71), (14, 1024, 121)):
                jobspecs.append(((name, power, steps, birth_seed, "off", 0), mark))
                for seed in (18431, 70163, 91283):
                    jobspecs.append(((name, power, steps, birth_seed, "vector", seed), mark))
    with tempfile.TemporaryDirectory(prefix="pole-collisions-") as temp:
        work = Path(temp)
        binary = build(study, work)
        controls = absolute_clock_control(work)

        def one(specs):
            (name, power, steps, birth_seed, law, seed), fraction = specs
            with tempfile.TemporaryDirectory(dir=work) as folder:
                folder = Path(folder)
                counts = folder/"counts"
                wrapper = folder/"wrapped"
                tail = [str(density if law != "off" else 0.), str(seed), str(int(law == "isotropic")), str(counts)]
                wrapper.write_text("#!/bin/sh\nexec "+shlex.quote(str(binary))+' "$@" '+shlex.join(tail)+"\n")
                wrapper.chmod(0o755)
                result, momenta = feedback["integrate"](snapshots[fraction], label, name, sources[name], power, steps,
                    orbits["SCRAMBLE"]+birth_seed, True, wrapper, folder)
                proposed, events, hit_fraction, dE, dp, peak, residual, variance = map(float,counts.read_text().split())
            if dE > 1e-10 or dp > 1e-12 or not np.isfinite(result["coarse_cell_total_variation"]):
                raise ValueError("Binary collision failed a conservation check")
            return dict(work_fraction=float(fraction), marker_time=snapshots[fraction][label][2]["time"], source=name, birth_pairs=2**power, steps=steps, birth_seed=birth_seed, law=law, rng_seed=seed,
                candidates=int(proposed), collisions=int(events), fraction_ever_hit=hit_fraction,
                encounters_per_particle=events/2**power, maximum_kernel_over_envelope=peak,
                energy_residual_per_event=dE, momentum_residual_per_event=dp,
                angular_second_moment_residual=residual, angular_second_moment_variance=variance,
                tv20=result["coarse_cell_total_variation"], sideways_pressure_over_target=result["local_lateral_pressure_over_thermal"],
                combined_principal_pressures_over_rest_energy=combined_stress(momenta), ledger=result["ledger"])

        with ThreadPoolExecutor(jobs) as pool:
            rows = list(pool.map(one, jobspecs))
    original = json.loads((study/"feedback/FEEDBACK.json").read_bytes())["runs"]
    for row in rows:
        if row["law"] != "off" or row["birth_pairs"] != 4096:
            continue
        name = "two_pressure_fit" if row["source"] == "pressure_fit" else "equal_early_and_late_4pct"
        prior = next(r for r in original if r["history"] == label and r["case"] == name and
                     r["steps"] == 512 and r["pairs"] == 4096 and r["feedback"])
        if abs(row["tv20"]-prior["coarse_cell_total_variation"]) > 2e-10:
            raise ValueError("Disabled pole does not reproduce the reference counts")
    all_events = [r for r in rows if r["law"] == "vector"]
    angular_z = sum(r["angular_second_moment_residual"] for r in all_events) / np.sqrt(sum(r["angular_second_moment_variance"] for r in all_events))
    return {"runs":[r for r in rows if r["work_fraction"] == 0.01],
            "other_work_markers":[r for r in rows if r["work_fraction"] != 0.01],
            "absolute_clock_controls": controls, "pooled_angular_z_diagnostic":float(angular_z)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", required=True, type=Path)
    parser.add_argument("--fields", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--jobs", default=2, type=int)
    options = parser.parse_args()
    result = simulate(options.study.resolve(), options.fields.resolve(), options.jobs)
    options.output.parent.mkdir(parents=True, exist_ok=True)
    options.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n")
