"""Reproduce the homogeneous, relativistic two-stream scattering experiment."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from scipy.optimize import brentq
from scipy.special import kve

HERE = Path(__file__).resolve().parent
MASS, WIDTH = 8300., .30
SEEDS = (1139, 2269, 3517, 4717, 5849, 6971)
TIMES = (0., .125, .25, .5, 1., 2., 4., 8., 16., 32., 64., 128.)
FIELDS = ("time", "charge_odd_momentum", "pressure_quadrupole", "number_legendre4",
          "ks_energy_pooled", "ks_energy_worst_sign", "mean_squared_energy_over_mean_squared",
          "never_scattered", "collisions_per_particle", "event_rate_over_initial_mono",
          "momentum_loss_over_initial_mono", "total_energy_residual", "total_momentum_residual",
          "max_pair_energy_residual", "c2_kernel_residual")


def mj_energy(theta):
    return float(3 * theta + kve(1, 1/theta) / kve(2, 1/theta))


def matched_temperature(energy):
    return brentq(lambda t: mj_energy(t)-energy, (energy-1)/4, energy)


def kernel(ea, eb):
    """Head-on monoenergetic wall-frame m² sigma v_M, for the initial beam pairs."""
    pa, pb = math.sqrt(ea**2-1), math.sqrt(eb**2-1)
    s = 2 + 2*(ea*eb + pa*pb)
    beta = math.sqrt(1-4/s)
    g = WIDTH * (s+2)
    return 6*math.pi*beta*g*g/((s-MASS**2)**2+(g*beta)**2)/(ea*eb)


def config(path, energy, times):
    theta = matched_temperature(energy)
    def pdf(x):
        return (1+theta*x) * np.sqrt(x*(2+theta*x)) * np.exp(-x)
    grid = np.linspace(0., 60., 24001)
    cdf = cumulative_trapezoid(pdf(grid), grid, initial=0.)
    cdf /= cdf[-1]
    norm = quad(pdf, 0, np.inf)[0]
    check = max(abs(cdf[j] - quad(pdf, 0, grid[j])[0]/norm)
                for j in np.linspace(1, len(grid)-1, 50, dtype=int))
    if check > 1e-5:
        raise ValueError("Maxwell–Jüttner energy CDF under-resolved")
    second = quad(lambda x: pdf(x)*((1+theta*x)/energy)**2, 0, np.inf)[0]/norm
    with path.open("w") as f:
        f.write(f"{theta:.17g} 60 {len(cdf)} {len(times)}\n")
        np.savetxt(f, cdf, fmt="%.17g")
        np.savetxt(f, times, fmt="%.17g")
    return dict(matched_theta=theta, equilibrium_energy_second_moment=second,
                sampled_cdf_interpolation_error=check, initial_mono_m2_sigma_vm=kernel(energy, energy))


def initial_controls(energy, shape):
    if shape == "thermal":
        return dict(exact_initial_event_rate=None, exact_initial_momentum_loss=0.)
    low, high = 1+.15*(energy-1), 1+1.85*(energy-1)
    energies = [energy] if shape == "mono" else [low, high]
    momenta = np.sqrt(np.square(energies)-1.)
    kernels = np.array([[kernel(a,b)/kernel(energy,energy) for a in energies] for b in energies])
    return dict(exact_initial_event_rate=float(kernels.mean()),
                exact_initial_momentum_loss=float((kernels*(momenta[:, None]+momenta[None, :])).mean()/(2*momenta.mean())))


def wall_column(prior, branch, energy):
    """Bound the hot reference optical-time coordinate from doubling to turnaround."""
    radius = branch["doubling_radius"]
    last_radius = branch["maximum_radius"]
    delta = branch["maximum_tau"]-branch["doubling_time"]
    coefficient = prior["action"]*kernel(energy,energy)/(6*prior["daughter_mass_L"]**3)
    frequency = coefficient*branch["doubling_rest"]/radius**3
    return dict(mono_reference_time_over_L_at_doubling=1/frequency,
                reference_times_per_expansion_time_at_doubling=frequency*radius/branch["doubling_aH"],
                interval_from_doubling_to_turnaround_over_L=delta,
                lower_reference_column=coefficient*branch["doubling_rest"]/last_radius**3*delta,
                upper_reference_column=coefficient*branch["end_rest"]/radius**3*delta)


def run(args):
    prior = json.loads((args.prior / "boost/EVOLUTION.json").read_text())
    branch = prior["runs"]["pole_cap_half"]
    energies = dict(drift_scale=branch["max_gamma"],
                    hot_energy=mj_energy(branch["doubling_theta"]))
    # time is dimensionless lab column density n_target,lab t K0/m²
    specs = {}
    for anchor in energies:
        for shape in ("mono", "split", "thermal"):
            specs[f"{anchor}_{shape}"] = (anchor, shape, 4096, 0., SEEDS)
    specs["drift_scale_mono_double_n"] = ("drift_scale","mono",8192,0.,SEEDS[:4])
    specs["hot_energy_mono_double_n"] = ("hot_energy","mono",8192,0.,SEEDS[:4])
    for step in (.1, .05, .025):
        specs[f"hot_energy_mono_batched_{step}"] = ("hot_energy","mono",4096,step,SEEDS[:4])
    if args.quick:
        specs = {key: (a, s, 512, dt, SEEDS[:2]) for key,(a,s,n,dt,seeds) in specs.items()
                 if key in ("drift_scale_mono", "hot_energy_mono", "hot_energy_split", "hot_energy_mono_batched_0.05")}
    times = TIMES[:8] if args.quick else TIMES
    with tempfile.TemporaryDirectory(prefix="hdblast-box-") as temp:
        working = Path(temp)
        ref = {key: config(working/key, e, times) for key,e in energies.items()}
        exe = working/"box"
        subprocess.run(["g++", "-std=c++17", "-O3", "-Wall", "-Wextra", "-pedantic",
                        str(HERE/"particles.cpp"), "-o", str(exe)], check=True)
        def simulate(task):
            key, seed = task
            anchor, shape, n, step, _ = specs[key]
            cmd = [str(exe), str(n), repr(energies[anchor]), shape, str(seed), str(step), str(working/anchor)]
            raw = subprocess.run(cmd, capture_output=True, text=True)
            if raw.returncode:
                raise RuntimeError(f"{key} {seed}: {raw.stderr}")
            data = np.fromstring(raw.stdout, sep=" ").reshape(len(times), len(FIELDS))
            if not np.isfinite(data).all() or not np.array_equal(data[:, 0], times):
                raise ValueError("Invalid particle-box result")
            err = data[:, 11:14].max()
            if err > 1e-8:
                raise ValueError(f"{key} conservation residual {err}")
            return key, data
        tasks = [(name, seed) for name,(_,_,_,_,seeds) in specs.items() for seed in seeds]
        groups = {name: [] for name in specs}
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for name, data in pool.map(simulate, tasks):
                groups[name].append(data)
        result = {}
        for name, runs in groups.items():
            anchor, shape, n, step, seeds = specs[name]
            arr = np.array(runs)
            columns = list(FIELDS[1:11])
            result[name] = dict(energy_anchor=anchor, shape=shape, particles_per_sign=n,
                                batch_dt=step, seeds=seeds, times=times,
                                **initial_controls(energies[anchor],shape),
                                mean={k:arr[:,:,j+1].mean(axis=0).tolist() for j,k in enumerate(columns)},
                                standard_error={k:(arr[:,:,j+1].std(axis=0,ddof=1)/math.sqrt(len(runs))).tolist()
                                                for j,k in enumerate(columns)},
                                maximum_residual={k:float(arr[:,:,j].max()) for j,k in enumerate(FIELDS) if 11<=j<=13},
                                final_mean_angular_c2_residual=float(arr[:,-1,14].mean()),
                                final_individual_energy_ks=arr[:,-1,5].tolist())
    return dict(quick=args.quick, source_run="pole_cap_half", initial_energies=energies, normalization=ref,
                late_wall_reference_column=wall_column(prior, branch, energies["hot_energy"]),
                pole_mass_over_m=MASS, fractional_width=WIDTH, runs=result)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--prior", type=Path, required=True)
    cli.add_argument("--output", type=Path, required=True)
    cli.add_argument("--jobs", type=int, default=2)
    cli.add_argument("--quick", action="store_true")
    args = cli.parse_args()
    args.output.write_text(json.dumps(run(args), indent=2, sort_keys=True, allow_nan=False)+"\n")
