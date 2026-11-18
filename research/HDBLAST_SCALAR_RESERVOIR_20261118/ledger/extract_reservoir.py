"""Reconstruct the B1 work rates. Requires NumPy, SciPy and the original B1 snapshots.

Usage: python extract_reservoir.py B1_derived_chi_source OUTPUT_DIRECTORY
The public CSVs contain only the eight weighted observables on the tested interval.
"""
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import roots_genlaguerre

PREFIX = "main_ps0.5_G100_y1_bcut_"
PINS = {
    "B1_RESULTS.json": "2a75112eaf2ea70822e8e889ba676f205163d9d667d1c29f47f9abae67719873",
    "evolve_b1.py": "a0ae79e12ff90f549a0d1ca674fa3285b8ad456b780880df405b9597d3434f0f",
    "shell_matter.py": "58c4c4cc18132f86aaa45ac9f0930490ef3debe160fe8f6ed59274451aadc7d1",
    "runs/main/" + PREFIX + "dzf1e-3_summary.json": "2f2b61e745c438ef057533d75329f998621983a3e565226c9d725e4925e34716",
    "runs/main/" + PREFIX + "dzf1e-3_timeseries.npz": "3c9e732322715bc33d09eff6212066b15448da70c7ea10a177edc39e9f6ac9b2",
    "runs/main/" + PREFIX + "dzf5e-4_summary.json": "bf162614a84fd207f2f2736b8b77cd1cb173c2db012469cf1d4fe1ccc50c8d0f",
    "runs/main/" + PREFIX + "dzf5e-4_timeseries.npz": "9eae8828655af2ec8cc91bb658ecf0e768ca9f54215b8927313937064f1b1646",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def work_on_filling(tr, mass2, event):
    """Reconstruct work per unit volume and s from the single cohort's saved progress."""
    nodes, weights = roots_genlaguerre(24, 0.5)
    q = event["q"]
    momentum = np.sqrt(q * nodes / np.pi)
    number = weights * (q / np.pi)**1.5 * event["number_factor"] / (4 * np.pi**2)
    rel_a = np.exp(tr["A_b"] - event["A_star"])
    energy = np.sqrt(mass2[:, None] + (momentum / rel_a[:, None])**2)
    # A stable, explicit reduction fixes the order of summation across the 24 nodes.
    target_energy = np.array([math.fsum(row*number) for row in energy]) / rel_a**3
    direction = np.sign(event["v_star"])
    v = tr["v_over_H0"]
    speed = np.clip(direction * v / (0.05 * abs(event["v_star"])), 0, 1)
    switch = speed**2 * (3 - 2*speed)
    u = tr["u_last"]
    ramp_derivative = np.where((u > 0) & (u < 1), 30*u**2*(1-u)**2, 0.0)
    return v * direction * switch * ramp_derivative * target_energy / event["ramp_dphi"]


def export(files, label, spacing, directory):
    tag = PREFIX + spacing
    path = "runs/main/" + tag
    summary = json.loads(files[path + "_summary.json"])
    report = json.loads(files["B1_RESULTS.json"])["runs"][tag]
    meta = summary["params"]["matter"]
    b, rb = meta["b"], summary["rho_b"]
    require((meta["G"], meta["phistar"], meta["y"], b, meta["mu0"], meta["nk"]) ==
            (100, 0.5, 1, 6.49e-6, 0, 24), "unsupported cell or particle physics")
    require(summary["params"]["delta"] == 0.1 and meta["mode"] == "ramp"
            and meta["sign_j"] == 1 and summary["params"]["Y"] == 0, "source branch")
    require(len(summary["matter"]["events"]) == 1, "expected a single cohort")
    event = summary["matter"]["events"][0]
    with np.load(io.BytesIO(files[path + "_timeseries.npz"]), allow_pickle=False) as npz:
        tr = dict(zip(npz["cols"].tolist(), npz["rec"].T))
    s, h, v = (tr[k] for k in ("H0tau", "H_over_H0", "v_over_H0"))
    require(s[0] == 0 and np.all(np.diff(s) > 0) and report["H0tau_end"] > 4
            and s[-1] >= 4, "proper clock coverage or reliability limit")
    a3 = np.exp(3*(tr["A_b"]-tr["A_b"][0]))
    dphi = tr["phi_b"]-meta["phistar"]
    mass2 = meta["G"]**2*dphi**2
    # The trace of the on-shell particle distribution reconstructs the mass force.
    adiabatic_force = np.divide(tr["rho_hat"]-3*tr["p_hat"], dphi,
                                out=np.zeros_like(s), where=dphi != 0)
    ramp_power = work_on_filling(tr, mass2, event)
    decay = meta["y"]**2 * mass2 * tr["n_hat"] / (8*np.pi)
    for computed, native, name in [
        (adiabatic_force, tr["j_hat"], "mass force"),
        (ramp_power, tr["P_hat"], "cohort production power"),
        (ramp_power, tr["jprod_hat"]*v, "production force times velocity"),
        (decay, tr["Q_hat"], "internal decay"),
        (rb*tr["X"], b*tr["rho_hat"], "particle energy"),
        (rb*tr["P"], b*tr["p_hat"]+rb*tr["R"]/3, "total pressure"),
        (rb*tr["J"]*v, b*(adiabatic_force*v+ramp_power), "junction scalar work")]:
        require(np.allclose(computed, native, rtol=2e-13, atol=5e-13), name)
    chi = a3*b*tr["rho_hat"]
    radiation = a3*rb*tr["R"]
    series = {
        "chi_energy": chi, "radiation_energy": radiation,
        "mass_work_rate": a3*b*adiabatic_force*v, "filling_work_rate": a3*b*ramp_power,
        "decay_rate": a3*b*decay, "chi_pressure_rate": 3*h*a3*b*tr["p_hat"],
        "radiation_pressure_rate": h*radiation,
    }
    points = np.r_[1.0, s[(s > 1) & (s < 4)], 4.0]
    cols = [points] + [np.interp(points,s,col) for col in series.values()]
    csv_path = directory / ("reservoir_" + label + ".csv")
    with csv_path.open("w", newline="") as out:
        writer = csv.writer(out, lineterminator="\n")
        writer.writerow(["s", *series])
        writer.writerows([format(float(x),".17g") for x in row] for row in zip(*cols))
    return {
        "timeseries": path+"_timeseries.npz", "summary": path+"_summary.json",
        "csv_sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest(),
        "samples": len(points), "rho_b_inverse_H0": rb,
        "reliable_end_s": report["H0tau_end"], "registered_plateau": report["plateau"],
        "maximum_reported_radiation_fraction": report["max_Omega_r"], "crossing_s": event["s_star"],
        "crossing_dimensionless_density_jump": event["kappa5sq_rho_jump_over_H0"],
        "approximate_jump_times_interpolated_crossing_a3":
            event["kappa5sq_rho_jump_over_H0"] * math.exp(3*(event["A_star"]-tr["A_b"][0])),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("b1_directory",type=Path)
    parser.add_argument("output_directory",type=Path)
    args = parser.parse_args()
    files = {name: (args.b1_directory/name).read_bytes() for name in PINS}
    for name, raw in files.items():
        require(hashlib.sha256(raw).hexdigest() == PINS[name], name)
    args.output_directory.mkdir(parents=True,exist_ok=True)
    rows = {label: export(files,label,spacing,args.output_directory)
            for label, spacing in [("coarse","dzf1e-3"), ("fine","dzf5e-4")]}
    result = {"source_sha256": PINS, "public_solver_revision": "a011dee2b490faa56ae81d3165b4b5fa5bbfe636",
              "archived_input_revision": "140a0e9c672b03e930400d5109213f662842bdf4",
              "cutoff": "b_cut, lambda_c=1", "b": 6.49e-6, "phi_star": 0.5, "G": 100, "y": 1,
              "bare_mass": 0, "detuning": 0.1, "clock_interval": [1,4], "runs": rows}
    (args.output_directory / "EXTRACTION.json").write_text(json.dumps(result,indent=2)+"\n")


if __name__ == "__main__":
    main()
