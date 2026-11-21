"""Give an optimistic upper limit to what collisions can change before the early work marker."""
import argparse
import json
from pathlib import Path
import runpy

import numpy as np
from scipy.optimize import brentq


# Conservative, upward-rounded bound on m^2 sigma v_M for the archived vector pole.
KERNEL_UPPER = 1.118569e-6
TOLERANCES = (0.10, 0.05, 0.01)


def optical_depth(history, births, end, prefactor):
    """Overestimate the exposure in each cell of the monotone recorded history."""
    t, a, rest = (history[key] for key in ("t", "a", "D"))
    last = np.searchsorted(t, end, side="left")
    if not 0 < last < len(t) or min(births) < t[0] or max(births) > end:
        raise ValueError("Marker or sampled birth outside the history")
    grid = np.append(t[:last], end)
    if np.any(np.diff(grid) <= 0) or np.any(np.diff(rest[:last+1]) < 0) or np.any(np.diff(a[:last+1]) < 0):
        raise ValueError("Time, daughter mass and radius must increase")
    rate = prefactor * rest[1:last+1] / a[:last]**3
    # Each birth uses only the occupied fraction of its first interval.
    exposure = np.append(np.cumsum((np.diff(grid) * rate)[::-1])[::-1], 0.0)
    return np.interp(births, grid, exposure)


def assign_cells(p, sign, theta, drift, energy_edges, cosine_edges):
    energy = np.hypot(1, np.linalg.norm(p, axis=1))
    gamma = np.hypot(1, drift)
    along = sign * np.sign(drift) * p[:, 2]
    parallel = gamma * along - abs(drift) * energy
    transverse = np.linalg.norm(p[:, :2], axis=1)
    magnitude = np.hypot(parallel, transverse)
    # The kinetic energy from the Lorentz-transformed three-momentum avoids
    # subtracting nearly equal terms in gamma'-1 for slow particles.
    kinetic = magnitude**2 / ((np.hypot(1, magnitude) + 1) * theta)
    cosine = np.divide(parallel, magnitude, out=np.zeros_like(parallel), where=magnitude > 0)
    rows = np.searchsorted(energy_edges[1:-1], kinetic, side="right")
    cols = np.searchsorted(cosine_edges[1:-1], cosine, side="right")
    return rows * (len(cosine_edges)-1) + cols


def tv_floor(labels, depths, target, multiplier):
    weights = np.exp(-multiplier * depths)
    unhit = np.bincount(labels, weights, minlength=target.size) / labels.size
    excess = np.maximum(unhit - target, 0)
    return float(excess.sum()), unhit


def required_multiplier(labels, depths, target, tolerance):
    def residual(x):
        return tv_floor(labels, depths, target, x)[0] - tolerance
    if residual(0) <= 0:
        return 0.0
    upper = 1.0
    while residual(upper) > 0:
        upper *= 2
        if upper > 1e8:
            raise ValueError("Unresolved required multiplier")
    return brentq(residual, 0, upper, xtol=1e-10)


def calculate(study, ensembles):
    history_dir = study / "parent/parent/history"
    parameters = json.loads((history_dir / "markers.json").read_bytes())
    coarse = runpy.run_path(str(study / "orbits/coarse.py"))
    result = []
    for label in ("dedicated_half", "shared_half"):
        config = parameters["runs"][label]
        prefactor = config["channel_fraction"] * parameters["action"] * KERNEL_UPPER / (6 * parameters["daughter_mass_L"]**3)
        with np.load(history_dir / f"{label}.npz") as data:
            history = {key: data[key] for key in ("t", "a", "D")}
        for case in ("pressure_fit", "small_pulse"):
            for power in (12, 13, 14):
                with np.load(ensembles / f"{label}_{case}_{power}.npz") as data:
                    born, p, sign = (data[key] for key in ("born", "p", "sign"))
                    end, theta, drift, reference = map(float, (data["time"], data["theta"], data["drift"], data["reference_tv"]))
                if len(born) != 2**(power+1) or np.any(np.diff(born) < 0) or not np.all(born[::2] == born[1::2]):
                    raise ValueError("Birth times must stay paired with the chronologically sorted trajectories")
                if not (np.all(sign[::2] == 1) and np.all(sign[1::2] == -1)):
                    raise ValueError("Unexpected charge ordering")
                depths = optical_depth(history, born, end, prefactor)
                labels = assign_cells(p, sign, theta, drift, coarse["ENERGY"], coarse["COSINE"])
                target = coarse["thermal_cells"](theta, drift).ravel()
                target /= target.sum()  # Remove the quadrature residual before using the probability simplex.
                original, observed = tv_floor(labels, depths, target, 0)
                if abs(original - reference) > 2e-10:
                    raise ValueError("Baseline differs from the 198 comparison")
                floor, unhit = tv_floor(labels, depths, target, 1)
                free = 1 - unhit.sum()
                generic = max(0.0, original - free)
                if floor < generic-1e-12 or floor > original+1e-12:
                    raise ValueError("Survival bound violates a counting inequality")
                # A constructive optimizer puts changed mass in target deficits,
                # approaching the protected excess regardless of the collision law.
                deficit = np.maximum(target-unhit, 0)
                optimistic = unhit + free * deficit / deficit.sum()
                witness = 0.5 * np.abs(optimistic-target).sum()
                if abs(witness-floor) > 1e-11:
                    raise ValueError("Optimal redistribution does not attain the bound")
                minimum = {str(tol): required_multiplier(labels, depths, target, tol) for tol in TOLERANCES}
                result.append(dict(history=label, case=case, birth_pairs=2**power, reference_tv=reference,
                    cap_prefactor=prefactor, max_depth=float(depths.max()), max_hit_fraction=float(free),
                    tv_bound_without_cells=generic, tv_bound=floor,
                    rate_multiplier_necessary_for_tolerance=minimum, unhit_lower_cells=unhit.reshape(5,4).tolist(),
                    observed_cells=observed.reshape(5,4).tolist(), target_cells=target.reshape(5,4).tolist()))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", type=Path, required=True, help="Unpacked #198 study")
    parser.add_argument("--particles", type=Path, required=True, help="Regenerated equal-weight particle ensembles")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows = calculate(args.study.resolve(), args.particles)
    args.output.write_text(json.dumps(dict(kernel_upper=KERNEL_UPPER, runs=rows), indent=2, sort_keys=True, allow_nan=False) + "\n")
