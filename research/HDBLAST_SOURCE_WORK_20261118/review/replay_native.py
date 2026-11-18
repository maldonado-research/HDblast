"""Propagate an independent Gaussian transfer matrix on the archived coupled wall histories."""

import argparse
import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from scipy.interpolate import BarycentricInterpolator, CubicHermiteSpline

ARMS = ("quiet", "plus", "minus", "far")
NAMES = ("N", "rad_hat", "rho_hat", "p_hat", "j_hat")
EPOCH = 2.10


def read_csv(archive, name):
    with archive.open(name + ".csv") as stream:
        data = np.genfromtxt(stream, delimiter=",", names=True, dtype=float)
    return {key: data[key] for key in data.dtype.names}


def trajectories(archive, previous, grid):
    summary = previous / "classical/runs" / grid / "quiet/rolling_summary.json"
    loga_ref = json.loads(summary.read_bytes())["matter"]["events"][0]["A_star"]
    histories = []
    for decay in (0, 1):
        for feedback in (0, 1):
            base = f"{grid}/y{decay}_b{feedback}/"
            pre = read_csv(archive, base + "pre")
            for arm in ARMS:
                post = read_csv(archive, base + arm)
                clock = np.r_[pre["H0tau"], post["H0tau"]]
                if min(np.diff(clock)) <= 0 or not clock[0] < EPOCH < clock[-1]:
                    raise ValueError(f"Invalid clock for {base}{arm}")

                def wall(key, dkey):
                    values = np.r_[pre[key], post[key]]
                    derivatives = np.r_[pre[dkey], post[dkey]]
                    return CubicHermiteSpline(clock, values, derivatives)

                # Pressure and force retain fast phases even after the number stabilizes.
                center = np.searchsorted(post["H0tau"], EPOCH)
                near = slice(center - 4, center + 4)
                sample_times = post["H0tau"][near] - EPOCH
                stored = {}
                for key in NAMES:
                    field = "comoving_number" if key == "N" else key
                    stored[key] = float(BarycentricInterpolator(
                        sample_times, post[field][near], rng=0
                    )(0.0))
                histories.append(dict(grid=grid, y=decay, feedback=feedback, arm=arm,
                                      start=clock[0], loga_ref=loga_ref, stored=stored,
                                      loga=wall("A_b", "H_over_H0"), phi=wall("phi_b", "v_over_H0")))
    return histories


def propagate(histories, steps):
    midpoints = (np.arange(steps) + .5) / steps
    times_fraction = np.r_[0., midpoints, 1.]
    a_series, mass_series, dt_series = [], [], []
    for history in histories:
        duration = EPOCH - history["start"]
        samples = history["start"] + duration * times_fraction
        a_series.append(np.exp(history["loga"](samples) - history["loga_ref"]))
        mass_series.append(100 * (history["phi"](samples) - .5))
        dt_series.append(duration / steps)
    a = np.array(a_series).T
    mass = np.array(mass_series).T
    dt = np.array(dt_series)[:, None]
    decay = np.array([item["y"] for item in histories])[:, None]
    nodes, weights = np.polynomial.legendre.leggauss(24)
    k = 4 * (1 + nodes)
    weights = 4 * weights * k**2 / (2 * np.pi**2)

    def frequencies(scale, mass_value):
        energy = np.hypot(k / scale[:, None], mass_value[:, None])
        return energy, scale[:, None]**3 * energy

    _, b = frequencies(a[0], mass[0])
    var_x, cov, var_p = .5 / b, np.zeros_like(b), .5 * b
    radiation = np.zeros_like(b)
    for scale, mu in zip(a[1:-1], mass[1:-1], strict=True):
        energy, b = frequencies(scale, mu)
        angle = energy * dt
        # Ground-state attenuation commutes with Hamiltonian rotation at fixed a and mass.
        co, si = np.cos(angle), np.sin(angle)
        gamma = (decay * mu[:, None])**2 / (8 * np.pi * energy)
        loss = -np.expm1(-gamma * dt)
        number = .5 * (b * var_x + var_p / b - 1)
        radiation += scale[:, None] * energy * number * loss
        var_x, cov, var_p = (
            (1-loss) * (co**2 * var_x + 2*co*si*cov/b + si**2 * var_p/b**2) + .5*loss/b,
            (1-loss) * (-co*si*b*var_x + (co**2-si**2)*cov + co*si*var_p/b),
            (1-loss) * (si**2*b**2*var_x - 2*co*si*b*cov + co**2*var_p) + .5*loss*b,
        )
        if np.min(var_x * var_p - cov**2) < .25 - 2e-10:
            raise ValueError("Transfer violated the uncertainty bound")

    scale = a[-1][:, None]
    mu = mass[-1][:, None]
    energy, b = frequencies(a[-1], mass[-1])
    dx = var_x - .5 / b
    dp = (var_p - .5*b) / scale**6
    number = .5 * (b * var_x + var_p/b - 1)
    pressure = .5 * (dp - (mu**2 + (k/scale)**2/3) * dx)
    return np.array([number @ weights, (radiation @ weights) / a[-1]**4,
                     (.5*(dp + energy**2*dx)) @ weights, pressure @ weights,
                     (100*mu*dx) @ weights]).T


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path, help="coupled/WALLS.zip")
    parser.add_argument("previous_capsule", type=Path, help="HDBLAST_PULSE_PARTICLE_RESPONSE_20261118")
    args = parser.parse_args()
    report = {}
    with ZipFile(args.archive) as archive:
        for grid in ("coarse", "fine"):
            histories = trajectories(archive, args.previous_capsule, grid)
            coarse, fine = (propagate(histories, n) for n in (2200, 4400))
            extrapolated = fine + (fine - coarse) / 3
            results, errors = [], []
            for history, row, ex in zip(histories, fine, extrapolated, strict=True):
                entry = {key: history[key] for key in ("grid", "y", "feedback", "arm", "stored")}
                entry["transfer"] = dict(zip(NAMES, row.tolist(), strict=True))
                entry["extrapolated"] = dict(zip(NAMES, ex.tolist(), strict=True))
                stored = np.array([entry["stored"][key] for key in NAMES])
                errors.append(np.abs(ex - stored) / np.maximum(1., np.abs(stored)))
                results.append(entry)
            maxima = dict(zip(NAMES, np.max(errors, axis=0).tolist(), strict=True))
            if max(maxima.values()) > 1.e-4:
                raise ValueError((grid, maxima))
            report[grid] = {"results": results, "max_scaled_errors": maxima}
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
