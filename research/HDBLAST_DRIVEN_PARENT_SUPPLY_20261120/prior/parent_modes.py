"""Neutral conformal scalar on HDblast's archived neutral precursor trajectory.

Usage: python parent_modes.py <path to the archived CHARGED.json>
Input: https://github.com/maldonado-research/HDblast/blob/add30482b0aba462518c062f119bfa991f3d038e/research/HDBLAST_NEUTRAL_PARENT_DECAY_20261120/prior/CHARGED.json
Requires numpy and scipy. All counts refer to a real neutral scalar on one
round sphere in the instantaneous oscillator basis. The prescribed background
already includes the fluid precursor. Decay and backreaction are switched off.
"""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq
from scipy.special import expit


DIGEST = "5dbed28e158584403a4cd238d392aab15c00480522cf489d1ee0a36f1d550f90"
MAX_SHELL = 128


def history(source):
    raw = source.read_bytes()
    if hashlib.sha256(raw).hexdigest() != DIGEST:
        raise ValueError("Expected the pinned HDblast neutral-precursor input")
    archive = json.loads(raw)
    run = archive["runs"]["neutral_precursor"]
    columns = dict(zip(archive["trace_columns"], np.asarray(run["trace"]).T))
    times = columns["t"]
    a = CubicHermiteSpline(times, columns["a"], columns["a"] * columns["H"])
    psi = CubicHermiteSpline(times, columns["psi"], columns["pi"])
    end = run["maximum_tau"]
    opening = brentq(lambda t: psi(t) - 0.77, 0.0, end)
    return a, psi, [opening, 10.5017, end]


def parent_masses(scalar):
    z = np.atleast_1d((scalar - 0.77) / 0.03 + 0.5)
    step = np.zeros_like(z)
    step[z >= 1.0] = 1.0
    inner = (z > 0) & (z < 1)
    step[inner] = expit(1 / (1 - z[inner]) - 1 / z[inner])
    # First row reproduces HDblast's parent mass law with daughter mass mL=1.
    # Second freezes it at the original mass on the identical geometry.
    return np.stack((1.9 + 0.2*step, np.full_like(step, 1.9)))


def propagate(a, psi, readouts, max_step):
    j = np.arange(1, MAX_SHELL + 2, dtype=float)
    initial_omega = np.hypot(j, a(0.) * parent_masses(psi(0.)))
    u = (1 / np.sqrt(2*initial_omega)).astype(complex)
    p = -1j*initial_omega*u  # p = du/d(conformal time)
    start = 0.0
    results = []
    for end in readouts:
        nsteps = math.ceil((end-start)/max_step)
        dt = (end-start)/nsteps
        mid = start + (np.arange(nsteps) + 0.5)*dt
        radii = a(mid)
        scaled_masses = parent_masses(psi(mid)) * radii
        for radius, masses in zip(radii, scaled_masses.T):
            omega = np.hypot(j, masses[:, None])
            phase = omega*(dt/radius)
            co, si = np.cos(phase), np.sin(phase)
            u, p = co*u + (si/omega)*p, co*p - (omega*si)*u
        omega = np.hypot(j, a(end)*parent_masses(psi(end)))
        occupation = abs((omega*u-1j*p)/np.sqrt(2*omega))**2
        if np.max(abs(2*(u*p.conjugate()).imag-1)) > 5e-10:
            raise ArithmeticError("Oscillator Wronskian drift")
        counts = np.cumsum(occupation*j**2, axis=1)
        results.append(counts[:, [0, 32, 64, MAX_SHELL]])
        start = end
    return np.asarray(results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    args = parser.parse_args()
    a, psi, readouts = history(args.source)
    coarse = propagate(a, psi, readouts, 0.00025)
    fine = propagate(a, psi, readouts, 0.000125)
    if np.max(abs(coarse-fine)) > 1e-6:
        raise ArithmeticError("Time refinement differs by more than 1e-6")
    print(json.dumps({
        "input_sha256": DIGEST, "proper_times_in_L": readouts,
        "highest_harmonic_in_each_column": [0, 32, 64, MAX_SHELL],
        "mass_laws_in_each_row": ["HDblast_gate_mu_1", "constant_parent_ML_1.9"],
        "count": fine.tolist(),
        "max_count_change_on_time_refinement": float(np.max(abs(coarse-fine))),
    }, indent=2))
