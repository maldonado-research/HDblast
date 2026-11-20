"""Check the scalar-parent mode counts with adaptive Bogoliubov transport.

Run: OPENBLAS_NUM_THREADS=1 python control.py /path/to/CHARGED.json
The control compares 24 finite-band counts to MODES.json beside this script.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import BPoly
from scipy.optimize import brentq
from scipy.special import expit


def transport(source):
    data = json.loads(source.read_bytes())
    reference = json.loads(Path(__file__).with_name("MODES.json").read_text())
    if hashlib.sha256(source.read_bytes()).hexdigest() != reference["input_sha256"]:
        raise ValueError("Mismatched input archive")
    run = data["runs"]["neutral_precursor"]
    col = dict(zip(data["trace_columns"], np.array(run["trace"]).T))
    t = col["t"]

    # Reconstruct each cubic independently in the Bernstein representation.
    def interpolate(value, slope):
        return BPoly.from_derivatives(t, np.column_stack((value, slope)))

    a = interpolate(col["a"], col["a"]*col["H"])
    psi = interpolate(col["psi"], col["pi"])
    da, dpsi = a.derivative(), psi.derivative()
    opening = brentq(lambda tau: psi(tau)-0.77, 0., run["maximum_tau"])
    readouts = np.array([opening, 10.5017, run["maximum_tau"]])
    j = np.arange(1, 130, dtype=float)
    state = np.zeros((3, 2, len(j)), dtype=complex)
    state[1] = 1.0  # theta=0, alpha=1, beta=0 for each mode

    def system(tau, y):
        radius, psi_value = float(a(tau)), float(psi(tau))
        z = (psi_value-0.77)/0.03 + 0.5
        if 0.0 < z < 1.0:
            x = 1/(1-z) - 1/z
            step = expit(x)
            slope = expit(-x)*step*(z**-2 + (1-z)**-2)*float(dpsi(tau))/0.03
        else:
            step, slope = float(z >= 1.), 0.
        mass = np.array([1.90 + 0.20*step, 1.90])[:, None]
        mass_prime = np.array([0.20*slope, 0.])[:, None]
        am = radius*mass
        omega2 = j*j + am*am
        frequency = np.sqrt(omega2)/radius
        mixing = am*(float(da(tau))*mass + radius*mass_prime)/(2*omega2)
        phase, alpha, beta = y.reshape(state.shape)
        rotor = np.exp(2j*phase.real)
        return np.array([frequency, mixing*rotor*beta,
                         mixing*rotor.conj()*alpha]).ravel()

    # A separate integration between knots avoids crossing a nonsmooth
    # higher derivative and supplies the exact readout slices.
    stops = np.unique(np.r_[t[(t > 0) & (t < readouts[-1])], readouts])
    start = 0.
    counts, residual = [], 0.
    for stop in stops:
        sol = solve_ivp(system, (start, stop), state.ravel(), method="DOP853",
                        rtol=2e-10, atol=3e-12, t_eval=[stop])
        if not sol.success:
            raise RuntimeError(sol.message)
        state = sol.y[:, -1].reshape(state.shape)
        start = stop
        if stop in readouts:
            occupation = abs(state[2])**2
            residual = max(residual, np.max(abs(abs(state[1])**2-occupation-1)))
            counts.append(np.cumsum(j*j*occupation, axis=1)[:, [0,32,64,128]])
    difference = float(np.max(abs(np.array(counts)-reference["count"])))
    if difference > 1e-6 or residual > 1e-8:
        raise ArithmeticError("Adaptive cross-check exceeds numerical tolerances")
    return dict(times=readouts.tolist(), counts=np.asarray(counts).tolist(),
                maximum_count_difference=difference, maximum_norm_error=float(residual))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    args = parser.parse_args()
    print(json.dumps(transport(args.source), indent=2))
