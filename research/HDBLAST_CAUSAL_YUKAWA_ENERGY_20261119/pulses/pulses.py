"""Paired perturbative energy and work differences on the four released walls."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_simpson, simpson, solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.signal import fftconvolve

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("budget", HERE.parent / "solver/budget.py")
budget = importlib.util.module_from_spec(spec)
spec.loader.exec_module(budget)


def independent_switch(q, weights, source, background, cutoff, steps=1024):
    """Integrate the retarded commutator and vacuum noise directly in conformal time."""
    start, loga, mass = background
    s = np.linspace(start, 2.10, 4 * steps + 1)
    clock = cumulative_simpson(np.exp(-loga(s)), x=s, initial=0.)
    eta = np.linspace(0., clock[-1], steps + 1)
    dt = eta[1]
    s = PchipInterpolator(clock, s)(eta)
    a = np.exp(loga(s))
    phase = np.pi * (s - 1.35) / .70
    active = (s > 1.35) & (s < 2.05)
    f = np.where(active, np.sin(phase)**4, 0.)
    df = np.where(active, 4 * np.pi / .70 * np.sin(phase)**3 * np.cos(phase), 0.)
    inertia = a[0]**3 * np.hypot(q / a[0], mass(start))
    mode0 = 1 / np.sqrt(2 * inertia)

    def rhs(t, z):
        at = np.exp(loga(t))
        u, p = np.split(z, 2)
        return np.r_[p / at**3, -at**3 * ((q / at)**2 + mass(t)**2) * u]

    free = solve_ivp(rhs, (start, s[-1]), np.r_[mode0, -1j * inertia * mode0],
                     t_eval=s, method="DOP853", max_step=.006, rtol=2e-11, atol=1e-13)
    if not free.success:
        raise RuntimeError(free.message)
    X, C, P = source
    occupation = (inertia * X + P / inertia - 1) / 2
    coherence = (inertia * X - P / inertia) / 2 + 1j * C
    nodes, weights_omega = np.polynomial.legendre.leggauss(96)
    total, correction = 0., 0.
    for k, wq, u, n, c in zip(q, weights, free.y[:len(q)], occupation, coherence, strict=True):
        omega = k + (cutoff - k) * (nodes + 1) / 2
        dw = weights_omega * (cutoff - k) / 2
        spectral = dw * (omega**2 - k**2) / (8 * np.pi**2)
        angles = omega[:, None] * eta
        response = 2 * (spectral @ np.sin(angles))
        noise = spectral @ np.cos(angles)
        drive = a * f * u

        # Composite trapezoid endpoint at the current time; the drive at the
        # initial time vanishes. The noise contact is real after projection.
        def causal(kernel):
            return dt * (fftconvolve(kernel, drive)[:len(eta)] - .5 * kernel[0] * drive)

        response_mode, noise_mode = causal(response), causal(noise)
        # -G_chi times the bath symmetric correlator plus the field covariance
        # times the bath response. G_chi(t,t') = i[u(t)u*(t')-u*(t)u(t')].
        mixed = (2 * (u * noise_mode.conj()).imag
                 - (2*n + 1) * (u * response_mode.conj()).real
                 - 2 * (c * u * response_mode).real)
        x_free = (2*n + 1) * abs(u)**2 + 2 * (c * u**2).real
        static = np.sum(2 * spectral / omega)
        total += wq * simpson(a * df * mixed, x=eta)
        correction += wq * static * simpson(a**2 * f * df * x_free, x=eta)
    return np.array([total, total + correction])


def summarize(row):
    e, w = row["energy"], row["work"]
    return {"daughter": e["daughter"], "parent": e["parent"],
            "switch": w["switch"] + w["switch_counterterm"], "mass": w["mass"],
            "scale": w["scale_parent"] + w["scale_daughter"] + w["scale_counterterm"],
            "closure": row["residual"]["total"]}


def calculate(native):
    output = {"units": "coefficients of y^2 H0^4, per comoving volume", "runs": {}}
    for grid in ("fine", "coarse"):
        for cutoff in (60, 120):
            arms = {}
            for arm in ("quiet", "plus", "minus", "far"):
                q, w, states, bg = budget.read_native(native, grid, arm)
                init = {"retain": states["retain"]}
                data = budget.run_budget(q, w, init, bg, cutoff, steps=1024)["retain"]
                coarse = independent_switch(q, w, init["retain"], bg, cutoff, steps=512)
                fine = independent_switch(q, w, init["retain"], bg, cutoff, steps=1024)
                independent = (4 * fine - coarse) / 3
                arms[arm] = {}
                for i, scheme in enumerate(("bare", "static")):
                    row = summarize(data[scheme]["checkpoints"]["cut"])
                    row["independent_switch"] = float(independent[i])
                    arms[arm][scheme] = row
            delta = {arm: {scheme: {key: value - arms["quiet"][scheme][key]
                                    for key, value in row.items()} for scheme, row in values.items()}
                     for arm, values in arms.items() if arm != "quiet"}
            output["runs"][f"{grid}_{cutoff}"] = {"quiet": arms["quiet"], "pulse_minus_quiet": delta}
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("native", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.native), indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
