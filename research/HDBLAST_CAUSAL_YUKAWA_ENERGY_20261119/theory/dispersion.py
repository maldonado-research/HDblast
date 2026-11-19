"""Independent static dispersion and integrated switch-work calibration checks."""
import argparse
import cmath
import json
import math
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import CubicHermiteSpline


def self_energy(q, cutoff, z):
    log_ratio = cmath.log(cutoff**2 - z**2) - cmath.log(q**2 - z**2)
    return (cutoff**2 - q**2 + (z**2 - q**2) * log_ratio) / (8 * math.pi**2)


def static_checks():
    q, mass, cutoff = 3.0, 4.0, 64.0
    z = 5.0 + 0.7j
    def integrand(omega):
        return omega * (omega**2 - q**2) / (4 * math.pi**2 * (omega**2 - z**2))
    dispersive = quad(integrand, q, cutoff, complex_func=True, epsabs=1e-12)[0]
    quadrature_error = abs(self_energy(q, cutoff, z) - dispersive)
    energy = math.hypot(q, mass)
    # Twice the negative imaginary part of the shifted pole is the number width.
    width = self_energy(q, cutoff, energy + 1e-9j).imag / energy
    reference = mass**2 / (8 * math.pi * energy)
    width_error = abs(width / reference - 1)
    # At imaginary z the subtracted integral is also numerically well behaved.
    # Subtract in the integrand to avoid cancellation of separate O(cutoff**2) terms.
    xi = 2.0
    def subtracted_integrand(omega):
        return -xi**2 * (omega**2 - q**2) / (4 * math.pi**2 * omega * (omega**2 + xi**2))
    cutoff_change = quad(subtracted_integrand, 1e4, 2e4, epsabs=1e-14)[0]
    leading_log = -xi**2 * math.log(4) / (8 * math.pi**2)
    if quadrature_error > 1e-10 or width_error > 1e-7 or abs(cutoff_change / leading_log - 1) > 1e-6:
        raise ArithmeticError("Static dispersion check failed")
    return dict(dispersion_abs_error=quadrature_error, moving_width_over_y2=width,
                width_relative_error=width_error, cutoff_doubling_change=cutoff_change,
                predicted_asymptotic_change=leading_log)


def switch_calibration(native, grid, cutoffs):
    fork = json.loads((native / "FORK.json").read_text())[grid]
    with ZipFile(native / "WALLS.zip") as archive:
        with archive.open(f"{grid}/retain/quiet/wall.csv") as stream:
            data = np.genfromtxt(stream, delimiter=",", names=True)
    times = data["H0tau"]
    loga = CubicHermiteSpline(times, data["A_b"] - fork["reference_A"], data["H_over_H0"])
    mass = CubicHermiteSpline(times, 100 * (data["phi_b"] - 0.5), 100 * data["v_over_H0"])
    q, weight = np.array(fork["q"]), np.array(fork["weight"])
    state = fork["states"]["retain"]
    initial = np.concatenate([state[key] for key in ("X", "C", "P")])
    def rhs(s, moments):
        a = np.exp(loga(s))
        X, C, P = np.split(moments, 3)
        stiffness = a * q**2 + a**3 * mass(s)**2
        return np.r_[2*C/a**3, P/a**3 - stiffness*X, -2*stiffness*C]
    lo, hi = 1.35, 2.05
    sol = solve_ivp(rhs, (fork["H0tau"], hi), initial, method="DOP853",
                    dense_output=True, rtol=2e-11, atol=2e-13, max_step=0.003)
    if not sol.success:
        raise ArithmeticError(sol.message)
    nodes, weights = np.polynomial.legendre.leggauss(256)
    s = lo + (nodes + 1) * (hi-lo)/2
    ds = weights * (hi-lo)/2
    X, C, _ = np.split(sol.sol(s), 3)
    a = np.exp(loga(s))
    phase = math.pi * (s-lo)/(hi-lo)
    f = np.sin(phase)**4
    df = 4 * np.sin(phase)**3 * np.cos(phase) * math.pi/(hi-lo)
    weighted = weight * np.sum(ds * a * f * df * X, axis=1)
    derivative = a * loga(s, 1) * X + 2 * C / a**2
    by_parts = -weight / 2 * np.sum(ds * f**2 * derivative, axis=1)
    error = float(np.max(abs(by_parts-weighted)))
    if error > 1e-8:
        raise ArithmeticError("Independent switch-work integrals disagree")
    # Shift = quadratic * L² + logarithmic * log(L) + constant, for all L > max(q).
    quadratic = np.sum(weighted) / (8 * math.pi**2)
    logarithmic = np.sum(-2 * q**2 * weighted) / (8 * math.pi**2)
    constant = np.sum(q**2 * (2 * np.log(q) - 1) * weighted) / (8 * math.pi**2)
    shifts = {}
    for cutoff in cutoffs:
        c_static = (cutoff**2 - q**2 - 2*q**2*np.log(cutoff/q)) / (8 * math.pi**2)
        value = np.sum(weighted * c_static)
        if not np.isclose(value, quadratic * cutoff**2 + logarithmic * np.log(cutoff) + constant,
                          rtol=1e-13, atol=1e-13):
            raise ArithmeticError("Cutoff dependence disagrees with spectral integral")
        shifts[str(cutoff)] = float(value)
    return dict(quadratic=float(quadratic), logarithmic=float(logarithmic), constant=float(constant),
                max_abs_mode_integration_by_parts_error=error, static_minus_bare_switch_work=shifts)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("native", type=Path, help="archived native directory with FORK.json and WALLS.zip")
    cli.add_argument("--cutoffs", type=int, nargs="+", default=[30, 40, 60, 80, 120, 160])
    args = cli.parse_args()
    if min(args.cutoffs) <= 8:
        cli.error("supply comoving energy cutoffs above all the native nodes (q<8)")
    result = dict(static=static_checks(), switch_calibration={
        grid: switch_calibration(args.native, grid, args.cutoffs) for grid in ("coarse", "fine")})
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
