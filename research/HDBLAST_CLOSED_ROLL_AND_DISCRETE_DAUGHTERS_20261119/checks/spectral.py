#!/usr/bin/env python3
"""Independent continuum reference and bounds for a compact sin^4 drive."""
import argparse
import json
from itertools import combinations
from math import cos, floor, fsum, log, pi, sin, sqrt
from pathlib import Path

from scipy.integrate import quad
from scipy.special import sici


def harmonics(cycles):
    """Cosine coefficients of f'', merging reflected and degenerate frequencies."""
    twice_cycles = round(2 * cycles)
    if cycles <= 0 or twice_cycles != 2 * cycles:
        raise ValueError("cycles must be a positive integer or half integer")
    weights = {}
    for shift, weight in zip(range(-2, 3), (1/16, -1/4, 3/8, -1/4, 1/16)):
        index = abs(twice_cycles + 2 * shift)
        if index:
            weights[index] = weights.get(index, 0) + weight
    return [(i / twice_cycles, -v * (i / twice_cycles)**2)
            for i, v in sorted(weights.items())]


def continuum(cycles, epsilon):
    duration = 2 * pi * cycles
    terms = harmonics(cycles)
    phase_sign = (-1)**round(2 * cycles)
    sine_integral = {}
    cosine_integral = {}
    for frequency, _ in terms:
        phase = frequency * duration
        si, ci = sici(phase)
        sine_integral[frequency] = si
        # The Euler constant cancels in the full result, because f''(0)=0.
        cosine_integral[frequency] = ci - log(phase)
    integrals = []
    for a, weight in terms:
        integrals.append(weight**2 * (cosine_integral[a] + a * duration / 2
                         * sine_integral[a] - (1 - phase_sign) / 2))
    for (a, wa), (b, wb) in combinations(terms, 2):
        divided_difference = ((a**2 * cosine_integral[a] - b**2 * cosine_integral[b])
                              / (a**2 - b**2))
        integrals.append(2 * wa * wb * divided_difference)
    energy = epsilon**2 / (16 * pi**2) * fsum(integrals)
    number = epsilon**2 * duration / (32 * pi) * fsum(w**2 for _, w in terms)
    width_energy = epsilon**2 * duration * (35/128) / (64 * pi)
    return dict(number=number, energy=energy, harmonic_width_energy=width_energy,
                energy_over_width=energy / width_energy if width_energy else None)


def tail_bounds(cycles, radius_mass, cutoff, epsilon):
    """Upper bounds on omitted leading-order energies at k/m > cutoff."""
    nu = 2 * cutoff
    b = 1 / cycles
    if nu <= 1 + 2 * b:
        raise ValueError("UV bound requires a cutoff above the last sideband")
    scale = (3 * b**4 / ((1 - 1/nu)**5 * (1 - (b/(nu - 1))**2)
                         * (1 - (2*b/(nu - 1))**2)))
    first_j = max(2, floor(sqrt(1 + (radius_mass * cutoff)**2)) + 1)
    continuum_tail = epsilon**2 * scale**2 / (128 * pi**2 * nu**4)
    sphere_tail = (epsilon**2 * scale**2 * radius_mass**4 / (512 * pi**2)
                   * (1 - 1/first_j**2)**(-3.5)
                   * (first_j**-5 + first_j**-4 / 4))
    return continuum_tail, sphere_tail


def time_integral_check(cycles, epsilon):
    """Evaluate the real-time expression without the harmonic autocorrelation."""
    duration = 2 * pi * cycles
    terms = harmonics(cycles)
    def derivative_integrand(t, h):
        g = sum(w * cos(a * t) for a, w in terms)
        gprime = -sum(w * a * sin(a * (t + h)) for a, w in terms)
        return -g * gprime / h
    def integrand(h):
        value, _ = quad(derivative_integrand, 0, duration - h, args=(h,),
                        epsabs=2e-12, epsrel=3e-12, limit=200)
        return value
    integrated, _ = quad(integrand, 0, duration, epsabs=2e-10, epsrel=2e-11,
                         limit=160)
    return epsilon**2 * integrated / (16 * pi**2)


def reference_results():
    epsilon = 1e-3
    cases = {}
    for cycles in (4, 8, 16):
        case = continuum(cycles, epsilon)
        if cycles == 4:
            case["direct_time_integrals_energy"] = time_integral_check(cycles, epsilon)
            assert abs(case["direct_time_integrals_energy"] / case["energy"] - 1) < 2e-11
        case["continuum_energy_tail_upper_above_k4"] = tail_bounds(cycles, 1, 4, epsilon)[0]
        cases[str(cycles)] = case
    return {"units": "m_parent=1, number per physical volume, excitation energy per physical volume",
            "epsilon": epsilon, "window": "sin(pi*t/T)^4 cos(t), T=2*pi*N", "cases": cases}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    options = parser.parse_args()
    options.output.write_text(json.dumps(reference_results(), indent=2, allow_nan=False) + "\n")
