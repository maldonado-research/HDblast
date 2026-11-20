"""Replay the reference scalar run with a tabulated bulk mass and the temporal junction."""

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq


def calculate(scalar):
    width, kick, height = scalar.DURATION, 0.485, 0.111

    def force(t):
        return kick * scalar.bump(t / width) / (width * scalar.BUMP_AREA)

    def source_rhs(t, state):
        # Reconstruct the source on the short absorption interval only.
        a, h, mass, phi, pi, u, deposited = state
        potential, slope = scalar.potential(phi, height)
        e, p = pi * pi / 2 + potential, pi * pi / 2 - potential
        j = force(t)
        q = pi * j
        b = 1 + e + h
        acceleration = (1 + e) * (1 - 2 * e - 3 * p) - 1 - mass / a**4 - q
        return [
            a * h,
            acceleration - h * h,
            -2 * a**4 * b * q,
            pi,
            -3 * h * pi - slope + j,
            1 / (a * b),
            a**3 * q,
        ]

    def integrate(rhs, interval, initial, step, events=None):
        sol = solve_ivp(
            rhs,
            interval,
            initial,
            method="DOP853",
            rtol=2e-13,
            atol=4e-15,
            max_step=step,
            dense_output=True,
            events=events,
        )
        if not sol.success:
            raise ArithmeticError(sol.message)
        return sol

    source = integrate(source_rhs, (0.0, width), [2.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0], width / 80)
    times = np.linspace(0.0, width, 6001)
    y = source.sol(times)
    slopes = np.array([source_rhs(t, z) for t, z in zip(times, y.T)])
    rates = slopes[:, 2] / slopes[:, 5]
    if max(rates) > 1e-12 or min(y[2]) < 0 or min(slopes[:, 5]) <= 0:
        raise ArithmeticError(
            "Source violates the outgoing null-energy condition or retarded clock"
        )
    profile = CubicHermiteSpline(y[5], y[2], rates)
    retarded_end, mass_final = y[5, -1], y[2, -1]

    def bulk_mass(u):
        if u <= 0:
            return 4.0, 0.0
        if u >= retarded_end:
            return mass_final, 0.0
        return float(profile(u)), float(profile(u, 1))

    def wall_rhs(t, state):
        # Once the mass is prescribed on null slices, the geometric evolution
        # no longer integrates the primary acceleration or mass-balance ODE.
        a, v, phi, pi, u, deposit, work = state
        mass, rate = bulk_mass(u)
        x = math.sqrt(1 + a * a - mass / a**2 + v * v)
        clock = 1 / (x + v)
        q_from_bulk = -rate * clock**2 / (2 * a**3)
        potential, slope = scalar.potential(phi, height)
        e, p = pi * pi / 2 + potential, pi * pi / 2 - potential
        acceleration = (1 - 2 * e - 3 * p) * x - a - mass / a**3 + rate * clock**2 / (2 * a**2)
        return [
            v,
            acceleration,
            pi,
            -3 * v * pi / a - slope + force(t),
            clock,
            a**3 * q_from_bulk,
            3 * a * a * v * p,
        ]

    def bounce(t, state):
        return state[1]

    bounce.direction = 1

    def plateau(t, state):
        return state[2] - scalar.HIGH

    plateau.direction = 1

    pulse = integrate(wall_rhs, (0.0, width), [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], width / 150)
    after = integrate(wall_rhs, (width, 1.2), pulse.y[:, -1], 0.003, events=(bounce, plateau))
    if len(after.t_events[0]) != 1 or len(after.t_events[1]) != 1:
        raise ArithmeticError("Independent wall failed to reach the plateau and bounce")
    bounce_a = after.y_events[0][0, 0]
    a_entry, _, _, pi_entry, *_ = after.y_events[1][0]
    momentum = a_entry**3 * pi_entry

    def speed_squared(a):
        return (
            height * (2 + height) * a * a
            - 1
            + mass_final / a**2
            + (1 + height) * momentum**2 / a**4
            + momentum**4 / (4 * a**10)
        )

    def speed_squared_derivative(a):
        return (
            2 * height * (2 + height) * a
            - 2 * mass_final / a**3
            - 4 * (1 + height) * momentum**2 / a**5
            - 2.5 * momentum**4 / a**11
        )

    minimum = brentq(speed_squared_derivative, 0.1, 4.0)
    small_root = brentq(speed_squared, 0.1, minimum)
    large_root = brentq(speed_squared, minimum, 4.0)
    source_radius = scalar.INNER

    def source_lapse(mass):
        return math.sqrt(1 + source_radius**2 - mass / source_radius**2)

    source_energy = source_radius**2 * (source_lapse(mass_final) - source_lapse(4.0))

    angular_error = first_law_error = 0.0
    for segment in (pulse, after):
        for state in segment.sol(np.linspace(segment.t[0], segment.t[-1], 800)).T:
            a, v, phi, pi, u, deposit, work = state
            mass, _ = bulk_mass(u)
            e = pi * pi / 2 + scalar.potential(phi, height)[0]
            angular_error = max(
                angular_error, abs(math.sqrt(1 + a * a - mass / a**2 + v * v) / a - 1 - e)
            )
            first_law_error = max(first_law_error, abs(a**3 * e + work - deposit))
    deposit = pulse.y[5, -1]
    errors = dict(
        angular_junction=angular_error,
        proper_energy_budget=first_law_error,
        plateau_root=abs(bounce_a - large_root),
        incident_flux_vs_scalar_work=abs(deposit - y[6, -1]),
    )
    if (
        max(errors.values()) > 5e-9
        or not small_root < source_radius < bounce_a < a_entry
        or not 0 < deposit < source_energy
    ):
        raise ArithmeticError("Independent metric, energy, or accessible-root test failed")
    return dict(
        retarded_width=float(retarded_end),
        mass_after=float(mass_final),
        emitted_proper_energy=float(source_energy),
        wall_deposited_energy=float(deposit),
        wall_pressure_work_during_pulse=float(pulse.y[6, -1]),
        bounce_time=float(after.t_events[0][0]),
        bounce_radius=float(bounce_a),
        plateau_momentum=float(momentum),
        lower_root=small_root,
        upper_root=large_root,
        errors={k: float(v) for k, v in errors.items()},
    )
