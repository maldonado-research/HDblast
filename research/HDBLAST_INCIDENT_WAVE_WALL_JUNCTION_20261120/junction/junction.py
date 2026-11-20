"""Absorption of an outward, positive-energy 5D null pulse by a closed radiation brane."""

import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

RADIUS = 2.0
INNER = 1.6


def mass_profile(s, initial, loss):
    if s <= 0:
        return initial, 0.0
    if s >= 1:
        return initial - loss, 0.0
    if s > 0.5:
        t = 1 - s
        before = t**4 * (35 + t * (-84 + t * (70 - 20 * t)))
        mass = initial - loss + loss * before
    else:
        after = s**4 * (35 + s * (-84 + s * (70 - 20 * s)))
        mass = initial - loss * after
    return mass, -loss * 140 * s**3 * (1 - s) ** 3


def scales(initial_h, fraction, initial_e=0.0):
    initial_m = RADIUS**4 * (initial_h**2 + 1 / RADIUS**2 - initial_e * (2 + initial_e))
    if initial_m < 0 or not 0 <= fraction <= 1:
        raise ValueError("Negative mass or invalid integrated pulse")
    return initial_m, initial_m * fraction


def matter(a, b, mass):
    return ((b - 1) ** 2 + 1 / a**2 - mass / a**4) / (2 * b)


def horizon(mass):
    return math.sqrt(2 * mass / (1 + math.sqrt(1 + 4 * mass)))


def snapshot(s, state, initial_m, loss):
    a, b, tau = state[:3]
    mass, dmu = mass_profile(s, initial_m, loss)
    e = matter(a, b, mass)
    h = b - 1 - e
    return dict(
        progress=float(s),
        tau=float(tau),
        a=float(a),
        B=float(b),
        H=float(h),
        e=float(e),
        mass=float(mass),
        comoving_radiation=float(a**4 * e),
        exterior_f=float(1 + a * a - mass / a**2),
    )


def reduced(initial_h, fraction, duration, initial_e=0.0, step=1 / 20):
    initial_m, loss = scales(initial_h, fraction, initial_e)

    def rhs(s, state):
        a, b, tau = state
        mass, _ = mass_profile(s, initial_m, loss)
        e = matter(a, b, mass)
        clock = duration * a * b
        return [clock * a * (b - 1 - e), -clock * (1 / a**2 + 2 * (b - 1) ** 2), clock]

    def turn(s, state):
        a, b, _ = state
        return b - 1 - matter(a, b, mass_profile(s, initial_m, loss)[0])

    turn.direction = -1
    sol = solve_ivp(
        rhs,
        (0.0, 1.0),
        [RADIUS, 1 + initial_h + initial_e, 0.0],
        events=turn,
        method="DOP853",
        rtol=3e-12,
        atol=2e-14,
        max_step=step,
        dense_output=True,
    )
    if not sol.success:
        raise ArithmeticError(sol.message)
    return sol, initial_m, loss


def independent(initial_h, fraction, duration, initial_e=0.0):
    # Proper time, unsquared spatial normal and temporal extrinsic curvature.
    initial_m, loss = scales(initial_h, fraction, initial_e)

    def rhs(tau, state):
        u, a, vel, e, heat = state
        mu, mu_s = mass_profile(u / duration, initial_m, loss)
        mu_u = mu_s / duration
        f = 1 + a * a - mu / a**2
        f_a = 2 * a + 2 * mu / a**3
        f_u = -mu_u / a**2
        x = math.sqrt(f + vel * vel)
        clock = 1 / (x + vel)
        source = -mu_u * clock**2 / (2 * a**3)
        return [
            clock,
            vel,
            (1 - 3 * e) * x - f_a / 2 - f_u * clock**2 / 2,
            source - 4 * vel * e / a,
            a**4 * source,
        ]

    def finish(tau, state):
        return state[0] - duration

    finish.terminal = True
    finish.direction = 1
    sol = solve_ivp(
        rhs,
        (0, max(1.0, 8 * duration)),
        [0.0, RADIUS, RADIUS * initial_h, initial_e, 0.0],
        events=finish,
        method="DOP853",
        rtol=4e-12,
        atol=3e-14,
        max_step=duration / 40,
    )
    if not sol.success or len(sol.t_events[0]) != 1:
        raise ArithmeticError("Outgoing null packet did not finish in proper time")
    u, a, v, e, heat = sol.y_events[0][0]
    return dict(
        H=float(v / a),
        a=float(a),
        e=float(e),
        tau=float(sol.t_events[0][0]),
        deposited_comoving=float(heat),
    )


def post_pulse(end):
    # On a free radiation plateau, (a')²=-1+(mu+2F)/a²+F²/a^6.
    f = max(0.0, end["comoving_radiation"])
    mu = end["mass"]

    def speed2(a):
        return -1 + (mu + 2 * f) / a**2 + f**2 / a**6

    a = end["a"]
    if end["H"] <= 0:
        return None
    upper = 2 * a
    while speed2(upper) > 0:
        upper *= 2
    maximum = brentq(speed2, a, upper, xtol=8e-15)

    # Factor the difference from the root to avoid cancellation at turnaround.
    def smooth(v):
        r = maximum - (maximum - a) * v * v
        p = (mu + 2 * f) / (r * r * maximum**2)
        p += f**2 * (r**4 + r * r * maximum**2 + maximum**4) / (r**6 * maximum**6)
        return 2 * math.sqrt(maximum - a) / math.sqrt((r + maximum) * p)

    time = quad(smooth, 0.0, 1.0, epsabs=2e-13)[0]
    return dict(
        a_max=maximum,
        proper_time_to_maximum=time,
        normalized_dark_fraction_at_maximum=mu / (mu + 2 * f + f**2 / maximum**4),
    )


def run(initial_h, fraction, duration, initial_e=0.0):
    sol, initial_m, loss = reduced(initial_h, fraction, duration, initial_e)
    initial = snapshot(0.0, sol.y[:, 0], initial_m, loss)
    end = snapshot(1.0, sol.y[:, -1], initial_m, loss)
    turning = [
        snapshot(s, state, initial_m, loss) for s, state in zip(sol.t_events[0], sol.y_events[0])
    ]
    summaries = [snapshot(s, sol.sol(s), initial_m, loss) for s in np.linspace(0.0, 1.0, 201)]
    inner_gap = min(row["a"] - INNER for row in summaries)
    inner_f = 1 + INNER**2 - initial_m / INNER**2
    if inner_gap < 0 or inner_f <= 0 or min(row["B"] for row in summaries) <= 0:
        raise ArithmeticError("Calculation left the untrapped annulus")
    # The maximum of the monotone radius segments is at a turnaround or an endpoint.
    if min(initial["a"], end["a"]) <= INNER:
        raise ArithmeticError("Wall crossed source tube")
    monotonic = max(np.diff([row["H"] for row in summaries]))
    positive = min(row["e"] for row in summaries)
    if monotonic > 3e-12 or positive < -3e-13:
        raise ArithmeticError("Raychaudhuri or positive matter failed")
    alt = independent(initial_h, fraction, duration, initial_e)
    errors = {key: abs(alt[key] - end[key]) for key in ("a", "H", "e", "tau")}
    deposited = end["comoving_radiation"] - initial["comoving_radiation"]
    errors["accumulated_heat"] = abs(alt["deposited_comoving"] - deposited)
    integral = quad(
        lambda s: -mass_profile(s, initial_m, loss)[1] / (2 * sol.sol(s)[1]), 0.0, 1.0, epsabs=2e-13
    )[0]
    errors["retarded_flux"] = abs(integral - deposited)
    if max(errors.values()) > 1e-8:
        raise ArithmeticError(errors)

    # Proper emitter and absorber energies, divided by tension times unit S³ volume.
    def fluence(s):
        a, b, tau = sol.sol(s)
        mass, derivative = mass_profile(s, initial_m, loss)
        e = matter(a, b, mass)
        h = b - 1 - e
        clock = duration * a * b
        f_source = 1 + INNER**2 - mass / INNER**2
        return (
            -derivative / (2 * math.sqrt(f_source)),
            -derivative / (2 * a * b),
            e * a**3 * h * clock,
        )

    source_energy, wall_heat, work = [
        quad(lambda s: fluence(s)[i], 0.0, 1.0, epsabs=1e-13)[0] for i in range(3)
    ]
    wall_internal = end["a"] ** 3 * end["e"] - RADIUS**3 * initial_e
    errors["proper_first_law"] = abs(wall_heat - wall_internal - work)
    if errors["proper_first_law"] > 2e-10:
        raise ArithmeticError("Proper brane heat budget differs")
    physical = dict(
        emitter_proper_energy=source_energy,
        wall_proper_heat=wall_heat,
        wall_internal_energy_change=wall_internal,
        wall_expansion_work=work,
        heat_over_source=wall_heat / source_energy if source_energy else None,
    )
    observed = 2 * deposited / loss if loss else None
    if loss and not (1 / initial["B"] - 1e-10 <= observed <= 1 / end["B"] + 1e-10):
        raise ArithmeticError("Deposited energy escaped redshift bounds")
    return dict(
        initial=initial,
        end=end,
        turns=turning,
        post_pulse=post_pulse(end),
        fraction_absorbed_of_initial_mass_parameter=fraction,
        retarded_duration=duration,
        normalized_comoving_transfer=observed,
        proper_energy=physical,
        independent=alt,
        errors=errors,
        minimum_wall_clearance_from_source=inner_gap,
        minimum_source_f=inner_f,
        initial_apparent_horizon=horizon(initial_m),
        minimum_sampled_matter=positive,
    )


def critical(initial_h, duration):
    def final_h(fraction):
        sol, m, l = reduced(initial_h, fraction, duration)
        return snapshot(1.0, sol.y[:, -1], m, l)["H"]

    if final_h(0.0) < 0 or final_h(1.0) > 0:
        return None
    return brentq(final_h, 0.0, 1.0, xtol=5e-13)


def load_sibling(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).with_name(name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def calculate():
    fast = math.sqrt(2) / 4
    slow = 0.04
    cases = {}
    for h, tag, choices in [
        (fast, "expanding", [(0.0, 0.2), (0.5, 0.2), (1.0, 0.2), (0.5, 0.05), (0.5, 0.4)]),
        (slow, "near_turnaround", [(0.0, 0.025), (0.15, 0.025), (0.5, 0.025), (1.0, 0.025)]),
        (0.0, "at_rest", [(0.0, 0.05), (0.5, 0.05)]),
    ]:
        for frac, dur in choices:
            cases[f"{tag}/{frac:g}/{dur:g}"] = run(h, frac, dur)
    impulse = {}
    b0 = 1 + slow
    m0 = scales(slow, 1.0)[0]
    impulse_limit = 2 * RADIUS**4 * b0 * slow / m0
    for dur in (0.04, 0.02, 0.01, 0.005, 0.001, 0.0002):
        threshold = critical(slow, dur)
        row = run(slow, 0.5, dur)
        delta = m0 * 0.5 / (2 * RADIUS**4 * b0)
        impulse[str(dur)] = dict(
            critical_mass_fraction=threshold,
            endpoint_H=row["end"]["H"],
            endpoint_e=row["end"]["e"],
            absolute_offset_from_kick=max(
                abs(row["end"]["H"] - (slow - delta)), abs(row["end"]["e"] - delta)
            ),
            width_in_wall_proper_time=row["end"]["tau"],
        )
    if abs(impulse["0.0002"]["critical_mass_fraction"] / impulse_limit - 1) > 0.005:
        raise ArithmeticError("Smooth pulse threshold differs from impulse limit")
    ref = cases["expanding/0/0.2"]
    a_analytic = math.sqrt(6 - (ref["end"]["tau"] - math.sqrt(2)) ** 2)
    if (
        abs(a_analytic - ref["end"]["a"]) > 1e-11
        or abs(ref["post_pulse"]["a_max"] - math.sqrt(6)) > 3e-13
    ):
        raise ArithmeticError("Pulse-free closed dark radiation solution failed")
    unequal = load_sibling("asymmetric").calculate(mass_profile)
    pressure = load_sibling("normal_pressure").calculate()
    prescribed = load_sibling("prescribed_source").calculate()
    radiation_exit = load_sibling("radiation_exit").calculate()
    scalar = load_sibling("scalar_capture")
    capture = scalar.calculate()
    validation = load_sibling("scalar_validation").calculate(scalar)
    primary = capture["cases"]["primary"]
    replay_checks = (
        (validation["bounce_radius"], primary["bounce"]["a"]),
        (validation["bounce_time"], primary["bounce"]["t"]),
        (validation["emitted_proper_energy"], primary["energy"]["two_side_emitter_energy"]),
        (validation["wall_deposited_energy"], primary["energy"]["two_side_wall_deposition"]),
    )
    if max(abs(actual - expected) for actual, expected in replay_checks) > 5e-9:
        raise ArithmeticError("Prescribing the reconstructed bulk mass changed the scalar result")
    symmetric = cases["expanding/0.5/0.2"]["end"]
    for left, right in [
        ("end_H", "H"),
        ("end_e", "e"),
        ("end_radius", "a"),
        ("end_proper_time", "tau"),
    ]:
        if abs(unequal["two_sided_half"][left] - symmetric[right]) > 2e-10:
            raise ArithmeticError("Two-sided independent systems differ")
    return dict(
        asymmetric=unequal,
        normal_pressure=pressure,
        prescribed_source=prescribed,
        radiation_exit=radiation_exit,
        scalar_capture=capture,
        scalar_validation=validation,
        units="AdS radius=1, 5D gravitational coupling absorbed into tuned tension, spherical k=+1",
        smooth_mass_profile="regularized beta(4,4) CDF in retarded time; X_i=a_i(1+e_i)",
        source_radius=INNER,
        launch_radius=RADIUS,
        assumed_wall_eos="Radiation p=rho/3 in cases, asymmetric and impulse. Other named sections specify their wall stresses.",
        cases=cases,
        impulse=dict(
            initial_H=slow,
            critical_fraction_limit=impulse_limit,
            half_mass_kick=m0 * 0.5 / (2 * RADIUS**4 * b0),
            smooth=impulse,
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(calculate(), indent=2, allow_nan=False) + "\n")
