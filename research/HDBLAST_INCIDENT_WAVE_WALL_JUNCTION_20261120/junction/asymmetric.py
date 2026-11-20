"""Compare independent irradiation of the ordinary interiors of a radiation brane."""

import math

import numpy as np
from scipy.integrate import solve_ivp


def solve(initial_h, fractions, width, profile, initial_e=0.0, comparison_time=None):
    a0 = 2.0
    m0 = a0**4 * (initial_h**2 + 1 / a0**2 - initial_e * (2 + initial_e))
    if m0 < 0 or initial_e < 0 or width <= 0 or max(fractions) > 1 or min(fractions) < 0:
        raise ValueError("Invalid mass, radiation, width or loss fraction")
    if 1 + 1.6**2 - m0 / 1.6**2 <= 0:
        raise ValueError("Emitter lies on or inside its initial horizon")
    losses = m0 * np.array(fractions)
    active = losses > 0
    if not any(active):
        raise ValueError("At least one incident stream is required")

    def bulk(clocks):
        values = np.array([profile(u / width, m0, float(loss)) for u, loss in zip(clocks, losses)])
        return values[:, 0], values[:, 1] / width

    def physical(a, v, clocks, e):
        masses, rates = bulk(clocks)
        f = 1 + a * a - masses / a**2
        x = np.sqrt(f + v * v)
        radial = x / a
        udot = 1 / (x + v)
        weights = (1 / radial) / sum(1 / radial)
        hm = 2 / sum(1 / radial)
        normal_pressure = masses / a**4 - rates * udot**2 / (2 * a**3)
        avg = float(weights @ normal_pressure)
        acc_per_a = (1 - 3 * e) * hm - 1 - avg
        source = -sum(rates * udot**2) / (4 * a**3)
        chi = (masses[0] - masses[1]) / (4 * a**4 * (1 + e))
        constraint = (v / a) ** 2 - (e * (2 + e) - 1 / a**2 + np.mean(masses) / a**4 + chi**2)
        return (
            udot,
            source,
            acc_per_a,
            dict(
                masses=masses,
                radial=radial,
                hm=hm,
                avg=avg,
                asym=chi,
                angular_residual=abs(np.mean(radial) - (1 + e)),
                friedmann_residual=abs(constraint),
                shifted_acceleration=acc_per_a + avg,
                upper_shifted=-2 * e - 3 * e * e if e <= 1 / 3 else -1.0,
            ),
        )

    def full_rhs(t, z):
        a, v, u1, u2, e, deposit = z
        du, src, acc, _ = physical(a, v, (u1, u2), e)
        return [v, a * acc, *du, src - 4 * (v / a) * e, src * a**4]

    def end(t, z):
        return min(np.asarray(z[2:4])[active]) - width

    end.terminal = True
    end.direction = 1

    def turn(t, z):
        return z[1]

    turn.direction = -1
    full = solve_ivp(
        full_rhs,
        (0.0, 8 * width),
        [a0, a0 * initial_h, 0.0, 0.0, initial_e, 0.0],
        rtol=3e-12,
        atol=3e-14,
        method="DOP853",
        max_step=width / 28,
        events=[end, turn],
        dense_output=True,
    )
    if not full.success or len(full.t_events[0]) != 1:
        raise ArithmeticError("Asymmetric pulse failed")
    a, v, u1, u2, e, deposit = full.y_events[0][0]
    stop = full.t_events[0][0]
    errs = dict(
        angular=0.0,
        friedmann=0.0,
        heat=float(abs(deposit - a**4 * e + a0**4 * initial_e)),
        harmonic_bound=0.0,
    )
    max_acc = -math.inf
    min_clearance = math.inf
    min_e = math.inf
    for t in np.linspace(0, stop, 201):
        r, vel, s1, s2, matter, _ = full.sol(t)
        _, _, acc, diag = physical(r, vel, (s1, s2), matter)
        errs["angular"] = max(errs["angular"], float(diag["angular_residual"]))
        errs["friedmann"] = max(errs["friedmann"], float(diag["friedmann_residual"]))
        errs["harmonic_bound"] = max(
            errs["harmonic_bound"], float(diag["shifted_acceleration"] - diag["upper_shifted"])
        )
        max_acc = max(max_acc, float(acc))
        min_e = min(min_e, float(matter))
        min_clearance = min(min_clearance, float(r - 1.6))
    if max(errs.values()) > 3e-10 or min_e < -1e-12 or min_clearance <= 0 or max_acc >= 0:
        raise ArithmeticError(errs)
    direct = {}
    if v > 0:
        # Independent expanding first-order equations, without the temporal Israel equation.
        def algebraic(t, z):
            r, sa, sb, matter = z
            mus, mp = bulk((sa, sb))
            chi = (mus[0] - mus[1]) / (4 * r**4 * (1 + matter))
            h = math.sqrt(matter * (2 + matter) + mus.mean() / r**4 + chi**2 - 1 / r**2)
            radial = np.array([1 + matter - chi, 1 + matter + chi])
            du = 1 / (r * (radial + h))
            return [r * h, *du, -float(mp @ du**2) / (4 * r**3) - 4 * h * matter]

        reduced = solve_ivp(
            algebraic,
            (0.0, stop),
            [a0, 0.0, 0.0, initial_e],
            method="DOP853",
            rtol=3e-12,
            atol=3e-14,
            max_step=width / 33,
        )
        if not reduced.success:
            raise ArithmeticError(reduced.message)
        discrepancy = abs(reduced.y[:, -1] - [a, u1, u2, e])
        direct = dict(zip(("a", "u_left", "u_right", "e"), map(float, discrepancy)))
        if max(direct.values()) > 1e-9:
            raise ArithmeticError(direct)
    comparison = None
    if comparison_time is not None:
        if comparison_time < stop:
            raise ValueError("Common-time readout must lie after every active stream")
        plateau = solve_ivp(
            full_rhs,
            (stop, comparison_time),
            full.y_events[0][0],
            method="DOP853",
            rtol=3e-12,
            atol=3e-14,
            max_step=width / 10,
        )
        if not plateau.success:
            raise ArithmeticError(plateau.message)
        rad, vel, sa, sb, density, _ = plateau.y[:, -1]
        _, _, acceleration, check = physical(rad, vel, (sa, sb), density)
        if max(check["angular_residual"], check["friedmann_residual"]) > 3e-10:
            raise ArithmeticError("Common-time plateau violates the spatial constraint")
        comparison = dict(
            tau=comparison_time,
            a=float(rad),
            H=float(vel / rad),
            e=float(density),
            acceleration_over_radius=float(acceleration),
        )
    _, _, acc, final = physical(a, v, (u1, u2), e)
    return dict(
        initial_H=initial_h,
        initial_e=initial_e,
        initial_mass_per_side=m0,
        width_per_side=width,
        left_mass_loss_fraction=fractions[0],
        right_mass_loss_fraction=fractions[1],
        end_proper_time=float(stop),
        end_radius=float(a),
        end_H=float(v / a),
        end_e=float(e),
        added_comoving_radiation=float(deposit),
        end_masses=list(map(float, final["masses"])),
        end_acceleration_over_radius=float(acc),
        end_radial_curvatures=list(map(float, final["radial"])),
        end_harmonic_mean=float(final["hm"]),
        end_asymmetry_squared=float(final["asym"] ** 2),
        maximum_sampled_acceleration_over_radius=max_acc,
        minimum_source_clearance=min_clearance,
        turning_proper_times=list(map(float, full.t_events[1])),
        common_time_plateau=comparison,
        errors=errs,
        independent_expanding_solution_errors=direct,
    )


def calculate(profile):
    fast = math.sqrt(2) / 4
    rows = {
        label: solve(h, frac, w, profile, comparison_time=0.6 if h == fast else None)
        for label, h, frac, w in [
            ("two_sided_half", fast, (0.5, 0.5), 0.2),
            ("one_sided_same_total", fast, (1.0, 0.0), 0.2),
            ("two_sided_unequal_same_total", fast, (0.75, 0.25), 0.2),
            ("one_sided_exchange", fast, (0.0, 1.0), 0.2),
            ("one_sided_slow_full", 0.04, (1.0, 0.0), 0.025),
            ("one_sided_at_rest", 0.0, (0.5, 0.0), 0.05),
        ]
    }
    rows["one_sided_dense"] = solve(1.0, (0.5, 0.0), 0.05, profile, initial_e=0.4)
    for field in ["end_H", "end_e", "end_radius", "end_acceleration_over_radius"]:
        if abs(rows["one_sided_same_total"][field] - rows["one_sided_exchange"][field]) > 2e-11:
            raise ArithmeticError("Swapping retained sides changed the result")
    return rows
