"""Constrained 5D Einstein-Maxwell triplet coupled to a moving timelike spherical sheet."""

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

INNER = 1.65
INITIAL_RADIUS = 2.0
INITIAL_MASS = 4.0
DURATION = 0.083
END = 0.50
K0 = 0.055
CAP = 20.0
SCALE = 0.017


def coefficient(psi, variable, cap):
    y = math.tanh(psi / SCALE)
    derivative = math.log(cap) / SCALE if variable else 0.0
    k = K0 * math.exp((math.log(cap) * y) if variable else 0.0)
    return k, k * derivative * (1 - y * y)


def pulse(t, strength):
    s = t / DURATION
    return strength * (math.sin(math.pi * s) ** 4) if 0 < s < 1 else 0.0


def derivative(values, h, upstream_left):
    # Third-order biased differences in the interior with second-order closures.
    out = np.empty_like(values)
    if upstream_left:
        out[2:] = (3 * values[2:] - 4 * values[1:-1] + values[:-2]) / (2 * h)
        out[2:-1] = (values[:-3] - 6 * values[1:-2] + 3 * values[2:-1] + 2 * values[3:]) / (6 * h)
        out[1] = (values[2] - values[0]) / (2 * h)
        out[0] = (-3 * values[0] + 4 * values[1] - values[2]) / (2 * h)
    else:
        out[:-2] = (-3 * values[:-2] + 4 * values[1:-1] - values[2:]) / (2 * h)
        out[1:-2] = (-2 * values[:-3] - 3 * values[1:-2] + 6 * values[2:-1] - values[3:]) / (6 * h)
        out[-2] = (values[-1] - values[-3]) / (2 * h)
        out[-1] = (3 * values[-1] - 4 * values[-2] + values[-3]) / (2 * h)
    return out


def potential(psi, ramp):
    low, high, height = ramp
    if psi <= low or height == 0.0:
        return 0.0, 0.0
    if psi >= high:
        return height, 0.0
    q = (psi - low) / (high - low)
    return height * q**3 * (10 + q * (-15 + 6 * q)), height * 30 * q * q * (1 - q) ** 2 / (
        high - low
    )


def run(
    cells,
    strength,
    variable,
    end=END,
    ramp=(0.035, 0.06, 0.0),
    stop_on_turn=False,
    initial_h=0.0,
    cap=CAP,
):
    # x=(r-r_source)/(a-r_source). Left boundary supplies the growing-r
    # characteristic. At the wall the sheet equation supplies the other.
    initial_mass = INITIAL_MASS + INITIAL_RADIUS**4 * initial_h**2
    if cells < 8 or cells % 4 or end <= 0 or cap < 1 or strength < 0:
        raise ValueError(
            "The arm mesh requires a positive multiple of four and valid boundary data"
        )
    x = np.linspace(0.0, 1.0, cells + 1)
    dx = 1.0 / cells

    def evaluate(t, state, diagnostics=False):
        a, v, ms, field, w, psi, pi, mflux, sheet_flux, sheet_work, scalar_work, scalar_pdilat = (
            state[:12]
        )
        out = np.empty(cells + 1)
        inward = np.empty(cells + 1)
        tube_clock = state[12]
        out[0] = pulse(tube_clock, strength)
        out[1:] = state[24 : 24 + cells]
        inward[:-1] = state[24 + cells :]
        k, kp = coefficient(psi, variable, cap)
        electric = w * w / a**2
        magnetic = 4 * field**2 / a**4
        e_g = k * (electric + magnetic) / 2
        pot, grad = potential(psi, ramp)
        e_s = 0.5 * pi * pi + pot
        p_s = 0.5 * pi * pi - pot
        e = e_g + e_s
        p = e_g / 3 + p_s
        # Using the spatial Israel jump to close the incoming characteristic
        # is algebraic. The independently integrated radial metric below
        # supplies its diagnostic residual and the temporal Israel source.
        xwall = a * (1 + e)
        incoming_boundary = (2 * w - (xwall - v) * out[-1]) / (xwall + v)
        inward[-1] = incoming_boundary
        P = (out + inward) * 0.5
        D = (inward - out) * 0.5
        r = INNER + x * (a - INNER)
        dr = (a - INNER) * dx

        def integ(z):
            # Cumulative Simpson quadrature on each endpoint pair and its first half.
            result = np.zeros_like(z)
            left, middle, right = z[:-2:2], z[1::2], z[2::2]
            result[2::2] = np.cumsum((dr / 3) * (left + 4 * middle + right))
            result[1::2] = result[:-2:2] + (dr / 12) * (5 * left + 8 * middle - right)
            return result

        primitive = integ(D)
        A = field + primitive - primitive[-1]
        logarithmic_slope = 2 * (D * D + P * P) / r
        log_factor = integ(logarithmic_slope)
        factor = np.exp(log_factor)
        # m' + (log N)' m = 2r(1+r²)(D²+P²) + 8A²/r.
        mass = (
            ms + integ(factor * (2 * r * (1 + r * r) * (D * D + P * P) + 8 * A * A / r))
        ) / factor
        f = 1 + r * r - mass / r**2
        X = math.sqrt(f[-1] + v * v)
        lapse = (X / f[-1]) * np.exp(log_factor - log_factor[-1])
        c = lapse * f
        # Cancellation removes the wave-gradient squares from (Nf)'.
        c_r = lapse * (2 * r + 2 * mass / r**3 - 8 * A * A / r**3)
        source = c * D / r - 4 * lapse * A / r**2
        u_rhs = -(c - v * x) * derivative(out, dr, True) - c_r * out + source
        d_rhs = (c + v * x) * derivative(inward, dr, False) + c_r * inward + source
        nA = w - (xwall - v) * out[-1]
        flux_in = -2 * w * nA / a**2
        normal_pressure = (w * w + nA * nA) / (2 * a * a) - 2 * field * field / a**4
        bulk_mass_change = 4 * a**3 * (X * w * nA / a**2 - v * normal_pressure)
        mass_inner_change = 4 * INNER * lapse[0] * f[0] ** 2 * P[0] * D[0]
        # Kinetic-only canonical scalar. No potential and no damping.
        scalar_force = 0.5 * kp * (electric - magnetic)
        accel = a * ((1 + e) * (1 - 2 * e - 3 * p) - 1 - mass[-1] / a**4 - 2 * normal_pressure)
        wall_rhs = [
            v,
            accel,
            mass_inner_change,
            w,
            -(v / a + kp * pi / k) * w - 4 * field / a**2 - 2 * nA / k,
            pi,
            -3 * v * pi / a + scalar_force - grad,
            bulk_mass_change,
            a**3 * flux_in,
            3 * a * a * v * p,
            a**3 * pi * scalar_force,
            3 * a * a * v * p_s,
            lapse[0] * math.sqrt(f[0]),
            0.5 * INNER * lapse[0] * f[0] ** 1.5 * out[0] ** 2,
            0.5 * INNER * lapse[0] * f[0] ** 1.5 * inward[0] ** 2,
            0.5 * a * (w - nA) ** 2,
            0.5 * a * (w + nA) ** 2,
        ]
        mid = np.array([cells // 4, cells // 2, 3 * cells // 4])
        mr = 2 * r * f * (P * P + D * D) + 8 * A * A / r
        mt = 4 * r * lapse * f * f * P * D
        middle_mass_rhs = (mt + x * v * mr)[mid]
        middle_field_rhs = (c * P + x * v * D)[mid]
        if not diagnostics:
            return np.r_[
                wall_rhs,
                middle_mass_rhs,
                middle_field_rhs,
                -(a**3) * pi * kp * magnetic / 2,
                u_rhs[1:],
                d_rhs[:-1],
            ]
        return dict(
            t=float(t),
            a=float(a),
            v=float(v),
            H=float(v / a),
            psi=float(psi),
            pi=float(pi),
            k=k,
            scalar_density=float(e_s),
            surface_em_density=float(e_g),
            energy=float(e),
            pressure=float(p),
            potential=float(pot),
            wall_field=float(field),
            wall_field_dot=float(w),
            wall_electric_squared=float(electric),
            wall_magnetic_squared=float(magnetic),
            sheet_proper_flux=float(sheet_flux),
            sheet_expansion_work=float(sheet_work),
            scalar_proper_work=float(scalar_work),
            tube_clock=float(tube_clock),
            source_input=float(state[13]),
            source_return=float(state[14]),
            wall_input=float(state[15]),
            wall_return=float(state[16]),
            scalar_magnetic_work=float(state[23]),
            scalar_electric_work=float(scalar_work - state[23]),
            scalar_expansion_work=float(scalar_pdilat),
            inner_mass=float(ms),
            wall_mass=float(mass[-1]),
            middle_mass_residual=float(max(abs(mass[mid] - initial_mass - state[17:20]))),
            middle_A_residual=float(max(abs(A[mid] - state[20:23]))),
            wall_trace_residual=float((c * P + v * x * D)[-1] - w),
            mass_ledger=float(mflux),
            spatial_jump=float(X / a - 1 - e),
            wall_mass_residual=float(mass[-1] - initial_mass - mflux),
            sheet_balance=float(a**3 * e + sheet_work - sheet_flux),
            scalar_balance=float(a**3 * e_s + scalar_pdilat - scalar_work),
            minimum_f=float(min(f)),
            minimum_radius_speed=float(min(c - abs(v * x))),
            maximum_radius_speed=float(max(c + abs(v * x))),
            Tnn=float(normal_pressure),
            wall_acceleration_over_radius=float(accel / a),
            inward_at_source=float(inward[0]),
            source_outgoing=float(out[0]),
            f_source=float(f[0]),
            N_source=float(lapse[0]),
            total_packet_energy_mass=float(mass[-1] - ms),
        )

    initial = np.zeros(24 + 2 * cells)
    initial[:3] = INITIAL_RADIUS, INITIAL_RADIUS * initial_h, initial_mass
    # Fixed fourth-order Runge-Kutta steps with a local characteristic CFL check.
    state = initial
    t = 0.0
    contracted = False
    rows = [evaluate(t, state, True)]
    for target in np.linspace(0.0, end, 101)[1:]:
        while t < target - 1e-15:
            max_speed = rows[-1]["maximum_radius_speed"]
            dt = min(target - t, 0.75 * (state[0] - INNER) / (cells * max_speed))
            k1 = evaluate(t, state)
            k2 = evaluate(t + dt / 2, state + dt * k1 / 2)
            k3 = evaluate(t + dt / 2, state + dt * k2 / 2)
            k4 = evaluate(t + dt, state + dt * k3)
            state = state + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
            t += dt
            contracted = contracted or state[1] < -1e-4
            if state[0] < INNER + 0.025 or (stop_on_turn and contracted and state[1] > 0):
                rows.append(evaluate(t, state, True))
                return rows
        rows.append(evaluate(t, state, True))
    return rows


def summarize(rows, cells, strength, varying, ramp, initial_h, cap):
    errors = {
        key: float(max(abs(row[key]) for row in rows))
        for key in (
            "spatial_jump",
            "wall_mass_residual",
            "middle_mass_residual",
            "middle_A_residual",
            "wall_trace_residual",
            "sheet_balance",
            "scalar_balance",
        )
    }
    end = rows[-1]
    result = dict(
        cells=cells,
        source_amplitude=strength,
        variable_coefficient=varying,
        coefficient_cap_ratio=cap,
        initial_H=initial_h,
        initial_mass=4 + 16 * initial_h**2,
        ramp=dict(start=ramp[0], end=ramp[1], height=ramp[2]),
        final=end,
        errors=errors,
        minimum_source_mass=min(row["inner_mass"] for row in rows),
        minimum_f=min(row["minimum_f"] for row in rows),
        minimum_characteristic_speed=min(row["minimum_radius_speed"] for row in rows),
        minimum_wall_radius=min(row["a"] for row in rows),
        maximum_acceleration_over_radius=max(row["wall_acceleration_over_radius"] for row in rows),
        maximum_potential=max(row["potential"] for row in rows),
    )
    diff = abs(end["wall_input"] - end["wall_return"] - end["sheet_proper_flux"])
    if (
        diff > 2e-8
        or min(
            result[k] for k in ("minimum_source_mass", "minimum_f", "minimum_characteristic_speed")
        )
        <= 0
    ):
        raise ArithmeticError("Invalid boundary flux budget or untrapped annulus")
    result["errors"]["incoming_minus_return"] = float(diff)
    result["scalar_work_over_wall_incident"] = (
        end["scalar_proper_work"] / end["wall_input"] if end["wall_input"] > 0 else None
    )
    events = []
    for left, right in zip(rows, rows[1:]):
        if left["v"] * right["v"] >= 0:
            continue
        dt = right["t"] - left["t"]

        def hermite(s, yl, yr, vl, vr):
            return (
                (1 - 3 * s * s + 2 * s**3) * yl
                + (3 * s * s - 2 * s**3) * yr
                + dt * (s * (1 - s) ** 2 * vl - s * s * (1 - s) * vr)
            )

        def velocity(s):
            return hermite(
                s,
                left["v"],
                right["v"],
                left["a"] * left["wall_acceleration_over_radius"],
                right["a"] * right["wall_acceleration_over_radius"],
            )

        s = brentq(velocity, 0, 1)
        events.append(
            dict(
                kind="minimum" if left["v"] < 0 else "maximum",
                time=left["t"] + s * dt,
                radius=hermite(s, left["a"], right["a"], left["v"], right["v"]),
                bracket=[left["t"], right["t"]],
            )
        )
    result["turning_points"] = events
    return result


def calculate():
    cases = {}
    zero = (0.035, 0.065, 0.0)
    rise = (0.035, 0.065, 0.10)

    def case(
        label, cells, strength, variable, end=0.5, ramp=zero, h=0.0, cap=20.0, save_series=False
    ):
        rows = run(cells, strength, variable, end=end, ramp=ramp, initial_h=h, cap=cap)
        data = summarize(rows, cells, strength, variable, ramp, h, cap)
        if save_series:
            fields = (
                "t",
                "a",
                "v",
                "H",
                "wall_mass",
                "inner_mass",
                "psi",
                "pi",
                "energy",
                "pressure",
                "scalar_density",
                "potential",
                "surface_em_density",
                "wall_acceleration_over_radius",
            )
            data["samples"] = [{key: row[key] for key in fields} for row in rows]
        cases[label] = data

    for n in (80, 160, 320):
        case(f"loading_{n}", n, 4.5, True, save_series=n == 320)
    case("constant_K", 160, 4.5, False)
    case("smaller_pulse", 160, 2.3, True)
    case("vacuum", 80, 0.0, True)
    for n in (80, 160, 240, 320, 640):
        case(f"rising_potential_{n}", n, 4.8, True, 2.2, rise, 0.12, 40.0, n == 320)
    case("rising_no_pulse", 80, 0.0, True, 2.2, rise, 0.12, 40.0)
    case("rising_no_potential", 160, 4.8, True, 2.2, zero, 0.12, 40.0)
    case("rising_constant_K", 160, 4.8, False, 2.2, rise, 0.12, 40.0)
    case("rising_at_rest", 160, 4.6, True, 1.3, rise, 0.0, 40.0)
    for n in (48, 160, 320):
        case(f"marginal_{n}", n, 4.8, True, 2.4, (0.035, 0.065, 0.09), 0.12)
    vacuum = cases["vacuum"]["final"]
    exact = math.sqrt(4 - vacuum["t"] ** 2)
    if abs(vacuum["a"] - exact) > 2e-10 or abs(vacuum["v"] + vacuum["t"] / exact) > 2e-10:
        raise ArithmeticError("Vacuum analytic solution differs")
    for label in ("vacuum", "constant_K", "rising_constant_K"):
        if cases[label]["final"]["scalar_density"] != 0:
            raise ArithmeticError("Scalar acquired energy without a force")
    if (
        cases["rising_no_pulse"]["maximum_potential"] != 0
        or cases["loading_320"]["maximum_potential"] != 0
    ):
        raise ArithmeticError("Unexpected potential energy")
    for branch in ("loading", "rising_potential"):
        ratios = {}
        for quantity in (
            "spatial_jump",
            "wall_mass_residual",
            "middle_mass_residual",
            "middle_A_residual",
        ):
            errors = [cases[f"{branch}_{n}"]["errors"][quantity] for n in (80, 160, 320)]
            ratios[quantity] = [errors[0] / errors[1], errors[1] / errors[2]]
            if min(ratios[quantity]) < 3.0:
                raise ArithmeticError(("Refinement failed", branch, quantity, errors))
        cases[branch + "_320"]["constraint_reduction_factors"] = ratios
    for n in (160, 240, 320, 640):
        model = cases[f"rising_potential_{n}"]
        if [event["kind"] for event in model["turning_points"]] != ["maximum", "minimum"] or model[
            "final"
        ]["v"] <= 0:
            raise ArithmeticError("Refined expanding branch failed to turn twice")
    for label in (
        "rising_no_pulse",
        "rising_no_potential",
        "rising_constant_K",
        "rising_at_rest",
        "marginal_160",
        "marginal_320",
    ):
        if (
            any(event["kind"] == "minimum" for event in cases[label]["turning_points"])
            or cases[label]["final"]["v"] >= 0
        ):
            raise ArithmeticError(
                "The control unexpectedly expanded before reaching the inner cutoff"
            )
    for label in ("loading_320", "rising_potential_320"):
        result = cases[label]
        if (
            not 0 < result["scalar_work_over_wall_incident"] < 1
            or result["final"]["scalar_magnetic_work"] >= 0
            or result["errors"]["spatial_jump"] > 3.5e-5
            or result["errors"]["middle_mass_residual"] > 7e-4
        ):
            raise ArithmeticError("Scalar work or refined gravitational convergence failed")
    return dict(
        units="AdS length and tuned wall tension = 1, kappa_5 squared = 6",
        source_radius=INNER,
        initial_radius=INITIAL_RADIUS,
        duration_in_source_proper_time=DURATION,
        coefficient_at_origin=K0,
        coefficient_field_scale=SCALE,
        annulus_stop_clearance=0.025,
        cases=cases,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(calculate(), indent=2, allow_nan=False) + "\n")
