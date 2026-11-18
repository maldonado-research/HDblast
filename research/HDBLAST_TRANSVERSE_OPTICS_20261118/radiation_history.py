"""Recover the closed FLRW radiation and vacuum terms from a proper-time volume radius."""
import json
from fractions import Fraction


def reconstruct(radius, velocity, acceleration, jerk):
    """Use c=1, curvature k=+1 and derivatives in a common proper-time unit."""
    radius, velocity, acceleration, jerk = map(Fraction, (radius, velocity, acceleration, jerk))
    if radius <= 0:
        raise ValueError("A regular positive radius is required")
    expansion = velocity / radius
    friedmann = expansion**2 + 1 / radius**2
    relative_acceleration = acceleration / radius
    friedmann_rate = 2 * expansion * (relative_acceleration - friedmann)
    acceleration_rate = jerk / radius - relative_acceleration * expansion
    vacuum = (friedmann + relative_acceleration) / 2
    vacuum_rate = (friedmann_rate + acceleration_rate) / 2
    radiation = radius**4 * (friedmann - relative_acceleration) / 2
    radiation_rate = (4 * expansion * radiation
                      + radius**4 * (friedmann_rate - acceleration_rate) / 2)
    if radiation_rate != -radius**4 * vacuum_rate:
        raise ArithmeticError("Inconsistent proper-time derivative of the inferred components")
    return {
        "radiation_fraction": radiation / (radius**4 * friedmann),
        "Lambda_over_3": vacuum,
        "eight_pi_G_rho_r_a4_over_3": radiation,
        "derivative_Lambda_over_3": vacuum_rate,
        "derivative_eight_pi_G_rho_r_a4_over_3": radiation_rate,
    }


def benchmarks():
    # Integers or rational strings specify (a, da/dt, d²a/dt², d³a/dt³).
    # Expected outputs are in the function's order above.
    cases = [
        ("de_Sitter_minimum", [2, 0, "1/2", 0], [0, "1/4", 0, 0, 0]),
        ("de_Sitter_expanding", ["5/2", "3/4", "5/8", "3/16"], [0, "1/4", 0, 0, 0]),
        ("closed_radiation_maximum", [5, 0, "-1/5", 0], [1, 0, 25, 0, 0]),
        ("closed_radiation_contracting", [4, "-3/4", "-25/64", "-225/1024"], [1, 0, 25, 0, 0]),
        ("radiation_and_vacuum", [2, 1, "1/4", "7/8"], ["3/8", "5/16", 3, 0, 0]),
        ("linear_radius", [2, 1, 0, 0], ["1/2", "1/4", 4, "-1/4", 4]),
    ]
    records = {}
    for name, derivatives, expected in cases:
        result = reconstruct(*derivatives)
        if list(result.values()) != list(map(Fraction, expected)):
            raise ArithmeticError(name)
        records[name] = {"proper_time_jet": list(map(str, derivatives)),
                         "inferred": {key: str(value) for key, value in result.items()}}
    return records


if __name__ == "__main__":
    print(json.dumps(benchmarks(), indent=2))
