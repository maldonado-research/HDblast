#!/usr/bin/env python3
"""Exact-rational evaluation of the characteristic observation majorant.

This is a conditional arithmetic tool, not a physical HDBLAST waveform solver.
Input values must already be rigorous rational enclosures for the true problem.
Run without arguments for meaningful exact self-checks, or pass one JSON file.
All JSON numbers must be integers or rational strings, never floating point.
"""

from __future__ import annotations

from fractions import Fraction as Q
import json
from pathlib import Path
import sys


def rat(value):
    if isinstance(value, float):
        raise TypeError("Use exact integers or rational strings, not float inputs")
    return Q(value)


def matvec(matrix, vector):
    return [sum((a * b for a, b in zip(row, vector)), Q(0)) for row in matrix]


def exp_action_upper(matrix, vector, tolerance=Q(1, 10**30)):
    """Enclose exp(matrix) vector for a nonnegative rational matrix/vector.

    Positive Taylor partial sums are lower bounds. After term n, the omitted
    first term has norm d; the future term norms have ratio at most r/(n+2),
    where r is the infinity matrix norm. Hence d/(1-r/(n+2)) bounds the entire
    omitted vector componentwise. Everything, including this tail, is rational.
    """
    size = len(vector)
    if len(matrix) != size or any(len(row) != size for row in matrix):
        raise ValueError("Square matrix with matching vector required")
    if tolerance <= 0 or any(x < 0 for x in vector):
        raise ValueError("Positive tolerance and nonnegative vector required")
    if any(x < 0 for row in matrix for x in row):
        raise ValueError("Nonnegative matrix required")
    radius = max(sum(row, Q(0)) for row in matrix)
    term = list(vector)
    total = list(vector)
    order = 0
    while order < 10000:
        next_term = [x / (order + 1) for x in matvec(matrix, term)]
        ratio = radius / (order + 2)
        if ratio < 1:
            tail = max(next_term) / (1 - ratio)
            if tail <= tolerance:
                return total, [x + tail for x in total], order, tail
        total = [x + y for x, y in zip(total, next_term)]
        term = next_term
        order += 1
    raise ArithmeticError("Taylor enclosure exceeded 10000 terms")


def certify(inputs):
    keys = (
        "T", "M", "B", "mu_lower", "A0", "b_plus", "b_minus", "V0",
        "F", "G", "gamma_abs", "gamma_error", "u_boundary_abs",
    )
    values = {key: rat(inputs[key]) for key in keys}
    if any(value < 0 for value in values.values()) or values["mu_lower"] <= 0:
        raise ValueError("Nonnegative bounds and a strictly positive mu_lower required")
    v = values
    a = 1 / v["mu_lower"]
    j0 = max(v["b_plus"], v["b_minus"])
    # Augment y'=Ay+c with a last component constantly equal to one.
    generator = [
        [Q(0), Q(1), Q(1), j0],
        [v["M"], Q(0), Q(0), v["F"]],
        [a * v["B"], a, Q(0), a * (v["b_plus"] + v["G"])],
        [Q(0), Q(0), Q(0), Q(0)],
    ]
    scaled = [[v["T"] * x for x in row] for row in generator]
    initial = [v["A0"], Q(0), v["V0"], Q(1)]
    tolerance = rat(inputs.get("tolerance", "1/1000000000000000000000000000000"))
    lower, upper, order, tail = exp_action_upper(scaled, initial, tolerance)
    e_upper, i_upper, v_upper = upper[:3]
    trace_upper = (
        v["b_plus"] + i_upper + v_upper + v["gamma_abs"] * e_upper
        + v["gamma_error"] * v["u_boundary_abs"]
    )
    eps2_upper = v["T"] * trace_upper**2
    return {
        "kind": "conditional exact-rational characteristic majorant",
        "physical_waveform_supplied": False,
        "majorant_lower_is_not_actual_error_lower_bound": True,
        "inputs": {key: str(value) for key, value in values.items()},
        "majorant_lower_at_T": dict(zip(("E", "I", "V"), map(str, lower[:3]))),
        "majorant_upper_at_T": dict(zip(("E", "I", "V"), map(str, upper[:3]))),
        "trace_error_sup_upper": str(trace_upper),
        "trace_error_L2_squared_upper": str(eps2_upper),
        "taylor_order": order,
        "taylor_component_tail_upper": str(tail),
    }


def self_checks():
    base = {
        "T": 1, "M": 0, "B": 0, "mu_lower": 1,
        "A0": 0, "b_plus": 0, "b_minus": 0, "V0": 0,
        "F": 0, "G": 0, "gamma_abs": 0, "gamma_error": 0,
        "u_boundary_abs": 0,
    }
    zero = certify(base)
    assert zero["trace_error_L2_squared_upper"] == "0"

    # Exact nonphysical operator control: e=(t-z)_+^3 is C^2, has zero bulk
    # residual and zero initial data. With mu=1, beta=0, boundary g=3t^2+6t,
    # so G=9 on [0,1]. Its exact observation is -3t^2 and norm squared 9/5.
    # The majorant is I=0, V=9t, E=9t^2/2; the nilpotent Taylor sum is exact.
    boundary = certify(dict(base, G=9))
    assert boundary["majorant_lower_at_T"] == {"E": "9/2", "I": "0", "V": "9"}
    assert boundary["majorant_upper_at_T"] == boundary["majorant_lower_at_T"]
    assert Q(boundary["trace_error_L2_squared_upper"]) == 81
    assert Q(9, 5) <= Q(boundary["trace_error_L2_squared_upper"])
    omitted_boundary = certify(base)
    assert Q(omitted_boundary["trace_error_L2_squared_upper"]) < Q(9, 5)

    # Polynomial forced control e=t^2 has f=g=2. The exact positive majorant
    # is I=2t, V=t^2+2t, E=t^3/3+2t^2, again captured without roundoff.
    forced = certify(dict(base, F=2, G=2))
    assert forced["majorant_lower_at_T"] == {"E": "7/3", "I": "2", "V": "3"}
    assert forced["majorant_upper_at_T"] == forced["majorant_lower_at_T"]
    # Dropping its boundary forcing would give V(1)=1 < |e_t(1,0)|=2.
    bad_forced = certify(dict(base, F=2, G=0))
    assert Q(bad_forced["majorant_upper_at_T"]["V"]) == 1 < 2

    nonnilpotent = certify(dict(
        base, T="1/2", M=2, B=3, mu_lower="2/3", A0="1/11",
        b_plus="1/13", b_minus="1/17", V0="1/19", F="1/7", G="1/9",
        gamma_abs="2/5", gamma_error="1/23", u_boundary_abs="4/3",
    ))
    assert 0 < Q(nonnilpotent["taylor_component_tail_upper"]) <= Q(1, 10**30)
    assert all(
        Q(nonnilpotent["majorant_lower_at_T"][key])
        <= Q(nonnilpotent["majorant_upper_at_T"][key])
        for key in ("E", "I", "V")
    )
    try:
        certify(dict(base, F=0.1))
    except TypeError:
        pass
    else:
        raise AssertionError("Uncertified floating input was accepted")
    return {
        "exact_checks": 11,
        "wrong_controls_rejected": [
            "omit nonzero boundary residual with zero bulk and initial residuals",
            "omit boundary residual in forced polynomial velocity estimate",
            "accept floating point as an exact input enclosure",
        ],
        "physical_waveform_supplied": False,
        "zero_control": zero,
        "boundary_control": boundary,
        "polynomial_control": forced,
        "positive_matrix_control": nonnilpotent,
    }


if __name__ == "__main__":
    result = certify(json.loads(Path(sys.argv[1]).read_text())) if len(sys.argv) > 1 else self_checks()
    print(json.dumps(result, indent=2))
