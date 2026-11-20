"""Open-wall hot expansion and finite charge-work diagnostics (standard library)."""
import argparse
import bisect
import hashlib
import json
import math
from pathlib import Path


def interpolate(l, r, t):
    f = (t - l["t"]) / (r["t"] - l["t"])
    return {k: l[k] + f * (r[k] - l[k]) for k in l}


def at(trace, t):
    i = bisect.bisect_left([r["t"] for r in trace], t)
    if i < len(trace) and trace[i]["t"] == t:
        return trace[i]
    if i == 0 or i == len(trace):
        raise ValueError("Requested event outside trace")
    return interpolate(trace[i - 1], trace[i], t)


def root(fn, lo, hi):
    fl = fn(lo)
    if fl == 0:
        return lo
    if fl * fn(hi) > 0:
        raise ValueError("Root has no sign change")
    for _ in range(80):
        mid = (lo + hi) / 2
        if fl * fn(mid) <= 0:
            hi = mid
        else:
            lo, fl = mid, fn(mid)
    return (lo + hi) / 2


def power_integral(b, delta):
    if delta == 1.0:
        return math.log(b)
    return math.expm1((1.0 - delta) * math.log(b)) / (1.0 - delta)


def geometry_margin(s, b, delta, comoving_rest_cap):
    return (1 + (s["a"] * s["H"]) ** 2
            + 2 * comoving_rest_cap / s["a"] * power_integral(b, delta)
            - b ** (2 - delta))


def growth_cap(s, delta, comoving_rest_cap):
    if delta >= 2:
        return None
    hi = 2.
    while geometry_margin(s, hi, delta, comoving_rest_cap) > 0:
        hi *= 2
    return root(lambda b: geometry_margin(s, b, delta, comoving_rest_cap), 1., hi)


def moment(s):
    return s["a"] ** 2 * (1 + (s["a"] * s["H"]) ** 2)


def describe(s):
    return dict(t=s["t"], a=s["a"], aH=s["a"] * s["H"], e=s["e"], K=s["K"],
                electric_share=s["X"] ** 2 / (s["B"] ** 2 + s["X"] ** 2),
                decayed_fraction=s["f"], internal_share=s["intrinsic_random"] / s["e"],
                C=moment(s))


def span(trace, start, stop):
    return [start, *(s for s in trace if start["t"] < s["t"] < stop["t"]), stop]


def trap(rows, f):
    return sum((r["t"] - l["t"]) * (f(l) + f(r)) / 2.
               for l, r in zip(rows, rows[1:]))


def stress_diagnostic(rows, fraction):
    result = {}
    for tol in (.01, .05, .10):
        start = None
        best = None
        for r in rows:
            share = r["intrinsic_random"] / r["e"]
            if share < fraction or abs(3 * r["p_total"] / r["e"] - 1) > tol:
                start = None
                continue
            if start is None:
                start = r
            growth = r["a"] / start["a"]
            if best is None or growth > best["radius_factor"]:
                best = dict(start=start["t"], stop=r["t"], radius_factor=growth)
        result[f"{tol:.2f}"] = best
    return dict(
        max_gas_relative_radiation_pressure_error=max(abs(3 * r["pgas"] / r["gas"] - 1)
                                                      for r in rows),
        min_total_relative_radiation_pressure_error=min(abs(3 * r["p_total"] / r["e"] - 1)
                                                        for r in rows),
        max_total_relative_radiation_pressure_error=max(abs(3 * r["p_total"] / r["e"] - 1)
                                                        for r in rows),
        longest_hot_saved_span_by_total_pressure_tolerance=result,
    )


def curvature_integrals(rows):
    """Signed contributions to C from resolved fields, integrated in ln(a)."""
    terms = {
        "stress_trace": lambda s: s["e"] - 3 * s["p_total"],
        "quadratic_density": lambda s: -s["e"] * (s["e"] + 3 * s["p_total"]),
        "boundary_phase": lambda s: s["B"] ** 2 - s["X"] ** 2,
    }
    sums = {}
    for key, term in terms.items():
        sums[key] = sum(math.log(r["a"] / l["a"]) *
                        (l["a"] ** 4 * term(l) + r["a"] ** 4 * term(r))
                        for l, r in zip(rows, rows[1:]))
    diff = moment(rows[-1]) - moment(rows[0])
    return dict(sums, C_change=diff, implied_nonnegative_radial_term=sum(sums.values()) - diff)


def finite_work(rows, rest_comoving, charge):
    start, end = rows[0], rows[-1]
    base = start["a"] ** 3 * (start["gas"] - start["rest"])
    cap = rest_comoving * abs(charge) * trap(rows, lambda r: abs(r["X"]) * r["f"])
    # Exact accumulated birth kinetic energy and pressure integral in the original ODE.
    birth = end["decay_kinetic"] - start["decay_kinetic"]
    pressure_work = end["p_work"] - start["p_work"]
    actual_work = end["charge_work"] - start["charge_work"]
    return dict(start_comoving_gas_above_rest=base, birth_kinetic=birth,
                actual_work=actual_work, speed_limited_work_quadrature=cap,
                pressure_work=pressure_work,
                comoving_internal=end["a"] ** 3 * end["intrinsic_random"],
                upper_internal_using_actual_work=base + birth + actual_work - pressure_work,
                upper_internal_using_speed_limit=base + birth + cap - pressure_work)


def diagnose(trace, exit_time, hot_time, turn_time, rest_comoving, rmax, charge, fraction=0.9):
    # Interpolated horizons avoid including the first contracting sample.

    # Use saved physical states for energy audits. Interpolating primitive
    # and conservative variables independently across fast charge creation
    # can violate the original identities even when all saved states obey them.
    raw_post = [s for s in trace if exit_time <= s["t"] <= turn_time]
    start = raw_post[0]
    hot = next(s for s in raw_post if s["t"] >= hot_time
               and s["intrinsic_random"] >= fraction * s["e"])
    turn = at(trace, turn_time)
    sampled = [*raw_post, turn]
    after = [s for s in sampled if s["t"] >= hot["t"]]

    kmin = min(s["K"] for s in raw_post)
    if min(s["p_total"] for s in raw_post) < 0 or min(s["mass"] for s in raw_post) < 0:
        raise ValueError("Positive pressure and bulk mass required for this bound")
    if max(abs(s["V"]) for s in raw_post) > 1e-10:
        raise ValueError("Zero late scalar potential required for this bound")
    delta = 2 * (1 - fraction) / kmin
    rest_cap = rest_comoving * max(1., rmax)
    capacity = lambda s: geometry_margin(s, 2., delta, rest_cap)
    crossings = []
    for l, r in zip(sampled, sampled[1:]):
        if capacity(l) >= 0 and capacity(r) < 0:
            t = root(lambda t: capacity(interpolate(l, r, t)), l["t"], r["t"])
            crossings.append(describe(interpolate(l, r, t)))
    predicted_trace = lambda s: (4 * s["V"] - s["pi"] ** 2 + s["parent"]
                                + s["rest"] / math.hypot(1., s["u"])
                                / (1.5 * s["theta"] + math.hypot(1., 1.5 * s["theta"])))
    matter_trace_error = max(abs(s["e"] - 3 * s["p_total"] - predicted_trace(s))
                             for s in raw_post)
    geom_error = max(abs(moment(s) - s["a"] ** 4 * (2 * s["e"] + s["e"] ** 2)
                         - s["mass"]) / moment(s) for s in raw_post)
    work_rows = span(trace, start, hot)
    work = finite_work(work_rows, rest_comoving, charge)
    sparse = [work_rows[0], *work_rows[2:-1:2], work_rows[-1]]
    return dict(exit_event_time=exit_time, first_hot_event_time=hot_time,
                first_postexit_sample=describe(start), first_confirmed_hot_sample=describe(hot),
                interpolated_turn=describe(turn),
                min_postexit_sampled_K=kmin, min_postexit_sampled_mass=min(s["mass"] for s in raw_post),
                delta=delta, comoving_rest_cap=rest_cap,
                doubling_margin_at_first_hot=capacity(hot),
                min_initial_aH_for_doubling_at_hot_radius=math.sqrt(max(0.,
                    (hot["a"] * hot["H"]) ** 2 - capacity(hot))),
                subsequent_radius_factor_upper_if_hot=growth_cap(hot, delta, rest_cap),
                subsequent_radius_factor_to_turn=turn["a"] / hot["a"],
                last_doubling_headroom_crossings=crossings, finite_work_between_saved_samples=work,
                speed_cap_spacing_sensitivity=abs(work["speed_limited_work_quadrature"] -
                    finite_work(sparse, rest_comoving, charge)["speed_limited_work_quadrature"]),
                postexit_trace_identity_max_error=matter_trace_error,
                stress_diagnostic=stress_diagnostic([s for s in raw_post if s["t"] >= hot["t"]], fraction),
                curvature_integrals_hot_to_turn=curvature_integrals(after),
                curvature_integrals_half_sampling=curvature_integrals([after[0], *after[2:-1:2], after[-1]]),
                sampled_geometric_constraint_max_relative_error=geom_error,
                minimum_sampled_interior_share_after_hot=min(s["intrinsic_random"] / s["e"] for s in after))


def analyze_archive(path, case, branch):
    raw = path.read_bytes()
    archive = json.loads(raw)
    r = archive["runs"][case]
    cfg = archive["parameters"][case]
    trace = [dict(zip(archive["trace_columns"], full)) for full in r["trace"]]
    if branch == "two_stage":
        d0, ri, rf, charge = (cfg[k] for k in ("parent_initial_density", "initial_mass_ratio",
                                              "final_mass_ratio", "charge_to_mass"))
    elif branch == "dual":
        # These masses are positional constants in axion/growth.py.
        d0, ri, rf, charge = cfg["D"], .95, 1.05, cfg["Q"]
    else:
        d0, ri, rf, charge = (cfg[k] for k in ("D", "R_low", "R_high", "Q"))
    result = diagnose(trace, r["exit_complete_tau"], r["first_ninety_intrinsic_time"],
                      r["maximum_tau"], trace[0]["a"] ** 3 * d0 / ri, rf, charge)
    if "first_ninety_intrinsic_aH" in r:
        a = r["first_ninety_intrinsic_radius"]
        event = {"a": a, "H": r["first_ninety_intrinsic_aH"] / a}
        result["integrator_hot_onset"] = {
            "a": a, "aH": r["first_ninety_intrinsic_aH"],
            "growth_to_turn": r["maximum_radius"] / a,
            "growth_upper_bound_if_hot": growth_cap(event, result["delta"], result["comoving_rest_cap"]),
        }
    if branch == "previous":
        result["archived_input_sha256"] = hashlib.sha256(raw).hexdigest()
    result["run"] = case
    return result


def dual_preparation(path):
    raw = path.read_bytes()
    data = json.loads(raw)
    states = [dict(zip(data["columns"], v)) for v in data["rows"]]
    start, end = states[0], states[-1]

    def b(s):
        return s["a"] ** 2 * s["B"]

    def z(s):
        return s["a"] ** 2 * (s["K"] * s["X"] + s["theta"] * s["B"])

    observed = end["X"] / end["B"]
    impulsive = (z(start) / b(start) - end["theta"]) / end["K"]
    return dict(points=len(states), duration=end["t"] - start["t"],
                elapsed_conformal=trap(states, lambda s: 1 / s["a"]),
                comoving_B_fractional_change=b(end) / b(start) - 1,
                displacement_change_over_final_odd_term=(z(end) - z(start)) / (end["theta"] * b(end)),
                predicted_X_over_B=impulsive, measured_X_over_B=observed,
                predicted_electric_share=impulsive ** 2 / (1 + impulsive ** 2),
                measured_electric_share=observed ** 2 / (1 + observed ** 2))


def main():
    base = Path(__file__).resolve().parents[1]
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prior", type=Path, default=base / "prior/DECAY.json")
    p.add_argument("--two-stage", type=Path, default=base / "exit/EXIT.json")
    p.add_argument("--dual", type=Path, default=base / "axion/GROWTH.json")
    p.add_argument("--handoff", type=Path, default=base / "axion/HANDOFF.json")
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    result = {key: analyze_archive(path, "primary", key) for key, path in
              (("previous", args.prior), ("two_stage", args.two_stage), ("dual", args.dual))}
    result["dual_handoff"] = dual_preparation(args.handoff)
    output = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(output)
    else:
        print(output, end="")


if __name__ == "__main__":
    main()
