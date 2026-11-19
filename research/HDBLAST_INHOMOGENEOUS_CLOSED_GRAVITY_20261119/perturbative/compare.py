"""Compare full Einstein–scalar evolution with the quadratic volume response and exact bound."""
import argparse
import json
from pathlib import Path


OBSERVABLES = {
    "radius_fractional_coefficient": "volume_radius_over_background_minus_one",
    "H_volume_difference_coefficient": "H_volume_minus_background",
    "acceleration_over_radius_difference_coefficient": "acceleration_difference_from_background",
    "fractional_growth_since_initial_coefficient": "fractional_growth_since_initial",
}


def observations(run):
    probe = run["at_probe"]
    initial = run["samples"]["0"]
    r0, rx = (d["volume_radius_over_background_minus_one"] for d in (initial, probe))
    return {**probe, "fractional_growth_since_initial": (rx-r0)/(1+r0)}


def compare(volume, coupled):
    x = coupled["probe_x"]
    result = {}
    for key, branch in volume["branches"].items():
        target = next(row for row in branch["times"] if row["x"] == x)
        source = coupled["branches"][key]["runs"]
        records = {}
        for A in (.001, .0005):
            positive, negative = (observations(source[str(a)]) for a in (A, -A))
            fields = {}
            for coefficient, observable in OBSERVABLES.items():
                even = (positive[observable] + negative[observable])/(2*A**2)
                odd = (positive[observable] - negative[observable])/(2*A**2)
                expected = target[coefficient]
                fields[coefficient] = dict(coupled_sign_average_per_A2=even,
                                           coupled_half_difference_per_A2=odd,
                                           perturbative=expected, relative_difference=even/expected-1)
            records[str(A)] = fields
        refinement = {}
        for name in OBSERVABLES:
            big, small = (records[str(A)][name]["relative_difference"] for A in (.001, .0005))
            ratio = big/small
            if abs(big) > .002 or abs(ratio/4-1) > .04:
                raise ArithmeticError(f"Quadratic amplitude check failed for {key} {name}")
            refinement[name] = ratio
        # Restrict this check to the five standard launches and their saved snapshots.
        snapshots = [s for run in source.values()
                     for s in [*run["samples"].values(), run["at_stop"]]]
        max_local = max(s["maximum_local_expansion_difference"] for s in snapshots)
        max_averaged = max(s["H_volume_minus_background"] for s in snapshots)
        peak = max(s["maximum_field"] for s in snapshots)
        cap = branch["potential_concavity"]["symmetric_field_radius"]
        if max(max_local, max_averaged) > 2e-10 or peak > cap + 2e-10:
            raise ArithmeticError("Focusing or field-window check failed")
        result[key] = dict(amplitude_pairs_at_probe=records, ratio_of_relative_errors_upon_halving=refinement,
                           max_sampled_local_expansion_difference=max_local,
                           max_sampled_H_volume_minus_background=max_averaged,
                           max_sampled_absolute_field=peak,
                           last_snapshot_time=max(s["x"] for s in snapshots))
    return dict(probe_x=x, branches=result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--volume", required=True, type=Path)
    parser.add_argument("--coupled", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output = compare(json.loads(args.volume.read_text()), json.loads(args.coupled.read_text()))
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
