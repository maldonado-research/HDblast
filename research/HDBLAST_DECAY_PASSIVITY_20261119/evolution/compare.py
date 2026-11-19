"""Compare full evolving-target reports on the same prescribed walls."""
import argparse
import json
from pathlib import Path

ARMS = ("quiet", "plus", "minus", "far")
MODES = ("physical", "sheared", "sheared_prepared")


def max_difference(left, right, field, channels):
    return max(abs(left[a][m]["end"][field][c] - right[a][m]["end"][field][c])
               for a in ARMS for m in MODES for c in channels)


def compare(left, right):
    x, y = left["arms"], right["arms"]
    return dict(
        nodes=[left["order"], right["order"]],
        max_abs_end_density=max_difference(x, y, "redshifted_density", ("H_loss", "count", "drive")),
        max_abs_contrast_percentage_points=max(
            abs(x[a]["contrasts"][m][v] - y[a]["contrasts"][m][v])
            for a in ARMS for m in MODES for v in ("H_loss", "count")),
        max_abs_min_signed_a4_difference=max(
            abs(x[a][m]["min_signed_a4_radiation"]["radiation_a4"]["H_loss"]
                - y[a][m]["min_signed_a4_radiation"]["radiation_a4"]["H_loss"])
            for a in ARMS for m in MODES))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("reports", type=Path, nargs="+")
    reports = [json.loads(path.read_text()) for path in p.parse_args().reports]
    by_order = {report["order"]: report for report in reports}
    if len(by_order) != len(reports) or 48 not in by_order:
        raise ValueError("Use unique orders and include the 48-node reference")
    runs = sorted(reports, key=lambda r: r["order"])
    arms = by_order[48]["arms"]
    quiet = arms["quiet"]
    checks = [compare(l, r) for l, r in zip(runs, runs[1:])]
    if max(row["max_abs_end_density"] for row in checks) > 2e-6:
        raise ValueError("Quadrature endpoint comparison failed")
    seed = {}
    for m in MODES:
        prep, minimum = quiet[m]["initial"], quiet[m]["min_signed_a4_radiation"]
        seed[m] = max(0.0, -minimum["radiation_a4"]["H_loss"])/prep["scale_factor"]**4
    summary = dict(quadrature=checks, fixed_quiet_wall_conditional_seed_initial_density=seed,
                   max_abs_native_control_density_disagreement=max(
                       abs(arms[a]["physical"]["end"]["redshifted_density"]["count"]
                           - arms[a]["archived_physical_density"]) for a in ARMS))
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
