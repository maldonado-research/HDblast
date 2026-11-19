"""Check identities and report refinements from the saved causal evolution."""
import argparse
import json
from pathlib import Path


def compare(left, right, point):
    errors = {}
    for state, first in left.items():
        for scheme, case in first.items():
            one = case["checkpoints"][point]
            two = right[state][scheme]["checkpoints"][point]
            for group in ("energy", "work", "transfer"):
                for k, v in one[group].items():
                    index = f"{group}.{k}"
                    errors[index] = max(errors.get(index, 0.), abs(v - two[group][k]))
    return errors


def check(record):
    runs = record["runs"]
    pairs = (("fine_wide_120", "fine_wide_120_N64_T2048"),
             ("fine_wide_120", "fine_wide_120_N96_T1024"),
             ("fine_wide_120", "fine_wide_120_N96_T2048"),
             ("fine_wide_240_N128_T2048", "fine_wide_240_N160_T4096"),
             ("fine_wide_480_N192_T4096", "fine_wide_480_N240_T8192"))
    errors = {}
    sector_max = {}
    primary_max = {}
    largest_identity_disagreement = 0.
    for key, states in runs.items():
        for cases in states.values():
            for mode, case in cases.items():
                for k, v in case["max_abs_residual"].items():
                    sector_max[k] = max(v, sector_max.get(k, 0.))
                    if key in ("fine_wide_120_N96_T2048", "fine_wide_480_N240_T8192"):
                        primary_max[k] = max(v, primary_max.get(k, 0.))
                common = cases["bare"]["checkpoints"]
                for label, row in case["checkpoints"].items():
                    energy, work, residual = (row[k] for k in ("energy", "work", "residual"))
                    point = f"{key}/{mode}/{label}"
                    if label in ("fork", "switch_on") and any(
                            v != 0 for group in (energy, work) for v in group.values()):
                        raise ValueError(f"{point}: nonzero correction before switch")
                    if label in ("switch_off", "cut") and (
                            energy["interaction"] != 0 or energy["counterterm"] != 0):
                        raise ValueError(f"{point}: interaction remaining after switch-off")
                    if energy["daughter"] < -1e-12:
                        raise ValueError(f"{point}: negative daughter energy")
                    for k in ("daughter", "interaction"):
                        if energy[k] != common[label]["energy"][k]:
                            raise ValueError(f"{point}: subtraction changed {k} at order y²")
                    if mode == "bare" and energy["counterterm"] != 0:
                        raise ValueError(f"{point}: unexpected bare counterterm")
                    if abs(sum(energy.values()) - sum(work.values()) - residual["total"]) >= 2e-11:
                        raise ValueError(f"{point}: inconsistency in recorded total residual")
                    sector_sum = sum(v for k, v in residual.items() if k != "total")
                    largest_identity_disagreement = max(largest_identity_disagreement,
                                                        abs(sector_sum - residual["total"]))
    for left, right in pairs:
        errors[left + "__vs__" + right] = {s: compare(runs[left], runs[right], s)
                                            for s in ("peak", "cut")}
    max_refinement = max(v for pair in errors.values() for point in pair.values()
                         for v in point.values())
    if max(sector_max.values()) > 5e-6 or max_refinement > 1e-3 or largest_identity_disagreement > 1e-10:
        raise ValueError((sector_max, max_refinement, largest_identity_disagreement))
    return dict(max_sampled_absolute_balance_residuals=sector_max,
                max_sampled_residuals_refined_fine_wide_120_and_480=primary_max,
                largest_sector_sum_vs_total_residual_roundoff=largest_identity_disagreement,
                component_absolute_refinement_differences=errors)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("result", type=Path)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    report = json.dumps(check(json.loads(args.result.read_bytes())), indent=2, allow_nan=False)
    if args.output:
        args.output.write_text(report + "\n")
    else:
        print(report)
