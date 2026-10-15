#!/usr/bin/env python3
"""Independent rational reconstruction of the exported waveform envelopes.

No waveform implementation is imported. The continuous comparator ranges are
separately assessed by analytic/source review; this audits the full exported
inventory and its arithmetic consequences.
"""
from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import csv
import json


def compute(root):
    report_path = root / "computed/WAVEFORM_ENCLOSURE.json"
    report = json.loads(report_path.read_text())
    checks = {}

    def check(name, condition):
        checks[name] = bool(condition)
        if not condition:
            raise ValueError(name)

    def ceil(q):
        return -((-q.numerator) // q.denominator)

    scale = 2**96

    def mul(a, b):
        products = [x*y for x in a for y in b]
        return min(products) // scale, ceil(F(max(products), scale))

    check("declared_value_denominator", report["cell_data"]["value_denominator"] == scale)
    step = report["cell_data"]["time_grid_denominator"]
    check("support_cells_match_grid", step == report["quadrature"]["support_cells"] == 4096)
    path = root / "computed/WAVEFORM_CELLS.csv"
    raw = path.read_bytes()
    check("csv_hash_matches_receipt", sha256(raw).hexdigest() == report["cell_data"]["sha256"])
    check("csv_bytes_match_receipt", len(raw) == report["cell_data"]["bytes"])
    with path.open(newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        rows = [[int(x) for x in row] for row in reader]
    check("csv_columns", header == ["left_grid", "right_grid", "comparison_lo", "comparison_hi",
                                    "H_lo", "H_hi", "ell_S_lo", "ell_S_hi"])
    check("csv_rows_match_receipt", len(rows) == report["cell_data"]["rows"] == 4098)
    check("full_closed_domain_coverage", rows[0][0] == 0 and rows[-1][1] == 3*step and
          all(a[1] == b[0] for a, b in zip(rows, rows[1:])))
    check("nonempty_ordered_cells", all(row[0] < row[1] for row in rows))
    check("one_support_cell_per_grid_step", all(r[1]-r[0] == 1 for r in rows[1:-1]))
    eta = F(report["error_constants"]["published_eta"])
    tail = F(report["error_constants"]["published_tail"])
    zeta = F(report["error_constants"]["norm_error_upper"])
    ell = F(report["error_constants"]["ell"])
    eta_units, tail_units = ceil(eta*scale), ceil(tail*scale)
    check("exact_prefront_zero", rows[0] == [0, step, 0, 0, 0, 0, 0, 0])
    check("tail_is_bounded_not_discarded", rows[-1][:6] ==
          [2*step, 3*step, 0, 0, -tail_units, tail_units])
    check("all_interval_orders", all(r[2] <= r[3] and r[4] <= r[5] and r[6] <= r[7]
                                     for r in rows))
    check("uniform_support_envelope", all(r[4] == r[2]-eta_units and r[5] == r[3]+eta_units
                                           for r in rows[1:-1]))
    gain = tuple(F(x) for x in report["physical_gain_interval"])
    gain_units = gain[0]*scale//1, ceil(gain[1]*scale)
    check("negative_physical_gain", gain[0] < gain[1] < 0)
    physical_equal = True
    for row in rows:
        xmax = F(row[1], step)
        factor = (1-(ell*xmax)**2)*scale//1, scale
        physical = mul(gain_units, mul(factor, (row[4], row[5])))
        physical_equal = physical_equal and tuple(row[6:8]) == physical
    check("all_physical_cells_exactly_reconstructed", physical_equal)
    low_sum = high_sum = 0
    for row in rows[1:-1]:
        lo, hi = row[2:4]
        low = 0 if lo <= 0 <= hi else min(lo*lo, hi*hi)
        low_sum += low//scale
        high_sum += ceil(F(max(lo*lo, hi*hi), scale))
    comparison = tuple(F(x) for x in report["quadrature"]["comparison_integral_interval"])
    check("quadrature_reconstructed_from_all_continuous_cells", comparison ==
          (F(low_sum, step*scale), F(high_sum, step*scale)))
    check("support_aware_norm_error_gate", F(report["error_constants"]["eta_upper"])**2 +
          F(report["error_constants"]["eta_tail_upper"])**2 < zeta*zeta)

    def sqrt_bracket(q):
        sc = 2**128
        k = isqrt(q.numerator*sc*sc//q.denominator)
        return F(k, sc), F(k+1, sc)

    real_lo = max(F(0), sqrt_bracket(comparison[0])[0]-zeta)**2
    real_hi = (sqrt_bracket(comparison[1])[1]+zeta)**2
    actual = tuple(F(x) for x in report["quadrature"]["actual_H_integral_interval"])
    check("reported_integral_outward_vs_128bit_recheck", actual[0] <= real_lo and real_hi <= actual[1])
    prefactor = tuple(F(x) for x in report["physical_norm_prefactor_interval"])
    coefficient = tuple(F(x) for x in report["quadrature"]["C_interval"])
    check("coefficient_jacobian_and_prefactor", coefficient ==
          (prefactor[0]*actual[0]/ell, prefactor[1]*actual[1]/ell))
    check("planned_outward_interval", F(606, 100)*10**15 < coefficient[0] <
          coefficient[1] < F(744, 100)*10**15)
    for rel, pin in report["source_pins"].items():
        check("source_pin_"+rel, sha256((root/rel).read_bytes()).hexdigest() == pin)
    check("gain_input_pin", sha256((root/report["gain_input"]["path"]).read_bytes()).hexdigest() ==
          report["gain_input"]["sha256"])
    return {
        "status": "PASS_INDEPENDENT_CSV_AND_NORM_AUDIT", "checks": checks,
        "finite_check_groups": len(checks),
        "input_pins": {"computed/WAVEFORM_ENCLOSURE.json": sha256(report_path.read_bytes()).hexdigest(),
                       "computed/WAVEFORM_CELLS.csv": sha256(raw).hexdigest(), **report["source_pins"]},
        "scope": "Independent exact-rational reconstruction of all 4098 exported cell envelopes, physical gain multiplication and response-norm arithmetic. Continuous comparator range validity is separately assessed by source-code and analytic review; this is not formal PDE certification."
    }


def main():
    parser = ArgumentParser()
    parser.add_argument("--source-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = compute(args.source_root)
    encoded = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(encoded)
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()
