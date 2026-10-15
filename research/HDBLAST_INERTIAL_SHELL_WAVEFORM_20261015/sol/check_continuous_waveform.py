#!/usr/bin/env python3
"""Independent CSV coverage, whole-cell squares, physical scaling and C audit.

This checker imports no scientific producer. All mathematical decisions use
integers and Fraction. A PASS audits the exported enclosures; the separately
reviewed analytic PDE and residual assumptions remain inputs.
"""
from __future__ import annotations

import argparse
import copy
import csv
from fractions import Fraction as F
from hashlib import sha256
import io
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
HEADERS = ["left_grid", "right_grid", "H0_lo", "H0_hi", "error_upper", "H_lo", "H_hi", "ell_S_lo", "ell_S_hi"]
INPUT_PINS = {
    "inputs/previous/archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json": "8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8",
    "inputs/previous/archive/PHYSICAL_GAIN_ENCLOSURE.json": "2a2dbc398ac2799418edcddeb917f4dcecce9b9a8049ac73b54541e617b616a4",
    "archive/BACKGROUND_DERIVATIVE_BOUNDS.json": "fd6fd277176033703a4e75a2f66ab94894318aac49c83264ef979e9a27f0a02c",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def ceil_div(n, d):
    require(d > 0, "positive arithmetic denominator")
    return -((-n) // d)


def product(a, b):
    products = [x * y for x in a for y in b]
    return min(products), max(products)


def square(a):
    return (0 if a[0] <= 0 <= a[1] else min(x*x for x in a), max(x*x for x in a))


def encoded(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {k: encoded(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encoded(v) for v in value]
    return value


def parse_csv(raw):
    reader = csv.reader(io.StringIO(raw.decode("utf-8"), newline=""))
    require(next(reader) == HEADERS, "CSV schema")
    rows = []
    for fields in reader:
        require(len(fields) == len(HEADERS), "CSV row width")
        rows.append([int(v) for v in fields])
    return rows


def audit_rows(report, rows):
    require(report["status"] == "PASS", "producer report PASS")
    state, cell = report["state"], report["cell_data"]
    ell = F(state["ell"])
    require(ell == F(1, 10**6) and F(state["T"]) == 3*ell and F(state["Z"]) == 8*ell, "new wider-state scale")
    require(state["domain"] == [0, 3], "full observation domain")
    N, S = cell["time_grid_denominator"], cell["value_denominator"]
    require(N == 4096 and S == 2**112, "registered grid denominators")
    require(S == report["arithmetic"]["value_denominator"], "report denominator agreement")
    require(cell["rows"] == len(rows) == 4353, "registered CSV row count")
    require(cell["support_cells"] == 4096 and cell["tail_cells"] == 256, "registered component counts")
    G = tuple(map(F, report["physical_gain_interval"]))
    K = tuple(map(F, report["physical_norm_prefactor_interval"]))
    require(G[0] < G[1] < 0 < K[0] < K[1], "gain and prefactor signs")
    totals = {name: [F(0), F(0)] for name in ("H0_raw", "H_raw", "H0_outward", "H_outward")}
    tail_norm = [F(0), F(0)]
    tail_H = [None, None]
    tail_physical = [None, None]
    tail_count = 0
    previous = 0
    physical_quantum = F(1, S)
    for index, row in enumerate(rows):
        ln, rn, a, b, error, h, j, p, q = row
        require(ln == previous and rn > ln, "contiguous positive-width CSV coverage")
        previous = rn
        require(a <= b and h <= j and p <= q and error >= 0, "ordered cell enclosures")
        require(h == a-error and j == b+error, "residual inflation matches exported error")
        if index == 0:
            require(row == [0, N, 0, 0, 0, 0, 0, 0, 0], "exact whole prefront cell")
        elif index <= N:
            require(ln == N+index-1 and rn == ln+1, "complete support cell sequence")
        else:
            require(ln >= 2*N and rn-ln == 16, "complete postpulse cell sequence")
        width = F(rn-ln, N)
        for name, interval in (("H0", (a, b)), ("H", (h, j))):
            lo, hi = square(interval)
            totals[name+"_raw"][0] += width * F(lo, S*S)
            totals[name+"_raw"][1] += width * F(hi, S*S)
            totals[name+"_outward"][0] += width * F(lo//S, S)
            totals[name+"_outward"][1] += width * F(ceil_div(hi, S), S)
        # Reconstruct the exact physical interval before any fixed-grid
        # outward rounding. It must be contained in the exported interval.
        factor = (1-(ell*F(rn, N))**2, F(1))
        require(factor[0] > 0, "positive proper-time clock factor")
        physical = product(G, product(factor, (F(h, S), F(j, S))))
        require(F(p, S) <= physical[0] <= physical[1] <= F(q, S), "physical coefficient contains exact gain-clock scaling")
        # Each multiplication incurs < one fixed-grid unit; the nested
        # implementation should not introduce a macroscopic surplus.
        surplus = (physical[0]-F(p, S)) + (F(q, S)-physical[1])
        require(surplus < 100*physical_quantum, "physical outward rounding surplus")
        if ln >= 2*N:
            tail_count += 1
            require(0 < h <= j and p <= q < 0, "signed postpulse H and physical coefficient")
            tail_H[0] = h if tail_H[0] is None else min(tail_H[0], h)
            tail_H[1] = j if tail_H[1] is None else max(tail_H[1], j)
            tail_physical[0] = p if tail_physical[0] is None else min(tail_physical[0], p)
            tail_physical[1] = q if tail_physical[1] is None else max(tail_physical[1], q)
            lo, hi = square((h, j))
            tail_norm[0] += width * F(lo//S, S)
            tail_norm[1] += width * F(ceil_div(hi, S), S)
    require(previous == 3*N and tail_count == 256, "coverage reaches complete endpoint including tail")
    quadrature = report["quadrature"]
    require(tuple(totals["H0_outward"]) == tuple(map(F, quadrature["comparison_H0_integral_interval"])), "comparison whole-cell squared integral")
    require(tuple(totals["H_outward"]) == tuple(map(F, quadrature["actual_H_integral_interval"])), "actual whole-cell squared integral")
    require(totals["H_outward"][0] <= totals["H_raw"][0] <= totals["H_raw"][1] <= totals["H_outward"][1], "outward square sum contains exact raw-square sum")
    C = K[0]*totals["H_outward"][0]/ell, K[1]*totals["H_outward"][1]/ell
    require(C == tuple(map(F, quadrature["C_interval"])), "response coefficient uses K integral divided by ell")
    require(F(606, 100)*10**11 < C[0] < C[1] < F(744, 100)*10**11, "public strict C bounds 6.06e11 and 7.44e11")
    require(tail_norm[0] > 0, "positive tail contribution retained in squared norm")
    theorem_tail = tuple(map(F, report["tail_theorem"]["actual_H_interval"]))
    require(theorem_tail[0] <= F(tail_H[0], S) <= F(tail_H[1], S) <= theorem_tail[1], "all CSV tail cells inside uniform tail theorem")
    return {
        "coverage": {"domain": [0, 3], "rows": len(rows), "prefront_cells": 1, "pulse_cells": N, "postpulse_cells": tail_count, "no_gaps_or_overlaps": True},
        "squared_integrals": totals, "postpulse_squared_integral_outward": tail_norm,
        "actual_C_interval": C, "public_strict_C_bounds": [F(606)*10**9, F(744)*10**9],
        "postpulse_H_from_csv": [F(x, S) for x in tail_H],
        "postpulse_ell_S_from_csv": [F(x, S) for x in tail_physical],
        "physical_scaling": "Every CSV physical interval contains G*[1-(ell*right)^2,1]*H; G<0 and the clock factor is positive",
    }


def run(source_root):
    source_root = Path(source_root)
    report_path = source_root/"computed/INERTIAL_WAVEFORM_ENCLOSURE.json"
    csv_path = source_root/"computed/INERTIAL_WAVEFORM_CELLS.csv"
    report_raw, csv_raw = report_path.read_bytes(), csv_path.read_bytes()
    report = json.loads(report_raw)
    require(report["cell_data"]["sha256"] == sha256(csv_raw).hexdigest(), "CSV digest matches completed report")
    require(report["cell_data"]["bytes"] == len(csv_raw), "CSV byte count")
    source_pins = {}
    for path, expected in report["source_pins"].items():
        digest = sha256((source_root/path).read_bytes()).hexdigest()
        require(digest == expected, "current source matches completed report: "+path)
        source_pins[path] = digest
    for path, expected in INPUT_PINS.items():
        require(sha256((source_root/path).read_bytes()).hexdigest() == expected, "independently pinned mathematical input: "+path)
    gain = json.loads((source_root/"inputs/previous/archive/PHYSICAL_GAIN_ENCLOSURE.json").read_text())
    retained = gain["correlated_retained_eta_box"]
    require(report["physical_gain_interval"] == retained["G"]["outward_rational_at_12_decimal_places"], "same retained shell gain")
    require(report["physical_norm_prefactor_interval"] == retained["K"]["outward_rational_at_12_decimal_places"], "same retained norm prefactor")
    rows = parse_csv(csv_raw)
    result = audit_rows(report, rows)
    controls = []
    for name, mutation, expected in (
        ("coverage_gap", lambda r, x: x[100].__setitem__(0, x[100][0]+1), "contiguous"),
        ("omit_postpulse", lambda r, x: x.pop(), "row count"),
        ("alter_squared_integral", lambda r, x: r["quadrature"]["actual_H_integral_interval"].__setitem__(0, "0"), "actual whole-cell squared"),
        ("wrong_C_scaling", lambda r, x: r["quadrature"]["C_interval"].__setitem__(0, str(F(r["quadrature"]["C_interval"][0])*F(r["state"]["ell"]))), "response coefficient"),
        ("reverse_physical_tail_sign", lambda r, x: x[-1].__setitem__(slice(7, 9), [-x[-1][8], -x[-1][7]]), "physical coefficient"),
    ):
        bad_report, bad_rows = copy.deepcopy(report), copy.deepcopy(rows)
        mutation(bad_report, bad_rows)
        try:
            audit_rows(bad_report, bad_rows)
        except RuntimeError as exc:
            require(expected in str(exc), "fault failed at unrelated gate: "+name+": "+str(exc))
            controls.append({"name": name, "rejected_at": str(exc)})
        else:
            raise RuntimeError("fault accepted: "+name)
    return encoded({
        "status": "PASS_INDEPENDENT_CONTINUOUS_WAVEFORM_REVIEW",
        "python_optimization": sys.flags.optimize,
        "arithmetic": "Integer and Fraction only; no floating certification and no producer import",
        "inputs": {
            "report": {"path": "computed/INERTIAL_WAVEFORM_ENCLOSURE.json", "bytes": len(report_raw), "sha256": sha256(report_raw).hexdigest()},
            "csv": {"path": "computed/INERTIAL_WAVEFORM_CELLS.csv", "bytes": len(csv_raw), "sha256": sha256(csv_raw).hexdigest()},
            "mathematical_input_pins": INPUT_PINS, "source_pins": source_pins,
        },
        "review": result, "fault_controls": controls, "fault_control_count": len(controls),
        "scope": "Independent reconstruction of continuous CSV coverage, residual inflation, physical interval containment, full-domain squared integrals and C. Analytic background, PDE, residual and nonlinear-family hypotheses are not proved by this finite audit.",
        "own_source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    })


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=HERE.parent)
    args = parser.parse_args()
    try:
        print(json.dumps(run(args.source_root), indent=2, sort_keys=True))
    except (RuntimeError, ValueError, OSError, KeyError) as exc:
        print("FAIL: "+str(exc), file=sys.stderr)
        raise SystemExit(1)
