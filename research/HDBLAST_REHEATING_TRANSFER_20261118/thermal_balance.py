"""Reproduce weighted radiation gains and source integrals on the published grids."""
import csv
import hashlib
import io
import json
import re
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
COLUMNS = ["x", "weighted_radiation", "weighted_source"]
WINDOWS = {
    "b1_coarse": (1, 4),
    "b1_fine": (1, 4),
    "b2_early_s1": (1, 2),
    "b2_control_s1": (1, 2),
    "b2_early_s2": (1, 2),
    "b2_control_s2": (1, 2),
    "b2_s1": (10, 11),
    "b2_s2": (10, 11),
}
NUMBER = re.compile(r"[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")


def decimal_text(number):
    """Round an exact rational once for display, uniformly across Python platforms."""
    with localcontext() as context:
        context.prec = 16
        return str(Decimal(number.numerator) / Decimal(number.denominator))


def parse_csv(raw, name):
    """Accept a finite, nonnegative radiation density and source on the declared grid."""
    reader = csv.reader(io.StringIO(raw.decode("ascii"), newline=""), strict=True)
    if next(reader, None) != COLUMNS:
        raise ValueError(f"{name}: unexpected CSV header")
    points = []
    for number, row in enumerate(reader, start=2):
        if len(row) != 3 or not all(NUMBER.fullmatch(value) for value in row):
            raise ValueError(f"{name}:{number}: require three finite decimal numbers")
        x, radiation, source = map(Fraction, row)
        if points and x <= points[-1][0]:
            raise ValueError(f"{name}:{number}: clock must increase strictly")
        if radiation < 0 or source < 0:
            raise ValueError(f"{name}:{number}: negative radiation or heating")
        if name.startswith("b2_control_") and (radiation or source):
            raise ValueError(f"{name}:{number}: nonzero source-off control")
        points.append((x, radiation, source))
    if len(points) < 2 or (points[0][0], points[-1][0]) != WINDOWS[name]:
        raise ValueError(f"{name}: incomplete clock interval")
    return points


def balance(points):
    radiation_start = points[0][1]
    gain = points[-1][1] - radiation_start
    integral = Fraction()
    largest_defect = Fraction()
    for (x0, _, source0), (x1, radiation1, source1) in zip(points, points[1:]):
        integral += (x1 - x0) * (source0 + source1) / 2
        defect = integral - (radiation1 - radiation_start)
        largest_defect = max(largest_defect, abs(defect))
    difference = integral - gain
    return {
        "rows": len(points),
        "clock_interval": [decimal_text(points[0][0]), decimal_text(points[-1][0])],
        "initial_weighted_radiation": decimal_text(radiation_start),
        "final_weighted_radiation": decimal_text(points[-1][1]),
        "radiation_gain": decimal_text(gain),
        "integrated_source": decimal_text(integral),
        "integrated_source_minus_gain": decimal_text(difference),
        "relative_to_gain_ppm": decimal_text(10**6 * difference / gain) if gain else None,
        "largest_absolute_cumulative_defect": decimal_text(largest_defect),
    }


def compute():
    information = json.loads((HERE / "PROVENANCE.json").read_text(encoding="utf-8"))
    if set(information["series"]) != set(WINDOWS):
        raise ValueError("source register and expected series differ")
    results = {}
    for name in WINDOWS:
        raw = (HERE / "data" / (name + ".csv")).read_bytes()
        metadata = information["series"][name]
        if hashlib.sha256(raw).hexdigest() != metadata["csv_sha256"]:
            raise ValueError(f"{name}: CSV checksum mismatch")
        points = parse_csv(raw, name)
        if len(points) != metadata["rows"]:
            raise ValueError(f"{name}: row count differs from source register")
        results[name] = balance(points)
    return {
        "method": "Exact composite trapezoid on decimal CSV nodes. Reported numbers rounded to 16 significant digits.",
        "residual_sign": "integrated weighted source minus endpoint weighted radiation gain",
        "series": results,
    }


if __name__ == "__main__":
    print(json.dumps(compute(), indent=2, allow_nan=False))
