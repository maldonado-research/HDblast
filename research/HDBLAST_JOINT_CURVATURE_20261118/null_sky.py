"""Positive angular weights and independent cellwise integration for the null Ricci signal."""
import csv
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GRID = 4096
DEN = 2**112
ELL = F(1, 10**6)
END = F(25, 16)
FRONT = F(1)
FLOOR = 1 - F(1, 10**11)
TRANSFER = F(3, 10**6)
WINDOWS = [(F(5,4), F(11,8)), (F(49,32), F(25,16))]
CSV_HASH = "2f76fd1ae2195b445582a29a352bd76658dea17c2233cbb9acc395b255172bf7"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def integrate_cells():
    path = ROOT / "inputs" / "FUTURE_CELLS.csv"
    require(hashlib.sha256(path.read_bytes()).hexdigest() == CSV_HASH, "waveform cells")
    # z'' = input for the auxiliary zero-potential transfer. Keep independent
    # lower and upper paths and integrate each constant-acceleration bound.
    z, speed = [F(0), F(0)], [F(0), F(0)]
    means = [[F(0), F(0)] for _ in WINDOWS]
    widths = [F(0) for _ in WINDOWS]
    peak = F(0)
    edge = FRONT
    with path.open() as stream:
        for row in csv.DictReader(stream):
            left = max(F(int(row["left_grid"]), GRID), FRONT)
            right = min(F(int(row["right_grid"]), GRID), END)
            if right <= left:
                continue
            require(edge == left, "cell gap")
            edge = right
            dt = right - left
            accel = [F(int(row["H_lo"]), DEN) - TRANSFER,
                     F(int(row["H_hi"]), DEN) + TRANSFER]
            peak = max(peak, *map(abs, accel))
            for i, (start, stop) in enumerate(WINDOWS):
                if left >= stop or right <= start:
                    continue
                require(start <= left < right <= stop, "unaligned window")
                widths[i] += dt
                for side in range(2):
                    means[i][side] += z[side]*dt + speed[side]*dt**2/2 + accel[side]*dt**3/6
            for side in range(2):
                z[side] += speed[side]*dt + accel[side]*dt**2/2
                speed[side] += accel[side]*dt
    require(edge == END, "waveform range")
    length = END - FRONT
    # The retarded integral of |z| has sup norm <= length² ||z|| / 2.
    theta = F(5, 4)*ELL**2*length**2
    bound = peak*length**2/(2*(1-theta))
    metric_error = (theta + 1-FLOOR)*bound
    null_error = (ELL*END)**2 * bound
    result = []
    for i,(start,stop) in enumerate(WINDOWS):
        require(widths[i] == stop-start, "exposure width")
        lo = means[i][0]/widths[i] - metric_error
        hi = means[i][1]/widths[i] + metric_error
        result.append([lo,hi])
    return result, null_error


def calculate():
    # Uniform orthonormal pair in R4. For a uniform S³ coordinate n:
    # E[n_i^4] = 3/(4*6), E[n_i^2 n_j^2] = 1/(4*6) for i != j.
    y2 = 2*(F(3, 24)-F(1, 24))
    yd = -y2/3   # Conditional second moment of v on the tangent two-sphere.
    zd = y2+yd/3
    multiplier = 4*zd
    require(multiplier == F(16, 27), "angular null contraction")
    metric, error = integrate_cells()
    registered = json.loads((ROOT/"METRIC.json").read_text())
    inherited_means = [list(map(F, pair)) for pair in registered["metric_means"]]
    require(metric == inherited_means, "independent retarded integration")
    # Output units are mu*rho_b^4/ell times R_kk first variation.
    shape_bytes = (ROOT / "inputs" / "SHAPE.json").read_bytes()
    require(hashlib.sha256(shape_bytes).hexdigest() == "247d0a9cd51aada5dac103732278429981aeade6d1f6943cb3ee179d629eec1a", "trace reference")
    scalar_means = json.loads(shape_bytes)["scaled_window_means"]
    rows = []
    for i,(lo,hi) in enumerate(metric):
        require(lo > error, "positive null response")
        negative = [-multiplier*(hi+error), -multiplier*(lo-error)]
        positive = [-negative[1], -negative[0]]
        denominators = list(map(F, scalar_means[i]))
        require(denominators[0]*denominators[1] > 0, "trace mean cannot vanish")
        ratios = [ELL**2*n/(6*d) for n in negative for d in denominators]
        lo_ratio, hi_ratio = min(ratios), max(ratios)
        rounded = [(-9677, -9602), (159167, 160653)][i]
        require(F(rounded[0], 10**19) < lo_ratio <= hi_ratio < F(rounded[1], 10**19),
                "physical curvature ratio")
        rows.append({"physical_plus_to_trace_ratio":list(map(str, (lo_ratio,hi_ratio))),
                     "window":list(map(str,WINDOWS[i])),
                     "mean_plus":list(map(str,negative)), "mean_minus":list(map(str,positive))})
    return {"csv_sha256":CSV_HASH, "angular_multiplier":str(multiplier),
            "units":"(mu*rho_b^4/ell) times the first variation of the null Ricci exposure",
            "null_factor_error":str(error), "full_sky_null_exposures":rows}


if __name__=="__main__":
    print(json.dumps(calculate(),indent=2))
