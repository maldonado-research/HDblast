"""Do the freely propagating daughters dominate the modeled expansion?"""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicHermiteSpline, PchipInterpolator

OPEN, CUT = 2.055, 2.10
LAMBDA = .01
TIMES = np.linspace(OPEN, CUT, 81)


def read(path):
    """Evaluate at equal shell proper time, on each run's own clock."""
    with np.load(path, allow_pickle=False) as run:
        s = run["s"]
        if not (s[0] < OPEN < CUT < s[-1]):
            raise ValueError(f"Late interval unavailable in {path}")
        field = lambda k: PchipInterpolator(s, run[k], extrapolate=False)(TIMES)
        a = np.exp(CubicHermiteSpline(s, run["A"], run["h"])(TIMES) - run["Aref"])
        return dict(a=a, h=field("h"), rho=field("rho"), rad=field("rad"), sigma=field("sigma"),
                    rb=float(run["rb"]), b=float(run["b"]), strength=float(run["lam"]))


def radiation_law(a, h):
    # For H=C/a², optimize the maximum fractional H error at these two
    # measured scale factors. The optimizer is the harmonic mean of a²H.
    q = a[[0,-1]]**2 * h[[0,-1]]
    if min(q) <= 0:
        raise ValueError("The two-point fit assumes an expanding interval")
    best = 2*q[0]*q[1]/(q[0]+q[1])
    error = max(abs(best/q - 1))
    if not np.isclose(error, abs(np.diff(q)[0])/sum(q), rtol=5e-13):
        raise ValueError("Two-point minimax identity failed")
    return dict(initial=float(q[0]), final=float(q[-1]), fractional_change=float(q[-1]/q[0]-1),
                minimum_two_point_fractional_error=float(error))


def measure(zero, plus, minus):
    h0, a0 = zero["h"], zero["a"]
    if (max(abs(zero["rad"])) > 1e-12 or (plus["strength"], minus["strength"]) != (LAMBDA, -LAMBDA)
            or np.any(plus["rad"] <= 0) or np.any(plus["rho"] <= 0)):
        raise ValueError("Radiation fraction or derivative is undefined")
    rad1 = (plus["rad"]-minus["rad"])/(2*LAMBDA)
    # The exact O(λ) product uses a₀: the daughter density vanishes at λ=0.
    redshift = a0**4 * rad1
    # Partial derivative of r_b² (σ+ρ_native)^2/36 holding every input
    # except the explicit physical daughter density fixed. ρ_native=b ρhat/r_b.
    mult = zero["b"]/zero["rb"]
    direct = zero["rb"]**2/18 * (zero["sigma"] + mult*zero["rho"]) * mult * rad1
    h2_response = (plus["h"]**2-minus["h"]**2)/(2*LAMBDA)
    if direct[-1] <= 0:
        raise ValueError("Direct daughter-density contribution must be positive in these runs")
    return dict(
        daughter_density_at_positive_strength=float(plus["rad"][-1]),
        all_modeled_matter_density_at_positive_strength=float(plus["rho"][-1]),
        daughter_fraction_of_modeled_matter=float(plus["rad"][-1]/plus["rho"][-1]),
        daughter_to_baseline_parent_per_lambda=float(rad1[-1]/zero["rho"][-1]),
        post_switch_max_fractional_drift_in_a4rho_daughter=float(max(abs(redshift/redshift[0]-1))),
        pure_radiation_fit=radiation_law(plus["a"],plus["h"]),
        h2_response_over_h02_per_lambda=float(h2_response[-1]/h0[-1]**2),
        h2_direct_daughter_over_h02_per_lambda=float(direct[-1]/h0[-1]**2),
        full_to_direct_h2_slope=float(h2_response[-1]/direct[-1]),
    )


def calculate(runs):
    result = dict(proper_interval=[OPEN,CUT], positive_lambda=LAMBDA, grids={})
    for grid in ("coarse","fine"):
        arms = result["grids"][grid] = {}
        for arm in ("quiet","plus","minus","far"):
            name = f"{grid}_{arm}_"
            zero = read(runs / (name+"zero.npz"))
            row = arms[arm] = dict(baseline_radiation_fit=radiation_law(zero["a"],zero["h"]), cases={})
            for cutoff in (60,120):
                for scheme in ("bare","static"):
                    prefix = runs/(name + f"{cutoff}_{scheme}_")
                    plus, minus = [read(Path(str(prefix)+side+"0.01.npz")) for side in ("p","m")]
                    row["cases"][f"{cutoff}_{scheme}"] = measure(zero, plus, minus)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.runs), indent=2, allow_nan=False)+"\n")
