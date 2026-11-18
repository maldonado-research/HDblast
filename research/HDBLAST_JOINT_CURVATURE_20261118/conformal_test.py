"""Certify joint-channel separation and the sign of the second-order cut test."""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json

SHAPE_HASH = "247d0a9cd51aada5dac103732278429981aeade6d1f6943cb3ee179d629eec1a"


def certify(source):
    raw = source.read_bytes()
    if hashlib.sha256(raw).hexdigest() != SHAPE_HASH:
        raise ValueError("The source differs from the archived waveform certificate.")
    shape = json.loads(raw)
    ell = F(1, 10**6)
    limits = {}
    for key, required in (("waveform", F(10246, 100000)), ("enlarged", F(9192, 100000))):
        first, second = shape[key]["mass_squared_intervals"]
        x = F(first[0]) + 4*ell**2
        y = F(second[1]) + 4*ell**2
        if not x > y > 0:
            raise ValueError("The transformed mass intervals must be positive and disjoint.")
        residual = (x-y)/(x+y)
        if residual <= required:
            raise ValueError("Rounded mismatch bound.")
        limits[key] = {"ratio_error_lower_bound": str(residual), "decimal": float(residual)}
    # Uniform measure on S^3: r=n1²+n2² is uniform on [0,1]; mean cos²(2θ)=1/2.
    avg_w_y2 = F(1,2) * (3*F(1,3) - 4*F(1,4))
    avg_w_grad_y2 = 4 * (3*F(1,2) - 4*F(1,3) - avg_w_y2)
    if avg_w_y2 != 0 or avg_w_grad_y2 != F(2,3):
        raise ValueError("Angular projection.")
    return {"archived_shape_sha256": SHAPE_HASH, "collinear_conformal_test": limits,
            "sphere_averages": {"W_times_Y_squared": str(avg_w_y2),
                                "W_times_gradient_norm_squared": str(avg_w_grad_y2)},
            "fixed_tower_strict_lower_bound_rho_b_squared_times_cut_coefficient_over_kappa4":
                str(6/ell**3)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shape", type=Path, nargs="?",
                        default=Path(__file__).resolve().parent / "inputs" / "SHAPE.json")
    print(json.dumps(certify(parser.parse_args().shape), indent=2))
