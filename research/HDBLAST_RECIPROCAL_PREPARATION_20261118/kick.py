"""Check the uniform modal bounds and transfer the published late-mean interval."""
import argparse
import hashlib
import json
from fractions import Fraction as F
from pathlib import Path

SHAPE_SHA256 = "247d0a9cd51aada5dac103732278429981aeade6d1f6943cb3ee179d629eec1a"


def require(condition, message):
    if not condition:
        raise ValueError(message)


def certificate(shape_path):
    raw = shape_path.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == SHAPE_SHA256, "registered window report")
    shape = json.loads(raw)
    late_lower, late_upper = map(F, shape["scaled_window_means"][0])
    require(0 < late_lower <= late_upper, "signed window endpoints")

    ell = F(1, 10**6)
    horizon = F(25, 16)
    late_end = F(11, 8)
    width = F(1, 8)
    # B=min(1, omega^-2). These constants apply to the two normalized exposures.
    positive_floor = F(19, 60)
    flat_ratios = [
        late_end / positive_floor,
        F(2) / positive_floor,
        F(4),
        F(10) / (2 * (1 - F(1, 4))),
        (2 / width) / (2 * (1 - F(1, 10))),
    ]
    flat_ratio = max(flat_ratios)
    require(flat_ratio == F(80, 9), "flat ratio over all five frequency bands")

    potential = F(35, 4) * ell**4 * horizon**2
    require(potential * horizon**2 < 1, "retarded integral bound")
    wave_error = potential * horizon**3 / (1 - potential * horizon**2)
    dilution = F(3, 4) * ell**2 * horizon**2
    dilution_derivative = F(3, 2) * ell**2 * horizon
    # Low frequencies use the bound on sin(omega*x)/omega. High frequencies
    # use the total variation of the two zero-extended exposure weights.
    exposure_error = wave_error + max(horizon*dilution, 16*dilution + dilution_derivative)
    require(exposure_error < F(1, 10**9), "uniform error per B")
    ratio = F(9)
    require((ratio - flat_ratio)*positive_floor > (1 + ratio)*exposure_error,
            "positive early kernel and curved-space ratio")
    min_residual = late_lower / (1 + ratio)
    noise = F(1, 2000)
    missing_sign_mass = (late_lower - (1 + ratio)*noise) / (2 * ratio)
    require(min_residual > F(4013, 10000), "residual bound")
    require(missing_sign_mass > F(2226, 10000), "both signs required")

    # A positive measure of impulse times in [0, 1/4] models an extended source.
    # Relative oscillation frequency is now omega*h with h=1-start_time.
    shortest = F(3, 4)
    pulse_ratios = [
        flat_ratios[0] / shortest**3,
        flat_ratios[1] / shortest**2,
        flat_ratios[2] / shortest**2,
        flat_ratios[3] / shortest**2,
        flat_ratios[4] / shortest,
    ]
    pulse_early_floor = shortest**3 * positive_floor
    pulse_ratio = F(12)
    require(max(pulse_ratios) == F(320, 27), "launch-time band limits")
    require((pulse_ratio - max(pulse_ratios))*pulse_early_floor >
            (1 + pulse_ratio)*2*exposure_error, "finite-duration drive")
    require((16 + 2*exposure_error) / (positive_floor - 2*exposure_error) < 51,
            "coefficient for arbitrary source-to-front gap")
    pulse_residual = late_lower / (1 + pulse_ratio)
    require(pulse_residual > F(3087, 10000), "finite-duration exposure error")

    gap_cubed = 51 * noise / (late_lower - noise)
    require(gap_cubed < F(1853, 10000)**3, "latest permitted source time before signal")

    return {
        "shape_sha256": SHAPE_SHA256,
        "early_triangle": ["0", "1"],
        "late_box": ["5/4", "11/8"],
        "physical_late_mean": [str(late_lower), str(late_upper)],
        "flat_ratio_limits": list(map(str, flat_ratios)),
        "flat_positive_floor": str(positive_floor),
        "curvature_and_dilution_error_per_min_1_inverse_omega_squared": str(exposure_error),
        "collocated_inequality_constant": str(ratio),
        "lower_bound_on_max_exposure_error": str(min_residual),
        "relative_to_target_late_mean": "1/10",
        "max_exposure_error_for_signed_extension": str(noise),
        "lower_bound_on_each_sign_early_contribution": str(missing_sign_mass),
        "one_signed_drive": {
            "source_support": ["0", "1/4"],
            "flat_ratio_limits": list(map(str, pulse_ratios)),
            "early_floor": str(pulse_early_floor),
            "curved_ratio": str(pulse_ratio),
            "lower_bound_on_max_exposure_error": str(pulse_residual),
            "relative_to_target_late_mean": "1/13",
            "latest_drive_time_strict_lower_at_error_1_over_2000": "8147/10000",
            "source_to_front_gap_cubed_upper_at_that_error": str(gap_cubed),
            "ratio_bound_when_drive_ends_by_1_minus_delta": "51/delta^3",
        },
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shape", type=Path, nargs="?", default=Path(__file__).resolve().parent / "inputs" / "SHAPE.json",
                        help="SHAPE.json from the 18 November comparison edition")
    args = parser.parse_args()
    print(json.dumps(certificate(args.shape), indent=2))
