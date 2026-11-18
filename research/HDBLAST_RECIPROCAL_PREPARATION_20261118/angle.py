"""Check the signed preparation limits against the archived physical first mean."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

DIGEST = "247d0a9cd51aada5dac103732278429981aeade6d1f6943cb3ee179d629eec1a"


def bounds(physical_shape):
    raw = physical_shape.read_bytes()
    if hashlib.sha256(raw).hexdigest() != DIGEST:
        raise ValueError("expected the published physical two-window certificate")
    late = Fraction(json.loads(raw)["scaled_window_means"][0][0])
    error = Fraction(1, 2000)
    cosine_upper = error / (late - error)
    smaller_sign = (late - 2 * error) / 2
    impulse_ratio = 9
    free_mixed_error = late / (2 * (1 + impulse_ratio))
    signed_mixed = late - 2 * (1 + impulse_ratio) * error
    if not (cosine_upper < Fraction(125, 10**6) and smaller_sign > Fraction(20061, 10000)):
        raise ValueError("displacement witnesses")
    if not (free_mixed_error > Fraction(20066, 100000) and signed_mixed > 4):
        raise ValueError("mixed-release witnesses")
    return {
        "physical_shape_sha256": DIGEST,
        "physical_late_mean_lower": str(late),
        "allowed_absolute_error": str(error),
        "absolute_cosine_upper_for_displacement_only": str(cosine_upper),
        "each_signed_total_lower_for_displacement_only": str(smaller_sign),
        "both_sign_coherent_mixed_error_floor": str(free_mixed_error),
        "floor_on_20_N_displacement_plus_18_N_velocity": str(signed_mixed),
        "impulse_comparison_constant": impulse_ratio,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("physical_shape", type=Path, nargs="?",
                        default=Path(__file__).resolve().parent / "inputs" / "SHAPE.json")
    print(json.dumps(bounds(parser.parse_args().physical_shape), indent=2))
