"""Crossover fraction for a prescribed family of source accelerations."""
import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

RED = "40b61502f6877d6c91614d2eded6db2c043fd38a69a84928ccc7309bbeb08a09"


def control(path):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != RED:
        raise ValueError("unexpected fixed-redshift report")
    report = json.loads(data)
    gain = list(map(F, report["doppler_to_abs_supported_term"]))
    supported = list(map(F,report["comparisons"]["supported"]["fractional_distance_cap_half_difference"]))
    geodesic = list(map(F,report["comparisons"]["geodesic"]["fractional_distance_cap_half_difference"]))
    if not supported[0] <= supported[1] < 0 < geodesic[0] <= geodesic[1] or not 1 < gain[0] <= gain[1]:
        raise ValueError("clock signs or Doppler ratio")
    critical = [1/gain[1], 1/gain[0]]
    if not F("0.0580844") < critical[0] <= critical[1] < F("0.0585507"):
        raise ValueError("critical fraction")
    return {"redshift_sha256":RED, "supported_fractional_slope":list(map(str,supported)),
            "geodesic_fractional_slope":list(map(str,geodesic)),
            "doppler_over_abs_supported":list(map(str,gain)),
            "critical_velocity_fraction":list(map(str,critical)),
            "critical_percent_display_only":[float(100*v) for v in critical]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("redshift_report", nargs="?", type=Path,
                        default=Path(__file__).resolve().with_name("REDSHIFT.json"))
    print(json.dumps(control(parser.parse_args().redshift_report), indent=2))
