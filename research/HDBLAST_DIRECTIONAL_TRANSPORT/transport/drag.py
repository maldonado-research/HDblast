"""Spin and flavor bounds on momentum drag by one neutral scalar resonance."""
import argparse
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path


INPUT_SHA = "0f3421ac6d30927ba32f3c8a8a431b591c27fb3307a161c7740c77cda0c0caad"


def cross_section_check():
    # At the pole divide the spin-averaged Breit-Wigner cross section by
    # the reference 4*pi/k_cm**2. The total width sums all open flavors.
    for couplings in ((0.1, 0.1, 0.1), (0.12, 0.03, 0.07)):
        for mass_ratio in (1.02, 1.05, 3.):
            parent_mass = 2 * mass_ratio
            s = parent_mass**2
            beta = math.sqrt(1-4/s)  # common Dirac mass = 1
            widths = [y*y*parent_mass*beta**3/(8*math.pi) for y in couplings]
            total = sum(widths)
            for i, yi in enumerate(couplings):
                for j, yj in enumerate(couplings):
                    cross_section = ((yi*yj)**2*(s-4)**2 /
                                     (16*math.pi*s*parent_mass**2*total**2))
                    reference = 16*math.pi/(s-4)
                    expected = widths[i]*widths[j]/(4*total**2)
                    if not math.isclose(cross_section/reference, expected, rel_tol=5e-13):
                        raise ArithmeticError("Spin trace and branching normalization disagree")


def bounds(path):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != INPUT_SHA:
        raise ValueError("Unexpected archived sweep receipt")
    cross_section_check()
    output = {}
    for key, run in json.loads(data, parse_float=F)["runs"].items():
        floor = run["reference_drag_floor"]
        cap = run["runs"]["half_step"]["continuum_pointwise_S0_upper"]
        required = floor/cap
        species = 3
        elastic = F(1, 4*species**2)
        inclusive = F(1, 4*species)
        if not 0 < elastic < inclusive < required:
            raise ArithmeticError("No incompatibility between these archived screens")
        output[key] = {
            "mass_mL": float(run["mu"]),
            "reference_drag_floor": float(floor),
            "angle_integrated_sweep_action_ceiling": float(cap),
            "required_fraction_at_published_action": float(floor/run["action"]),
            "required_fraction_at_sweep_ceiling": float(required),
            "elastic_fraction_upper": float(elastic),
            "inclusive_fraction_upper": float(inclusive),
            "elastic_action_floor": float(floor/elastic),
            "inclusive_action_floor": float(floor/inclusive),
            "elastic_minimum_shortfall_factor": float(required/elastic),
            "inclusive_minimum_shortfall_factor": float(required/inclusive),
        }
    return {"source_sha256": INPUT_SHA, "runs": output}


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sweep_json", type=Path)
    args=parser.parse_args()
    print(json.dumps(bounds(args.sweep_json), indent=2))
