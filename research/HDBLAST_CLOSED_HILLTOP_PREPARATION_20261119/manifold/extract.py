"""Extract the local phi and conformal-factor spline, with matching data from the prior study."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np


def extract(prior):
    sources = ["background/background.py", "hilltop/HILLTOP.json",
               "archive/blast_common.py", "archive/BLAST_REHEATING_RESULTS.json"]
    hashes = {p: hashlib.sha256((prior / p).read_bytes()).hexdigest() for p in sources}
    spec = importlib.util.spec_from_file_location("prior_background", prior / sources[0])
    background = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(background)
    module, histories = background.source(prior / "archive")
    tables = module.build_tables()
    thresholds = json.loads((prior / sources[1]).read_text())["branches"]
    branches = {}
    for key, row in thresholds.items():
        pm = float(key.split("=")[1])
        potential = background.Potential(module, tables, pm, histories[key])
        keep = np.flatnonzero(abs(potential.theta) < .030)
        first, last = keep[0] - 1, keep[-1] + 1
        critical = row["detunings"]["0.001"]["primary"]
        branches[key] = dict(minimum=pm, b=potential.b, v0=potential.v0,
                             top=row["stationary_theta"], critical_x=critical["x_match"],
                             critical_log_d=critical["log_abs_bounce_seed"],
                             source_mu_squared=row["negative_mass_squared_over_Htop_squared"],
                             spline_knots=potential.theta[first:last+1].tolist(),
                             phi_f_coefficients=potential.coeff[:,first:last].tolist())
    return dict(prior_sha256=hashes, branches=branches)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--prior", required=True, type=Path)
    cli.add_argument("--output", required=True, type=Path)
    args = cli.parse_args()
    args.output.write_text(json.dumps(extract(args.prior), separators=(",", ":")) + "\n")
