"""Randomly scrambled Sobol convergence at the most discrepant small-sample probe."""
import argparse
import json
import math
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import sys

import numpy as np
from scipy.stats import qmc, t

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "boost"))
import sampling
import thermal

THETA = 10.0
RAPIDITY = 4.0
SEEDS = (431, 522, 613, 704, 795, 886, 977, 1068, 1159, 1250, 1341, 1432)
MIN_POWER = 16
MAX_POWER = 21


def scramble(seed):
    """Use nested prefixes and equal-sized blocks to limit memory per worker."""
    sampler = qmc.Sobol(d=5, scramble=True, seed=seed)
    running_sum = 0.0
    estimates = []
    for index in range(1, 2**(MAX_POWER - MIN_POWER) + 1):
        pair_samples = sampler.random(2**MIN_POWER)
        running_sum += sampling.rate(THETA, RAPIDITY, pair_samples)
        if index & (index - 1) == 0:
            estimates.append(running_sum / index)
    return estimates


def run(jobs):
    if jobs == 1:
        rates = list(map(scramble, SEEDS))
    else:
        with ProcessPoolExecutor(max_workers=jobs) as pool:
            rates = list(pool.map(scramble, SEEDS))
    reference = thermal.rate(THETA, RAPIDITY, order=32)
    result = []
    for power, values in enumerate(np.array(rates).T, MIN_POWER):
        mean = float(values.mean())
        sem = float(values.std(ddof=1) / math.sqrt(len(values)))
        margin = sem * t.ppf(.975, len(values) - 1)
        result.append(dict(draws_per_scramble=2**power, mean_rate=mean,
                           relative_difference=mean / reference - 1,
                           standard_error=sem, relative_standard_error=sem / mean,
                           approximate_95pct_interval=[mean - margin, mean + margin],
                           rates_by_seed=values.tolist()))
    return dict(theta=THETA, rapidity=RAPIDITY, seeds=SEEDS, quadrature_rate=reference,
                rate_units="wall drag multiplied by mass squared / proper target density",
                sampling="nested prefixes per scramble; randomizations independent across seeds",
                convergence=result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args()
    args.output.write_text(json.dumps(run(args.jobs), allow_nan=False, indent=2) + "\n")
