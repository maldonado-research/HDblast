"""Regenerate one table from its original NPZ and public run summary. Requires NumPy."""
import argparse
import csv
import hashlib
import io
import json
import sys
from pathlib import Path

import numpy as np


def authenticated(path, digest):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != digest:
        raise ValueError(f"Input checksum mismatch: {path.name}")
    return data


def derive(name, meta, arrays, params):
    x = arrays["H0tau"]
    if name.startswith("b1_"):
        a4 = np.exp(4 * (arrays["A_b"] - arrays["A_b"][0]))
        radiation = params["rho_b"] * a4 * arrays["R"]
        matter = params["params"]["matter"]
        mass2 = matter["mu0"]**2 + matter["G"]**2 * (arrays["phi_b"] - matter["phistar"])**2
        decay = matter["y"]**2 * mass2 * arrays["n_hat"] / (8 * np.pi)
        source = matter["b"] * a4 * decay
    else:
        a4 = np.exp(4 * arrays["ln_a"])
        radiation = a4 * arrays["R"]
        source = (params["params"]["Y"] / params["rho_b"]) * a4 * arrays["v_over_H0"]**2
    start, end = meta["clock_interval"]
    if not (np.all(np.diff(x) > 0) and x[0] <= start < end <= meta["reported_reliable_end"] <= x[-1]):
        raise ValueError("Invalid recorded grid or requested interval")
    times = np.r_[start, x[(x > start) & (x < end)], end]
    rows = zip(times, np.interp(times, x, radiation), np.interp(times, x, source))
    table = io.StringIO(newline="")
    writer = csv.writer(table, lineterminator="\n")
    writer.writerow(["x", "weighted_radiation", "weighted_source"])
    writer.writerows([format(float(value), ".17g") for value in row] for row in rows)
    return table.getvalue().encode("ascii")


def main():
    registry = json.loads(Path(__file__).with_name("PROVENANCE.json").read_text())["series"]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("series", choices=sorted(registry))
    parser.add_argument("source_npz", type=Path)
    parser.add_argument("run_summary", type=Path)
    args = parser.parse_args()
    meta = registry[args.series]
    source = authenticated(args.source_npz, meta["history_npz_sha256"])
    params = json.loads(authenticated(args.run_summary, meta["run_summary"]["sha256"]))
    with np.load(io.BytesIO(source), allow_pickle=False) as stored:
        arrays = dict(zip(stored["cols"].tolist(), stored["rec"].T, strict=True))
        table = derive(args.series, meta, arrays, params)
    if hashlib.sha256(table).hexdigest() != meta["csv_sha256"]:
        raise ValueError("Derived table differs from the registered bytes")
    sys.stdout.buffer.write(table)


if __name__ == "__main__":
    main()
