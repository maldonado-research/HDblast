"""Extract B1 clock, energy-density and geometric diagnostics from authenticated runs."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np


def load(path, sha):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != sha:
        raise ValueError(f"Unexpected input bytes for {path.name}")
    return data


def export(archive, grid):
    spec = json.loads(Path(__file__).with_name("INPUTS.json").read_text())[grid]
    run = archive / "runs" / "main" / spec["tag"]
    summary = json.loads(load(run.with_name(run.name + "_summary.json"), spec["summary_sha256"]))
    with np.load(io.BytesIO(load(run.with_name(run.name + "_timeseries.npz"), spec["npz_sha256"])),
                 allow_pickle=False) as raw:
        values = dict(zip(raw["cols"].tolist(), raw["rec"].T, strict=True))
    clock = values["H0tau"]
    if not (np.all(np.diff(clock) > 0) and clock[0] <= 1 < 6 < spec["reliable_end"] <= clock[-1]):
        raise ValueError("Clock or reliability interval")
    if (summary["rho_b"],summary["params"]["matter"]["b"]) != (spec["rho_b"],spec["b"]):
        raise ValueError("Model normalization")
    sample = np.unique(np.r_[1,4,6,clock[(clock > 1) & (clock < 6)]])
    keys = ("H_over_H0", "R", "ln_a", "rad", "Wy", "vac")
    stream = io.StringIO(newline="")
    output = csv.writer(stream, lineterminator="\n")
    output.writerow(["s", *keys])
    with np.errstate(invalid="raise"):
        columns = [np.interp(sample,clock,values[key]) for key in keys]
    for row in zip(sample,*columns,strict=True):
        output.writerow([format(float(v),".17g") for v in row])
    return stream.getvalue()


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive",type=Path,help="B1_derived_chi_source directory containing original runs")
    parser.add_argument("grid",choices=("coarse","fine"))
    args=parser.parse_args()
    print(export(args.archive,args.grid),end="")
