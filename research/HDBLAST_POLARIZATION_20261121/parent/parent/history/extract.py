"""Extract all accepted native steps spanning the first 5% of whole-run absolute electrical work."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid


def extract(native, destination):
    metadata = json.loads((native / "metadata.json").read_text())
    destination.mkdir(parents=True, exist_ok=True)
    summary = {key: metadata[key] for key in ("action", "daughter_mass_L", "charge_to_mass")}
    summary["runs"] = {}
    columns = metadata["columns"]
    for label, run in metadata["runs"].items():
        raw = np.memmap(native / (label + ".bin"), dtype="=f8", mode="r").reshape(-1, len(columns))
        if len(raw) != run["rows"]:
            raise ValueError(f"Incomplete native history for {label}")
        data = {k: raw[:, i] for i, k in enumerate(columns)}
        first = np.searchsorted(data["D"], 0., side="right") - 1
        if first < 0 or data["D"][first] != 0.:
            raise ValueError("Positive population at startup")
        limit = np.searchsorted(data["t"], run["first_hot_time"]) + 1
        startup = slice(first, limit)
        t, x, u, rest = (data[k][startup] for k in ("t", "X", "u", "D"))
        abs_rate = np.abs(summary["charge_to_mass"] * x * rest * u / np.hypot(1., u))
        work = cumulative_trapezoid(abs_rate, t, initial=0.)
        fractions = (.001, .01, .05)
        marks = {str(f): {"time": float(np.interp(f * run["abs_electric"], work, t))}
                 for f in fractions}
        if work[-1] < fractions[-1] * run["abs_electric"]:
            raise ValueError("Insufficient work window in input")
        last = first + np.searchsorted(t, marks["0.05"]["time"]) + 1
        keep = slice(first, last)
        # The native initial state has C=0, a=2, X=0 and d(aC)/dt=aX.
        ts = np.concatenate(([0.], data["t"][:last]))
        ax = np.concatenate(([0.], data["a"][:last] * data["X"][:last]))
        integral = cumulative_trapezoid(ax, ts, initial=0.)[1:]
        arrays = {k: data[k][keep].copy() for k in columns}
        arrays["B"] = 2. * integral[keep] / arrays["a"]**2
        arrays["I"] = summary["charge_to_mass"] * cumulative_trapezoid(
            arrays["a"] * arrays["X"], arrays["t"], initial=0.)
        np.savez_compressed(destination / (label + ".npz"), **arrays)
        for point in marks.values():
            for k in ("a", "X", "D", "theta", "R", "B"):
                point[k] = float(np.interp(point["time"], arrays["t"], arrays[k]))
        summary["runs"][label] = dict(markers=marks, channel_fraction=run["channel_fraction"],
                                      native_rows=run["rows"], absolute_work=run["abs_electric"],
                                      first_birth_left_knot=float(arrays["t"][0]),
                                      retained_steps=len(arrays["t"]))
    (destination / "markers.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    extract(args.native, args.output_dir)
