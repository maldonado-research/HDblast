#!/usr/bin/env python3
"""Draw completed continuous physical-coefficient enclosures; certify nothing.

Every plotted band comes directly from the audited whole-cell CSV. Floating
values are created only after input authentication, solely for rendering.
"""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction as F
from hashlib import sha256
import io
import json
import os
from pathlib import Path
import tempfile

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def plot(source_root, output):
    source_root, output = Path(source_root), Path(output)
    report_raw = (source_root/"computed/INERTIAL_WAVEFORM_ENCLOSURE.json").read_bytes()
    csv_raw = (source_root/"computed/INERTIAL_WAVEFORM_CELLS.csv").read_bytes()
    review = json.loads((source_root/"sol/CONTINUOUS_WAVEFORM_REVIEW.json").read_text())
    report = json.loads(report_raw)
    require(report["status"] == "PASS", "completed producer PASS required")
    require(review["status"] == "PASS_INDEPENDENT_CONTINUOUS_WAVEFORM_REVIEW", "independent review PASS required")
    require(sha256(report_raw).hexdigest() == review["inputs"]["report"]["sha256"], "reviewed report authentication")
    require(sha256(csv_raw).hexdigest() == review["inputs"]["csv"]["sha256"] == report["cell_data"]["sha256"], "reviewed CSV authentication")
    data = list(csv.DictReader(io.StringIO(csv_raw.decode("utf-8"), newline="")))
    N, S = report["cell_data"]["time_grid_denominator"], report["cell_data"]["value_denominator"]
    # All conversions below this point serve display only. Their values do
    # not establish signs, error bounds, integral endpoints or PASS status.
    x = [float(F(int(row["left_grid"]), N)) for row in data]
    x.append(float(F(int(data[-1]["right_grid"]), N)))
    lower = [float(F(int(row["ell_S_lo"]), S)) for row in data]
    upper = [float(F(int(row["ell_S_hi"]), S)) for row in data]
    lower.append(lower[-1]); upper.append(upper[-1])
    center = [(a+b)/2 for a,b in zip(lower, upper)]
    with tempfile.TemporaryDirectory(prefix="inertial-mpl-") as cache:
        os.environ["MPLCONFIGDIR"] = cache
        os.environ["XDG_CACHE_HOME"] = cache
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import MultipleLocator

        plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                             "axes.spines.top": False, "axes.spines.right": False,
                             "axes.labelcolor": "#213547", "text.color": "#213547",
                             "xtick.color": "#526371", "ytick.color": "#526371"})
        fig, ax = plt.subplots(figsize=(13.4, 6.8), dpi=200)
        fig.subplots_adjust(left=.085, right=.97, top=.82, bottom=.18)
        fig.text(.085, .935, "The shell's linear first variation", fontsize=23, weight="semibold")
        fig.text(.085, .885, r"Physical quadrupole coefficient $\ell S$  |  retained branch, $\delta=10^{-3}$, $\ell=10^{-6}$", fontsize=12, color="#526371")
        ax.axvspan(2, 3, color="#dceef0", alpha=.42, zorder=0)
        ax.axhline(0, color="#a8b6c0", linewidth=.8, linestyle=(0,(4,4)), zorder=1)
        ax.fill_between(x, lower, upper, step="post", color="#2a8e96", alpha=.32, linewidth=0, label="Continuous enclosure", zorder=2)
        ax.step(x, center, where="post", color="#17636e", linewidth=1.2, label="Enclosure midpoint guide", zorder=3)
        ax.set_xlim(0, 3)
        ax.set_ylim(-.052, .052)
        ax.set_xlabel(r"Rescaled proper time  $x=\tau/(\rho_b\ell)$", labelpad=12)
        ax.set_ylabel(r"Physical coefficient  $\ell S$", labelpad=12)
        ax.xaxis.set_major_locator(MultipleLocator(.5))
        ax.yaxis.set_major_locator(MultipleLocator(.02))
        ax.grid(axis="y", color="#e5eaed", linewidth=.7)
        ax.text(.43, .006, "Exact prefront zero", color="#71818d", fontsize=10, ha="center")
        ax.text(1.27, -.050, "Negative physical lobe", color="#17636e", fontsize=10, ha="center")
        ax.text(1.72, .049, "Positive physical lobe", color="#17636e", fontsize=10, ha="center")
        ax.legend(loc="lower right", frameon=False, fontsize=9.5, bbox_to_anchor=(1.0,-.005))

        inset = ax.inset_axes([.665, .57, .315, .34])
        tail_start = next(i for i, row in enumerate(data) if int(row["left_grid"]) >= 2*N)
        tx = x[tail_start:]
        tl = [v*1e10 for v in lower[tail_start:]]
        tu = [v*1e10 for v in upper[tail_start:]]
        inset.fill_between(tx, tl, tu, step="post", color="#c66a45", alpha=.30, linewidth=0)
        inset.step(tx, [(a+b)/2 for a,b in zip(tl,tu)], where="post", color="#a94d2b", linewidth=1)
        inset.set_xlim(2,3)
        inset.set_ylim(-6.3,-2.7)
        inset.set_xticks([2,2.5,3])
        inset.set_yticks([-6,-5,-4,-3])
        inset.set_title("Signed post-pulse tail", fontsize=11, pad=9, loc="left", weight="semibold")
        inset.set_ylabel(r"$\ell S\ \times 10^{10}$", fontsize=10, labelpad=4)
        inset.tick_params(labelsize=9)
        inset.grid(axis="y", color="#e9e3de", linewidth=.6)
        inset.set_facecolor("#fffdfb")
        fig.text(.085, .065, "Bands enclose every time in each exported cell. Midpoints are visual guides; floats are used only to draw this figure.", fontsize=9.5, color="#62727e")
        fig.text(.085, .035, "Conditional linear result. No finite-amplitude nonlinear prediction or observational time calibration.", fontsize=9.5, color="#62727e")
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, dpi=200, facecolor="white", metadata={
            "Title": "Inertial shell linear first-variation waveform",
            "SourceCSV_SHA256": sha256(csv_raw).hexdigest(),
            "SourceReport_SHA256": sha256(report_raw).hexdigest(),
            "Description": "Visualization of audited continuous physical coefficient bands; floating values are for display only.",
        })
        plt.close(fig)
    return {"path": str(output), "bytes": output.stat().st_size, "sha256": sha256(output.read_bytes()).hexdigest()}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=HERE.parent)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or args.source_root/"computed/INERTIAL_SHELL_WAVEFORM.png"
    print(json.dumps(plot(args.source_root, output), indent=2))
