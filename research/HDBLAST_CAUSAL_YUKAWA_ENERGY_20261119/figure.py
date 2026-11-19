"""Plot the different energy accounts and the outgoing pulse signal."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "energy.png")
    args = parser.parse_args()
    runs = json.loads((HERE / "pulses/PULSES.json").read_text())["runs"]
    results = runs["fine_120"]["quiet"]
    clocks = json.loads((HERE / "clock/CLOCK.json").read_text())["grids"]["fine"]["arms"]
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.8, 4.15),
                                     gridspec_kw={"width_ratios": [1.35, 1.]}, layout="constrained")
    names = ["parent", "daughter", "switch", "mass", "scale"]
    labels = ["Parent", "Pairs", "Switch", "Mass", "Expansion"]
    x = np.arange(len(names))
    colors = ["#245777", "#bf622a"]
    for i, (label, scheme) in enumerate([("Bare", "bare"), ("Static calibration", "static")]):
        left.bar(x + (i - .5)*.36, [results[scheme][key] for key in names], .35,
                 color=colors[i], label=label, zorder=2)
    left.axvline(1.5, color="#b6b6b6", ls="--", lw=1)
    left.set_xticks(x, labels)
    left.set_ylabel(r"Coefficient per $y^2 H_0^4$ (comoving volume)")
    left.set_title("Who supplies the pair energy?", loc="left", weight="bold")
    left.legend(frameon=False, ncols=2, loc="lower right", fontsize=9)
    left.set_ylim(-50, 29)
    arms = ["plus", "minus", "far"]
    x = np.arange(3)
    for i, (label, clock, color) in enumerate([("Proper-time clock", "proper_time", "#246b68"),
                                                ("Scalar clock", "scalar", "#96bbab")]):
        delta = [clocks[arm]["energy"]["120"][clock + "_pulse_percent_of_quiet"] for arm in arms]
        right.bar(x + (i - .5)*.35, delta, .33, label=label, color=color, zorder=2)
    right.set_xticks(x, ["Plus pulse", "Minus pulse", "Far pulse"])
    right.set_title("What sets the pulse signal?", loc="left", weight="bold")
    right.set_ylabel("Pair-energy change from quiet (%)")
    right.legend(frameon=False, loc="upper right", fontsize=9)
    right.set_ylim(-.49, 2.0)
    for ax in (left, right):
        ax.axhline(0, color="#38444a", lw=1)
        ax.grid(axis="y", color="#dddddd", lw=.65, zorder=0)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(length=0, pad=5)
    fig.savefig(args.output, dpi=180, facecolor="white")

if __name__ == "__main__":
    main()
