"""Plot the certified shell trace and the two-exposure impulse inequality."""
import csv
import json
from fractions import Fraction
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parent
TARGET = json.loads((ROOT / "KICK.json").read_text())
lower, upper = map(lambda value: float(Fraction(value)), TARGET["physical_late_mean"])
time, lo, hi = [], [], []
with (ROOT / "inputs" / "FUTURE_CELLS.csv").open() as stream:
    for cell in csv.DictReader(stream):
        left = int(cell["left_grid"]) / 4096
        right = int(cell["right_grid"]) / 4096
        if left >= 25 / 16:
            break
        lval = int(cell["H_lo"]) / (2**112) - 3e-6
        hval = int(cell["H_hi"]) / (2**112) + 3e-6
        lower_a = min(lval, lval * (1 - 1e-11))
        upper_a = max(hval, hval * (1 - 1e-11))
        if right <= 1:  # The analytic support result gives exact silence here.
            lower_a = upper_a = 0
        time.extend((left, min(right, 25 / 16)))
        lo.extend((lower_a, lower_a))
        hi.extend((upper_a, upper_a))

fig, (ax, bx) = plt.subplots(1, 2, figsize=(10, 4.25), layout="constrained")
navy, aqua, purple, orange = "#163c52", "#317d79", "#634795", "#b44721"
mid = (np.asarray(lo) + np.asarray(hi)) / 2
ax.axvspan(0, 1, color=aqua, alpha=0.13)
ax.axvspan(1.25, 1.375, color=purple, alpha=0.14)
ax.axhline(0, color="0.78", lw=0.8)
ax.fill_between(time, lo, hi, color=navy)
ax.plot(time, mid, lw=1.15, color=navy)
ax.text(0.48, 3.7, "Early triangle", color=aqua, ha="center", fontsize=10)
ax.annotate("Late mean", xy=(1.32, 3.9), xytext=(0.92, 4.9),
            arrowprops=dict(arrowstyle="-", color=purple), color=purple, fontsize=10)
ax.set(xlim=(0, 25 / 16), ylim=(-2.4, 5.4), xlabel="Elapsed model time  x",
       ylabel="Normalized shell quadrupole  A", title="Two exposures of the predicted signal")
x = np.linspace(-0.63, 0.63, 500)
bx.fill_between(x, -9*np.abs(x), 9*np.abs(x), color=aqua, alpha=0.19)
bx.plot(x, 9*np.abs(x), color=aqua, lw=1.2)
bx.plot(x, -9*np.abs(x), color=aqua, lw=1.2)
bx.vlines(0, lower, upper, color=orange, lw=3)
bx.plot(0, (lower + upper) / 2, "o", color=orange, markersize=6)
err = lower / 10
bx.add_patch(Rectangle((-err, lower-err), 2*err, 2*err, fill=False, ls="--",
                       color=orange, lw=1.0))
bx.text(0, 4.68, "Shell target", ha="center", color=orange, fontsize=10)
bx.text(0.42, 0.62, "Satisfies the\nkick bound", ha="center", color=aqua, fontsize=10)
bx.set(xlim=(-0.62, 0.62), ylim=(-0.65, 5.35), xlabel="Early triangular exposure  T",
       ylabel="Late mean  L", title="A matched four-dimensional kick")
for axis in (ax, bx):
    axis.spines[["top", "right"]].set_visible(False)
    axis.grid(axis="y", lw=0.6, color="0.87")
    axis.set_axisbelow(True)
fig.savefig(ROOT / "exposures.png", dpi=175, facecolor="white")
