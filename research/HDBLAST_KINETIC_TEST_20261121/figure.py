"""Render the collision-budget and momentum-distribution diagnostics."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

ROOT = Path(__file__).resolve().parent
reversals = json.loads((ROOT / "ballistic/REVERSAL.json").read_text())["runs"]
oriented = json.loads((ROOT / "driven/ORIENTED.json").read_text())["groups"]
fig, (left, right) = plt.subplots(1, 2, figsize=(11.8, 5.25))
fig.subplots_adjust(left=.08, right=.97, top=.78, bottom=.25, wspace=.32)
fig.suptitle("How quickly can the wall particles spread out?", fontsize=20, y=.96)
blue, gold, plum = "#176985", "#BC6715", "#754D86"
for pos, (channel, label) in enumerate([("dedicated", "Dedicated"), ("shared", "Shared")]):
    bounds = reversals[channel + "_half"]["cell_enclosure"]
    for offset, key, color, series in [
        (-.19, "upper_collision_incidences", blue, "Maximum share reached by encounters"),
        (.19, "minimum_hit_fraction_for_longitudinal_stress", gold, "Minimum participation required"),
    ]:
        value = bounds[key]
        left.bar(pos + offset, value, .34, color=color, label=series if pos == 0 else None)
        left.text(pos + offset, value + .012, f"{100 * value:.1f}%", ha="center", fontsize=10)
left.set_title("A. Can the early fluid state be reached?", loc="left", fontsize=13)
left.set_ylabel("Share of particles present")
left.set_xticks([0, 1], ["Dedicated channel", "Shared channel"])
left.set_ylim(0, .8)
left.legend(loc="upper center", bbox_to_anchor=(.48, -.14), frameon=False, fontsize=9)
markers = ["half_work", "reference_first_hot", "radius_doubling"]
labels = ["Half early\nwork", "Fluid first\npasses 90%", "Reference\nradius doubles"]
for channel, label, color, marker in [("dedicated", "Dedicated channel", blue, "o"),
                                      ("shared", "Shared channel", plum, "s")]:
    rows = oriented["more_particles_" + channel + "_magnetic"]
    pressure = [rows[m]["mean"]["rest_longitudinal_share"] for m in markers]
    right.plot(range(3), pressure, marker=marker, color=color, lw=2, label=label)
right.axhline(1 / 3, color="#444444", ls="--", lw=1)
right.text(-.03, .25, "Equal pressure in all directions (⅓)", fontsize=9)
right.set_title("B. What do individual particles do?", loc="left", fontsize=13)
right.set_ylabel("Axial share of stream pressure")
right.set_xticks(range(3), labels)
right.set_xlim(-.12, 2.12)
right.set_ylim(0, 1.04)
right.legend(loc="upper center", bbox_to_anchor=(.50, -.19), ncol=2, frameon=False, fontsize=9)
for ax in (left, right):
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, alpha=.18)
    ax.yaxis.set_major_formatter(PercentFormatter(1))
    ax.spines[["right", "top"]].set_visible(False)
fig.text(.08,.03, "A tests each fluid at its first current reversal. B follows particles on the original expansion history.", fontsize=9)
fig.savefig(ROOT / "kinetic-tests.png", dpi=175, facecolor="white")
