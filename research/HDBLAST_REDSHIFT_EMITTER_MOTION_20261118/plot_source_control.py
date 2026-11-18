"""Draw the fixed-redshift cap-distance derivative for prescribed emitter motion."""
from fractions import Fraction
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
import numpy as np

ROOT = Path(__file__).resolve().parent
report = json.loads((ROOT/"SOURCE_CONTROL.json").read_text())


def pair(name):
    return np.array([float(Fraction(s)) for s in report[name]])


static = pair("supported_fractional_slope") * 1e10
free = pair("geodesic_fractional_slope") * 1e10
gain = pair("doppler_over_abs_supported")
cross = pair("critical_velocity_fraction")


def envelope(alpha):
    # Convex endpoint interpolation and H*(alpha*K-1) enclose the same physical
    # line. The latter retains the common overall scale through the zero.
    convex_low = (1-alpha)*static[0] + alpha*free[0]
    convex_high = (1-alpha)*static[1] + alpha*free[1]
    corners = np.array([(-s)*(alpha*k-1) for s in static for k in gain])
    return np.maximum(convex_low, corners.min(axis=0)), np.minimum(convex_high, corners.max(axis=0))


fig, axes = plt.subplots(1, 2, figsize=(10.1, 4.5), layout="constrained",
                         gridspec_kw={"width_ratios": [1.35, 1]})
line, band, threshold = "#134a65", "#9bbecf", "#96602e"
for ax, (left, right) in zip(axes, [(0, 1), (0.055, 0.062)]):
    x = np.linspace(left, right, 1401)
    low, high = envelope(x)
    ax.fill_between(x, low, high, color=band)
    ax.plot(x, low, color=line, lw=1.15)
    ax.plot(x, high, color=line, lw=1.15)
    ax.axvspan(*cross, color="#d8ad6c", alpha=0.30, zorder=0)
    ax.axhline(0, color="#373c3b", lw=0.95, ls="--")
    ax.set_xlim(left, right)
    ax.grid(axis="y", color="#e4e7e8", lw=0.8)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlabel("Fraction of free-fall velocity retained, α")
axes[0].set(ylabel="Fractional distance derivative  (×10⁻¹⁰)", ylim=(-1.05, 8.9),
            title="Full range of emitter motions", xticks=[0,.2,.4,.6,.8,1])
axes[1].set(ylim=(-.04,.04), title="Detail near the sign change",
            xticks=[.055,.057,.059,.061], yticks=[-.04,-.02,0,.02,.04])
axes[0].text(.025,-.80,"Supported", color=line, fontsize=9)
axes[0].text(.98,8.45,"Free fall", ha="right", color=line, fontsize=9)
axes[1].text(.05833,.035,"5.808–5.856%", ha="center", va="center", color=threshold, fontsize=9)
axes[1].legend(handles=[Patch(facecolor=band, edgecolor=line, label="Enclosed model response"),
                         Patch(facecolor="#f0dfc5", label="Zero-crossing interval")],
               loc="lower right", frameon=False, fontsize=8.5)
fig.suptitle("Two sky caps at the same observed frequency ratio", fontsize=14)
fig.savefig(ROOT/"source_control.png", dpi=180, facecolor="white")
