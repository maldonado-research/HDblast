"""Plot the expanding primary and the matched-grid decay-rate comparison in JOINED.json."""
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter


HERE = Path(__file__).resolve().parent
INK, TEAL, MUTED, RULE = "#172b36", "#007870", "#4e6069", "#d5dfe0"
AMBER = "#a44d12"


def expanding_curve(data, run):
    rows = [dict(zip(data["trace_columns"], v)) for v in run["trace"]]
    selected = [r for r in rows if run["bounce_tau"] < r["t"] < run["maximum_tau"] and r["H"] > 0]
    a = [r["a"] for r in selected]
    fraction = [100 * r["intrinsic"] / r["e"] for r in selected]
    # Interpolate the energy at the native H=0 endpoint without tracing contraction.
    before = selected[-1]
    after = next(r for r in rows if r["t"] >= run["maximum_tau"])
    weight = (run["maximum_tau"] - before["t"]) / (after["t"] - before["t"])
    energy = before["e"] + weight * (after["e"] - before["e"])
    heat = before["intrinsic"] + weight * (after["intrinsic"] - before["intrinsic"])
    return a + [run["maximum_radius"]], fraction + [100 * heat / energy]


def main():
    data = json.loads((HERE / "coupled" / "JOINED.json").read_bytes())
    run = data["runs"]["primary"]
    a, fraction = expanding_curve(data, run)
    hot, double = run["first_hot_radius"], run["doubling_radius"]
    maximum = run["maximum_radius"]
    at_double = 100 * run["doubling_hot_fraction"]

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11,
                         "axes.labelcolor": INK, "text.color": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED})
    fig, (left, right) = plt.subplots(1, 2, figsize=(12.3, 5.65),
                                      gridspec_kw={"width_ratios": [1.12, 1]})
    fig.subplots_adjust(left=.073, right=.966, top=.785, bottom=.213, wspace=.34)
    fig.suptitle("Particle heating during wall expansion", x=.073, y=.97,
                 ha="left", fontsize=23, weight="semibold")
    fig.text(.073, .89, "The hot phase begins when random particle motion carries 90% of the added wall energy.",
             fontsize=11.5, color=MUTED)

    left.set_title("Energy carried by particle heat", loc="left", fontsize=14, weight="semibold", pad=14)
    left.axvspan(hot, maximum, color="#e6f3f0", linewidth=0)
    left.axhline(90, color=AMBER, linestyle=(0, (4, 3)), linewidth=1.1)
    left.plot(a, fraction, color=TEAL, linewidth=2.5)
    left.scatter([hot, double, maximum], [90, at_double, fraction[-1]],
                 s=[36, 43, 36], facecolors=[TEAL, TEAL, INK], edgecolors="white", linewidth=.9, zorder=5)
    notes = [("Hot phase begins", (hot, 90), (12.0, 53), "left"),
             (f"{at_double:.1f}% at double\nthe starting hot radius", (double, at_double), (23.6, 69), "center"),
             ("Expansion stops", (maximum, fraction[-1]), (35.8, 42), "right")]
    for label, point, position, align in notes:
        left.annotate(label, xy=point, xytext=position, ha=align, va="center", fontsize=9.5,
                      color=INK, arrowprops={"arrowstyle": "-", "color": MUTED, "lw": .8,
                                              "shrinkA": 4, "shrinkB": 4})
    left.text(2.0, 92.5, "90% threshold", color=AMBER, fontsize=9)
    left.set(xlim=(1, 36.5), ylim=(-2.5, 104), xlabel="Wall radius (model units)",
             ylabel="Share of added wall energy")
    left.set_xticks([2, 10, 20, 30, 35])
    left.set_yticks([0, 25, 50, 75, 100])
    left.yaxis.set_major_formatter(PercentFormatter(100, decimals=0))
    left.grid(axis="y", color=RULE, linewidth=.55, zorder=0)
    left.set_axisbelow(True)

    right.set_title("Effect of changing the decay rate", loc="left", fontsize=14, weight="semibold", pad=14)
    cases = ["slower004", "slower04", "weak_primary", "reference18", "faster_primary", "rapid_primary"]
    kappas = [data["parameters"][k]["kappa"] for k in cases]
    counts = [data["runs"][k].get("hot_efolds", 0.) / math.log(2) for k in cases]
    right.axvline(1, color=AMBER, linestyle=(0, (4, 3)), linewidth=1.1, zorder=0)
    for i, value in enumerate(counts):
        color = TEAL if value >= 1 else MUTED
        right.hlines(i, 0, value, color="#b8cdcb", linewidth=2.5)
        right.scatter(value, i, s=42, color=color, zorder=3)
        # Label the nearly sufficient case on the left of the criterion line.
        if .8 < value < 1.:
            right.text(value-.045, i-.18, f"{value:.2f}", ha="right", va="center", fontsize=9.5, color=color)
        else:
            right.text(value+.042, i, f"{value:.2f}" if value else "0", va="center", fontsize=9.5, color=color)
    right.set(xlim=(-.03, 1.72), ylim=(len(cases)-.6, -.6),
              xlabel="Radius doublings during the hot phase", ylabel="Decay-rate coefficient (κ)")
    right.set_xticks([0, .5, 1, 1.5])
    right.set_yticks(range(len(cases)), [f"{k:g}" for k in kappas])
    right.text(1.012, -.48, "One doubling", ha="left", va="center", color=AMBER, fontsize=9)
    for ax in (left, right):
        ax.spines[["right", "top"]].set_visible(False)
        ax.spines[["left", "bottom"]].set_color("#8c9a9d")
        ax.tick_params(length=3, width=.7)
        ax.xaxis.labelpad = 8
    fig.text(.073, .068, "Reference model, 48 early and 6 late radial elements.", color=MUTED, fontsize=9.3)
    fig.text(.574, .068, "Rate comparison, 18 early and 3 late radial elements.", color=MUTED, fontsize=9.3)
    fig.savefig(HERE / "overview.png", dpi=190, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
