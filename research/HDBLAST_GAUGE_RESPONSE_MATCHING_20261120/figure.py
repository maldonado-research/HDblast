"""Plot the stored first-minimum photon fractions. Run with python figure.py."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter

ROOT = Path(__file__).resolve().parent
BLUE = "#16638A"
ORANGE = "#AD4515"
INK = "#22323D"


def series(records, delta, highest):
    """Preserve uncompleted trajectories as gaps in photon share."""
    data = []
    for key, row in records.items():
        harmonic, contrast = key.split("/")
        if float(contrast) != delta or int(harmonic) > highest:
            continue
        if row["status"] == "kinetic_reserve":
            share = np.nan
        elif row["status"] == "first_minimum":
            share = 100 * row["end"]["photon_fraction"]
        else:
            raise ValueError(f"Unexpected termination at {key}: {row['status']}")
        data.append((int(harmonic), share))
    return np.array(sorted(data)).T


def main():
    matched = json.loads((ROOT / "matched/MATCHED.json").read_text())
    reach = json.loads((ROOT / "reach/REACH.json").read_text())
    assert matched["reserve"] == reach["reserve"] == 0.05
    branches = ("phi_m=0.5", "phi_m=-0.5")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "text.color": INK, "axes.labelcolor": INK,
                         "xtick.color": INK, "ytick.color": INK,
                         "axes.spines.right": False, "axes.spines.top": False})
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.95))
    fig.subplots_adjust(left=.06, right=.987, bottom=.285, top=.785, wspace=.20)
    fig.text(.06, .926, "How much energy reaches the photon field?",
             fontsize=20, fontweight="bold")
    fig.text(.06, .864, "Closed-universe model  ·  Photon share at the scalar field’s first potential minimum",
             fontsize=12)
    for ax in axes:
        ax.spines[["left", "bottom"]].set_color("#94A2AA")
        ax.tick_params(length=3)
        ax.grid(axis="y", color="#D9E0E4", linewidth=.7)
        ax.set_axisbelow(True)
        ax.set_xlabel("Highest harmonic included, $K$")

    ax = axes[0]
    for branch, color, sign in zip(branches, [BLUE, ORANGE], ["+0.5", "−0.5"]):
        for key, line, marker, label in [
            ("runs", "-", "o", "Matched"), ("uncalibrated", "--", "s", "Bare")
        ]:
            x, y = series(matched["branches"][branch][key], .1, 512)
            ax.plot(x, y, line, color=color, marker=marker, ms=3.8, lw=1.5,
                    mfc=color if key == "runs" else "white", label=f"{label}, {sign}")
    ax.set_title(r"Weak coupling   $\delta=0.1$", loc="left", fontweight="bold", fontsize=12)
    ax.set_ylabel("Photon share (%)")
    ax.set_xscale("log", base=2)
    ax.set_xticks([48, 96, 192, 384, 512], ["48", "96", "192", "384", "512"])
    ax.set_ylim(0, 1.28)
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax.legend(loc="lower left", ncols=2, fontsize=8.5, frameon=False,
              columnspacing=1, handlelength=1.8)

    for ax, branch, sign in zip(axes[1:], branches, ["+0.5", "−0.5"]):
        ax.set_title(r"Strong   $\delta=3$   ·   $\phi_m=$" + sign,
                     loc="left", fontweight="bold", fontsize=12)
        ax.set_xlim(47, 73)
        ax.set_xticks([48, 52, 56, 60, 64, 68, 72])
        ax.set_ylim(0, 100)
        ax.axhline(50, color="#84929C", lw=.9, ls=":", zorder=1)
        for runs, color, marker, label in [
            (matched["branches"][branch]["runs"], BLUE, "o", "Peaked"),
            (reach["branches"][branch]["tapered_profile_runs"], ORANGE, "s", "Tapered")
        ]:
            x, y = series(runs, 3, 72)
            ax.plot(x, y, color=color, marker=marker, ms=4, lw=1.5, label=label)
            first_stop = x[np.isnan(y)][0]
            ax.plot(first_stop, .966, marker="v", color=color, ms=6.5,
                    transform=ax.get_xaxis_transform())
        ax.legend(loc="lower left", ncols=2, frameon=False, fontsize=9)
    fig.text(.06, .128, "Weak panel colors identify the scalar minimum. Each point gives the photon percentage of total scalar-plus-photon energy.",
             fontsize=9.6)
    fig.text(.06, .072, "Strong panels show the two switch profiles. Triangles mark the first tested cutoff stopped by the 0.05 kinetic reserve. "
             "Those runs have no endpoint value.", fontsize=9.6)
    fig.savefig(ROOT / "response-matching.png", dpi=240, bbox_inches="tight", pad_inches=.17,
                metadata={"Title": "Short-wave response matching in two closed-universe model branches"})
    plt.close(fig)


if __name__ == "__main__":
    main()
