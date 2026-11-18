"""Plot particle-number contrasts on the rolling shell from the archived fine grid."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D


BLUE = "#246088"
RUST = "#97492E"
INK = "#202B34"
MUTED = "#4B5862"
GRID = "#D6DDE2"


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def populations(report):
    return (report["event_comoving_N_8"], report["event_leading_comoving_N_8"],
            report["ends"]["2.1"]["comoving_N_8"])


def differences(pulse, quiet):
    return [100 * (value / baseline - 1)
            for value, baseline in zip(populations(pulse), populations(quiet))]


def signed(value, digits):
    return f"{value:+.{digits}f}%".replace("-", "−")


def style(axis):
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.spines["bottom"].set_color(GRID)
    axis.tick_params(colors=MUTED, length=0, pad=6)
    axis.set_axisbelow(True)
    axis.grid(axis="x", color=GRID, linewidth=0.65)
    axis.axvline(0, color=INK, lw=0.9, zorder=1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capsule", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--classical", type=Path, help="Directory containing near-control modes/")
    parser.add_argument("--remote", type=Path, help="Directory containing far-control modes/")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    near = (args.classical or args.capsule / "classical") / "modes"
    far = (args.remote or args.capsule / "far") / "modes"
    quiet = read(near / "base-fine-quiet.json")
    positive = differences(read(near / "base-fine-0.001.json"), quiet)
    negative = differences(read(near / "base-fine--0.001.json"), quiet)
    isolated = differences(read(far / "plus-fine.json"), quiet)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5,
                         "axes.titleweight": "bold", "axes.titlesize": 13,
                         "axes.labelcolor": INK, "text.color": INK})
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.6, 4.25),
                                          gridspec_kw={"width_ratios": [1.4, 1]})
    figure.subplots_adjust(left=0.17, right=0.98, bottom=0.28, top=0.80, wspace=0.52)
    for axis in (left, right):
        style(axis)
    left.set_title("Particle estimates with nearby compensation", loc="left", pad=16)
    left.set_xlim(-5.5, 5.6)
    left.set_xticks([-4, -2, 0, 2, 4])
    left.set_ylim(2.52, -0.52)
    left.set_yticks([0, 1, 2], ["Event rule with\nderivative correction",
                                "Leading event rule", "Free-wave evolution"])
    left.set_xlabel("Change in comoving number from quiet (%)", labelpad=9)
    for values, shift, color, marker in ((positive, -0.13, BLUE, "o"),
                                          (negative, 0.13, RUST, "s")):
        for i, value in enumerate(values):
            left.scatter(value, i + shift, color=color, marker=marker, s=30, zorder=3)
            towards_left = value < 0
            left.annotate(signed(value, 2), (value, i + shift),
                          xytext=(-6 if towards_left else 6, 0), textcoords="offset points",
                          ha="right" if towards_left else "left", va="center",
                          fontsize=9.5, weight="semibold", color=color)
    keys = [Line2D([], [], color=BLUE, marker="o", lw=0, ms=5, label="Positive pulse"),
            Line2D([], [], color=RUST, marker="s", lw=0, ms=5, label="Negative pulse")]
    figure.legend(handles=keys, loc="upper left", bbox_to_anchor=(0.162, 1.008),
                  ncol=2, frameon=False, fontsize=9.7, columnspacing=1.5, handletextpad=0.25)

    right.set_title("Free waves, positive pulse", loc="left", pad=16)
    right.set_xlim(-0.32, 0.149)
    right.set_xticks([-0.3, -0.2, -0.1, 0, 0.1], ["−0.3", "−0.2", "−0.1", "0", "0.1"])
    right.set_ylim(1.65, -0.65)
    right.set_yticks([0, 1], ["Nearby\ncompensation", "Distant\ncompensation"])
    right.set_xlabel("Change in comoving number from quiet (%)", labelpad=9)
    for row, value in enumerate([positive[2], isolated[2]]):
        right.barh(row, value, height=0.27, color=BLUE, zorder=2)
        towards_left = value < 0
        right.annotate(signed(value, 3), (value, row), xytext=(-5 if towards_left else 5, 0),
                       textcoords="offset points", ha="right" if towards_left else "left",
                       va="center", fontsize=10, weight="semibold", color=BLUE)

    figure.text(0.17, 0.095,
                r"All estimates use $0\leq k/H_0\leq8$. The event rules use the crossing spectrum. Free modes are read at $s=2.10$.",
                fontsize=9.5, color=MUTED)
    figure.text(0.17, 0.032,
                "The nearby compensating scalar pulse can affect the late count. The distant pulse lies outside the readout’s causal past.",
                fontsize=9.5, color=MUTED)
    figure.savefig(args.output or args.capsule / "rolling_response.png", dpi=200, facecolor="white")
    plt.close(figure)


if __name__ == "__main__":
    main()
