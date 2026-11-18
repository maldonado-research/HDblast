"""Draw the detector and stationary-pulse results from capsule JSON files."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter


INK = "#202B34"
MUTED = "#4B5862"
GRID = "#D6DDE2"
FULL = "#246088"
HALF = "#97492E"


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def format_axis(axis):
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.spines["bottom"].set_color(GRID)
    axis.tick_params(colors=MUTED, length=0, pad=6)
    axis.set_axisbelow(True)


def plot_detector(axis, replay):
    entries = [
        ("plus", "Full positive", FULL),
        ("half_plus", "Half positive", HALF),
        ("half_minus", "Half negative", HALF),
        ("minus", "Full negative", FULL),
    ]
    rows = list(range(len(entries)))
    for row, (key, _label, color) in enumerate(entries):
        value = 100 * replay[key]["relative_difference"]
        axis.barh(row, value, height=0.39, color=color, zorder=2)
        offset = 0.012 if value >= 0 else -0.012
        label = f"{value:+.3f}%".replace("-", "−")
        axis.text(value + offset, row, label, ha="left" if value >= 0 else "right",
                  va="center", color=color, fontsize=10, weight="semibold")
    axis.set_yticks(rows, [label for _key, label, _color in entries])
    axis.set_ylim(len(rows) - 0.5, -0.5)
    axis.set_xlim(-0.47, 0.47)
    axis.set_xticks([-0.4, -0.2, 0, 0.2, 0.4], ["−0.4", "−0.2", "0", "0.2", "0.4"])
    axis.set_xlabel("Change from quiet response (%)", labelpad=11)
    axis.grid(axis="x", color=GRID, linewidth=0.65)
    axis.axvline(0, color=INK, linewidth=0.85, zorder=3)
    axis.set_title("Switched detector on the shell", loc="left", pad=13)


def plot_scattering(axis, protocols, bounds):
    for amplitude, color in [("full", FULL), ("half", HALF)]:
        plus_key = f"g8_{amplitude}_positive"
        minus_key = f"g8_{amplitude}_negative"
        plus = protocols[plus_key]["spectrum"]
        minus = protocols[minus_key]["spectrum"]
        # The union gives a conservative enclosure for either sign at any k.
        lo = min(float(bounds[key]["occupation_lower"]) for key in (plus_key, minus_key))
        hi = max(float(bounds[key]["occupation_upper"]) for key in (plus_key, minus_key))
        wave_numbers = [k for k, _value in plus]
        axis.fill_between([min(wave_numbers), max(wave_numbers)], 1e6 * lo, 1e6 * hi,
                          color=color, alpha=0.11, linewidth=0)
        axis.plot(wave_numbers, [1e6 * value for _k, value in plus], color=color,
                  linewidth=1.9)
        # Sampled markers keep both signs visible where their values nearly coincide.
        samples = minus[::4]
        axis.plot([k for k, _value in samples], [1e6 * value for _k, value in samples],
                  linestyle="none", marker="o", markersize=3.7, markerfacecolor="white",
                  markeredgecolor=color, markeredgewidth=1.0, clip_on=False)
        axis.text(5.4, 1e6 * hi / 1.24, f"{amplitude.capitalize()} amplitude", color=color,
                  fontsize=10, weight="semibold", va="center")

    axis.set_yscale("log")
    axis.set_ylim(0.47, 8.5)
    axis.set_xlim(4.9, 18.15)
    axis.set_xticks([5, 8, 11, 14, 18])
    axis.yaxis.set_major_locator(FixedLocator([0.5, 1, 2, 4, 8]))
    axis.yaxis.set_major_formatter(FuncFormatter(lambda x, _pos: f"{x:g}"))
    axis.yaxis.set_minor_formatter(NullFormatter())
    axis.tick_params(axis="y", which="minor", length=0)
    axis.grid(axis="y", color=GRID, linewidth=0.6)
    axis.set_xlabel(r"Comoving wave number $k/H_0$", labelpad=11)
    axis.set_ylabel(r"Out occupation $|\beta_k|^2$ ($\times 10^{-6}$)", labelpad=8)
    axis.set_title("Stationary mass pulse", loc="left", pad=13)
    keys = [Line2D([], [], color=INK, lw=1.7, label="Positive pulse"),
            Line2D([], [], marker="o", color=INK, markerfacecolor="white",
                   ms=4, lw=0, label="Negative pulse")]
    axis.legend(handles=keys, loc="upper right", bbox_to_anchor=(0.997, 0.961),
                frameon=False, ncol=2, fontsize=9, columnspacing=1.0,
                handletextpad=0.55, handlelength=1.5)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capsule", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output", type=Path)
    options = parser.parse_args()
    root = options.capsule
    replay = read_json(root / "detector" / "REPLAY.json")["s1200_k64"]
    protocols = read_json(root / "scattering" / "SCATTER.json")["protocols"]
    bounds = read_json(root / "scattering" / "BOUNDS.json")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5,
                         "axes.titleweight": "bold", "axes.titlesize": 14,
                         "axes.labelcolor": INK, "text.color": INK})
    figure, (left, right) = plt.subplots(1, 2, figsize=(11.8, 4.35))
    figure.subplots_adjust(left=0.116, right=0.985, bottom=0.23, top=0.82, wspace=0.34)
    for axis in (left, right):
        format_axis(axis)
    plot_detector(left, replay)
    plot_scattering(right, protocols, bounds)
    figure.text(left.get_position().x0, 0.05, r"Fixed detector and switch, $5 \leq k/H_0 \leq 18$.",
                fontsize=9.5, color=MUTED)
    figure.text(right.get_position().x0, 0.05, "Shading encloses the analytic bounds for both signs.",
                fontsize=9.5, color=MUTED)
    output = options.output or root / "particle_response.png"
    figure.savefig(output, dpi=200, facecolor="white")
    plt.close(figure)


if __name__ == "__main__":
    main()
