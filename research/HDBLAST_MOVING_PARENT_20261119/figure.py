"""Plot the common-phase pair-count response using the archived refined results."""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FixedLocator, FuncFormatter, MaxNLocator
import numpy as np


ROOT = Path(__file__).resolve().parent
CASES = [
    ("plus", "Nearby positive pulse", "#215b82"),
    ("minus", "Nearby negative pulse", "#a44c27"),
    ("far", "Distant positive pulse", "#266a5b"),
]
EXAMPLE = -0.5
SHADE = "#e8e0f4"


def phase_curve(data):
    """Reconstruct the exact cosine from the extremal phase and allowed range."""
    low = data["constant_phase_min_percent"]
    high = data["constant_phase_max_percent"]
    center = data["incoherent"]["difference_percent_of_incoherent_quiet"]
    radius = (high - low) / 2
    peak = data["maximizing_constant_phase_radians"]
    if not np.isclose((high + low) / 2, center, rtol=0, atol=1e-11):
        raise ValueError("Coherent extrema disagree with the incoherent baseline")

    def evaluate(phase):
        return center + radius * np.cos(phase - peak)

    for phase, record in data["constant_phases"].items():
        observed = record["difference_percent_of_incoherent_quiet"]
        if not np.isclose(evaluate(float(phase)), observed, rtol=0, atol=1e-10):
            raise ValueError("Reconstructed curve disagrees with a recorded phase")

    # The two roots closest to the example on either side delimit its reversal.
    if not abs(center) < radius or evaluate(EXAMPLE) * center >= 0:
        raise ValueError("Expected a sign reversal at the marked phase")
    angle = np.arccos(-center / radius)
    roots = [peak + sign * angle + period * 2 * np.pi
             for sign in (-1, 1) for period in (-1, 0, 1)]
    interval = (max(t for t in roots if t < EXAMPLE),
                min(t for t in roots if t > EXAMPLE))
    return evaluate, center, interval


def percent_tick(value, _position):
    return f"{value:.3f}".rstrip("0").rstrip(".").replace("-", "−")


def main():
    with (ROOT / "coherence" / "REFINED.json").open() as stream:
        fine = json.load(stream)["results"]["fine"]
    curves = [phase_curve(fine[key]["wide"]["number"]) for key, _, _ in CASES]
    left = max(interval[0] for _, _, interval in curves)
    right = min(interval[1] for _, _, interval in curves)
    if not left < EXAMPLE < right:
        raise ValueError("No common reversal interval at the marked phase")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.spines.left": False, "axes.edgecolor": "#929aa3",
                         "text.color": "#202b35", "axes.labelcolor": "#202b35",
                         "xtick.color": "#394651", "ytick.color": "#394651"})
    fig, axs = plt.subplots(3, 2, figsize=(9.5, 6.9), sharex="col",
                            gridspec_kw={"width_ratios": (1.46, 1)})
    fig.subplots_adjust(left=0.115, right=0.969, top=0.790, bottom=0.145,
                        hspace=0.53, wspace=0.26)
    fig.text(0.115, 0.964, "A common squeezing phase can reverse all three responses",
             fontsize=15, weight="bold", va="top")
    fig.text(0.115, 0.918, "Pulse minus quiet pair yield, relative to the incoherent quiet yield (%)",
             fontsize=11, va="top")
    legend = [Line2D([], [], color="#555b65", lw=1.4, ls=(0, (4, 3)),
                     label="Incoherent preparation"),
              Patch(facecolor=SHADE, edgecolor="none", label="All three signs reversed"),
              Line2D([], [], color="#303945", ls="none", marker="o", ms=4.5,
                     label="Example at −0.5 rad")]
    fig.legend(handles=legend, frameon=False, ncol=3, loc="upper left",
               bbox_to_anchor=(0.106, 0.893), columnspacing=2.0, fontsize=9)
    fig.text(0.115, 0.821, "One full phase cycle", fontsize=11, weight="bold")
    fig.text(axs[0, 1].get_position().x0, 0.821, "Detail near the reversal",
             fontsize=11, weight="bold")
    limits = [(-np.pi, np.pi), (-1.1, -0.1)]
    for row, ((_, name, color), (curve, baseline, _)) in enumerate(zip(CASES, curves)):
        for col, (lower, upper) in enumerate(limits):
            ax = axs[row, col]
            phase = np.linspace(lower, upper, 801)
            values = curve(phase)
            ax.axvspan(left, right, color=SHADE, lw=0)
            ax.axhline(0, lw=0.9, color="#515b64")
            if col == 0:
                ax.axhline(baseline, lw=1.3, color="#555b65", ls=(0, (4, 3)))
                ax.set_title(name, loc="left", fontsize=10.5, pad=4, color=color)
            ax.plot(phase, values, lw=2.15, color=color)
            value = curve(EXAMPLE)
            ax.plot(EXAMPLE, value, "o", color=color, ms=5.0,
                    markeredgecolor="white", markeredgewidth=0.8, zorder=5)
            if col == 1:
                ax.text(0.04, 0.94, f"At −0.5 rad, {value:+.4f}%".replace("-", "−"),
                        transform=ax.transAxes, ha="left", va="top", fontsize=8.5,
                        color=color, bbox={"facecolor": "white", "edgecolor": "none",
                                           "alpha": 0.92, "pad": 1.5})
            span = values.max() - values.min()
            ax.set_ylim(values.min() - 0.15 * span, values.max() + 0.18 * span)
            ax.set_xlim(lower, upper)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4, steps=(1, 2, 2.5, 5, 10)))
            ax.yaxis.set_major_formatter(FuncFormatter(percent_tick))
            ax.tick_params(length=0, labelsize=9, pad=5)
            ax.grid(axis="y", color="#a7b0b7", alpha=0.30, lw=0.65)
            ax.set_axisbelow(True)
    axs[-1, 0].set_xticks(np.arange(-2, 3) * np.pi / 2,
                           ["−π", "−π/2", "0", "π/2", "π"])
    axs[-1, 1].xaxis.set_major_locator(FixedLocator([-1.0, -0.75, -0.5, -0.25]))
    axs[-1, 1].xaxis.set_major_formatter(FuncFormatter(percent_tick))
    for ax in axs[-1]:
        ax.set_xlabel("Common squeezing phase (rad)", labelpad=8)
    fig.text(0.115, 0.04,
             "Fine shell mesh, wide interaction window, full daughter spectrum. "
             "Unit occupation for 0 ≤ q ≤ 8.", fontsize=9, color="#47535d")
    fig.savefig(ROOT / "phase_response.png", dpi=200, facecolor="white",
                metadata={"Title": "Moving-parent pair yield versus squeezing phase"})
    plt.close(fig)


if __name__ == "__main__":
    main()
