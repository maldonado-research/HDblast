"""Draw the pulse-by-pulse bounds after matching the initial shell sources."""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parent
CASES = (("plus", "Nearby positive"), ("minus", "Nearby negative"),
         ("far", "Distant positive"))
NAVY = "#225b80"
RUST = "#a74b25"
INK = "#1e2c36"


def main():
    with (ROOT / "source_check" / "SOURCE_TRANSFER.json").open() as stream:
        results = json.load(stream)["128"]
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "text.color": INK, "axes.labelcolor": INK,
                         "xtick.color": INK, "ytick.color": INK})
    fig, ax = plt.subplots(figsize=(8.4, 5.3))
    fig.subplots_adjust(left=0.105, right=0.95, bottom=0.215, top=0.695)
    fig.text(0.105, 0.944, "Can coherence reverse the pulse response?",
             fontsize=17, weight="bold", va="top")
    fig.text(0.105, 0.868, "Largest coherence correction opposing the incoherent "
             "pulse–quiet pair-count difference,", fontsize=10.4)
    fig.text(0.105, 0.826, "as a fraction of that difference. Values above 1 permit reversal.",
             fontsize=10.4)
    handles = [Line2D([], [], ls="none", marker="o", color=NAVY, ms=8,
                      label="Vacuum-subtracted excess"),
               Line2D([], [], ls="none", marker="D", color=RUST, ms=7,
                      label="Before vacuum subtraction")]
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.093, 0.80),
               frameon=False, fontsize=10, ncol=2, columnspacing=2)
    ax.axhline(1.0, lw=1.4, color="#36424c", ls=(0, (5, 3)), zorder=1)
    for x, (arm, _) in enumerate(CASES):
        result = results[arm]
        excess = result["quantum_bound"]
        absolute = result["quantum_absolute_bound"]
        if abs(excess - result["quantum_witness"]) > 1e-7:
            raise ValueError(f"Unresolved optimality gap for {arm}")
        ax.plot([x, x], [absolute, excess], lw=1.8, color="#c0c8cc", zorder=1)
        ax.plot(x, excess, marker="o", ms=8.5, color=NAVY, zorder=3)
        ax.plot(x, absolute, marker="D", ms=7.0, color=RUST, zorder=3)
        if arm == "far":
            ax.annotate(f"{excess:.4f}\n{1 - excess:.4f} below threshold",
                        (x, excess), xytext=(0, -12), textcoords="offset points",
                        ha="center", va="top", color=NAVY, fontsize=9.5,
                        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.8})
        else:
            ax.annotate(f"{excess:.4f}", (x, excess), xytext=(0, 8),
                        textcoords="offset points", ha="center", color=NAVY)
        ax.annotate(f"{absolute:.4f}", (x, absolute), xytext=(0, -11),
                    textcoords="offset points", va="top", ha="center", color=RUST)
    ax.set_ylim(0.70, 1.164)
    ax.set_xlim(-0.53, 2.55)
    ax.set_xticks(range(len(CASES)), [label for _, label in CASES])
    ax.set_yticks([0.75, 0.85, 0.95, 1.00, 1.05, 1.15],
                  ["0.75", "0.85", "0.95", "1.00", "1.05", "1.15"])
    ax.set_ylabel("Reversal factor (dimensionless)", fontsize=10)
    ax.spines[["top", "left", "right"]].set_visible(False)
    ax.spines["bottom"].set_color("#a7afb5")
    ax.grid(axis="y", color="#e1e5e7", lw=0.7)
    ax.set_axisbelow(True)
    ax.tick_params(length=0, pad=6, labelsize=9.5)
    fig.text(0.105, 0.095, "Each pulse is optimized separately at matched initial energy "
             "density, pressure and scalar source.", fontsize=9.1)
    fig.text(0.105, 0.05, "Wide interaction window, full daughter spectrum, 128 parent "
             "momentum intervals.", fontsize=9.1, color="#525f65")
    fig.savefig(ROOT / "matched_bounds.png", dpi=200, facecolor="white",
                metadata={"Title": "Pair-count reversal bounds with matched initial sources"})
    plt.close(fig)


if __name__ == "__main__":
    main()
