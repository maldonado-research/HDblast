"""Render the signed and counted energy histories near the first mass crossing."""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

HERE = Path(__file__).resolve().parent


def main():
    with (HERE/"evolution/QUIET48.csv").open() as f:
        rows = [row for row in csv.DictReader(f) if 1.30 <= float(row["s"]) <= 1.56]
    s = [float(row["s"]) for row in rows]
    fig, ax = plt.subplots(figsize=(7.5, 4.3), layout="constrained")
    ax.axhspan(-.125, 0, color="#f9eae3", zorder=0)
    ax.axhline(0, lw=.9, color="#414246")
    for column, label, color, style in (
        ("sheared_prepared_count_a4", "Radiation assigned from the particle count", "#14687a", "-"),
        ("sheared_prepared_signed_a4", "Energy removed from the physical field", "#a14228", "-"),
        ("physical_signed_a4", "Original decay model", "#696d76", "--"),
    ):
        ax.plot(s, [float(r[column]) for r in rows], color=color, ls=style, lw=2.1, label=label)
    ax.annotate("−0.0983", xy=(1.541,-.0983), xytext=(1.493,-.105),
                arrowprops=dict(arrowstyle="->", color="#76391f"), color="#76391f")
    ax.set(xlim=(1.30, 1.56), ylim=(-.125,.50), xlabel=r"Shell proper time $H_0\tau$",
           ylabel=r"Comoving radiation energy $a^4\rho_r/H_0^4$",
           )
    ax.set_title("Energy near the mass crossing on the quiet shell", loc="left", fontweight="bold", pad=12)
    ax.yaxis.set_major_locator(MultipleLocator(.1))
    ax.grid(axis="y", alpha=.18)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    fig.savefig(HERE/"radiation_near_crossing.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
