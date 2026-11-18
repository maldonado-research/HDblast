"""Plot cumulative particle-source radiation and the numerical balance residual."""
from pathlib import Path
import csv
import itertools

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent


def curve(name):
    with (ROOT / "data" / f"b1_{name}.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    times = [float(row["x"]) for row in rows]
    radiation = [float(row["weighted_radiation"]) for row in rows]
    source = [float(row["weighted_source"]) for row in rows]
    areas = [(t1-t0)*(q0+q1)/2 for (t0,t1),(q0,q1)
             in zip(itertools.pairwise(times), itertools.pairwise(source))]
    gain = [r - radiation[0] for r in radiation]
    integrated = list(itertools.accumulate(areas, initial=0.0))
    residual = [1e5 * (value - r) for value, r in zip(integrated, gain)]
    return times, gain, integrated, residual


def plot():
    fine = curve("fine")
    coarse = curve("coarse")
    fig, (signal, balance) = plt.subplots(2, 1, sharex=True, figsize=(7.2, 5.3),
                                         height_ratios=[1.8, 1.0], layout="constrained")
    signal.plot(fine[0], fine[2], lw=2.4, color="#235e73", label="Integrated particle source")
    signal.scatter(fine[0][::7], fine[1][::7], facecolors="white", edgecolors="#a04e2a",
                   s=21, linewidths=1.1, label="Evolved radiation gain", zorder=3)
    signal.set(ylabel="Comoving radiation gain\n(model units)", title="Radiation from the modeled particle decay")
    signal.legend(loc="upper left", frameon=False)
    balance.axhline(0.0, lw=0.8, color="0.4")
    for label, data, style, color in (("Fine radial grid",fine,"-","#235e73"),
                                      ("Coarse radial grid",coarse,"--","#a04e2a")):
        balance.plot(data[0], data[3], style, color=color, lw=1.6, label=label)
    balance.set(ylabel="Source minus gain\n($10^{-5}$ model units)",
                xlabel="Elapsed model time   $s=H_0\\tau$", xlim=(1,4))
    balance.legend(loc="upper right", frameon=False, fontsize=9)
    for ax in (signal, balance):
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(axis="y", color="0.88", linewidth=0.7)
    fig.savefig(ROOT / "b1_transfer.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    plot()
