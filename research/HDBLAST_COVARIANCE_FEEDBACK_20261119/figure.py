"""Plot the evolving quiet shell, weak-probe pair count and fixed-wall force sensitivity."""
import json
from pathlib import Path
from zipfile import ZipFile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent
STATES = ("erase", "reflect", "kick")
LABELS = ("Erase", "Reflect", "Impulse")
COLORS = ("#136278", "#a53e1b", "#745695")
INK = "#1a303b"
GRID = "#d8e0e1"


def axes_style(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["bottom", "left"]].set_color("#8d9da1")
    ax.tick_params(colors=INK, labelsize=10, length=3)
    ax.yaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)


def expansion(ax):
    with ZipFile(ROOT / "native" / "WALLS.zip") as archive:
        data = {}
        for state in ("retain", *STATES):
            with archive.open(f"fine/{state}/quiet/wall.csv") as source:
                data[state] = np.genfromtxt(source, names=True, delimiter=",")
    clock = np.linspace(max(row["H0tau"][0] for row in data.values()), 2.10, 1201)

    def interpolate(state):
        row = data[state]
        if row["H0tau"][-1] < clock[-1]:
            raise ValueError("Readout exceeds the archived evolution")
        return np.interp(clock, row["H0tau"], row["H_over_H0"])

    baseline = interpolate("retain")
    for state, label, color in zip(STATES, LABELS, COLORS):
        ax.plot(clock, 1e6*(interpolate(state)-baseline), color=color, label=label, lw=1.9)
    ax.set(xlim=(clock[0], clock[-1]), xticks=[1.4, 1.6, 1.8, 2.0],
           xlabel=r"Proper time $H_0\tau$", ylabel=r"Change in $H/H_0$  ($\times 10^{-6}$)")
    ax.legend(frameon=False, loc="upper left", ncol=1, fontsize=9.5, handlelength=1.8)


def pair_count(ax):
    report = json.loads((ROOT/"emission"/"EMISSION.json").read_bytes())
    walls = report["grids"]["fine"]["walls"]

    def total(geometry, state):
        return walls[geometry]["quiet"]["states"][state]["wide"]["full"]["number"]["total"]

    reference = total("retain", "retain")
    direct, geometry = [], []
    for state in STATES:
        fixed = total("retain", state)
        returned = total(state, state)
        direct.append(100*(fixed-reference)/reference)
        geometry.append(100*(returned-fixed)/reference)
    if min(direct) <= 0 or max(geometry) >= 0:
        raise ValueError("Figure legend sign differs from the calculated decomposition")
    x = np.arange(len(STATES))
    ax.vlines(x, -np.array(geometry), direct, color="#b9c5c8", lw=1.5)
    ax.scatter(x, direct, label="Initial state  (+)", marker="o", s=50, color="#2a6378", zorder=3)
    ax.scatter(x, -np.array(geometry), label="Wall response  (−)", marker="D", s=40,
               color="#90542a", zorder=3)
    ax.set_yscale("log")
    ax.set(xticks=x, xticklabels=LABELS, xlim=(-.45, 2.45), ylim=(0.00055, 2.7),
           xlabel="Scalar preparation", ylabel="Absolute change (% of retained count)")
    ax.set_yticks([.001, .01, .1, 1], ["0.001", "0.01", "0.1", "1"])
    ax.tick_params(axis="y", which="minor", length=0)
    ax.legend(loc="center", bbox_to_anchor=(.5, .52), frameon=True, framealpha=1,
              facecolor="white", edgecolor="none", fontsize=9.3, handletextpad=.35)


def conditioning(ax):
    report = json.loads((ROOT/"moments"/"CONDITIONING.json").read_bytes())
    finest = max(report["fine"], key=lambda row: row["momentum_nodes"])
    gains = np.array(finest["first_twelve_ratios"])
    x = np.arange(1, gains.size+1)
    if np.sum(gains > 1e-2) != finest["above_one_percent"] or not np.all(np.diff(gains)<0):
        raise ValueError("Unexpected singular spectrum")
    ax.plot(x, gains, marker="o", ms=4.0, color="#276b6c", lw=1.7)
    ax.axhline(1e-2, color="#7d7062", ls=(0, (3, 3)), lw=1.15)
    ax.text(7.0, 2.5e-2, "1% of leading gain", fontsize=9.1, ha="center", color="#695c50")
    ax.set_yscale("log")
    ax.set(xlim=(.6, 12.4), ylim=(5e-15, 3), xticks=[1, 3, 6, 9, 12],
           xlabel="Singular-value index", ylabel="Force gain / leading gain")
    ax.set_yticks([1, 1e-4, 1e-8, 1e-12])
    ax.tick_params(axis="y", which="minor", length=0)


def main():
    plt.rcParams.update({"font.size":11, "font.family":"DejaVu Sans", "text.color":INK,
                         "axes.labelcolor":INK, "axes.labelsize":10.5})
    fig, axes = plt.subplots(1,3,figsize=(13.0,5.8))
    fig.subplots_adjust(left=.075, right=.985, top=.735, bottom=.275, wspace=.45)
    fig.text(.055,.94,"How scalar correlations affect the rolling shell",fontsize=22.5,
             fontweight="bold",va="top")
    fig.text(.055,.85,"The preparations share the same initial energy density, pressure and scalar force.",
             fontsize=12.0,va="top")
    for ax,title,draw in zip(axes,["Expansion on the quiet shell","Change in quiet pair yield",
                                   "Sensitivity of the force record"],
                             [expansion,pair_count,conditioning]):
        axes_style(ax)
        ax.set_title(title,loc="left",fontsize=12.1,fontweight="bold",pad=12)
        draw(ax)
    fig.text(.055,.145,"The fine-mesh histories are compared with the quiet shell retaining its original state. Counts include the switched vacuum.",fontsize=10.0)
    fig.text(.055,.105,r"The probe uses $1.35\leq H_0\tau\leq 2.05$ and total pair momentum $q\leq8$. "
             "Each diamond holds the preparation fixed and changes its wall.",fontsize=10.0)
    fig.text(.055,.045,"The force sensitivity uses the common reference geometry, 96 momentum nodes and 1,025 proper-time samples. "
             "Only three directions exceed 1% of peak sensitivity.",fontsize=10.0)
    fig.savefig(ROOT/"feedback.png",dpi=175,facecolor="white",metadata={"Software":"HDblast"})
    plt.close(fig)


if __name__ == "__main__":
    main()
