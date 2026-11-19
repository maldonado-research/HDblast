"""Plot the gravitational response at a common proper time across the finite calibrations."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline

HERE = Path(__file__).resolve().parent


def curve(runs, key, grid):
    with np.load(runs / f"{grid}_quiet_{key}.npz", allow_pickle=False) as record:
        return CubicSpline(record["s"],record["h"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE/"gravity.png")
    args = parser.parse_args()
    s = np.linspace(1.34, 2.10, 501)
    fig, (ax, bars) = plt.subplots(1,2,figsize=(10.6,4.15),layout="constrained",
                                   gridspec_kw={"width_ratios":[1.4,1]})
    base = curve(args.runs, "zero","fine")(s)
    colors = ["#245d82","#b0572c","#317063"]
    ax.axvspan(1.35,2.05,color="#eeede8",zorder=0)
    for i,(scheme,label) in enumerate([("bare","Original interaction (κ = 0)"),
                                       ("static","Static calibration (κ = 1)")]):
        plus,minus = [curve(args.runs, f"120_{scheme}_{side}0.01","fine")(s) for side in ("p","m")]
        ax.plot(s, 1e6*(plus-minus)/(2*base), lw=2.3, label=label, color=colors[i])
    ax.set_title("Response while the interaction acts",loc="left",weight="bold")
    ax.set_xlabel(r"Shell proper time $H_0\tau$")
    ax.set_ylabel("Change in expansion rate (parts per million)")
    ax.legend(frameon=False, loc="lower left", fontsize=8)
    with (HERE / "bulk/READOUTS.csv").open() as f:
        rows=list(csv.DictReader(f))
    result=json.loads((HERE/"calibration/CALIBRATION.json").read_text())
    for i,(grid,color) in enumerate([("fine","#245d82"),("coarse","#b4ccd5")]):
        h0=float(next(row["zero_h"] for row in rows if row["grid"]==grid and row["arm"]=="quiet"))
        response=np.array(result["grids"][grid]["response"]["0.01"]["dh_dlambda"])*1e4/h0
        bars.bar(np.arange(3)+(i-.5)*.36,response,.35,color=color,label=grid.capitalize()+" grid",zorder=2)
    bars.set_xticks(range(3),["κ = 0","κ = 1","κ = 3"])
    bars.set_title("Response at readout",loc="left",weight="bold")
    bars.set_xlabel("Strength of the static compensation term")
    bars.set_ylabel("Change in expansion rate (parts per million)")
    bars.legend(frameon=False, loc="lower right",fontsize=8)
    for panel in (ax,bars):
        panel.axhline(0,color="#444",lw=.8)
        panel.grid(axis="y",color="#ddd",lw=.6)
        panel.set_axisbelow(True)
        panel.spines[["top","right"]].set_visible(False)
    fig.savefig(args.output,dpi=170,facecolor="white")


if __name__=="__main__":
    main()
