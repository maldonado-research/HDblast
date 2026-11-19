"""Plot the matched shell rates and the independently calculated delayed response."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--repo",type=Path,default=HERE.parents[1])
    cli.add_argument("--output",type=Path,default=HERE/"matching.png")
    args=cli.parse_args()
    gravity=json.loads((HERE/"gravity/GRAVITY.json").read_text())
    pulse=json.loads((HERE/"memory/PULSES.json").read_text())
    historical=args.repo/"research/HDBLAST_SCALAR_CLOCK_GRAVITY_20261119/bulk/READOUTS.csv"
    with historical.open() as stream:
        base={r["grid"]:float(r["zero_h"]) for r in csv.DictReader(stream) if r["arm"]=="quiet"}
    fig,(left,right)=plt.subplots(1,2,figsize=(9.2,3.9),layout="constrained")
    x=np.arange(2)
    for i,(grid,color) in enumerate([("coarse","#b5cddb"),("fine","#245d82")]):
        result=gravity["grids"][grid]["steps"]["0.005"]["observed"]
        values=[1e4*result[k]["h"]/base[grid] for k in ("60","120")]
        left.bar(x+(i-.5)*.32,values,.32,color=color,label=grid.capitalize(),zorder=3)
    left.set_title("Expansion with matched static response",loc="left",size=12,weight="bold")
    left.set_xticks(x,["60","120"])
    left.axhline(0,color="#555",lw=.7)
    left.set_xlabel("Upper energy cutoff (comoving units)")
    left.set_ylabel("Expansion-rate change (parts per million)")
    left.legend(frameon=False,loc="lower right",fontsize=9)
    right.set_title("Delayed force after two positive pulses",loc="left",size=12,weight="bold")
    rows=pulse["rows"]
    phases=[r["n"] for r in rows]
    right.loglog(phases,[r["half_to_full"] for r in rows],"o-",color="#8c4825",label="Calculated force ratio",lw=2)
    right.loglog(phases,[pulse["separation_half"]/r["v"] for r in rows],"--",color="#32594d",label="Late-time prediction")
    right.set_xticks(phases,phases)
    ticks=[.0078125,.015625,.03125,.0625,.125]
    right.set_yticks(ticks,["1/128","1/64","1/32","1/16","1/8"])
    right.minorticks_off()
    right.set_xlabel("Number of threshold periods (phase index)")
    right.set_ylabel("Half-period / full-period response")
    right.legend(frameon=False,loc="upper right",fontsize=8)
    for ax in (left,right):
        ax.grid(axis="y",color="#ddd",lw=.6)
        ax.set_axisbelow(True)
        ax.spines[["top","right"]].set_visible(False)
    fig.savefig(args.output,dpi=170,facecolor="white")


if __name__=="__main__":
    main()
