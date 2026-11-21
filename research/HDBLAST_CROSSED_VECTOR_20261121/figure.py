"""Plot the diagonal probability bound and the energy split of the extra stress derivative."""
import argparse
import json
from pathlib import Path
import runpy

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def draw(root, destination):
    wave = runpy.run_path(str(root/"unitarity/partial_wave.py"))
    bound = json.loads((root/"unitarity/PARTIAL.json").read_bytes())["root_s_over_M"]
    rates = json.loads((root/"crossed-stress/STRESS.json").read_bytes())["runs"]
    fine = sorted((r for r in rates if "-p14-" in r["snapshot"]),key=lambda r:r["snapshot"])
    if len(fine)!=2: raise ValueError("Expected two fine snapshots")
    colors = ("#1e6382","#a9462c")
    fig, (left,right) = plt.subplots(1,2,figsize=(9.8,4.1),layout="constrained")
    energies = np.geomspace(.30,8.,400)
    magnitude = wave["partial_magnitude"]((wave["MASS"]*energies)**2)
    left.plot(energies,magnitude,color=colors[0],lw=2.3)
    left.axhline(.5,c=colors[1],ls="--",lw=1.4,label="Allowed real part ≤ 0.5 in magnitude")
    left.axvline(bound,c="#666666",ls=":",lw=1.1)
    left.scatter([bound],[.5],color=colors[1],s=28,zorder=3)
    left.annotate(f"{bound:.3f} × mediator mass",(bound,.5),xytext=(9,-22),textcoords="offset points")
    left.set(xlim=(.3,8.),ylim=(0,2.),xlabel="Pair energy / mediator mass",ylabel="Magnitude of real J = 0 amplitude",
             title="Where the probability check fails")
    left.legend(loc="upper left",frameon=False,fontsize=8.5)
    x = np.arange(len(fine))
    for index,(field,label) in enumerate((("below_cut","Below the threshold"),("above_or_equal_cut","Above the threshold"))):
        y = [-r[field]["crossed_minus_pole_derivative"]["mean"] for r in fine]
        bars = right.bar(x+(-.17 if index==0 else .17),y,width=.32,color=colors[index],label=label)
        for bar,val in zip(bars,y):
            right.annotate(f"{val:.3g}",(bar.get_x()+bar.get_width()/2,val),xytext=(0,4),textcoords="offset points",ha="center",fontsize=8.5)
    right.set_xticks(x,["6.67","7.05"])
    right.set_yscale("log")
    right.set(ylim=(.008,160),xlabel="Shell time",ylabel="Added rate of stress loss (per model time)",
              title="Which pairs drive the apparent gain?")
    right.legend(loc="upper center",fontsize=8.2,ncol=2,frameon=False)
    for ax in (left,right):
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(axis="y",lw=.5,alpha=.2);ax.set_axisbelow(True)
    fig.savefig(destination,dpi=175)
    plt.close(fig)


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root",type=Path,default=Path(__file__).resolve().parent)
    parser.add_argument("--output",type=Path,default=Path(__file__).with_name("crossed-transport.png"))
    args = parser.parse_args()
    draw(args.root,args.output)
