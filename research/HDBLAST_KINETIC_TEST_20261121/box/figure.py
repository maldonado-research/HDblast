"""Plot kinetic relaxation measured in the event-driven two-stream experiments."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np


def draw(archive, output):
    rows = json.loads(archive.read_text())["runs"]
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 10, "axes.spines.right": False,
                         "axes.spines.top": False, "axes.edgecolor": "#8d969a", "axes.labelcolor": "#283841"}):
        fig, axes = plt.subplots(1, 2, figsize=(9.1, 4.8), sharey=True)
        fig.subplots_adjust(left=.082, right=.985, bottom=.24, top=.70, wspace=.17)
        fig.text(.082, .94, "Do collisions erase the beam before its particles thermalize?",
                 fontsize=15, fontweight="bold", color="#19343a")
        fig.text(.082, .885, "Particles start with one energy in two opposite beams. Each curve averages six simulated boxes.",
                 fontsize=9.5, color="#405258")
        for ax, case, title in zip(axes, ("drift_scale","hot_energy"),
                                   ("Matched to the peak stream speed", "Matched to the late gas energy")):
            r=rows[case+"_mono"]
            t=np.array(r["times"])
            x=np.log2(1.+t)
            for key, label, color, style in [
                    ("charge_odd_momentum","Directed momentum remaining","#00666a","-"),
                    ("ks_energy_pooled","Largest gap in energy distribution","#ad5428","-")]:
                y=100*np.array(r["mean"][key])
                err=100*np.array(r["standard_error"][key])
                ax.plot(x,y, color=color,ls=style,label=label,lw=2.2)
                ax.fill_between(x,y-2*err,y+2*err,color=color,alpha=.18,lw=0)
            thermal=100*np.array(rows[case+"_thermal"]["mean"]["ks_energy_pooled"])
            ax.plot(x, thermal, "--", label="Gap for initially thermal gas", lw=1.15, color="#626774")
            ax.set_title(title,loc="left",fontsize=10,fontweight="bold",pad=9,color="#283841")
            ticks=np.array([0,1,4,8,16,32,64,128])
            ax.set_xticks(np.log2(1+ticks),[str(i) for i in ticks])
            ax.set_xlim(-.04,7.14)
            ax.set_ylim(-5,104)
            ax.grid(axis="y",alpha=.15)
        axes[0].set_ylabel("Percent")
        h,l=axes[0].get_legend_handles_labels()
        fig.legend(h,l,loc="center left",bbox_to_anchor=(.072,.81),ncol=3,
                   frameon=False,fontsize=8.5,columnspacing=1.25)
        fig.text(.54,.135,"Elapsed time in initial monoenergetic-beam collision times",ha="center",fontsize=10)
        fig.text(.082,.056,"The distribution gap is the largest difference in the fraction of particles below any energy. Shading spans two standard errors.",
                 fontsize=8.1,color="#405258")
        fig.savefig(output,dpi=180,facecolor="white",metadata={"Software":"Matplotlib"})
        plt.close(fig)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args=p.parse_args()
    draw(args.input,args.output)
