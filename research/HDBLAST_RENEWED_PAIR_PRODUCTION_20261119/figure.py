"""Draw the quiet-shell energy and number comparison with its waveform-independent bound."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent


def render(path):
    q=json.loads((ROOT/"emission/EMISSION.json").read_text())["arms"]["quiet"]["cutoffs"]
    limit=json.loads((ROOT/"criteria/CRITERIA.json").read_text())["histories"]["quiet"]["bounds"]["whole_post_gate"]["continuum"]["lowest_upper"]
    cut=np.array(sorted(int(s) for s in q))
    factor=100*np.sqrt(.1)
    fig,(a,b)=plt.subplots(1,2,figsize=(10.1,4.0),layout="constrained",width_ratios=[1.1,1])
    for key,label,color,mrk in (("pair_count","Number of pairs","#2b617e","o"),
             ("both_sign_conformal_energy","Energy in pairs","#b35235","s")):
        y=[]
        for k in cut:
            p=q[str(k)]["profiles"]
            y.append(p["A+B"][key]/p["A"][key])
        a.plot(cut,y,color=color,marker=mrk,label=label,lw=2)
    a.set_xscale("log",base=2)
    a.set_xticks(cut,labels=[str(k) for k in cut])
    a.set(xlabel="Highest modeled pair frequency",ylabel="Two pulses / original pulse",title="What the second pulse adds")
    a.axhline(1,color="0.6",ls=":",lw=1)
    a.legend(frameon=False)
    vals=[q["120"]["profiles"]["A"]["equilibrium_number_ratio_per_sqrt_y"]*factor,
          q["120"]["profiles"]["A+B"]["equilibrium_number_ratio_per_sqrt_y"]*factor,
          q["1920"]["profiles"]["A+B"]["equilibrium_number_ratio_per_sqrt_y"]*factor,
          100*limit]
    labels=["Original · 120","Two pulses · 120","Two pulses · 1920","Any allowed late pulse"]
    b.barh(np.arange(3),vals[:3],color=["#638394","#2b617e","#2b617e"],height=.58)
    b.hlines(3,0,vals[-1],color="#aa552d",lw=2)
    b.plot(vals[-1],3,"|",color="#aa552d",ms=15,mew=2)
    for i,v in enumerate(vals):
        b.text(v+0.6,i, ("≤ " if i==3 else "")+f"{v:.2f}%",va="center",fontsize=9)
    b.set_yticks(range(4),labels=labels)
    b.invert_yaxis()
    b.set(xlim=(0,32),xlabel="Percent of the equilibrium particle count",title="Number at the energy actually produced")
    for ax in (a,b):
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(axis="y" if ax is a else "x",alpha=.18)
        ax.set_axisbelow(True)
    fig.savefig(path,dpi=170)
    plt.close(fig)


if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output",type=Path,required=True)
    render(cli.parse_args().output)
