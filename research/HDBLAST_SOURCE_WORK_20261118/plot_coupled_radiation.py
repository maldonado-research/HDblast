"""Draw the matched daughter-radiation contrasts from the native B1 shell model."""
import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.ticker import FuncFormatter


INK = "#21323B"
MUTED = "#485C64"
GREEN = "#286B60"
PALE = "#ACCDC2"
GRID = "#D9E1E1"
CASES = ("plus", "minus", "far")
LABELS = ("Near, positive", "Near, negative", "Distant, positive")


def fmt(v, digits=4):
    return f"{v:+.{digits}f}".replace("-", "−")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capsule", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    d = json.loads((args.capsule/"coupled/COUPLED.json").read_bytes())
    if d["cut_s"]!=2.10:
        raise ValueError("Expected s=2.10")
    rows = []
    for case in CASES:
        contrast = d["contrasts"][f"fine/y1/{case}"]["rad_hat"]
        values = [contrast["coupled_percent"],contrast["passive_percent"]]
        for b, value in zip((1,0), values, strict=True):
            quiet = d["cases"][f"fine/y1_b{b}/quiet"]["rad_hat"]
            ref = d["cases"][f"fine/y1_b{b}/{case}"]["rad_hat"]
            if not math.isclose(100*(ref-quiet)/quiet, value, abs_tol=1e-12):
                raise ValueError("Unexpected native quiet control")
        frac = 100*d["cases"][f"fine/y1_b1/{case}"]["rad_fraction"]
        rows.append((*values,frac))
    fig, (bars,table) = plt.subplots(1,2,figsize=(10.25,4.45),gridspec_kw={"width_ratios":(2.5,1.05)})
    fig.subplots_adjust(left=.143, right=.99, top=.70, bottom=.281, wspace=.085)
    fig.text(.143,.91,"Daughter radiation on the rolling shell",fontsize=17,weight="bold",color=INK)
    fig.text(.143,.831,r"Pulse response at $H_0\tau=2.10$ in the band $0\leq k/(a_{*,0}H_0)\leq8$",
             fontsize=11,color=MUTED)
    bars.set_axisbelow(True)
    bars.set_ylim(2.50,-.55)
    bars.set_xlim(-.32,1.15)
    bars.set_yticks(range(3),LABELS)
    bars.set_xticks([-.25,0,.25,.5,.75,1])
    bars.xaxis.set_major_formatter(FuncFormatter(lambda v,p: f"{v:g}".replace("-", "−")))
    bars.tick_params(length=0,pad=6, labelcolor=INK)
    bars.set_xlabel("Change in daughter density relative to quiet (%)",color=INK,labelpad=8,fontsize=10)
    bars.grid(axis="x",color=GRID,linewidth=.75)
    bars.axvline(0,color=MUTED,lw=1.1)
    for side in ("top","left","right"):
        bars.spines[side].set_visible(False)
    bars.spines["bottom"].set_color(GRID)
    table.set_xlim(0,1)
    table.set_ylim(bars.get_ylim())
    table.axis("off")
    bars.text(0,1.06,"With feedback",transform=bars.transAxes,fontsize=10.5,weight="bold",color=GREEN)
    xcols = (.41,.93)
    for x,title,color in zip(xcols,("With feedback","Passive"),(GREEN,MUTED)):
        table.text(x,1.06,title,transform=table.transAxes,ha="right",fontsize=10.3,weight="bold",color=color)
    for i,(coupled,passive,_) in enumerate(rows):
        bars.barh(i,coupled,height=.31,facecolor=PALE,edgecolor=GREEN,linewidth=.9,zorder=2)
        bars.scatter(coupled,i,s=28,color=GREEN,zorder=3)
        for x,value,color in zip(xcols,(coupled,passive),(GREEN,MUTED)):
            table.text(x,i,fmt(value)+"%",ha="right",va="center",fontsize=11.7,color=color,
                       fontfamily="DejaVu Sans Mono")
    fractions = [r[2] for r in rows]
    fig.text(.143,.116,f"The daughter term accounts for {min(fractions):.4f}%–{max(fractions):.4f}% of the model’s "
             + r"$H^2$ on these three coupled histories.", fontsize=10,color=INK)
    fig.text(.143,.054,"Near and distant locate the compensating pulse. Passive runs omit matter feedback. Each contrast uses a matched quiet control.",
             fontsize=9.5,color=MUTED)
    path=args.output or args.capsule/"coupled_radiation.png"
    fig.savefig(path,dpi=210,facecolor="white",metadata={"Software":"Matplotlib"})
    plt.close(fig)


if __name__=="__main__":
    main()
