"""Draw the three homogeneous expansion histories from the released traces."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter

ROOT = Path(__file__).resolve().parent

def read(path):
    product=json.loads((ROOT/path).read_bytes())
    run=product["runs"]["primary"]
    start=run.get("hot_from_radius",run["first_ninety_intrinsic_radius"])
    peak=run["maximum_tau"]
    rows=[dict(zip(product["trace_columns"],row,strict=True)) for row in run["trace"]]
    expanding=[r for r in rows if run["exit_complete_tau"]<=r["t"]<=peak and r["e"]>0]
    return ([r["a"]/start for r in expanding],
            [r["intrinsic_random"]/r["e"] for r in expanding],
            run["maximum_radius"]/start)

if __name__ == "__main__":
    runs=[("Previous release","prior/DECAY.json","#77716B"),
          ("Two-stage scalar release","exit/EXIT.json","#087C75"),
          ("Parity-odd scalar release","axion/GROWTH.json","#9E4F28")]
    fig,(left,right)=plt.subplots(1,2,figsize=(11.8,4.4),gridspec_kw={"width_ratios":[1.45,1]})
    for rank,(label,path,color) in enumerate(runs):
        x,y,growth=read(path)
        left.plot(x,y,label=label,color=color,lw=2.2)
        left.scatter(x[-1],y[-1],color=color,s=20,zorder=3)
        position=2-rank
        right.plot([1,growth],[position]*2,color=color,lw=5,solid_capstyle="round")
        right.scatter([1,growth],[position]*2,s=36,color=color)
        right.text(growth+.07,position,f"{growth:.3f}×",color=color,va="center",fontweight="bold")
    left.axhline(.9,color="#262626",linestyle="--",lw=.9);left.set_ylim(-.015,1.035)
    left.set_xlabel("Radius / radius at first 90% crossing");left.set_ylabel("Share of added energy in particle internal motion")
    left.yaxis.set_major_formatter(PercentFormatter(1))
    left.legend(loc="lower right",frameon=False,fontsize=9)
    right.axvline(2,color="#595959",linestyle="--",lw=.95)
    right.set_xlabel("Radius / radius at first 90% crossing")
    right.set_yticks([2,1,0],[r[0] for r in runs]);right.set_xlim(.9,3.95);right.set_ylim(-.55,2.55)
    left.set_title("Heating during expansion",loc="left",fontweight="bold")
    right.set_title("Growth during the hot interval",loc="left",fontweight="bold")
    for ax in (left,right):
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(axis="x",color="#DFE1DD",lw=.5);ax.set_axisbelow(True)
    fig.tight_layout(pad=1.4)
    fig.savefig(ROOT/"hot-expansion.png",dpi=180,bbox_inches="tight")
