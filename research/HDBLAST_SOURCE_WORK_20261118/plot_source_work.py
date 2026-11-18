"""Compare pulse-minus-quiet free-field work and change in modeled daughter radiation."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter


INK = "#1E303B"
MID = "#465965"
GRID = "#D7DEE1"
FIELD = "#205A7B"
WORK = "#B47731"
LIGHT_WORK = "#F4DEC5"
RADIATION = "#397364"
CASES = ["near_plus", "near_minus", "far_plus"]
NAMES = ["Near, positive", "Near, negative", "Distant, positive"]


def read(path):
    return json.loads(path.read_bytes())


def signed(value, decimals):
    return f"{value:+.{decimals}f}".replace("-", "−")


def style(ax, title, sub, xlim, ticks, xlabel):
    ax.set_axisbelow(True)
    ax.grid(axis="x", color=GRID, lw=.7)
    ax.axvline(0, color=MID, lw=.95, zorder=1.5)
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.set_ylim(2.49, -.52)
    ax.tick_params(axis="both", length=0, labelcolor=MID, pad=5)
    ax.set_xlim(*xlim)
    ax.set_xticks(ticks)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, pos: f"{v:g}".replace("-", "−")))
    ax.set_yticks(range(3), NAMES)
    ax.text(0, 1.28, title, transform=ax.transAxes, ha="left", va="baseline",
            size=14.5, weight="bold", color=INK)
    ax.text(0, 1.14, sub, transform=ax.transAxes, ha="left", va="baseline", size=9.5, color=MID)
    ax.set_xlabel(xlabel, color=INK, fontsize=10.7, labelpad=8)


def endpoint(ax, y, x, suffix="", digits=3):
    label = signed(x, digits)+suffix
    off = 5 if x>=0 else -5
    ax.annotate(label, (x,y), xytext=(off,0), textcoords="offset points",
                ha="left" if x>=0 else "right", va="center", fontsize=9.4, color=INK)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capsule", type=Path, default=Path(__file__).resolve().parent,
                        help="Capsule root containing ledger/ and decay/")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    ledger = read(args.capsule/"ledger/REFINED.json")
    radiation = {r["case"]:r for r in read(args.capsule/"decay/CHECK.json")}
    if (ledger["start"],ledger["end"],ledger["band"],ledger["grid"]) != (.8,2.10,[0,8],"fine"):
        raise ValueError("Expected the published fine-grid rolling work experiment")
    changes = ledger["pulse_minus_quiet"]
    fig, (left, right) = plt.subplots(1,2,figsize=(11.5,4.8),gridspec_kw={"width_ratios":[1.25,1]})
    fig.subplots_adjust(left=.113,right=.966,bottom=.257,top=.738,wspace=.4)
    style(left, "Free field", "Second-order adiabatic initial state", (-.82,2.31),
          [-.5,0,.5,1,1.5,2], r"Change from quiet per comoving volume ($H_0^4$)")
    style(right, "Field with decay", "Instantaneous ground state and empty radiation reservoir",
          (-.45,1.34),[-.25,0,.25,.5,.75,1,1.25], "Change in physical daughter density (%)")
    left.legend(handles=[Patch(facecolor=FIELD,label="Final field energy"),
                         Patch(facecolor=LIGHT_WORK,edgecolor=WORK,label="Work from varying mass")],
                loc="lower left",bbox_to_anchor=(-.009,.975),ncol=2,frameon=False,labelcolor=INK,
                fontsize=8.6,handlelength=1.15,handleheight=.7,handletextpad=.45,
                columnspacing=1.25,borderaxespad=0)
    for i, arm in enumerate(CASES):
        energy = changes[arm]["energies"]["full"]["final"]
        work = changes[arm]["work"]["full"]["mass_work"]
        delta = radiation[arm]["change_percent"]
        quiet = radiation["quiet"]["physical_radiation"]
        if abs(100*(radiation[arm]["physical_radiation"]/quiet-1)-delta)>1e-9:
            raise ValueError("Inconsistent daughter-radiation baseline")
        left.barh(i-.15,energy,height=.22,color=FIELD,zorder=2)
        left.barh(i+.15,work,height=.22,facecolor=LIGHT_WORK,edgecolor=WORK,linewidth=1.05,zorder=2)
        endpoint(left,i-.15,energy)
        endpoint(left,i+.15,work)
        right.barh(i,delta,height=.36,color=RADIATION,zorder=2)
        endpoint(right,i,delta, digits=3)
    fig.text(.113,.092, r"Both calculations use $0\leq k/(a_{*,0}H_0)\leq8$ and end at $H_0\tau=2.10$ on the same shell histories.",
             fontsize=10,color=INK)
    fig.text(.113,.039,"Near and distant refer to the compensating disturbance. Each change uses the quiet control for that panel’s preparation.",
             fontsize=9.5,color=MID)
    target=args.output or args.capsule/"source_work.png"
    fig.savefig(target,dpi=200,facecolor="white",metadata={"Software":"Matplotlib"})
    plt.close(fig)


if __name__=="__main__":
    main()
