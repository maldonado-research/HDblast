"""Plot proper-volume growth and the generated field from the nonlinear solutions."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import numpy as np


def main(output):
    data=json.loads((Path(__file__).resolve().parent/"evolution/results.json").read_text())
    fig, axes=plt.subplots(1,2,figsize=(9.1,4.1))
    colors=["#145a77","#a34d22"]
    for (branch, model), color in zip(data["branches"].items(),colors):
        pair=[model["runs"][code]["samples"] for code in ("0.001","-0.001")]
        times=sorted(pair[0],key=float)
        changes=[[points[t]["mean_change"] for t in times] for points in pair]
        growth=[]
        for pts in pair:
            v0=pts["0"]["volume_over_background"]
            growth.append([1e6*(pts[t]["volume_over_background"]/v0-1) for t in times])
        x=list(map(float,times))
        label="Positive minimum" if branch=="phi_m=0.5" else "Negative minimum"
        axes[0].plot(x,np.mean(growth,axis=0),marker="o",ms=3.5,color=color,label=label)
        axes[1].plot(x,1e6*np.mean(changes,axis=0),marker="o",ms=3.5,color=color)
    axes[0].set_title("Change in volume growth")
    axes[1].set_title("Shift in the average field")
    axes[0].set_ylabel("Difference (parts per million)")
    axes[1].set_ylabel("Change from initial mean ($10^{-6}M_{\\rm Pl}$)")
    for ax in axes:
        ax.axhline(0,color="#6c7477",lw=.8)
        ax.grid(color="#e3e7e9",lw=.6)
        ax.spines[["right","top"]].set_visible(False)
        ax.set_xlabel("Proper time from the bounce (model units)")
        ax.set_xticks([0,.03,.06,.09,.12])
    fig.subplots_adjust(left=.10,right=.98,bottom=.24,top=.83,wspace=.34)
    fig.legend(*axes[0].get_legend_handles_labels(),loc="lower center",ncol=2,frameon=False)
    fig.savefig(output,dpi=165,bbox_inches="tight")
    plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True)
    main(p.parse_args().output)
