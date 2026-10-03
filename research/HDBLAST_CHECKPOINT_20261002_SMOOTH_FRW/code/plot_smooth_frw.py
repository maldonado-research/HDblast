"""Render only executed, finite saved arrays; do not recompute mode evolution."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run",type=Path)
    args=parser.parse_args()
    summary=json.loads((args.run/"summary.json").read_text())
    output=args.run/"figures"
    output.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,
                        "figure.dpi":140,"savefig.dpi":200})
    fig,axes=plt.subplots(3,2,figsize=(10,8),sharex=True,constrained_layout=True)
    ledger,ledger_axes=plt.subplots(2,2,figsize=(10,6),sharex=True,constrained_layout=True)
    trace,trace_axes=plt.subplots(2,2,figsize=(10,6),sharex=True,constrained_layout=True)
    for column,A in enumerate(("0.0","0.2")):
        record=summary["amplitudes"][A]
        selected=record["selected"]
        K=selected.split("_K")[-1]
        prefix=f"K{K}_"
        with np.load(args.run/selected/"arrays.npz",allow_pickle=False) as arrays:
            t=arrays["t"]
            for row,quantity in enumerate(("rho","p","Q")):
                ax=axes[row,column]
                ax.plot(t,arrays[prefix+quantity],color="#16697a",lw=1.5,label="full subtracted")
                for smaller in (24,48,96):
                    if str(smaller)!=K and f"K{smaller}_{quantity}" in arrays:
                        ax.plot(t,arrays[f"K{smaller}_{quantity}"],lw=.75,alpha=.5,label=f"K={smaller}")
                ax.axvline(.5,color="gray",lw=.7,ls=":")
                ax.axvline(1,color="gray",lw=.7,ls="--")
                ax.set_ylabel({"rho":r"$\rho$","p":r"$p$","Q":r"$Q$"}[quantity])
                if row==0:
                    ax.set_title(f"A={A}; K={K}; matrix {record['status']}")
                    ax.legend(fontsize=7)
                if row==2:
                    ax.set_xlabel(r"conformal time $\eta$")
            E=arrays[prefix+"conformal_energy"]
            ledger_axes[0,column].plot(t,E-E[0],label=r"$a^4\rho-a(0)^4\rho(0)$")
            ledger_axes[0,column].plot(t,arrays[prefix+"ledger"],ls="--",label="integrated paired work")
            ledger_axes[0,column].set_title(f"A={A}")
            ledger_axes[0,column].set_ylabel("conformal energy change")
            ledger_axes[0,column].legend(fontsize=8)
            ledger_axes[1,column].plot(t,arrays[prefix+"ledger_residual"],label="energy minus work")
            ledger_axes[1,column].set_ylabel("integrated exchange residual")
            ledger_axes[1,column].set_xlabel(r"$\eta$")
            trace_axes[0,column].plot(t,arrays[prefix+"trace"],label=r"$-\rho+3p$")
            trace_axes[0,column].plot(t,arrays[prefix+"trace_sampled_801"],ls="--",label="sampled Q derivatives + finite-K remainder")
            trace_axes[0,column].set_title(f"A={A}")
            trace_axes[0,column].set_ylabel("trace")
            trace_axes[0,column].legend(fontsize=7)
            trace_axes[1,column].plot(t,arrays[prefix+"trace"]-arrays[prefix+"trace_sampled_801"],label="sampled trace residual")
            trace_axes[1,column].set_ylabel("trace residual")
            trace_axes[1,column].set_xlabel(r"$\eta$")
    fig.suptitle("Prescribed smooth FRW control: complete stress and paired current")
    ledger.suptitle("Exchange ledger including mass-squared source work")
    trace.suptitle("Pressure validation using the finite-cutoff trace relation")
    for figure,name in ((fig,"observables"),(ledger,"exchange"),(trace,"trace")):
        figure.savefig(output/f"{name}.png")
        figure.savefig(output/f"{name}.pdf")
        plt.close(figure)
    print(output)


if __name__=="__main__":
    main()
