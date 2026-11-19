"""Plot the change from the classical energy estimate and the accuracy of the reduced quantum bath."""
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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE/"memory.png")
    args = parser.parse_args()
    values = json.loads((HERE/"quantum/SPECTRUM.json").read_text())["arms"]["quiet"]["all_preparations_conformal"]["retain"]
    with (HERE/"source/COMPARISON.csv").open() as stream:
        rows = [r for r in csv.DictReader(stream) if r["scheme"]=="static" and r["field"]=="rad"]
    fig,(energy,accuracy)=plt.subplots(1,2,figsize=(9.2,4.05),layout="constrained",
                                      gridspec_kw={"width_ratios":[.9,1.3]})
    y=[values["symmetrized"],values["commutator"],values["emission"]]
    energy.bar([0,1,2], [y[0],-y[1],y[2]],bottom=[0,y[2],0],width=.63,
               color=["#47789a","#c79565","#346e61"],zorder=3)
    energy.plot([.32,.68],[y[0]]*2,color="#777",lw=.8)
    energy.plot([1.32,1.68],[y[2]]*2,color="#777",lw=.8)
    for x,(v,cy) in enumerate(zip(y,[y[0]/2,y[2]-y[1]/2,y[2]/2])):
        energy.text(x,cy,f"{v:+.3f}" if x==1 else f"{v:.3f}",ha="center",va="center",
                    color="white" if x!=1 else "#2c2119",fontweight="bold")
    energy.set_xticks([0,1,2],["Driven\nestimate","Quantum\ncorrection","Emitted\nenergy"])
    energy.set_ylabel("High-band energy (conformal units)")
    energy.set_title("What the quiet shell emits",loc="left",fontweight="bold",size=12)
    for key,label,color,style in [("quiet","Quiet","#376d8b","o-"),("minus","Near negative","#aa5d35","s-"),
                                  ("plus","Near positive","#648740","^-."),("far","Distant pulse","#74658e","D:")]:
        orders=[1,2,4,8,16]
        family={r["bath"]:abs(float(r["readout_high_relative_error"]))*100 for r in rows if r["case"]=="fine_"+key}
        accuracy.loglog(orders,[family["N"+str(n)] for n in orders],style,color=color,
                        label=label,ms=4,linewidth=1.6)
    accuracy.axhline(1,color="#666",lw=.8,ls="--")
    accuracy.set_xticks([1,2,4,8,16],[1,2,4,8,16])
    accuracy.set_ylabel("Upper-band energy error at readout (%)")
    accuracy.set_xlabel("Oscillator frequencies per momentum in the upper band")
    accuracy.set_title("How many frequencies are needed?",loc="left",fontweight="bold",size=12)
    accuracy.legend(frameon=False,fontsize=8,loc="lower left")
    for ax in (energy,accuracy):
        ax.grid(axis="y",color="#ddd",lw=.6)
        ax.set_axisbelow(True)
        ax.spines[["top","right"]].set_visible(False)
    fig.savefig(args.output,dpi=170,facecolor="white")


if __name__=="__main__":
    main()
