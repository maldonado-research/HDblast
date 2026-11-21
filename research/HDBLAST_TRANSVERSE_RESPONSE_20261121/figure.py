"""Plot two tests of the origin of transverse particle momentum at 1% electrical work."""
import json
from pathlib import Path
import runpy

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent
caps=json.loads((ROOT/"causal/CURRENT_CAP.json").read_bytes())["runs"]
runs=json.loads((ROOT/"particles/RESULTS.json").read_bytes())["runs"]
thermal=runpy.run_path(str(ROOT/"global/clumping.py"))["ThermalCost"]
channels=["dedicated","shared"]
colors=["#146879","#a2492c"]
fig,(left,right)=plt.subplots(1,2,figsize=(11,5.2))
fig.subplots_adjust(left=.085,right=.975,bottom=.23,top=.79,wspace=.34)
fig.suptitle("What could turn the newly born particles sideways?",fontsize=19,y=.96)
for i,ch in enumerate(channels):
    row=caps[ch+"_half"]["0.01"]
    given=row["population_rms_impulse_upper"]
    needed=thermal(row["temperature"],24576).independent_integral(row["unhit_probability_floor_under_homogeneous_hazard"])
    for x,num,color,label in [(i-.18,given,"#226e83","Most one new wave can supply"),(i+.18,needed,"#b75e24","Least the thermal target requires")]:
        left.bar(x,num,.32,color=color,label=label if i==0 else None)
        left.text(x,num+14,f"{num:.0f}",ha="center",fontsize=10)
left.set_ylim(0,810);left.set_xticks([0,1],["Dedicated","Shared"])
left.set_ylabel("RMS transverse impulse (mc)")
left.set_title("Particle currents, no initial wave",loc="left",fontsize=13)
left.legend(loc="upper center",bbox_to_anchor=(.5,-.14),frameon=False,fontsize=9)
for ch,color,m in zip(channels,colors,["o","s"]):
    matches={r["initial_peak_field"]:r["result"] for r in runs
             if r["history"]==ch+"_native" and r["resolution"]=="base" and r["feedback"]
             and r["polarization"]=="By" and r["k_physical_opening"]==50}
    values=[100*matches[s]["transverse_cdf_distance"] for s in [0.,.005,.04]]
    right.plot(range(3),values,color=color,marker=m,lw=2,label=ch.capitalize())
right.axhline(1,color="#40474c",ls="--",lw=1)
right.text(.03,4,"One-percentage-point target",fontsize=9)
right.set_xlim(-.08,2.12);right.set_ylim(0,107)
right.set_xticks(range(3),["0","0.005","0.04"])
right.set_xlabel("Initial peak magnetic field (wall units)")
right.set_ylabel("Largest gap in cumulative fraction (points)")
right.set_title("A supplied wave, followed particle by particle",loc="left",fontsize=13)
right.legend(frameon=False,loc="upper right",fontsize=9)
for ax in [left,right]:
    ax.spines[["right","top"]].set_visible(False)
    ax.set_axisbelow(True);ax.grid(axis="y",alpha=.18)
fig.text(.085,.03,"Both panels use the 1% electrical-work marker. The supplied standing wave has initial physical wavenumber 50/L.",fontsize=9)
fig.savefig(ROOT/"transverse-response.png",dpi=165,facecolor="white")
