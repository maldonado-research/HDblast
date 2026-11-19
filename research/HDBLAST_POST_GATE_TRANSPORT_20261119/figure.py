"""Generate the particle-spectrum and first-collision figure for this study."""
import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import gammainc

ROOT = Path(__file__).resolve().parent


def render(output):
    with (ROOT / "daughters/SPECTRA.csv").open() as f:
        rows = [r for r in csv.DictReader(f) if r["arm"]=="quiet" and r["band"]=="full"]
    edge = np.r_[[float(r["k_low"]) for r in rows], float(rows[-1]["k_high"])]
    counts = np.array([float(r["one_sign_comoving_number"]) for r in rows])
    counts /= counts.sum()
    tests = json.loads((ROOT / "survival/SURVIVAL.json").read_text())["histories"]
    quiet = tests["quiet"]
    z, theta = quiet["thermal_fugacity"], quiet["thermal_comoving_scale"]
    j=np.arange(1,13)
    weight=(-1.)**(j-1)*z**j/j**3
    eq_cdf=np.sum(weight[:,None]*gammainc(3,j[:,None]*edge/theta),axis=0)/sum(weight)
    env=json.loads((ROOT/"envelope/ENVELOPE.json").read_text())["histories"]
    arms=("quiet","plus","minus","far")
    yields=[env[arm]["readouts"]["2.16"]["evolving_targets"]["first_collision_fraction_number_upper"]*1.e6 for arm in arms]
    fig, (left,right)=plt.subplots(1,2,figsize=(10,4.25),layout="constrained",width_ratios=[1.35,1])
    left.stairs(counts/np.diff(edge),edge,color="#265a8d",lw=2,label="Simulated particles")
    left.stairs(np.diff(eq_cdf)/np.diff(edge),edge,color="#ad512e",lw=2,label="Thermal reference")
    left.set(xlabel=r"Comoving momentum $k/H_0$",ylabel="Fraction per unit momentum",title="Quiet shell at equal particle count and energy",xlim=(0,64),ylim=(0,None))
    left.legend(frameon=False)
    right.bar(["Quiet","Near +","Near −","Far"],yields,color=["#365472","#5c8a9e","#477d96","#79a6aa"],width=.63)
    for k,val in enumerate(yields):right.text(k,val+1,f"{val:.1f}",ha="center",color="#273e50",fontsize=10)
    right.set(title="Upper bound on first collisions",ylabel="Particles per million",ylim=(0,83))
    for ax in (left,right):
        ax.spines[["right","top"]].set_visible(False)
        ax.grid(axis="y",alpha=.19)
        ax.set_axisbelow(True)
    fig.savefig(output,dpi=170)
    plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True)
    render(p.parse_args().output)
