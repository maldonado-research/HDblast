"""Compare the finite-energy spectrum with the soft limit and show formation windows."""
import argparse
import json
from pathlib import Path
import runpy

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(base, output):
    emission = runpy.run_path(str(base/"hard/emission.py"))
    model_mass = emission["MASS"]
    nodes, weights = np.polynomial.legendre.leggauss(40)
    angle, weights = (nodes+1.)/2., weights/2.
    fig, (left, right) = plt.subplots(1,2,figsize=(10.0,4.25),layout="constrained")
    colors = ("#196584", "#aa512b")
    for multiple, color in zip((1.,2.), colors):
        s = (multiple*model_mass)**2
        xi = np.linspace(.001,.99,161)
        x = (1.-4./s)*xi[:,None]
        rate = sum(emission["pair_rate"](s,x,angle,weights,leg) for leg in (False,True))
        ratio = x[:,0]*rate/(2.*emission["soft_coefficient"](s))
        left.plot(xi,ratio,c=color,lw=2.2,label=f"{multiple:g} × resonance mass")
    left.axhline(1,ls=":",lw=1,c="#555555")
    left.axvline(.75,ls="--",lw=.9,c="#aa512b",alpha=.6)
    left.set(xlim=(0,1), xlabel="Photon energy / kinematic maximum", ylabel="Spectrum / leading soft-photon extrapolation",
             title="Energy radiated in an isolated collision")
    left.legend(title="Pair energy before collision",loc="upper left",fontsize=9,title_fontsize=9,frameon=False)
    result = json.loads((base/"formation/FORMATION.json").read_bytes())
    run = next(r for r in result["runs"] if r["label"]=="t7.05000-p14-n8192-r91283-particle-pole")
    rows = run["band_energy_fractions"]
    x = np.arange(len(rows)+1)
    for idx,(key,label) in enumerate([("ten_radian_time","Spacing between collisions"),
                                      ("integrated_force_envelope_and_ten_radians","Spacing and field-bending margin")]):
        hard_key = "ten_radian_spacing_fraction" if idx==0 else "joint_both_sides_fraction"
        heights = [100*r[key] for r in rows]+[100*run["above_quarter_hard_band"][hard_key]]
        right.bar(x+(-.18 if idx==0 else .18),heights,width=.35,label=label,color=colors[idx])
    right.set_xticks(x,["0.001–0.01","0.01–0.1","0.1–1","1–10","25–100\n(recoil)"])
    right.set(ylim=(0,122),ylabel="Estimated band energy passing check (%)",
              xlabel="Photon band (% of kinematic maximum)",title="Separation on the shell at time 7.05")
    right.set_yticks(range(0,101,20))
    right.legend(loc="upper left",fontsize=8.5,frameon=False)
    for ax in (left,right):
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(axis="y",lw=.5,alpha=.2)
        ax.set_axisbelow(True)
    fig.savefig(output,dpi=175)
    plt.close(fig)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=Path(__file__).with_name("collision-light.png"))
    args=parser.parse_args()
    plot(Path(__file__).resolve().parent,args.output)
