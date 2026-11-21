"""Draw the angular crowding example and both history-integrated loss estimates."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

HERE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit",type=Path,default=HERE/"AUDIT.json")
    parser.add_argument("--directional",type=Path,default=HERE.parent/"science/DIRECTIONAL.json")
    parser.add_argument("--output",type=Path,default=HERE/"directions.png")
    opt=parser.parse_args()
    check=json.loads(opt.audit.read_text())["histories"]
    main_results=json.loads(opt.directional.read_text())["runs"]
    shot=check["mass_4e12"]["half_step"]["peak_directional_snapshot"]
    c=shot["shell_speed_over_requested"]
    d=shot["force_speed_magnitude_over_requested"]
    u=np.linspace(-1.,1.,1200)
    flux=np.abs(c-d*u)
    blue, orange, black = "#205575", "#bd4b32", "#1c2c33"
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":11,"text.color":black,
                         "axes.labelcolor":black,"xtick.color":black,"ytick.color":black,
                         "axes.spines.top":False,"axes.spines.right":False,
                         "axes.spines.left":False,"axes.edgecolor":"#9ba5a5"})
    fig, (ax, data_ax)=plt.subplots(1,2,figsize=(11.6,4.55),gridspec_kw={"width_ratios":[1.29,1.]})
    fig.subplots_adjust(.075,.21,.97,.73,wspace=.33)
    fig.text(.075,.945,"Particle states can run short in particular directions",
             fontsize=19,fontweight="semibold",ha="left",va="top")
    fig.text(.075,.858,"Sharp, collisionless pair creation on two simulated wall histories",
             fontsize=11.5,ha="left",va="top",color="#445761")
    pooled=(c*c+d*d)/(2*d) if d>abs(c) else abs(c)
    ax.axhline(1.,color=black,ls="--",lw=1.)
    ax.fill_between(u,flux,1.,where=flux<1.,interpolate=True,color="#f1b5a4",alpha=.85)
    ax.plot(u,flux,color=blue,lw=2.6)
    ax.set(xlim=(-1.,1.),ylim=(0,3.25),yticks=[0,1,2,3],xticks=[-1,0,1],
           xticklabels=["Opposite","Across","Along"],xlabel="Emission direction relative to the electric force",
           ylabel="Available rate / requested rate")
    ax.tick_params(axis="y",length=0)
    ax.grid(axis="y",alpha=.17,zorder=-5)
    ax.set_title("One instant near a field reversal",loc="left",fontweight="semibold",fontsize=13,pad=12)
    ax.annotate(f"Angular average   {pooled:.2f}",xy=(-.94,1.5),fontsize=10,color=blue)
    ax.annotate(f"{100*shot['loss']:.1f}% of requested births\nmissed at this instant",
                xy=(.64,.42),xytext=(.13,1.85),fontsize=10,color=orange,
                arrowprops=dict(arrowstyle="->",lw=1.1,color=orange))
    data_ax.set_title("Births missed over time",loc="left",fontweight="semibold",fontsize=13,pad=12)
    data_ax.set_xscale("symlog",linthresh=.001,linscale=.45)
    data_ax.set_xlim(-.00024,1.91);data_ax.set_ylim(-.46,1.55)
    data_ax.set_xticks([0,.003,.01,.1,1.])
    data_ax.set_xticklabels(["0",".003",".01",".1","1"])
    data_ax.set_yticks([])
    data_ax.set_xlabel("Percent of all requested births")
    for case,n,row,off in [("mass_4e10",10,1,.14),("mass_4e12",12,0,.14)]:
        values=main_results[case]["runs"]["half_step"]["action"]["extended"]
        orig=100*values["pooled_denied_over_all_requested"]
        directional=100*check[case]["half_step"]["loss_per_total_requested"]
        data_ax.axhline(row,color="#d6dddd",lw=.8,zorder=-1)
        data_ax.plot([orig,directional],[row,row],color="#8e9d9f",lw=2.)
        data_ax.plot(orig,row,'o',ms=8,mec=blue,mew=1.7,mfc="white")
        data_ax.plot(directional,row,'s',ms=7,color=orange)
        data_ax.text(.025,row+.25,rf"$mL=4\times10^{{{n}}}$",transform=data_ax.get_yaxis_transform(),fontsize=10.5)
        data_ax.annotate("0%" if orig==0 else f"{orig:.3f}%",(orig,row),xytext=(0,-18),textcoords="offset points",ha="center",color=blue)
        data_ax.annotate(f"{directional:.3f}%" if directional>.01 else f"{directional:.4f}%",(directional,row),
                         xytext=(2,8),textcoords="offset points",ha="center",color=orange)
    data_ax.legend(handles=[Line2D([0],[0],ls="",marker="o",mfc="white",mec=blue,mew=1.5,label="Directions pooled"),
                            Line2D([0],[0],ls="",marker="s",color=orange,label="Directions checked separately")],
                   loc="upper right",frameon=False,fontsize=8.6,bbox_to_anchor=(1.025,1.05),borderaxespad=0)
    fig.text(.075,.075,r"Left: $mL=4\times10^{12}$.  Right: half-step trajectories covering 99.998% of births. Previously occupied states can further reduce supply.",
             fontsize=9,color="#445761",va="center")
    fig.savefig(opt.output,dpi=200,facecolor="white")
    plt.close(fig)


if __name__=="__main__":
    main()
