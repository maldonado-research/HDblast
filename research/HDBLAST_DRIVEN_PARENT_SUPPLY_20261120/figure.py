"""Plot energy transfer to an initially empty coherent parent field and its later action."""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent

if __name__=="__main__":
    data=json.loads((HERE/"coupled/PARENT.json").read_bytes())
    run=data["runs"]["stable"]
    trace=[dict(zip(data["trace_columns"],row,strict=True)) for row in run["trace"]]
    trace=[r for r in trace if r["t"]<=run["maximum_tau"]]
    early=[r for r in trace if r["t"]<=2.45]
    target=.0008/1.9
    fig,(left,right)=plt.subplots(1,2,figsize=(11,4.2),gridspec_kw={"width_ratios":[1,1.45]})
    left.plot([r["t"] for r in early],[r["v"] for r in early],label="Moving potential minimum",color="#85622F",lw=2)
    left.plot([r["t"] for r in early],[r["chi"] for r in early],label="Neutral field",color="#126F70",lw=2)
    left.axhline(0,color="#999999",lw=.65)
    left.set_xlabel("Wall time (model units)")
    left.set_ylabel("Field displacement")
    left.set_title("A finite push starts an oscillation",loc="left",fontweight="bold")
    left.legend(frameon=False,loc="upper center",bbox_to_anchor=(.5,-.18),fontsize=9,ncol=2)
    post=[r for r in trace if r["t"]>=run["drive_end_time"]]
    right.plot([r["t"] for r in post],[r["a"]**3*r["rho_chi"]/r["M"]/target for r in post],color="#126F70",lw=1.9)
    right.axhline(1,color="#555555",ls="--",lw=1,label="Earlier inventory target")
    right.axvline(run["opening_time"],color="#A45E35",lw=1.2,ls=":",label="Decay becomes possible")
    right.axvline(run["maximum_tau"],color="#333333",lw=.9,ls="-.",label="Expansion turns around")
    right.set_xlabel("Wall time (model units)")
    right.set_ylabel("Parent oscillation action / target")
    right.set_title("Stored action with decay switched off",loc="left",fontweight="bold")
    right.set_ylim(.85,2.35)
    right.legend(frameon=False,loc="upper center",bbox_to_anchor=(.5,-.18),fontsize=8,ncol=2)
    for ax in (left,right):
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(axis="y",lw=.45,color="#E0E3E0")
    fig.tight_layout(pad=1.4)
    fig.savefig(HERE/"parent-supply.png",dpi=165,bbox_inches="tight")
