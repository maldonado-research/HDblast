"""Plot actual audit arrays and exact/verified oscillator controls."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
out=Path("outputs");dest=out/"figures";dest.mkdir(exist_ok=True)
report=json.loads((out/"geometry_audit.json").read_text())
z=np.load(out/"plot_reconstructions.npz",allow_pickle=False)
cols=list(z["columns"]);fig,axes=plt.subplots(2,1,figsize=(8,7),sharex=True)
for label,curve in zip(z["labels"],z["curves"]):
    name=str(label).split(":")[-1]
    axes[0].plot(curve[:,0],curve[:,cols.index("U")],label=name,lw=1.1)
    axes[1].plot(curve[:,0],curve[:,cols.index("p_R2_regular")],label=name,lw=1.1)
axes[0].set_ylabel("Conformal potential U");axes[0].set_yscale("symlog",linthresh=10)
axes[1].set_ylabel("Unit R² pressure: regular part");axes[1].set_yscale("symlog",linthresh=10)
axes[1].set_xlabel("s = H0 proper time")
axes[0].legend(fontsize=8);axes[0].set_title("Archived fine background: reconstruction sensitivity")
for ax in axes:ax.grid(alpha=.2)
fig.text(.5,.01,"R² coefficients are not fitted; knot-supported distributions are excluded from these curves.",ha="center",fontsize=8)
fig.tight_layout(rect=[0,.03,1,1]);fig.savefig(dest/"reconstruction_sensitivity.svg");plt.close(fig)
s=json.loads((out/"smooth_spectra.json").read_text());c=json.loads((out/"smooth_controls.json").read_text())
fig,ax=plt.subplots(figsize=(8,5))
ax.loglog(s["k"],s["sudden"],"k--",label="Abrupt step (exact)")
colors={}
for i,(eps,n) in enumerate(s["smooth"].items()):
    color="C"+str(i);colors[float(eps)]=color
    ax.loglog(s["k"],np.maximum(n,1e-300),color=color,label="Smooth ε="+eps+" (exact)")
for case in c["cases"]:
    ax.scatter(case["k"],case["occupation_numeric"],color=colors[case["eps"]],marker="o",s=22)
ax.set_ylim(1e-28,1);ax.set_xlabel("Momentum k");ax.set_ylabel("Occupation |β|²")
ax.set_title("Independent toy control: smoothing removes the abrupt UV tail")
ax.legend(fontsize=8);ax.grid(alpha=.2);fig.tight_layout();fig.savefig(dest/"smooth_vs_sudden.svg");plt.close(fig)
rr=[x for x in report["reconstructions"] if x["kind"]=="hermite"]
fig,ax=plt.subplots(figsize=(7,4.5));positions=np.arange(len(rr))
ax.bar(positions-.18,[x["energy_tail"]["diagonal_log_coefficient"] for x in rr],.36,label="Interior Hermite knots")
ax.bar(positions+.18,[x["state_boundaries"]["initial_physical_WKB_log_coefficient"] for x in rr],.36,label="Initial WKB state mismatch")
ax.set_xticks(positions,["Fine archive","Coarse archive"]);ax.set_yscale("log")
ax.set_ylabel("Separate formal coefficients of ln K");ax.set_title("Two distinct obstructions to an unchanged UV extension")
ax.legend(fontsize=8);ax.grid(axis="y",alpha=.2);fig.tight_layout();fig.savefig(dest/"uv_obstructions.svg");plt.close(fig)
print("HDBLAST_FRW_FIGURES="+json.dumps([str(p) for p in sorted(dest.glob("*.svg"))]))
for p in sorted(dest.glob("*.svg")):
    print("HDBLAST_FRW_SVG_JSON="+json.dumps({"path":str(p),"content":p.read_text()},separators=(",",":")))

