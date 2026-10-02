"""Export scientific figures from completed runs; no image synthesis."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
base=Path("outputs/pv_actual")
report=json.loads((base/"pv_actual.json").read_text())
figdir=Path("outputs/figures");figdir.mkdir(parents=True,exist_ok=True)
primary=np.load(base/"primary_curves.npz",allow_pickle=False)["curves"]
t=primary[:,0]
fig,axes=plt.subplots(3,1,figsize=(8,8),sharex=True)
for ax,label,col,direct in zip(axes,[r"$\rho$",r"$p$",r"$S=\langle\chi^2\rangle/2$"],[2,3,4],[9,10,11]):
    ax.plot(t,primary[:,col],label=r"matched PV, $\Lambda=16$")
    ax.plot(t,primary[:,direct],"--",label="direct physical-mode prescription")
    ax.axvline(0,color=".6",lw=.8);ax.set_ylabel(label);ax.grid(alpha=.2)
axes[0].legend(fontsize=8)
axes[-1].set_xlabel(r"$t/T$")
fig.suptitle("Smooth external mass crossing: matched quantum source\nFlat prescribed background; finite cutoffs, specified matching")
fig.tight_layout();fig.savefig(figdir/"matched_sources.svg");plt.close(fig)
fig,axes=plt.subplots(1,2,figsize=(11,4),sharey=True)
for name,label in [("lambda4",r"$\Lambda=4$"),("lambda8",r"$\Lambda=8$"),("primary",r"$\Lambda=16$")]:
    path=base/(name+"_curves.npz")
    if not path.exists():continue
    z=np.load(path,allow_pickle=False)["curves"]
    axes[0].plot(z[:,0],z[:,12],label=label)
    axes[1].plot(z[:,0],z[:,3],label=label)
for ax in axes:
    ax.set_xlim(-2,2);ax.set_xlabel(r"$t/T$");ax.grid(alpha=.2);ax.legend(fontsize=8)
axes[0].set_title("Potential matching only: negative control")
axes[1].set_title("Potential and curvature matching")
axes[0].set_ylabel("Pressure")
fig.suptitle("Energy conservation alone does not test flat-space pressure")
fig.tight_layout();fig.savefig(figdir/"pressure_matching.svg");plt.close(fig)
fig,ax=plt.subplots(figsize=(8,4))
for case in report["cases"]:
    vals=np.load(base/(case["name"]+"_curves.npz"),allow_pickle=False)
    k=vals["k"];nn=vals["n_out"][0]
    if case["name"] in ("primary","quadrature48","tight"):
        mask=k<8
        ax.semilogy(k[mask],np.maximum(nn[mask],1e-30),".",ms=3,label=case["name"])
data=np.load(base/"primary_curves.npz",allow_pickle=False)
k=data["k"];mask=k<8
ax.semilogy(k[mask],np.maximum(data["n_exact"][0][mask],1e-30),"k-",lw=1,label="exact scattering occupation")
ax.set_xlabel(r"$kT$");ax.set_ylabel(r"$n_k$");ax.set_title("Exact-state scattering control, physical field")
ax.grid(alpha=.2);ax.legend(fontsize=8);fig.tight_layout()
fig.savefig(figdir/"scattering_control.svg");plt.close(fig)
for p in sorted(figdir.glob("*.svg")):
    print("HDBLAST_PV_FIGURE="+json.dumps(dict(name=p.name,svg=p.read_text()),separators=(",",":")))
