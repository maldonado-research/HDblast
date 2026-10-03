"""Export compact permanent numerical curves and scientific SVG figures."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

base=Path("outputs/actual_modes")
report=json.loads((base/"actual_modes.json").read_text())
primary=next((x for x in report["cases"] if x["name"]=="primary"),None)
if primary is None:
    print("HDBLAST_PUBLIC_EXPORT="+json.dumps({"status":"primary case unavailable"}))
    raise SystemExit(0)
figdir=Path("outputs/figures")
figdir.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({"font.size":10,"svg.fonttype":"none"})
fig,ax=plt.subplots(figsize=(7,4.4))
for name in ("primary","tight","quadrature64","cubic","hamiltonian_in"):
    row=next((x for x in report["cases"] if x["name"]==name),None)
    if row:
        ax.semilogy(row["kappa"],np.maximum(row["spectrum_n"],1e-16),label=name.replace("_"," "),lw=1.3)
kap=np.linspace(.05,3,200)
ax.semilogy(kap,np.exp(-np.pi*kap**2),"k--",label="single-crossing Gaussian",lw=1.2)
ax.set(xlabel="Momentum at crossing / sqrt(q)",ylabel="Specified endpoint occupation",
       title="HDBLAST Y=0 prescribed geometry: finite-band scalar modes")
ax.set_ylim(bottom=1e-14)
ax.legend(fontsize=8,ncol=2)
ax.grid(alpha=.2)
fig.tight_layout()
fig.savefig(figdir/"endpoint_spectrum.svg")
plt.close(fig)
with np.load(base/"primary_diagnostics.npz",allow_pickle=False) as z:
    cols=["t","a","H","phi","v","delta_rho","delta_p","delta_J",
          "paired_scalar_work","paired_pressure_work","paired_balance",
          "particle_gas_ward_defect"]
    rows=np.column_stack([z[c] for c in cols])
    curves={"columns":cols,"rows":rows.tolist(),
            "meaning":"Exact same-background state difference on a finite fixed comoving band; not the absolute renormalized source."}
    (base/"primary_curves.json").write_text(json.dumps(curves,separators=(",",":"))+"\n")
    t=z["t"];energy=z["a"]**3*z["delta_rho"];energy-=energy[0]
    fig,axes=plt.subplots(2,1,figsize=(7,6),sharex=True)
    axes[0].plot(t,energy,label="Change in a^3 Delta rho",lw=1.7)
    axes[0].plot(t,z["paired_scalar_work"],label="Scalar work",lw=1.2)
    axes[0].plot(t,z["paired_pressure_work"],label="Pressure work",lw=1.2)
    axes[0].plot(t,z["paired_scalar_work"]+z["paired_pressure_work"],"--",label="Sum of work",lw=1.2)
    axes[0].set_ylabel("Comoving energy / H0^4")
    axes[0].legend(fontsize=8,ncol=2)
    axes[0].set_title("Relative stress and its independently integrated work")
    axes[1].plot(t,z["paired_balance"],label="Energy minus work residual")
    axes[1].set(xlabel="Proper time H0 tau",ylabel="Residual / H0^4")
    for ax in axes:ax.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(figdir/"relative_work.svg")
    plt.close(fig)
print("HDBLAST_CURVES_JSON="+json.dumps(curves,separators=(",",":"),allow_nan=False),flush=True)
for path in sorted(figdir.glob("*.svg")):
    print("HDBLAST_FIGURE_JSON="+json.dumps({"name":path.name,"svg":path.read_text()},separators=(",",":")),flush=True)
