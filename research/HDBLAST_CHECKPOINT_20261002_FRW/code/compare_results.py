"""Post hoc descriptive comparisons; no registered acceptance criterion changed."""
import json
from pathlib import Path
import numpy as np
out=Path("outputs")
g=json.loads((out/"geometry_audit.json").read_text())
d=json.loads((out/"smooth_controls.json").read_text());r=json.loads((out/"independent_smooth.json").read_text())
ra={(x["eps"],x["k"]):x for x in r["cases"]}
modes=[]
for x in d["cases"]:
    y=ra[(x["eps"],x["k"])]
    modes.append(dict(eps=x["eps"],k=x["k"],occupation_abs_difference=abs(x["occupation_numeric"]-y["occupation_numeric"])))
z=np.load(out/"plot_reconstructions.npz",allow_pickle=False);cols=list(z["columns"])
curves={str(k).split(":")[-1]:v for k,v in zip(z["labels"],z["curves"])}
a,b=curves["quintic_C4"],curves["septic_C6"]
comp={}
for i,n in enumerate(cols):
    if i==0:continue
    delta=np.abs(a[:,i]-b[:,i]);j=int(np.argmax(delta))
    comp[n]=dict(max_abs=float(delta[j]),time_of_max=float(a[j,0]),
                 relative_to_C4_sampled_max=float(delta[j]/max(np.max(np.abs(a[:,i])),1e-30)))
hermites=[x for x in g["reconstructions"] if x["kind"]=="hermite"]
outd=dict(status="EXECUTED_POST_HOC_DESCRIPTIVE_COMPARISONS",source_ref=g["source_ref"],
    DOP853_vs_Radau=modes,quintic_vs_septic_fine_full_uniform_grid=comp,
    initial_to_knot_log_coefficient_ratio={x["tag"]:x["state_boundaries"]["initial_physical_WKB_log_coefficient"]/x["energy_tail"]["diagonal_log_coefficient"] for x in hermites},
    caution="Sampled maxima only; R2 curves use unit coefficient and omit knot distributions. No physical FRW source is computed.")
(out/"posthoc_comparisons.json").write_text(json.dumps(outd,indent=2,allow_nan=False)+"\n")
print("HDBLAST_FRW_COMPARISONS_JSON="+json.dumps(outd,separators=(",",":"),allow_nan=False))

